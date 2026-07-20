"""Shared test infrastructure for the golden command-stream tests (Part 2).

The goldens freeze the exact sequence of LAMMPS commands each driver path emits.
They are captured on the *current* (pylammpsmpi) code by substituting a
``RecordingRunner`` for the live LAMMPS object: it records every ``command(str)``
but executes nothing.  After the ExecutableRunner refactor the same drivers must
emit the same logical command stream (Part 5.9), so these ``tests/golden/*.txt``
files are an immutable contract -- if one changes, the code is wrong.
"""
import os
import copy
import difflib
import shutil

import numpy as np
import pytest
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(HERE, "fixtures")
GOLDEN_DIR = os.path.join(HERE, "golden")
BASELINE_INPUTS = os.path.join(HERE, "baselines", "inputs")


def normalize_command(s, simfolder=None):
    """Normalize a command string for golden comparison.

    Collapses whitespace (so cosmetic padding never breaks diffs) and rewrites
    the absolute sim-folder path (which some drivers embed in ``read_data`` /
    ``write_data`` via ``os.path.join(self.simfolder, ...)``) to a bare relative
    name -- otherwise goldens would carry a machine- and run-specific tmp path.
    """
    s = " ".join(s.split())
    if simfolder:
        s = s.replace(simfolder + "/", "").replace(simfolder, ".")
    return s


class RecordingRunner:
    """Signature-compatible stand-in for the LAMMPS driver object.

    Records normalized command strings; every lifecycle method is a harmless
    no-op so the interactive driver code runs to completion without a binary.
    """

    def __init__(self):
        self.commands = []
        self.simfolder = None

    def command(self, s):
        self.commands.append(normalize_command(s, self.simfolder))

    def close(self):
        pass

    def clear(self):
        pass

    def write(self, path):
        # Old script-mode branches call this; keep it harmless.
        with open(path, "w") as fh:
            fh.write("\n".join(self.commands))


def assert_golden(commands, name):
    """Compare ``commands`` against ``tests/golden/<name>.txt``.

    On first run (file missing) the golden is written and the check passes -- the
    file is then committed.  Afterwards it is a read-only comparison with a
    readable unified diff on mismatch.
    """
    os.makedirs(GOLDEN_DIR, exist_ok=True)
    path = os.path.join(GOLDEN_DIR, name + ".txt")
    actual = "\n".join(commands) + "\n"
    if not os.path.exists(path):
        with open(path, "w") as fh:
            fh.write(actual)
        print(f"[golden] wrote new golden {name}.txt ({len(commands)} commands)")
        return
    with open(path) as fh:
        expected = fh.read()
    if actual != expected:
        diff = "\n".join(
            difflib.unified_diff(
                expected.splitlines(),
                actual.splitlines(),
                fromfile=f"golden/{name}.txt",
                tofile="recorded",
                lineterm="",
            )
        )
        raise AssertionError(f"golden stream mismatch for {name}:\n{diff}")


def _load_calc_data(scenario, overrides):
    with open(os.path.join(BASELINE_INPUTS, scenario + ".yaml")) as fh:
        data = yaml.safe_load(fh)
    calc_data = data["calculations"][0]
    # Portable, deterministic goldens: keep pair_coeff paths verbatim (no abspath).
    calc_data["fix_potential_path"] = False
    calc_data.update(copy.deepcopy(overrides))
    return calc_data


@pytest.fixture
def make_calc(tmp_path):
    """Factory: build a validated Calculation from a baseline input yaml.

    Structure generation during validation writes a lammps-data file into cwd;
    we contain that in a tmp dir and then pin ``calc.lattice`` to a fixed
    sentinel so the recorded ``read_data`` command is machine-independent.
    """
    from calphy.input import Calculation

    def _build(scenario, lattice_sentinel="structure.data", **overrides):
        calc_data = _load_calc_data(scenario, overrides)
        build_dir = tmp_path / ("build_" + scenario)
        build_dir.mkdir(exist_ok=True)
        prev = os.getcwd()
        os.chdir(build_dir)
        try:
            calc = Calculation(**calc_data)
        finally:
            os.chdir(prev)
        # natoms/composition are already captured on calc; the file path is only
        # ever emitted as a `read_data` string, so pin it to a stable sentinel.
        if lattice_sentinel is not None:
            calc.lattice = lattice_sentinel
        return calc

    return _build


class _Recorded:
    """Adapter exposing ``.commands`` (whitespace- and simfolder-normalized)
    over an ExecutableRunner's ``logical_commands``, so the golden tests read the
    same logical stream whether the runner records or (dry-run) drives LAMMPS."""

    def __init__(self, holder):
        self._holder = holder

    @property
    def commands(self):
        runner = self._holder["runner"]
        return [normalize_command(c, runner.directory) for c in runner.logical_commands]


@pytest.fixture
def recorded_job(tmp_path, monkeypatch):
    """Factory: stage fixtures, patch the runner + solid-fraction, build a job.

    Returns ``(job, recorder)``.  ``recorder.commands`` holds the logical command
    stream after the driver method is called.  create_object is monkeypatched to
    return an ``ExecutableRunner(dry_run=True)`` (Part 5.9): the driver's added
    ``sync()`` calls are exercised but no LAMMPS binary runs, and
    ``logical_commands`` excludes the runner-injected replay/restart lines.

    ``_make(..., runner="library")`` substitutes a ``LibraryRunner`` over a fake
    in-memory pylammpsmpi instead (Part 10), so the same goldens can assert both
    backends emit an identical logical command stream.
    """
    import sys
    import types

    import calphy.helpers as ph
    from calphy.runner import ExecutableRunner

    np.random.seed(42)
    holder = {}

    class _FakeLammpsLibrary:
        """Command-swallowing stand-in for pylammpsmpi.LammpsLibrary."""

        def __init__(self, cores=1, working_directory=".", cmdargs=None):
            self.working_directory = working_directory
            # inner concurrent handle; MLIAP activation (if the test env has
            # lammps.mliap importable) must be a harmless no-op here
            self.lmp = types.SimpleNamespace(
                activate_mliappy=lambda: None,
                activate_mliappy_kokkos=lambda: None,
            )

        def command(self, s):
            if s.startswith("log "):
                open(os.path.join(self.working_directory, s.split()[1]), "w").close()

        def close(self):
            pass

    def fake_create_object(calc, directory):
        if holder.get("mode") == "library":
            from calphy.library_runner import LibraryRunner

            lmp = LibraryRunner(cores=1, cmdargs="", directory=directory)
        else:
            lmp = ExecutableRunner(
                binary="lmp", mpi_command=None, cores=1, cmdargs="",
                directory=directory, dry_run=True,
            )
        ph.emit_init_commands(lmp, calc)
        holder["runner"] = lmp
        return lmp

    monkeypatch.setattr(ph, "create_object", fake_create_object)

    def _enable_fake_pylammpsmpi():
        mod = types.ModuleType("pylammpsmpi")
        mod.LammpsLibrary = _FakeLammpsLibrary
        monkeypatch.setitem(sys.modules, "pylammpsmpi", mod)

    # Controllable solid-fraction: an optional finite sequence (for the melt
    # loop), else a constant default.  Default keeps checks passing.
    state = {"seq": None, "default": None}

    def fake_find_solid_fraction(file):
        if state["seq"]:
            return state["seq"].pop(0)
        return state["default"]

    monkeypatch.setattr(ph, "find_solid_fraction", fake_find_solid_fraction)

    def _make(JobClass, calc, fixtures=(), solid_fraction_seq=None,
              default_solid_fraction=None, runner="executable"):
        holder["mode"] = runner
        if runner == "library":
            _enable_fake_pylammpsmpi()
        # Pin the master seed to the value the goldens were recorded with:
        # Phase.__init__ now does np.random.seed(md.seed), which must
        # reproduce the historical np.random.seed(42) stream exactly.
        calc.md.seed = 42
        simfolder = str(tmp_path)
        for fx in fixtures:
            shutil.copy(os.path.join(FIXTURE_DIR, fx), os.path.join(simfolder, fx))
        state["seq"] = list(solid_fraction_seq) if solid_fraction_seq else None
        state["default"] = (
            default_solid_fraction
            if default_solid_fraction is not None
            else calc._natoms
        )
        job = JobClass(calculation=calc, simfolder=simfolder)
        return job, _Recorded(holder)

    return _make
