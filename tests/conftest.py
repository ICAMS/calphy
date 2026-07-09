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


@pytest.fixture
def recorded_job(tmp_path, monkeypatch):
    """Factory: stage fixtures, patch the runner + solid-fraction, build a job.

    Returns ``(job, recorder)``.  ``recorder.commands`` holds the stream after
    the driver method is called.
    """
    import calphy.helpers as ph

    np.random.seed(42)
    recorder = RecordingRunner()
    original_create_object = ph.create_object

    def fake_create_object(*args, **kwargs):
        # Reuse the real init-command emission (units/boundary/atom_style/
        # timestep/box + init_commands merge) but onto the recorder, and force
        # the library (non-script) path.
        kwargs["lmp"] = recorder
        kwargs["script_mode"] = False
        return original_create_object(*args, **kwargs)

    monkeypatch.setattr(ph, "create_object", fake_create_object)

    # Controllable solid-fraction: an optional finite sequence (for the melt
    # loop), else a constant default.  Default keeps checks passing.
    state = {"seq": None, "default": None}

    def fake_find_solid_fraction(file):
        if state["seq"]:
            return state["seq"].pop(0)
        return state["default"]

    monkeypatch.setattr(ph, "find_solid_fraction", fake_find_solid_fraction)

    def _make(JobClass, calc, fixtures=(), solid_fraction_seq=None,
              default_solid_fraction=None):
        simfolder = str(tmp_path)
        recorder.simfolder = simfolder
        for fx in fixtures:
            shutil.copy(os.path.join(FIXTURE_DIR, fx), os.path.join(simfolder, fx))
        state["seq"] = list(solid_fraction_seq) if solid_fraction_seq else None
        state["default"] = (
            default_solid_fraction
            if default_solid_fraction is not None
            else calc._natoms
        )
        job = JobClass(calculation=calc, simfolder=simfolder)
        return job, recorder

    return _make
