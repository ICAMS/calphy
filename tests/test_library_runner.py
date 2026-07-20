"""Unit tests for calphy.library_runner (Part 9).  No pylammpsmpi required:
a fake module is injected so the lazy import inside LibraryRunner resolves."""
import os
import sys
import types

import numpy as np
import pytest

from calphy.errors import RunnerStateError

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(HERE, "fixtures")


class _FakeConcurrent:
    """The inner pylammpsmpi handle (LammpsLibrary.lmp); records MLIAP activation."""

    def __init__(self):
        self.activated = []

    def activate_mliappy(self):
        self.activated.append("mliappy")

    def activate_mliappy_kokkos(self):
        self.activated.append("mliappy_kokkos")


class FakeLammpsLibrary:
    """Mimics the pylammpsmpi.LammpsLibrary surface LibraryRunner touches."""

    def __init__(self, cores=1, working_directory=".", cmdargs=None):
        self.cores = cores
        self.working_directory = working_directory
        self.cmdargs = list(cmdargs or [])
        self.commands = []
        self.closed = False
        self.lmp = _FakeConcurrent()

    def command(self, s):
        self.commands.append(s)
        # LAMMPS opens the new log file as soon as the log command runs.
        if s.startswith("log "):
            open(os.path.join(self.working_directory, s.split()[1]), "w").close()

    def close(self):
        self.closed = True


@pytest.fixture
def fake_pylammpsmpi(monkeypatch):
    mod = types.ModuleType("pylammpsmpi")
    mod.LammpsLibrary = FakeLammpsLibrary
    monkeypatch.setitem(sys.modules, "pylammpsmpi", mod)
    # deterministic MLIAP state: absent unless a test injects the fake below
    monkeypatch.setitem(sys.modules, "lammps", None)
    monkeypatch.setitem(sys.modules, "lammps.mliap", None)
    return mod


@pytest.fixture
def fake_lammps_mliap(monkeypatch):
    """Make `import lammps.mliap` succeed (as when LAMMPS ships ML-IAP python)."""
    mliap = types.ModuleType("lammps.mliap")
    lammps_mod = types.ModuleType("lammps")
    lammps_mod.mliap = mliap
    monkeypatch.setitem(sys.modules, "lammps", lammps_mod)
    monkeypatch.setitem(sys.modules, "lammps.mliap", mliap)
    return mliap


def make_runner(directory, cores=1, cmdargs=""):
    from calphy.library_runner import LibraryRunner

    return LibraryRunner(cores=cores, cmdargs=cmdargs, directory=str(directory))


# --------------------------------------------------------------------------- #
# Construction
# --------------------------------------------------------------------------- #
def test_constructor_wires_library(tmp_path, fake_pylammpsmpi):
    run = make_runner(tmp_path, cores=4, cmdargs="-sf omp")
    assert run.lmp.cores == 4
    assert run.lmp.working_directory == str(tmp_path)
    assert run.lmp.cmdargs[:2] == ["-sf", "omp"]
    assert run.lmp.cmdargs[-4:] == ["-screen", "none", "-log", "calphy.live.0.log"]


def test_constructor_respects_existing_screen_and_log(tmp_path, fake_pylammpsmpi):
    run = make_runner(tmp_path, cmdargs=["-screen", "out", "-log", "my.log"])
    assert run.lmp.cmdargs == ["-screen", "out", "-log", "my.log"]


def test_missing_pylammpsmpi_raises_helpful_error(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "pylammpsmpi", None)  # import -> ImportError
    with pytest.raises(ImportError, match=r"calphy\[library\]"):
        make_runner(tmp_path)


# --------------------------------------------------------------------------- #
# Command dispatch
# --------------------------------------------------------------------------- #
def test_commands_forwarded_immediately_and_recorded(tmp_path, fake_pylammpsmpi):
    run = make_runner(tmp_path)
    run.command("fix    1  all   nve")
    assert run.lmp.commands == ["fix 1 all nve"]      # normalized, no buffering
    assert run.logical_commands == ["fix 1 all nve"]


def test_unknown_command_raises_and_is_not_forwarded(tmp_path, fake_pylammpsmpi):
    run = make_runner(tmp_path)
    with pytest.raises(RunnerStateError, match="unknown command"):
        run.command("min_style cg")
    assert run.lmp.commands == []


def test_sync_is_noop(tmp_path, fake_pylammpsmpi):
    run = make_runner(tmp_path)
    run.command("run 100")
    run.sync()
    assert run.lmp.commands == ["run 100"]            # no restart/replay traffic


def test_close_closes_library(tmp_path, fake_pylammpsmpi):
    run = make_runner(tmp_path)
    run.close()
    assert run.lmp.closed


# --------------------------------------------------------------------------- #
# MLIAP activation (upstream ICAMS/calphy#260)
# --------------------------------------------------------------------------- #
def test_no_mliap_activation_without_module(tmp_path, fake_pylammpsmpi):
    run = make_runner(tmp_path)
    assert run.lmp.lmp.activated == []


def test_mliap_activated_when_available(tmp_path, fake_pylammpsmpi, fake_lammps_mliap):
    run = make_runner(tmp_path)
    assert run.lmp.lmp.activated == ["mliappy"]


def test_mliap_kokkos_variant_with_k_flag(tmp_path, fake_pylammpsmpi, fake_lammps_mliap):
    run = make_runner(tmp_path, cmdargs="-k on g 1 -sf kk")
    assert run.lmp.lmp.activated == ["mliappy_kokkos"]


def test_mliap_skipped_on_old_pylammpsmpi(
    tmp_path, fake_pylammpsmpi, fake_lammps_mliap, monkeypatch, caplog
):
    # pylammpsmpi <= 0.4.1 has no activation proxies: warn and skip, never crash
    monkeypatch.delattr(_FakeConcurrent, "activate_mliappy")
    monkeypatch.delattr(_FakeConcurrent, "activate_mliappy_kokkos")
    run = make_runner(tmp_path)
    assert run.lmp.lmp.activated == []
    assert any("pylammpsmpi" in r.message for r in caplog.records)


# --------------------------------------------------------------------------- #
# rotate_logs
# --------------------------------------------------------------------------- #
def test_rotate_logs_renames_scratch_log(tmp_path, fake_pylammpsmpi):
    run = make_runner(tmp_path)
    with open(tmp_path / "calphy.live.0.log", "w") as fh:
        fh.write("averaging output\n")
    run.rotate_logs("averaging")
    assert run.lmp.commands == ["log calphy.live.1.log"]
    assert (tmp_path / "averaging.log.lammps").read_text() == "averaging output\n"
    assert not (tmp_path / "calphy.live.0.log").exists()
    # second rotation collects the next scratch (created by the log command)
    run.rotate_logs("integration")
    assert (tmp_path / "integration.log.lammps").exists()


def test_rotate_logs_without_scratch_writes_empty_file(tmp_path, fake_pylammpsmpi):
    run = make_runner(tmp_path)
    run.rotate_logs("averaging")   # fake created live.1 but live.0 never existed
    assert (tmp_path / "averaging.log.lammps").read_text() == ""


def test_rotate_logs_after_close_renames_without_commanding(tmp_path, fake_pylammpsmpi):
    # the drivers always close the session before rotating; the closed session
    # must not be sent a log-switch command, but the scratch log still rotates.
    run = make_runner(tmp_path)
    with open(tmp_path / "calphy.live.0.log", "w") as fh:
        fh.write("stage output\n")
    run.close()
    run.rotate_logs("averaging")
    assert run.lmp.commands == []
    assert (tmp_path / "averaging.log.lammps").read_text() == "stage output\n"


# --------------------------------------------------------------------------- #
# read_timeseries (inherited file-based accessor; no .seg files in library mode)
# --------------------------------------------------------------------------- #
def test_read_timeseries_single_file(tmp_path, fake_pylammpsmpi):
    import shutil

    shutil.copy(os.path.join(FIXTURE_DIR, "avg.dat"), tmp_path / "avg.dat")
    run = make_runner(tmp_path)
    got = run.read_timeseries("avg.dat", usecols=(1, 2, 3, 4))
    exp = np.loadtxt(os.path.join(FIXTURE_DIR, "avg.dat"), usecols=(1, 2, 3, 4))
    assert np.array_equal(got, exp)


# --------------------------------------------------------------------------- #
# Backend equivalence: full driver paths must emit the exact command stream the
# executable backend froze in tests/golden/ -- the "unified backends" contract.
# --------------------------------------------------------------------------- #
def test_solid_fe_averaging_matches_executable_golden(make_calc, recorded_job):
    from calphy.solid import Solid
    from conftest import assert_golden

    calc = make_calc(
        "B1", tolerance={"pressure": 1e12, "spring_constant": 1e12}
    )
    job, rec = recorded_job(
        Solid, calc, fixtures=["avg.dat", "msd.dat"], runner="library"
    )
    job.run_averaging()
    assert_golden(rec.commands, "solid_fe_averaging")


def test_liquid_fe_averaging_meltcycle_matches_executable_golden(
    make_calc, recorded_job
):
    from calphy.liquid import Liquid
    from conftest import assert_golden

    calc = make_calc(
        "B4", tolerance={"pressure": 1e12, "spring_constant": 1e12}
    )
    n = calc._natoms
    job, rec = recorded_job(
        Liquid, calc, fixtures=["avg.dat"], solid_fraction_seq=[n, n, 0],
        runner="library",
    )
    job.run_averaging()
    assert_golden(rec.commands, "liquid_fe_averaging_meltcycle")


# --------------------------------------------------------------------------- #
# create_object backend selection
# --------------------------------------------------------------------------- #
def test_create_object_library_mode(tmp_path, fake_pylammpsmpi, make_calc):
    import calphy.helpers as ph
    from calphy.library_runner import LibraryRunner

    calc = make_calc("B1", execution_mode="library")
    lmp = ph.create_object(calc, str(tmp_path))
    assert isinstance(lmp, LibraryRunner)
    # primed with the usual four init commands
    assert lmp.logical_commands[0] == "units metal"
    assert len(lmp.logical_commands) == 4


def test_execution_mode_validation(make_calc):
    with pytest.raises(Exception, match="execution_mode"):
        make_calc("B1", execution_mode="banana")


# --------------------------------------------------------------------------- #
# Import hygiene: calphy must import without pylammpsmpi installed
# --------------------------------------------------------------------------- #
def test_module_imports_without_pylammpsmpi():
    import subprocess

    code = (
        "import sys, importlib.abc\n"
        "class B(importlib.abc.MetaPathFinder):\n"
        "    def find_spec(self, name, path, target=None):\n"
        "        if name == 'pylammpsmpi' or name.startswith('pylammpsmpi.'):\n"
        "            raise ModuleNotFoundError('blocked ' + name)\n"
        "        return None\n"
        "sys.meta_path.insert(0, B())\n"
        "import calphy.library_runner\n"
        "assert 'pylammpsmpi' not in sys.modules, 'pylammpsmpi leaked in'\n"
        "print('OK')\n"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True,
        cwd=os.path.dirname(HERE),
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "OK" in proc.stdout
