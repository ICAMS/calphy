"""Unit tests for calphy.runner (Part 3).  No LAMMPS binary required."""
import os
import subprocess
import sys
import types

import numpy as np
import pytest

import calphy.runner as R
from calphy.errors import LammpsExecutionError, RunnerStateError

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(HERE, "fixtures")
REPO_ROOT = os.path.dirname(HERE)


def make_runner(directory, dry_run=True, cores=1, mpi_command=None,
                cmdargs="", timeout=None, binary="lmp"):
    return R.ExecutableRunner(
        binary=binary, mpi_command=mpi_command, cores=cores, cmdargs=cmdargs,
        directory=str(directory), dry_run=dry_run, timeout=timeout,
    )


def feed(runner, commands):
    for c in commands:
        runner.command(c)


def seg_lines(directory, k):
    with open(os.path.join(str(directory), "calphy.seg%d.lmp" % k)) as fh:
        return fh.read().splitlines()


# --------------------------------------------------------------------------- #
# Classification / command()
# --------------------------------------------------------------------------- #
ALL_TOKENS = sorted(R.KNOWN_TOKENS)


@pytest.mark.parametrize("token", ALL_TOKENS)
def test_every_vocabulary_token_classifies(tmp_path, token):
    run = make_runner(tmp_path)
    # a bare first token (plus a dummy id) must not raise on ingestion
    run.command("%s x y z" % token)
    assert run.logical_commands[-1].startswith(token)


def test_unknown_command_raises(tmp_path):
    run = make_runner(tmp_path)
    with pytest.raises(RunnerStateError, match="unknown command"):
        run.command("min_style cg")


def test_command_normalizes_whitespace(tmp_path):
    run = make_runner(tmp_path)
    run.command("fix              1   all    nve")
    assert run.logical_commands == ["fix 1 all nve"]


def test_empty_command_ignored(tmp_path):
    run = make_runner(tmp_path)
    run.command("   ")
    assert run.logical_commands == []


# --------------------------------------------------------------------------- #
# Segment 0 vs replay header
# --------------------------------------------------------------------------- #
BOOT = [
    "units metal", "boundary p p p", "atom_style atomic", "timestep 0.001",
    "box tilt large", "pair_style eam/alloy", "read_data conf.data",
    "pair_coeff * * pot Cu", "mass 1 63.5",
]


def test_segment0_verbatim_and_restart(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["fix 1 all nve", "run 100"])
    run.sync()
    lines = seg_lines(tmp_path, 0)
    assert "read_restart" not in "\n".join(lines)          # seg 0 has no restart read
    assert lines[:len(BOOT)] == BOOT                        # verbatim, in order
    assert lines[-1] == "write_restart calphy.restart"      # ends with restart write


def test_segment1_header_matches_spec_exactly(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + [
        "group g1 type 1", "compute c1 g1 msd com yes",
        "variable v1 equal 1.0", "fix 1 all nve",
        "fix 2 all ave/time 10 10 100 v_v1 file avg.dat",
        "thermo_style custom step pe", "thermo 10", "run 100",
    ])
    run.sync()                       # segment 0
    feed(run, ["run 200"])
    run.sync()                       # segment 1
    assert seg_lines(tmp_path, 1) == [
        "units metal",
        "atom_style atomic",
        "timestep 0.001",
        "read_restart calphy.restart",
        "pair_style eam/alloy",
        "pair_coeff * * pot Cu",
        "mass 1 63.5",
        "group g1 type 1",
        "compute c1 g1 msd com yes",
        "variable v1 equal 1.0",
        "fix 1 all nve",
        "fix 2 all ave/time 10 10 100 v_v1 file avg.dat.seg1",
        "thermo_style custom step pe",
        "thermo 10",
        "run 200",
        "write_restart calphy.restart",
    ]


def test_empty_buffer_sync_is_noop(tmp_path):
    run = make_runner(tmp_path)
    run.sync()
    assert not os.path.exists(os.path.join(str(tmp_path), "calphy.seg0.lmp"))


# --------------------------------------------------------------------------- #
# fix / unfix / fix_modify / redefinition
# --------------------------------------------------------------------------- #
def _replay_after(run, tmp_path, extra="run 1"):
    """sync once (seg0), issue `extra`, sync again, return seg1 lines."""
    run.sync()
    run.command(extra)
    run.sync()
    return seg_lines(tmp_path, 1)


def test_fix_appears_in_replay_header(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["fix 1 all nve"])
    run.sync()                                  # seg 0 (verbatim)
    run.command("run 1")
    run.sync()                                  # seg 1 replays fix 1
    assert "fix 1 all nve" in seg_lines(tmp_path, 1)


def test_unfix_removes_from_subsequent_replays(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["fix 1 all nve"])
    run.sync()                                  # seg 0
    # seg 1: header re-emits live fix 1, then the buffer unfixes it (correct
    # continuity -- read_restart does not restore fixes).
    run.command("unfix 1")
    run.command("run 1")
    run.sync()
    assert "fix 1 all nve" in seg_lines(tmp_path, 1)   # re-defined...
    assert "unfix 1" in seg_lines(tmp_path, 1)         # ...then removed
    assert "1" not in run.state.fixes                   # gone from state
    # seg 2: no trace of fix 1 anywhere
    run.command("run 1")
    run.sync()
    assert not any(l.startswith("fix 1 ") for l in seg_lines(tmp_path, 2))


def test_redefined_fix_keeps_position(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["fix a all nve", "fix b all nvt temp 1 1 1",
                      "fix a all langevin 1 1 1 12"])
    lines = _replay_after(run, tmp_path)
    fixes = [l for l in lines if l.startswith("fix ")]
    assert fixes == ["fix a all langevin 1 1 1 12", "fix b all nvt temp 1 1 1"]


def test_fix_modify_attaches_and_replays_after_fix(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["fix 1 all nve", "fix_modify 1 temp mytemp", "fix 2 all nvt temp 1 1 1"])
    lines = _replay_after(run, tmp_path)
    i = lines.index("fix 1 all nve")
    assert lines[i + 1] == "fix_modify 1 temp mytemp"


def test_fix_modify_unknown_id_raises(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["fix_modify 99 temp x"])
    with pytest.raises(RunnerStateError, match="unknown fix id"):
        run.sync()


# --------------------------------------------------------------------------- #
# variables / pair_block / dumps
# --------------------------------------------------------------------------- #
def test_variable_redefinition_keeps_position(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["variable x equal 1", "variable y equal 2", "variable x equal 3"])
    lines = _replay_after(run, tmp_path)
    vs = [l for l in lines if l.startswith("variable ")]
    assert vs == ["variable x equal 3", "variable y equal 2"]


def test_immediate_variable_crossing_boundary_raises(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["variable step0 equal $(step)"])
    run.sync()                          # seg 0 fine (verbatim)
    run.command("run 1")
    with pytest.raises(RunnerStateError, match="immediate-evaluation"):
        run.sync()                      # replay would capture a stale value


def test_new_pair_style_clears_old_block(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["pair_style meam", "pair_coeff * * lib SiC Si"])
    lines = _replay_after(run, tmp_path)
    pair = [l for l in lines if l.startswith(("pair_style", "pair_coeff", "mass"))]
    assert pair == ["pair_style meam", "pair_coeff * * lib SiC Si"]  # old eam/alloy gone


def test_overlay_pair_block_replays_completely(tmp_path):
    run = make_runner(tmp_path)
    feed(run, ["units metal", "atom_style atomic", "read_data conf.data",
               "pair_style hybrid/overlay eam/alloy eam/alloy",
               "pair_coeff * * eam/alloy 1 potA Cu",
               "pair_coeff * * eam/alloy 2 potB Cu", "mass 1 63.5"])
    lines = _replay_after(run, tmp_path)
    pair = [l for l in lines if l.startswith(("pair_style", "pair_coeff", "mass"))]
    assert pair == [
        "pair_style hybrid/overlay eam/alloy eam/alloy",
        "pair_coeff * * eam/alloy 1 potA Cu",
        "pair_coeff * * eam/alloy 2 potB Cu",
        "mass 1 63.5",
    ]


def test_live_dump_at_sync_raises(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["dump d1 all custom 1 t.dat id x", "run 1"])
    run.sync()                          # seg 0 verbatim (dump still live in state)
    run.command("run 1")
    with pytest.raises(RunnerStateError, match="dump"):
        run.sync()


def test_dump_undump_within_segment_is_fine(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["dump d1 all custom 1 t.dat id x", "run 0", "undump d1"])
    run.sync()
    run.command("run 1")
    run.sync()                          # no live dump -> OK
    assert not any("dump" in l for l in seg_lines(tmp_path, 1))


def test_uncompute_removes(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["compute c1 all pe", "uncompute c1"])
    lines = _replay_after(run, tmp_path)
    assert not any(l.startswith("compute ") for l in lines)


# --------------------------------------------------------------------------- #
# file redirection across segments
# --------------------------------------------------------------------------- #
def test_file_replay_rewritten_per_segment(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ["fix 2 all ave/time 10 10 100 v_x file avg.dat", "run 1"])
    run.sync()
    assert "fix 2 all ave/time 10 10 100 v_x file avg.dat" in seg_lines(tmp_path, 0)
    run.command("run 1"); run.sync()
    assert "fix 2 all ave/time 10 10 100 v_x file avg.dat.seg1" in seg_lines(tmp_path, 1)
    run.command("run 1"); run.sync()
    assert "fix 2 all ave/time 10 10 100 v_x file avg.dat.seg2" in seg_lines(tmp_path, 2)


def test_replayed_fix_print_raises(tmp_path):
    run = make_runner(tmp_path)
    feed(run, BOOT + ['fix p all print 1 "${x}" file out.dat', "run 1"])
    run.sync()
    run.command("run 1")
    with pytest.raises(RunnerStateError, match="fix print"):
        run.sync()


# --------------------------------------------------------------------------- #
# SessionState helpers direct
# --------------------------------------------------------------------------- #
def test_rewrite_replayed_fix_append_raises():
    with pytest.raises(RunnerStateError, match="append"):
        R._rewrite_replayed_fix("fix 2 all ave/time 1 1 1 v_x file a.dat append", 1)


def test_rewrite_replayed_fix_no_file_unchanged():
    cmd = "fix 1 all nvt temp 1 1 1"
    assert R._rewrite_replayed_fix(cmd, 3) == cmd


# --------------------------------------------------------------------------- #
# Execution layer (monkeypatched subprocess)
# --------------------------------------------------------------------------- #
class FakeRun:
    def __init__(self, returncode=0, log_content=None, stderr="", timeout_exc=False):
        self.returncode = returncode
        self.log_content = log_content
        self.stderr = stderr
        self.timeout_exc = timeout_exc
        self.calls = []

    def __call__(self, argv, cwd=None, capture_output=None, text=None, timeout=None):
        self.calls.append({"argv": list(argv), "cwd": cwd, "timeout": timeout})
        if self.timeout_exc:
            raise subprocess.TimeoutExpired(argv, timeout)
        if self.log_content is not None:
            log_name = argv[argv.index("-log") + 1]
            with open(os.path.join(cwd, log_name), "w") as fh:
                fh.write(self.log_content)
        return types.SimpleNamespace(returncode=self.returncode, stdout="", stderr=self.stderr)


def drive_one_segment(runner):
    feed(runner, BOOT + ["run 1"])
    runner.sync()


def test_argv_serial(tmp_path, monkeypatch):
    fake = FakeRun()
    monkeypatch.setattr(R.subprocess, "run", fake)
    run = make_runner(tmp_path, dry_run=False, cores=1)
    drive_one_segment(run)
    argv = fake.calls[0]["argv"]
    assert argv == ["lmp", "-in", "calphy.seg0.lmp", "-log", "calphy.seg0.log",
                    "-screen", "none"]
    assert fake.calls[0]["cwd"] == str(tmp_path)


def test_argv_mpi_and_cmdargs_str(tmp_path, monkeypatch):
    fake = FakeRun()
    monkeypatch.setattr(R.subprocess, "run", fake)
    run = make_runner(tmp_path, dry_run=False, cores=4, mpi_command="mpirun",
                      cmdargs="-sf gpu -pk gpu 1")
    drive_one_segment(run)
    argv = fake.calls[0]["argv"]
    assert argv[:3] == ["mpirun", "-np", "4"]
    assert argv[3] == "lmp" and "-screen" in argv
    assert argv[-5:] == ["-sf", "gpu", "-pk", "gpu", "1"]   # str cmdargs split on spaces


def test_cmdargs_list(tmp_path, monkeypatch):
    fake = FakeRun()
    monkeypatch.setattr(R.subprocess, "run", fake)
    run = make_runner(tmp_path, dry_run=False, cmdargs=["-sf", "omp"])
    drive_one_segment(run)
    assert fake.calls[0]["argv"][-2:] == ["-sf", "omp"]


def test_nonzero_rc_raises_with_excerpt(tmp_path, monkeypatch):
    fake = FakeRun(returncode=1, log_content="line1\nline2\nfatal boom\n")
    monkeypatch.setattr(R.subprocess, "run", fake)
    run = make_runner(tmp_path, dry_run=False)
    with pytest.raises(LammpsExecutionError) as ei:
        drive_one_segment(run)
    msg = str(ei.value)
    assert "return code 1" in msg and "fatal boom" in msg and "calphy.seg0.lmp" in msg


def test_error_in_log_with_rc0_raises(tmp_path, monkeypatch):
    fake = FakeRun(returncode=0, log_content="setup...\nERROR: bad pair_coeff\n")
    monkeypatch.setattr(R.subprocess, "run", fake)
    run = make_runner(tmp_path, dry_run=False)
    with pytest.raises(LammpsExecutionError, match="ERROR"):
        drive_one_segment(run)


def test_timeout_raises(tmp_path, monkeypatch):
    fake = FakeRun(timeout_exc=True)
    monkeypatch.setattr(R.subprocess, "run", fake)
    run = make_runner(tmp_path, dry_run=False, timeout=5)
    with pytest.raises(LammpsExecutionError, match="timed out"):
        drive_one_segment(run)


def test_failure_message_uses_stderr_when_log_missing(tmp_path, monkeypatch):
    fake = FakeRun(returncode=1, log_content=None, stderr="stderr detail here")
    monkeypatch.setattr(R.subprocess, "run", fake)
    run = make_runner(tmp_path, dry_run=False)
    with pytest.raises(LammpsExecutionError, match="stderr detail here"):
        drive_one_segment(run)


def test_dry_run_writes_files_but_does_not_execute(tmp_path, monkeypatch):
    fake = FakeRun()
    monkeypatch.setattr(R.subprocess, "run", fake)
    run = make_runner(tmp_path, dry_run=True)
    drive_one_segment(run)
    assert fake.calls == []
    assert os.path.exists(os.path.join(str(tmp_path), "calphy.seg0.lmp"))


# --------------------------------------------------------------------------- #
# close() and rotate_logs()
# --------------------------------------------------------------------------- #
def test_close_flushes_and_removes_restart(tmp_path, monkeypatch):
    monkeypatch.setattr(R.subprocess, "run", FakeRun())
    run = make_runner(tmp_path, dry_run=False)
    feed(run, BOOT + ["run 1"])
    # simulate LAMMPS leaving a restart file behind
    open(os.path.join(str(tmp_path), "calphy.restart"), "w").close()
    run.close()
    assert os.path.exists(os.path.join(str(tmp_path), "calphy.seg0.lmp"))
    assert not os.path.exists(os.path.join(str(tmp_path), "calphy.restart"))


def test_rotate_logs_concatenates_segment_logs(tmp_path, monkeypatch):
    monkeypatch.setattr(R.subprocess, "run", FakeRun(log_content="seg log body\n"))
    run = make_runner(tmp_path, dry_run=False)
    feed(run, BOOT + ["run 1"]); run.sync()
    run.command("run 1"); run.sync()
    run.rotate_logs("averaging")
    out = os.path.join(str(tmp_path), "averaging.log.lammps")
    text = open(out).read()
    assert "# --- segment 0: calphy.seg0.lmp ---" in text
    assert "# --- segment 1: calphy.seg1.lmp ---" in text
    assert text.count("seg log body") == 2
    # list reset: a second rotation writes an empty file
    run.rotate_logs("integration")
    assert open(os.path.join(str(tmp_path), "integration.log.lammps")).read() == ""


# --------------------------------------------------------------------------- #
# Binary resolution
# --------------------------------------------------------------------------- #
def test_resolution_explicit_wins(tmp_path, monkeypatch):
    monkeypatch.setattr(R.shutil, "which", lambda v: "/resolved/" + v)
    monkeypatch.setenv("CALPHY_LAMMPS_EXECUTABLE", "envlmp")
    assert R.resolve_lammps_executable("explicitlmp") == "/resolved/explicitlmp"


def test_resolution_env_second(tmp_path, monkeypatch):
    monkeypatch.setattr(R.shutil, "which", lambda v: "/resolved/" + v if v else None)
    monkeypatch.setenv("CALPHY_LAMMPS_EXECUTABLE", "envlmp")
    assert R.resolve_lammps_executable(None) == "/resolved/envlmp"


def test_resolution_path_third(monkeypatch):
    monkeypatch.delenv("CALPHY_LAMMPS_EXECUTABLE", raising=False)
    monkeypatch.setattr(R.shutil, "which", lambda v: "/usr/bin/lmp" if v == "lmp" else None)
    assert R.resolve_lammps_executable(None) == "/usr/bin/lmp"


def test_resolution_absolute_path_executable(tmp_path, monkeypatch):
    fake = tmp_path / "mylmp"
    fake.write_text("#!/bin/sh\n")
    fake.chmod(0o755)
    monkeypatch.setattr(R.shutil, "which", lambda v: None)  # not on PATH
    assert R.resolve_lammps_executable(str(fake)) == str(fake)


def test_resolution_failure_names_all_steps(monkeypatch):
    monkeypatch.delenv("CALPHY_LAMMPS_EXECUTABLE", raising=False)
    monkeypatch.setattr(R.shutil, "which", lambda v: None)
    with pytest.raises(ValueError) as ei:
        R.resolve_lammps_executable(None)
    msg = str(ei.value)
    assert "input key" in msg and "CALPHY_LAMMPS_EXECUTABLE" in msg and "PATH" in msg


def test_resolve_mpi_uses_mpirun(monkeypatch):
    monkeypatch.delenv("CALPHY_MPI_EXECUTABLE", raising=False)
    monkeypatch.setattr(R.shutil, "which", lambda v: "/usr/bin/mpirun" if v == "mpirun" else None)
    assert R.resolve_mpi_executable(None) == "/usr/bin/mpirun"


# --------------------------------------------------------------------------- #
# read_timeseries
# --------------------------------------------------------------------------- #
def test_read_timeseries_single_file_matches_loadtxt():
    got = R.read_timeseries(FIXTURE_DIR, "avg.dat", usecols=(1, 2, 3, 4))
    exp = np.loadtxt(os.path.join(FIXTURE_DIR, "avg.dat"), usecols=(1, 2, 3, 4))
    assert np.array_equal(got, exp)


def test_read_timeseries_single_column():
    got = R.read_timeseries(FIXTURE_DIR, "msd.dat", usecols=(1,))
    exp = np.loadtxt(os.path.join(FIXTURE_DIR, "msd.dat"), usecols=(1,))
    assert np.array_equal(got, exp)


def test_read_timeseries_multifile_concat(tmp_path):
    base = np.loadtxt(os.path.join(FIXTURE_DIR, "avg.dat"), usecols=(1, 2, 3, 4))
    # stage avg.dat + two synthetic segment files with comment headers
    import shutil as _sh
    _sh.copy(os.path.join(FIXTURE_DIR, "avg.dat"), tmp_path / "avg.dat")
    seg1 = base[:2] + 100.0
    seg2 = base[:1] + 200.0
    for name, arr in [("avg.dat.seg1", seg1), ("avg.dat.seg2", seg2)]:
        with open(tmp_path / name, "w") as fh:
            fh.write("# TimeStep a b c d\n")
            for row in np.atleast_2d(arr):
                fh.write("0 " + " ".join("%.10g" % v for v in row) + "\n")
    got = R.read_timeseries(str(tmp_path), "avg.dat", usecols=(1, 2, 3, 4))
    assert got.shape[0] == base.shape[0] + 2 + 1
    assert np.allclose(got[:base.shape[0]], base)
    assert np.allclose(got[base.shape[0]:base.shape[0] + 2], seg1)


def test_read_timeseries_segment_ordering(tmp_path):
    # ensure .seg10 sorts after .seg2 (numeric, not lexicographic)
    with open(tmp_path / "x.dat", "w") as fh:
        fh.write("0 0\n")
    for k, v in [(2, 2.0), (10, 10.0)]:
        with open(tmp_path / ("x.dat.seg%d" % k), "w") as fh:
            fh.write("0 %g\n" % v)
    got = R.read_timeseries(str(tmp_path), "x.dat", usecols=(1,))
    assert list(got) == [0.0, 2.0, 10.0]


def test_read_timeseries_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        R.read_timeseries(str(tmp_path), "nope.dat", usecols=(1,))


def test_read_timeseries_all_columns_multifile_mixed_rows(tmp_path):
    # regression: usecols=None with a multi-row base and a single-row seg file
    # (np.loadtxt squeezes the single row to 1D) must still concatenate.
    with open(tmp_path / "m.dat", "w") as fh:
        fh.write("# TimeStep v\n0 1.0\n100 2.0\n200 3.0\n")
    with open(tmp_path / "m.dat.seg1", "w") as fh:
        fh.write("# TimeStep v\n300 4.0\n")            # single row -> 1D from loadtxt
    got = R.read_timeseries(str(tmp_path), "m.dat")     # usecols=None -> all columns
    assert got.shape == (4, 2)
    assert list(got[:, 0]) == [0.0, 100.0, 200.0, 300.0]
    assert list(got[:, 1]) == [1.0, 2.0, 3.0, 4.0]


# --------------------------------------------------------------------------- #
# Import without pylammpsmpi
# --------------------------------------------------------------------------- #
def test_imports_without_pylammpsmpi():
    code = (
        "import sys, importlib.abc\n"
        "class B(importlib.abc.MetaPathFinder):\n"
        "    def find_spec(self, name, path, target=None):\n"
        "        if name == 'pylammpsmpi' or name.startswith('pylammpsmpi.'):\n"
        "            raise ModuleNotFoundError('blocked ' + name)\n"
        "        return None\n"
        "sys.meta_path.insert(0, B())\n"
        "import calphy.runner\n"
        "assert 'pylammpsmpi' not in sys.modules, 'pylammpsmpi leaked in'\n"
        "assert calphy.runner.KNOWN_TOKENS\n"
        "print('OK')\n"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True,
                          text=True, cwd=REPO_ROOT)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "OK" in proc.stdout
