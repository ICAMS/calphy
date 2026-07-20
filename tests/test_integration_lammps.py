"""Integration tests that run the real LAMMPS binary (Part 6).

These exercise the ExecutableRunner end to end: restart continuity across
segment boundaries, cumulative file redirection, and sticky-state replay.  They
are marked ``lammps`` and skipped when no ``lmp`` binary resolves.
"""
import os
import re
import glob

import numpy as np
import pytest

from calphy.runner import (
    ExecutableRunner,
    resolve_lammps_executable,
    resolve_mpi_executable,
    read_timeseries,
)
from calphy.errors import LammpsExecutionError

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(HERE, "fixtures")
POT = os.path.join(HERE, "Cu01.eam.alloy")
CONF = os.path.join(FIXTURE_DIR, "conf.equilibration.data")   # equilibrated Cu fcc, 500 atoms

try:
    BINARY = resolve_lammps_executable(None)
except ValueError:
    BINARY = None

pytestmark = [
    pytest.mark.lammps,
    pytest.mark.skipif(BINARY is None, reason="no lmp binary resolves"),
]


def make_runner(directory, cores=1, mpi_command=None):
    return ExecutableRunner(
        binary=BINARY, mpi_command=mpi_command, cores=cores, cmdargs="",
        directory=str(directory), dry_run=False, timeout=300,
    )


def boot(run, temp=500.0):
    """Init header + Cu structure + eam/alloy potential (calphy's ordering)."""
    for c in [
        "units metal", "boundary p p p", "atom_style atomic", "timestep 0.001",
        "box tilt large",
        "pair_style eam/alloy",
        "read_data %s" % CONF,
        "pair_coeff * * %s Cu" % POT,
        "mass 1 63.546",
    ]:
        run.command(c)


def log_values(directory, tag):
    """Return the floats printed via `print "<tag> $(...)"` across all seg logs."""
    vals = []
    for logpath in seg_logs(directory):
        with open(logpath) as fh:
            for line in fh:
                if line.startswith(tag + " "):
                    vals.append([float(x) for x in line.split()[1:]])
    return vals


def seg_logs(directory):
    return sorted(
        glob.glob(os.path.join(str(directory), "calphy.seg*.log")),
        key=lambda p: int(re.search(r"seg(\d+)", p).group(1)),
    )


def parse_thermo(logpath):
    """Parse LAMMPS thermo tables in a log into {column: np.array}."""
    cols = None
    rows = []
    with open(logpath) as fh:
        for line in fh:
            s = line.split()
            if not s:
                continue
            if s[0] == "Step":
                cols = s
                continue
            if cols is not None:
                try:
                    vals = [float(x) for x in s]
                except ValueError:
                    cols = None
                    continue
                if len(vals) == len(cols):
                    rows.append(vals)
                else:
                    cols = None
    arr = np.array(rows)
    return {c: arr[:, i] for i, c in enumerate(cols_all(logpath))} if arr.size else {}


def cols_all(logpath):
    with open(logpath) as fh:
        for line in fh:
            s = line.split()
            if s and s[0] == "Step":
                return s
    return []


# --------------------------------------------------------------------------- #
# T1 — exact-restart energy continuity (the strongest single test)
# --------------------------------------------------------------------------- #
def test_t1_restart_energy_continuity(tmp_path):
    run = make_runner(tmp_path)
    boot(run)
    run.command("fix 1 all nvt temp 500 500 0.1")
    run.command("velocity all create 500 12345 mom yes rot yes")
    run.command("run 100")
    run.command('print "PECHK $(pe)"')
    run.sync()                          # segment 0: run + write_restart
    run.command("run 0")
    run.command('print "PECHK $(pe)"')
    run.sync()                          # segment 1: read_restart + replay + run 0
    run.close()

    pe = [v[0] for v in log_values(tmp_path, "PECHK")]
    assert len(pe) == 2, "expected one PE per segment, got %r" % pe
    end_k, start_k1 = pe
    rel = abs(start_k1 - end_k) / max(abs(end_k), 1e-12)
    assert rel <= 1e-6, (
        "PE discontinuity across the restart boundary: "
        "end seg0=%.10g, start seg1=%.10g, rel=%.2e" % (end_k, start_k1, rel)
    )


# --------------------------------------------------------------------------- #
# T2 — Nose-Hoover state restoration (no discontinuity at the boundary)
# --------------------------------------------------------------------------- #
def test_t2_nose_hoover_restoration(tmp_path):
    run = make_runner(tmp_path)
    boot(run)
    run.command("thermo_style custom step temp press vol")
    run.command("thermo 1")
    run.command("fix 1 all npt temp 500 500 0.1 iso 0 0 0.1")
    run.command("velocity all create 500 12345 mom yes rot yes")
    run.command("run 1000")
    run.sync()                          # boundary at step 1000
    run.command("run 1000")
    run.sync()
    run.close()

    logs = seg_logs(tmp_path)
    assert len(logs) == 2
    t0 = parse_thermo(logs[0])
    t1 = parse_thermo(logs[1])
    for key in ("Volume", "Press", "Step"):
        assert key in t0 and key in t1, "missing thermo column %r" % key

    # exact restart: the restored state (first row of seg1) matches seg0's last row
    assert abs(t1["Volume"][0] - t0["Volume"][-1]) <= 1e-6 * abs(t0["Volume"][-1])
    assert abs(t1["Step"][0] - t0["Step"][-1]) < 1e-9

    # the first dynamics step after the restart must not kick volume/pressure
    # beyond a few times the normal step-to-step fluctuation of the continuous run
    for key, factor in (("Volume", 3.0), ("Press", 4.0)):
        normal = np.std(np.diff(t0[key]))
        jump = abs(t1[key][1] - t1[key][0])
        assert jump < factor * normal + 1e-9, (
            "%s jump %.4g at the boundary exceeds %g x normal fluctuation %.4g"
            % (key, jump, factor, normal)
        )


# --------------------------------------------------------------------------- #
# T3 — ave/time file redirection across a segment boundary
# --------------------------------------------------------------------------- #
def test_t3_avetime_redirection(tmp_path):
    run = make_runner(tmp_path)
    boot(run)
    run.command("fix 1 all nvt temp 500 500 0.1")
    run.command("velocity all create 500 12345 mom yes rot yes")
    run.command("variable mpress equal press")
    run.command("fix 2 all ave/time 10 10 100 v_mpress file avg.dat")
    run.command("run 200")
    run.sync()                          # seg0 -> avg.dat
    run.command("run 200")
    run.sync()                          # seg1 -> avg.dat.seg1 (replayed fix rewritten)
    run.close()

    assert os.path.exists(os.path.join(str(tmp_path), "avg.dat"))
    assert glob.glob(os.path.join(str(tmp_path), "avg.dat.seg*")), "no .seg redirection file"

    steps = read_timeseries(str(tmp_path), "avg.dat", usecols=(0,))
    assert steps.ndim == 1 and len(steps) >= 3
    assert np.all(np.diff(steps) > 0), "step column not strictly increasing: %r" % steps
    press = read_timeseries(str(tmp_path), "avg.dat", usecols=(1,))
    assert np.all(np.isfinite(press))


# --------------------------------------------------------------------------- #
# T4 — group / compute / variable replay (MSD) across a boundary
# --------------------------------------------------------------------------- #
def test_t4_compute_replay_msd(tmp_path):
    run = make_runner(tmp_path)
    boot(run)
    run.command("fix 1 all nvt temp 500 500 0.1")
    run.command("velocity all create 500 12345 mom yes rot yes")
    run.command("group g1 type 1")
    run.command("compute c1 g1 msd com yes")
    run.command("variable msd1 equal c_c1[4]")
    run.command("fix 4 all ave/time 10 10 100 v_msd1 file msd.dat")
    run.command("run 100")
    run.sync()
    run.command("run 200")
    run.sync()
    run.close()

    msd = read_timeseries(str(tmp_path), "msd.dat")
    assert msd.ndim == 2 and msd.shape[1] >= 2
    col = msd[:, 1]
    assert np.all(np.isfinite(col)), "MSD column has non-finite values"
    assert np.any(col > 0), "MSD never became non-zero -- compute/group replay failed"


# --------------------------------------------------------------------------- #
# 6.2 — End-to-end matrix: re-run the Part 0 baselines through the new code
# --------------------------------------------------------------------------- #
import subprocess
import sys
import shutil
import yaml

REPO = os.path.dirname(HERE)
BASELINE_INPUTS = os.path.join(HERE, "baselines", "inputs")
BASELINE_DIR = os.path.join(HERE, "baselines")


def _abs_pair_coeff(pc):
    def fix(one):
        t = one.split()
        if len(t) > 2:
            cand = t[2] if os.path.isabs(t[2]) else os.path.join(REPO, t[2])
            if os.path.exists(cand):
                t[2] = os.path.abspath(cand)
        return " ".join(t)
    return [fix(x) for x in pc] if isinstance(pc, list) else fix(pc)


def run_scenario(scenario, rundir):
    """Run `calphy_kernel -i input.yaml -k 0` for a baseline scenario; return report dict."""
    with open(os.path.join(BASELINE_INPUTS, scenario + ".yaml")) as fh:
        data = yaml.safe_load(fh)
    calc = data["calculations"][0]
    if "pair_coeff" in calc:
        calc["pair_coeff"] = _abs_pair_coeff(calc["pair_coeff"])
    os.makedirs(rundir, exist_ok=True)
    with open(os.path.join(rundir, "input.yaml"), "w") as fh:
        yaml.safe_dump(data, fh, sort_keys=False)

    kernel_exe = os.path.join(os.path.dirname(sys.executable), "calphy_kernel")
    proc = subprocess.run(
        [kernel_exe, "-i", "input.yaml", "-k", "0"],
        cwd=rundir, capture_output=True, text=True, timeout=900,
    )
    assert proc.returncode == 0, (
        "calphy_kernel failed for %s:\n%s\n%s" % (scenario, proc.stdout[-2000:], proc.stderr[-3000:])
    )
    sims = [d for d in glob.glob(os.path.join(rundir, "*")) if os.path.isdir(d)]
    assert sims, "no sim folder produced for %s" % scenario
    with open(os.path.join(sims[0], "report.yaml")) as fh:
        return yaml.safe_load(fh)


def baseline_report(scenario):
    with open(os.path.join(BASELINE_DIR, scenario, "report.yaml")) as fh:
        return yaml.safe_load(fh)


def _check_free_energy(scenario, rundir):
    new = run_scenario(scenario, rundir)
    base = baseline_report(scenario)
    fe_new = new["results"]["free_energy"]
    fe_base = base["results"]["free_energy"]
    err_new = float(new["results"].get("error", 0.0))
    err_base = float(base["results"].get("error", 0.0))
    assert np.isfinite(fe_new), "%s free energy is not finite" % scenario
    # Floor of 5e-3 eV/atom (0.15%): the baselines use deliberately tiny
    # settings (2500 equil / 5000 switching) whose run-to-run fe scatter is
    # ~3-5e-3 -- larger than the reported switching-hysteresis error bar, and
    # inherent to the settings (not the refactor; verified by re-runs). A real
    # regression would be far larger than this floor.
    tol = max(3 * (abs(err_new) + abs(err_base)), 5e-3)
    delta = abs(fe_new - fe_base)
    assert delta < tol, (
        "%s free energy delta %.2e eV/atom exceeds tol %.2e "
        "(new=%.6f base=%.6f)" % (scenario, delta, tol, fe_new, fe_base)
    )
    # Spring constant is informational only: calphy rounds the mean MSD to 2
    # decimals and does not pin the RNG seed, so k = 3 kB T / MSD is bistable
    # run-to-run (e.g. B1 flips between ~1.85 and ~2.15) in the OLD code too.
    # The free energy is insensitive to the exact k (Frenkel-Ladd), so we assert
    # on fe only and just report the spring delta.
    k_new = new.get("average", {}).get("spring_constant")
    k_base = base.get("average", {}).get("spring_constant")
    kinfo = ""
    if k_new is not None and k_base is not None:
        kinfo = " k_new=%.4f k_base=%.4f" % (float(k_new), float(k_base))
    print("[%s] fe_new=%.6f fe_base=%.6f |d|=%.2e tol=%.2e%s"
          % (scenario, fe_new, fe_base, delta, tol, kinfo))


@pytest.mark.parametrize("scenario", ["B1", "B4", "B5"])
def test_e2e_baseline_fast(scenario, tmp_path):
    _check_free_energy(scenario, str(tmp_path))


@pytest.mark.slow
@pytest.mark.parametrize("scenario", ["B2", "B3", "B6", "B7", "B8", "B9", "B10"])
def test_e2e_baseline_slow(scenario, tmp_path):
    _check_free_energy(scenario, str(tmp_path))


# --------------------------------------------------------------------------- #
# 6.3 — Error paths
# --------------------------------------------------------------------------- #
def test_error_broken_pair_coeff(tmp_path):
    """A bad potential path makes LAMMPS ERROR -> LammpsExecutionError w/ excerpt."""
    run = make_runner(tmp_path)
    for c in [
        "units metal", "boundary p p p", "atom_style atomic", "timestep 0.001",
        "box tilt large", "pair_style eam/alloy",
        "read_data %s" % CONF,
        "pair_coeff * * /nonexistent/does_not_exist.eam.alloy Cu",
        "mass 1 63.546", "run 0",
    ]:
        run.command(c)
    with pytest.raises(LammpsExecutionError) as ei:
        run.sync()
    msg = str(ei.value)
    assert "ERROR" in msg
    assert "calphy.seg0.lmp" in msg          # names the segment script
    assert "calphy.seg0.log" in msg          # and the log


def test_error_nonexistent_binary_resolution():
    """resolution names all three lookup steps when nothing resolves."""
    import calphy.runner as R
    saved = os.environ.pop("CALPHY_LAMMPS_EXECUTABLE", None)
    orig_which = R.shutil.which
    try:
        R.shutil.which = lambda v: None
        with pytest.raises(ValueError) as ei:
            R.resolve_lammps_executable("/no/such/lmp")
        m = str(ei.value)
        assert "input key" in m and "CALPHY_LAMMPS_EXECUTABLE" in m and "PATH" in m
    finally:
        R.shutil.which = orig_which
        if saved is not None:
            os.environ["CALPHY_LAMMPS_EXECUTABLE"] = saved


@pytest.mark.slow
def test_error_melted_solid(tmp_path):
    """fe-solid at 5000 K with melt detection on -> MeltedError + error log."""
    with open(os.path.join(BASELINE_INPUTS, "B1.yaml")) as fh:
        data = yaml.safe_load(fh)
    calc = data["calculations"][0]
    calc["pair_coeff"] = _abs_pair_coeff(calc["pair_coeff"])
    calc["temperature"] = 5000.0
    calc["tolerance"] = {"solid_fraction": 0.7}     # re-enable melt detection
    with open(os.path.join(str(tmp_path), "input.yaml"), "w") as fh:
        yaml.safe_dump(data, fh, sort_keys=False)

    kernel_exe = os.path.join(os.path.dirname(sys.executable), "calphy_kernel")
    proc = subprocess.run([kernel_exe, "-i", "input.yaml", "-k", "0"],
                          cwd=str(tmp_path), capture_output=True, text=True, timeout=900)
    assert proc.returncode != 0, "expected a failure at 5000 K"
    assert "MeltedError" in proc.stderr, proc.stderr[-2000:]
    sims = [d for d in glob.glob(os.path.join(str(tmp_path), "*")) if os.path.isdir(d)]
    assert sims and os.path.exists(os.path.join(sims[0], "melted_error.log.lammps"))


# --------------------------------------------------------------------------- #
# 6.4 — MPI variant: restart continuity under -np 2
# --------------------------------------------------------------------------- #
def _mpi_binary():
    try:
        return resolve_mpi_executable(None)
    except ValueError:
        return None


@pytest.mark.skipif(_mpi_binary() is None, reason="no mpirun resolves")
def test_t1_restart_energy_continuity_mpi(tmp_path):
    run = make_runner(tmp_path, cores=2, mpi_command=_mpi_binary())
    boot(run)
    run.command("fix 1 all nvt temp 500 500 0.1")
    run.command("velocity all create 500 12345 mom yes rot yes")
    run.command("run 100")
    run.command('print "PECHK $(pe)"')
    run.sync()
    run.command("run 0")
    run.command('print "PECHK $(pe)"')
    run.sync()
    run.close()

    pe = [v[0] for v in log_values(tmp_path, "PECHK")]
    assert len(pe) == 2
    rel = abs(pe[1] - pe[0]) / max(abs(pe[0]), 1e-12)
    assert rel <= 1e-6, "PE discontinuity under MPI: %r rel=%.2e" % (pe, rel)
