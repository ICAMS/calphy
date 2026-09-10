"""LAMMPS semantics the dynamic Clausius-Clapeyron driver relies on (DCCI_PLAN Part 0).

The d-CCI sweep is a Python block loop over two LAMMPS sessions.  It needs three
things from LAMMPS that no other calphy mode exercises, checked here against the
real binary through the ExecutableRunner (one block = one segment):

1. a single ``fix npt`` with a linear pressure ramp keeps following the *global*
   sweep progress when the sweep is chopped into ``run n start 0 stop N`` blocks;
2. the ``hybrid/scaled`` scale variable can be redefined between runs and takes
   effect immediately, also across a segment boundary;
3. ``fix ave/time 1 n n`` writes exactly one row per block after a restart.

Marked ``lammps``; skipped when no ``lmp`` binary resolves.
"""
import os
import re
import glob

import numpy as np
import pytest

from calphy.runner import ExecutableRunner, resolve_lammps_executable, read_timeseries

HERE = os.path.dirname(os.path.abspath(__file__))
POT = os.path.join(HERE, "Cu01.eam.alloy")
CONF = os.path.join(HERE, "fixtures", "conf.equilibration.data")  # Cu fcc, 500 atoms

try:
    BINARY = resolve_lammps_executable(None)
except ValueError:
    BINARY = None

pytestmark = [
    pytest.mark.lammps,
    pytest.mark.skipif(BINARY is None, reason="no lmp binary resolves"),
]

T = 500.0
P0, P1 = 0.0, 100000.0  # bar
N_BLOCK, N_BLOCKS = 1000, 3
N_TOTAL = N_BLOCK * N_BLOCKS


def make_runner(directory):
    os.makedirs(str(directory), exist_ok=True)
    return ExecutableRunner(
        binary=BINARY, mpi_command=None, cores=1, cmdargs="",
        directory=str(directory), dry_run=False, timeout=600,
    )


def boot(run, pair_style="pair_style eam/alloy", pair_coeff="pair_coeff * * %s Cu" % POT):
    for c in [
        "units metal", "boundary p p p", "atom_style atomic", "timestep 0.001",
        pair_style, "read_data %s" % CONF, pair_coeff, "mass 1 63.546",
    ]:
        run.command(c)


def seg_logs(directory):
    return sorted(
        glob.glob(os.path.join(str(directory), "calphy.seg*.log")),
        key=lambda p: int(re.search(r"seg(\d+)", p).group(1)),
    )


def thermo_rows(directory):
    """(step, col1, col2, ...) rows of every thermo table in every segment log."""
    rows = []
    for logpath in seg_logs(directory):
        cols = None
        with open(logpath) as fh:
            for line in fh:
                s = line.split()
                if not s:
                    continue
                if s[0] == "Step":
                    cols = s
                    continue
                if cols is None:
                    continue
                try:
                    rows.append([float(x) for x in s])
                except ValueError:
                    cols = None
    return np.array(rows)


def printed(directory, tag):
    vals = []
    for logpath in seg_logs(directory):
        with open(logpath) as fh:
            for line in fh:
                if line.startswith(tag + " "):
                    vals.append([float(x) for x in line.split()[1:]])
    return vals


def _window_mean(rows, lo, hi, col=1):
    sel = (rows[:, 0] > lo) & (rows[:, 0] <= hi)
    assert sel.any(), "no thermo rows in (%d, %d]" % (lo, hi)
    return rows[sel, col].mean()


def _ramped_npt(run):
    run.command("velocity all create %f 4321 mom yes rot yes" % T)
    run.command("fix 1 all npt temp %f %f 0.1 iso %f %f 0.1" % (T, T, P0, P1))
    run.command("thermo_style custom step press vol")
    run.command("thermo 10")


# --------------------------------------------------------------------------- #
# 1. fix npt ramp under `run n start 0 stop N` blocks
# --------------------------------------------------------------------------- #
def test_npt_ramp_follows_global_progress_across_blocks(tmp_path):
    # reference: the whole ramp in one run
    ref = make_runner(tmp_path / "single")
    boot(ref)
    _ramped_npt(ref)
    ref.command("run %d" % N_TOTAL)
    ref.sync()
    ref.close()
    single = thermo_rows(tmp_path / "single")

    # blocks with start/stop: the target must keep ramping across segments
    blk = make_runner(tmp_path / "blocks")
    boot(blk)
    _ramped_npt(blk)
    for _ in range(N_BLOCKS):
        blk.command("run %d start 0 stop %d" % (N_BLOCK, N_TOTAL))
        blk.sync()
    blk.close()
    blocks = thermo_rows(tmp_path / "blocks")

    # control without start/stop: every segment restarts the ramp from P0
    ctl = make_runner(tmp_path / "control")
    boot(ctl)
    _ramped_npt(ctl)
    for _ in range(N_BLOCKS):
        ctl.command("run %d" % N_BLOCK)
        ctl.sync()
    ctl.close()
    control = thermo_rows(tmp_path / "control")

    assert blocks[:, 0].max() == N_TOTAL and blocks.shape[0] >= N_TOTAL // 10

    # Right after the first boundary the target is ~1/3 of the way up the
    # ramp; the control has just been reset to P0 and its pressure collapses.
    lo, hi = N_BLOCK, N_BLOCK + 300
    p_single = _window_mean(single, lo, hi)
    p_blocks = _window_mean(blocks, lo, hi)
    p_control = _window_mean(control, lo, hi)
    print("press (%d,%d]: single %.0f  blocks %.0f  control %.0f" % (lo, hi, p_single, p_blocks, p_control))
    assert abs(p_blocks - p_single) < 0.15 * P1
    assert p_control < 0.5 * p_single

    # and the end of the sweep reaches P1 in both
    assert abs(_window_mean(blocks, N_TOTAL - 300, N_TOTAL) - P1) < 0.15 * P1
    assert abs(_window_mean(single, N_TOTAL - 300, N_TOTAL) - P1) < 0.15 * P1


# --------------------------------------------------------------------------- #
# 2. hybrid/scaled picks up a redefined scale variable
# --------------------------------------------------------------------------- #
def test_hybrid_scaled_variable_redefinition(tmp_path):
    run = make_runner(tmp_path)
    run.command("variable lam equal 1.0")
    boot(
        run,
        pair_style="pair_style hybrid/scaled v_lam eam/alloy",
        pair_coeff="pair_coeff * * eam/alloy %s Cu" % POT,
    )
    run.command("compute pr all pair eam/alloy")  # unscaled sub-style energy
    run.command("velocity all create %f 4321 mom yes rot yes" % T)
    run.command("fix 1 all nvt temp %f %f 0.1" % (T, T))
    run.command("thermo_style custom step pe c_pr")
    run.command("thermo 10")
    run.command("run 20")
    run.command('print "RATIO $(v_lam) $(pe/c_pr)"')
    # redefinition inside one segment, with pre no / post no as the block loop
    run.command("variable lam equal 0.5")
    run.command("run 20 pre no post no")
    run.command('print "RATIO $(v_lam) $(pe/c_pr)"')
    run.sync()
    # redefinition across a segment boundary (variable replayed, then replaced)
    run.command("variable lam equal 0.25")
    run.command("run 20")
    run.command('print "RATIO $(v_lam) $(pe/c_pr)"')
    run.sync()
    run.close()

    ratios = printed(tmp_path, "RATIO")
    assert [r[0] for r in ratios] == [1.0, 0.5, 0.25]
    for lam, ratio in ratios:
        assert ratio == pytest.approx(lam, rel=1e-6), ratios


# --------------------------------------------------------------------------- #
# 3. fix ave/time: one row per block, across restarts
# --------------------------------------------------------------------------- #
def test_ave_time_one_row_per_block(tmp_path):
    run = make_runner(tmp_path)
    boot(run)
    run.command("velocity all create %f 4321 mom yes rot yes" % T)
    run.command("fix 1 all npt temp %f %f 0.1 iso %f %f 0.1" % (T, T, P0, P0))
    run.command("variable pe_atom equal pe/atoms")
    run.command("variable vol_atom equal vol/atoms")
    run.command(
        "fix fav all ave/time 1 %d %d v_pe_atom v_vol_atom file blocks.dat"
        % (N_BLOCK, N_BLOCK)
    )
    for _ in range(N_BLOCKS):
        run.command("run %d start 0 stop %d" % (N_BLOCK, N_TOTAL))
        run.sync()
    run.close()

    data = read_timeseries(str(tmp_path), "blocks.dat")
    assert data.shape == (N_BLOCKS, 3), data
    assert list(data[:, 0]) == [N_BLOCK * (k + 1) for k in range(N_BLOCKS)]
    assert np.all(data[:, 1] < 0) and np.all(data[:, 2] > 10)
