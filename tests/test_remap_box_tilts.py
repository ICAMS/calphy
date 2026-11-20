"""The reference cell includes its tilt factors (pmrv/calphy#13).

The backward reversible-scaling leg reads the hot end of the forward leg and
resets the box to the cold reference cell.  ``change_box ... y final`` changes
the lengths only, so the hot tilts used to stay on the cold lengths: a shear
that ``fix npt iso/aniso`` (which rescales the tilts with the lengths) then
kept for the whole leg.  The same holds for the isotropic ``change_box ...
scale`` of the pressure convergence.

Recorded-stream tests check what is emitted; the ``lammps`` tests check, with
the real binary, that the emitted commands restore the cell -- and that the
old command did not, so the check can fail.
"""
import os
import shutil

import numpy as np
import pytest

import calphy.helpers as ph
from calphy.solid import Solid
from calphy.runner import ExecutableRunner, resolve_lammps_executable

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(HERE, "fixtures")
POT = os.path.join(HERE, "Cu01.eam.alloy")
CONF = os.path.join(FIXTURE_DIR, "conf.equilibration.data")

LOOSE_TOL = {"tolerance": {"pressure": 1e12, "spring_constant": 1e12}}


class _Rec:
    def __init__(self):
        self.commands = []

    def command(self, s):
        self.commands.append(s)


# --------------------------------------------------------------------------- #
# command strings
# --------------------------------------------------------------------------- #
ORTHO = "change_box     all x final 0.0 10.000000 y final 0.0 11.000000 z final 0.0 12.000000 remap units box"


@pytest.mark.parametrize("tilts", [(), (None, None, None), (0.0, 0.0, 0.0)])
def test_remap_box_orthogonal_is_unchanged(tilts):
    rec = _Rec()
    ph.remap_box(rec, 10.0, 11.0, 12.0, *tilts)
    assert rec.commands == ["run 0", ORTHO]


def test_remap_box_sets_tilts_of_a_triclinic_cell():
    rec = _Rec()
    ph.remap_box(rec, 10.0, 11.0, 12.0, -5.0, 0.0, 1.25)
    assert rec.commands[-1] == (
        "change_box     all x final 0.0 10.000000 y final 0.0 11.000000 z final 0.0 12.000000"
        " xy final -5.000000 xz final 0.000000 yz final 1.250000 remap units box"
    )


def test_scale_box_command():
    assert ph.scale_box_command(1.01) == (
        "change_box       all x scale 1.010000 y scale 1.010000 z scale 1.010000 remap"
    )
    assert ph.scale_box_command(1.01, triclinic=True) == (
        "change_box       all x scale 1.010000 y scale 1.010000 z scale 1.010000"
        " xy final $(xy*1.010000) xz final $(xz*1.010000) yz final $(yz*1.010000) remap"
    )


def test_avg_columns_keep_the_historical_order():
    # columns 1-7 are read by position elsewhere (usecols=(1, 2, 3, 4) etc.)
    assert [ph.AVG_COLUMNS[n] for n in ("mlx", "mly", "mlz", "mpress", "mpe", "metotal", "mtemp")] \
        == [1, 2, 3, 4, 5, 6, 7]
    assert [ph.AVG_COLUMNS[n] for n in ("mxy", "mxz", "myz")] == [8, 9, 10]


# --------------------------------------------------------------------------- #
# the averaging stage stores the tilts, the backward leg restores them
# --------------------------------------------------------------------------- #
def _triclinic_avg_dat(path, xy=-9.0, xz=0.5, yz=-0.25):
    """The orthogonal fixture with constant non-zero tilt columns."""
    rows = []
    with open(os.path.join(FIXTURE_DIR, "avg.dat")) as fh:
        for line in fh:
            if line.startswith("#"):
                rows.append(line.rstrip("\n"))
                continue
            f = line.split()
            f[8], f[9], f[10] = repr(xy), repr(xz), repr(yz)
            rows.append(" ".join(f))
    with open(path, "w") as fh:
        fh.write("\n".join(rows) + "\n")


def test_averaging_stores_zero_tilts_for_an_orthogonal_cell(make_calc, recorded_job):
    calc = make_calc("B1", **LOOSE_TOL)
    job, _ = recorded_job(Solid, calc, fixtures=["avg.dat", "msd.dat"])
    job.run_averaging()
    assert (job.xy, job.xz, job.yz) == (0.0, 0.0, 0.0)


def test_averaging_stores_tilts_and_forward_leg_restores_them(make_calc, recorded_job):
    calc = make_calc("B6", **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc, fixtures=["avg.dat", "msd.dat"])
    _triclinic_avg_dat(os.path.join(job.simfolder, "avg.dat"))
    job.run_averaging()
    assert (job.xy, job.xz, job.yz) == (-9.0, 0.5, -0.25)

    # the forward leg opens a runner of its own, which the recorder now follows
    job._reversible_scaling_forward(iteration=1)
    remaps = [c for c in rec.commands if c.startswith("change_box")]
    assert len(remaps) == 1
    assert remaps[0].endswith(
        "xy final -9.000000 xz final 0.500000 yz final -0.250000 remap units box"
    )


def test_backward_leg_continues_from_the_forward_end_state(make_calc, recorded_job):
    """The backward sweep starts from conf.ts.forward_<i>.data as written: no
    reset of the cell to the T0 reference box at the turning point."""
    calc = make_calc("B6", **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc)
    job.lx = job.ly = job.lz = 18.075
    job.xy, job.xz, job.yz = -9.0, 0.5, -0.25
    job._reversible_scaling_backward(iteration=1)
    assert any(c.startswith("read_data") and "conf.ts.forward_1.data" in c for c in rec.commands)
    assert not [c for c in rec.commands if c.startswith("change_box")]


# --------------------------------------------------------------------------- #
# real LAMMPS
# --------------------------------------------------------------------------- #
try:
    BINARY = resolve_lammps_executable(None)
except ValueError:
    BINARY = None

def needs_lammps(f):
    return pytest.mark.lammps(
        pytest.mark.skipif(BINARY is None, reason="no lmp binary resolves")(f)
    )


def _runner(directory):
    run = ExecutableRunner(
        binary=BINARY, mpi_command=None, cores=1, cmdargs="",
        directory=str(directory), dry_run=False, timeout=300,
    )
    for c in [
        "units metal", "boundary p p p", "atom_style atomic", "timestep 0.001",
        "box tilt large",
        "pair_style eam/alloy",
        "read_data %s" % CONF,
        "pair_coeff * * %s Cu" % POT,
        "mass 1 63.546",
    ]:
        run.command(c)
    run.command("change_box all triclinic")
    return run


def _cell(directory, tag):
    vals = []
    for name in sorted(os.listdir(str(directory))):
        if name.startswith("calphy.seg") and name.endswith(".log"):
            with open(os.path.join(str(directory), name)) as fh:
                vals += [
                    [float(x) for x in line.split()[1:]]
                    for line in fh if line.startswith(tag + " ")
                ]
    assert len(vals) == 1, vals
    return dict(zip(("lx", "ly", "lz", "xy", "xz", "yz"), vals[0]))


PRINT_CELL = 'print "%s $(lx) $(ly) $(lz) $(xy) $(xz) $(yz)"'


def _hot_cell(run, ref):
    """Set the reference cell, then a 'hot' one: every vector 3 % longer."""
    run.command(
        "change_box all x final 0.0 %f y final 0.0 %f z final 0.0 %f "
        "xy final %f xz final %f yz final %f remap units box"
        % tuple(ref[k] for k in ("lx", "ly", "lz", "xy", "xz", "yz"))
    )
    run.command("change_box all x scale 1.03 y scale 1.03 z scale 1.03 "
                "xy final $(xy*1.03) xz final $(xz*1.03) yz final $(yz*1.03) remap")


@needs_lammps
@pytest.mark.parametrize("restore_tilts", [False, True])
def test_remap_box_restores_the_reference_cell(tmp_path, restore_tilts):
    ref = dict(lx=18.0, ly=17.0, lz=16.0, xy=-6.0, xz=2.0, yz=-1.5)
    run = _runner(tmp_path)
    _hot_cell(run, ref)
    tilts = (ref["xy"], ref["xz"], ref["yz"]) if restore_tilts else ()
    ph.remap_box(run, ref["lx"], ref["ly"], ref["lz"], *tilts)
    run.command(PRINT_CELL % "CELL")
    run.command("run 0")
    run.sync()
    run.close()
    got = _cell(tmp_path, "CELL")
    for k in ("lx", "ly", "lz"):
        assert got[k] == pytest.approx(ref[k], abs=1e-5)
    restored = all(got[k] == pytest.approx(ref[k], abs=1e-5) for k in ("xy", "xz", "yz"))
    # without the tilts the hot ones survive: that is the bug this guards against
    assert restored is restore_tilts
    if not restore_tilts:
        assert got["xy"] == pytest.approx(1.03 * ref["xy"], abs=1e-5)


@needs_lammps
def test_scale_box_command_keeps_the_cell_shape(tmp_path):
    run = _runner(tmp_path)
    _hot_cell(run, dict(lx=18.0, ly=17.0, lz=16.0, xy=-6.0, xz=2.0, yz=-1.5))
    run.command(PRINT_CELL % "BEFORE")
    run.command(ph.scale_box_command(0.98, triclinic=True))
    run.command(PRINT_CELL % "AFTER")
    run.command("run 0")
    run.sync()
    run.close()
    a, b = _cell(tmp_path, "BEFORE"), _cell(tmp_path, "AFTER")
    for k in a:
        assert b[k] == pytest.approx(0.98 * a[k], rel=1e-6), k


# --------------------------------------------------------------------------- #
# the P-V fit rescale of the pressure convergence keeps a triclinic cell's shape
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("xy", [0.0, -4.0])
def test_pv_fit_rescale_scales_tilts_of_a_triclinic_cell(make_calc, recorded_job, xy):
    from test_pressure_convergence_pv_fit import _PVStub

    class _TiltedPVStub(_PVStub):
        def command(self, cmd):
            super().command(cmd)
            for row in self.rows:
                row[8] = xy

    calc = make_calc("B1")
    calc.tolerance.pressure = 1.0
    job, _ = recorded_job(Solid, calc)
    ncount = int(calc.md.n_small_steps) // int(
        calc.md.n_every_steps * calc.md.n_repeat_steps
    )
    stub = _TiltedPVStub(natoms=job.natoms, ncount=ncount)
    job.run_pressure_convergence(stub)

    rescales = [c for c in stub.commands if c.startswith("change_box")]
    assert rescales, "P-V fit path was never exercised"
    for c in rescales:
        assert ("xy final $(xy*" in c) is (xy != 0.0), c
    assert job.xy == xy
