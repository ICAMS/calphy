"""Golden command-stream tests (Part 2).

Each test drives one calphy method with a RecordingRunner (no LAMMPS binary) and
snapshots the emitted command stream against tests/golden/<name>.txt.

Loops that converge on data (pressure, spring constant) are forced to converge on
the first cycle by loosening the convergence tolerances -- the tolerances are a
Python-side check and never appear in the command stream, so this does not change
the recorded commands, only the number of cycles (the multi-cycle repetition is
exercised with real LAMMPS in Part 6).  The pre-set box/spring values used for the
integration / sweep methods stand in for what run_averaging would have measured.
"""
import pytest

from calphy.solid import Solid
from calphy.liquid import Liquid
from calphy.alchemy import Alchemy

from conftest import assert_golden

# Loosen convergence so the static avg.dat/msd.dat fixtures converge on cycle 1.
LOOSE_TOL = {"tolerance": {"pressure": 1e12, "spring_constant": 1e12}}

# Deterministic stand-ins for values run_averaging would set on the job.
BOX = 18.075
KSPRING = [2.0]


def _set_state(job, box=BOX, k=None):
    job.lx = job.ly = job.lz = box
    if k is not None:
        job.k = k


# --------------------------------------------------------------------------- #
# Solid fe: averaging (three equilibration paths) + integration
# --------------------------------------------------------------------------- #
def test_solid_fe_averaging(make_calc, recorded_job):
    calc = make_calc("B1", **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc, fixtures=["avg.dat", "msd.dat"])
    job.run_averaging()
    assert_golden(rec.commands, "solid_fe_averaging")


def test_solid_fe_averaging_finite_p(make_calc, recorded_job):
    calc = make_calc("B2", **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc, fixtures=["avg.dat", "msd.dat"])
    job.run_averaging()
    assert_golden(rec.commands, "solid_fe_averaging_finite_p")


def test_solid_fe_averaging_fixlattice(make_calc, recorded_job):
    calc = make_calc("B3", **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc, fixtures=["avg.dat", "msd.dat"])
    job.run_averaging()
    assert_golden(rec.commands, "solid_fe_averaging_fixlattice")


def test_solid_fe_integration(make_calc, recorded_job):
    calc = make_calc("B1", **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc)
    _set_state(job, k=KSPRING)
    job.run_integration(iteration=1)
    assert_golden(rec.commands, "solid_fe_integration")


# --------------------------------------------------------------------------- #
# Liquid fe: averaging (melt cycle / rattle) + integration
# --------------------------------------------------------------------------- #
def test_liquid_fe_averaging_meltcycle(make_calc, recorded_job):
    calc = make_calc("B4", **LOOSE_TOL)
    # solid-fraction sequence forces the melt loop to iterate 3 times (2 "still
    # solid" then a melted reading that breaks the loop).
    n = calc._natoms
    job, rec = recorded_job(
        Liquid, calc, fixtures=["avg.dat"], solid_fraction_seq=[n, n, 0]
    )
    job.run_averaging()
    assert_golden(rec.commands, "liquid_fe_averaging_meltcycle")


def test_liquid_fe_averaging_rattle(make_calc, recorded_job):
    calc = make_calc("B5", **LOOSE_TOL)
    job, rec = recorded_job(Liquid, calc, fixtures=["avg.dat"])
    job.run_averaging()
    assert_golden(rec.commands, "liquid_fe_averaging_rattle")


def test_liquid_fe_integration(make_calc, recorded_job):
    calc = make_calc("B4", **LOOSE_TOL)
    job, rec = recorded_job(Liquid, calc)
    _set_state(job)
    job.run_integration(iteration=1)
    assert_golden(rec.commands, "liquid_fe_integration")


# --------------------------------------------------------------------------- #
# ts (reversible scaling): forward / backward / uniform-temperature schedule
# --------------------------------------------------------------------------- #
def test_ts_forward_solid(make_calc, recorded_job):
    calc = make_calc("B6", **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc)
    _set_state(job)
    job._reversible_scaling_forward(iteration=1)
    assert_golden(rec.commands, "ts_forward_solid")


def test_ts_backward_solid(make_calc, recorded_job):
    calc = make_calc("B6", **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc)
    _set_state(job)
    job._reversible_scaling_backward(iteration=1)
    assert_golden(rec.commands, "ts_backward_solid")


def test_ts_uniform_temperature(make_calc, recorded_job):
    calc = make_calc("B6", lambda_schedule="uniform_temperature", **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc)
    _set_state(job)
    job._reversible_scaling_forward(iteration=1)
    assert_golden(rec.commands, "ts_uniform_temperature")


# --------------------------------------------------------------------------- #
# tscale / pscale
# --------------------------------------------------------------------------- #
def test_tscale(make_calc, recorded_job):
    calc = make_calc("B8", **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc)
    _set_state(job)
    job.temperature_scaling(iteration=1)
    assert_golden(rec.commands, "tscale")


def test_pscale(make_calc, recorded_job):
    calc = make_calc("B9", **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc)
    _set_state(job)
    job.pressure_scaling(iteration=1)
    assert_golden(rec.commands, "pscale")


# --------------------------------------------------------------------------- #
# alchemy
# --------------------------------------------------------------------------- #
def test_alchemy_averaging(make_calc, recorded_job):
    calc = make_calc("B10", **LOOSE_TOL)
    job, rec = recorded_job(Alchemy, calc, fixtures=["avg.dat"])
    job.run_averaging()
    assert_golden(rec.commands, "alchemy_averaging")


def test_alchemy_integration(make_calc, recorded_job):
    calc = make_calc("B10", **LOOSE_TOL)
    job, rec = recorded_job(Alchemy, calc)
    _set_state(job)
    job.run_integration(iteration=1)
    assert_golden(rec.commands, "alchemy_integration")


# --------------------------------------------------------------------------- #
# QTB (fe-qtb): averaging + integration.  Command recording needs no QTB binary.
# --------------------------------------------------------------------------- #
def test_qtb_averaging(make_calc, recorded_job):
    calc = make_calc("B11", **LOOSE_TOL)
    assert calc._qtb is True and calc.mode == "fe"
    job, rec = recorded_job(Solid, calc, fixtures=["avg.dat", "msd.dat"])
    job.run_averaging()
    assert_golden(rec.commands, "qtb_averaging")


def test_qtb_integration(make_calc, recorded_job):
    calc = make_calc("B11", **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc)
    _set_state(job, k=KSPRING)
    job.run_integration(iteration=1)
    assert_golden(rec.commands, "qtb_integration")


# --------------------------------------------------------------------------- #
# Overlay potential (hybrid/overlay of two component pair styles)
# --------------------------------------------------------------------------- #
OVERLAY = {
    "pair_mode": "overlay",
    "pair_style": ["eam/alloy", "eam/alloy"],
    "pair_coeff": [
        "* * tests/Cu01.eam.alloy Cu",
        "* * tests/Cu01.eam.alloy Cu",
    ],
}


def test_overlay_potential(make_calc, recorded_job):
    calc = make_calc("B1", **OVERLAY, **LOOSE_TOL)
    job, rec = recorded_job(Solid, calc, fixtures=["avg.dat", "msd.dat"])
    job.run_averaging()
    assert_golden(rec.commands, "overlay_potential")


# --------------------------------------------------------------------------- #
# ts with Monte Carlo swaps (forward integration, two swap types)
# --------------------------------------------------------------------------- #
def test_ts_mc_swaps(make_calc, recorded_job):
    calc = make_calc(
        "B6",
        monte_carlo={
            "n_swaps": 5,
            "n_steps": 100,
            "forward_swap_types": [1, 2],
            "reverse_swap_types": [1, 2],
        },
        **LOOSE_TOL,
    )
    job, rec = recorded_job(Solid, calc)
    _set_state(job)
    job._reversible_scaling_forward(iteration=1)
    assert_golden(rec.commands, "ts_mc_swaps")


# --------------------------------------------------------------------------- #
# Prescan (phase_transition_detection) -- DEFERRED.
# calphy.range_scan is missing from the codebase (pre-existing bug): any
# mode != none raises ModuleNotFoundError.  See tests/baselines/BASELINE_INFO.md
# (B12 deferral).  Kept as a skipped placeholder so the gap is visible.
# --------------------------------------------------------------------------- #
@pytest.mark.skip(
    reason="prescan blocked: calphy.range_scan module missing (pre-existing bug), "
    "phase_transition_detection.mode != none raises ModuleNotFoundError; see "
    "tests/baselines/BASELINE_INFO.md (B12 deferral)."
)
def test_prescan(make_calc, recorded_job):
    calc = make_calc("B12")
    job, rec = recorded_job(Solid, calc)
    _set_state(job)
    job.scan_temperature_range()
    assert_golden(rec.commands, "prescan")
