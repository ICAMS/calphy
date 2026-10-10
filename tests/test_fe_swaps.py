"""Solid fe integration with Monte Carlo swaps (configurational entropy).

Input validation, the reference terms (ideal mixing, CM correction with the
simulated mass) and the assembly of the free energy.  The LAMMPS command
streams are covered by the B14 goldens in test_golden_streams.py, the real
run by test_fe_swaps_lammps.py (needs the thermoatoms fork).
"""
import os

import numpy as np
import pytest

from calphy.solid import Solid
from calphy.integrators import (
    get_einstein_crystal_fe,
    get_ideal_mixing_fe,
    kb,
)
import calphy.helpers as ph

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
ZRCU = os.path.join(REPO, "examples", "example_10", "ZrCu.data")

LOOSE_TOL = {"tolerance": {"pressure": 1e12, "spring_constant": 1e12}}


def binary_calc(make_calc, **overrides):
    calc = make_calc("B14", lattice_sentinel=None, lattice=ZRCU, **LOOSE_TOL, **overrides)
    calc.lattice = "structure.data"
    return calc


# --------------------------------------------------------------------------- #
# input validation
# --------------------------------------------------------------------------- #
def test_fe_swaps_flag_set(make_calc):
    calc = binary_calc(make_calc)
    assert calc._fe_swaps is True
    assert calc.tolerance.swap_energy == pytest.approx(0.001)


def test_fe_swaps_off_without_n_swaps(make_calc):
    calc = binary_calc(make_calc, monte_carlo={"n_swaps": 0})
    assert calc._fe_swaps is False


def test_fe_swaps_need_fork(make_calc):
    with pytest.raises(ValueError, match="thermoatoms LAMMPS fork"):
        binary_calc(make_calc, monte_carlo={"n_swaps": 5, "use_custom_lammps": False})


def test_fe_swaps_need_two_elements(make_calc):
    with pytest.raises(ValueError, match="at least two elements"):
        make_calc("B1", monte_carlo={"n_swaps": 5, "use_custom_lammps": True})


def test_fe_swaps_solid_only(make_calc):
    with pytest.raises(ValueError, match="reference_phase solid"):
        binary_calc(make_calc, reference_phase="liquid")


def test_fe_swaps_not_with_qtb(make_calc):
    with pytest.raises(ValueError, match="quantum thermal bath"):
        binary_calc(make_calc, mode="fe-qtb")


def test_fe_swaps_not_for_other_modes(make_calc):
    # the flag is specific to mode fe; ts keeps its own swap handling
    calc = make_calc(
        "B13", lattice_sentinel=None, lattice=ZRCU, **LOOSE_TOL,
        monte_carlo={"n_swaps": 5, "use_custom_lammps": True},
    )
    assert calc._fe_swaps is False


# --------------------------------------------------------------------------- #
# reference terms
# --------------------------------------------------------------------------- #
def test_mean_mass_is_composition_weighted(make_calc):
    calc = binary_calc(make_calc)
    x = {e: calc._element_dict[e]["composition"] for e in calc.element}
    assert x["Zr"] == pytest.approx(0.5) and x["Cu"] == pytest.approx(0.5)
    assert ph.mean_mass(calc) == pytest.approx(0.5 * (91.224 + 63.546))


def test_ideal_mixing_binary(make_calc):
    calc = binary_calc(make_calc)
    f = get_ideal_mixing_fe(calc)
    assert f == pytest.approx(kb * 800.0 * np.log(0.5))
    assert f < 0


def test_ideal_mixing_unary_is_zero(make_calc):
    calc = make_calc("B1")
    assert get_ideal_mixing_fe(calc) == 0.0


def test_cm_mass_changes_only_the_cm_correction(make_calc):
    calc = binary_calc(make_calc)
    vol, k = 18.0**3, [2.0, 3.0]
    fe, fcm = get_einstein_crystal_fe(calc, vol, k, return_contributions=True)
    fe_u, fcm_u = get_einstein_crystal_fe(
        calc, vol, k, return_contributions=True, cm_mass=ph.mean_mass(calc)
    )
    assert fe_u == fe
    assert fcm_u != fcm
    # the correction is a 1/N effect: both are small and of the same order
    assert abs(fcm_u - fcm) < 0.1 * abs(fcm)


def test_cm_mass_equal_to_real_mass_is_identity(make_calc):
    calc = make_calc("B1")
    fe, fcm = get_einstein_crystal_fe(calc, 18.0**3, [2.0], return_contributions=True)
    fe_u, fcm_u = get_einstein_crystal_fe(
        calc, 18.0**3, [2.0], return_contributions=True, cm_mass=63.546
    )
    assert (fe_u, fcm_u) == pytest.approx((fe, fcm))


# --------------------------------------------------------------------------- #
# free energy assembly
# --------------------------------------------------------------------------- #
def _write_switching_files(folder, n_elements, n_iterations, w_target):
    """Linear lambda ramps whose forward and backward work both equal w_target."""
    lam = np.linspace(0, 1, 201)
    for it in range(1, n_iterations + 1):
        for name, lam_leg in (("forward", lam), ("backward", lam[::-1])):
            # dU_sys - sum(dU_ref) = -w_target at every lambda, so
            # int (dU_ref - dU_sys) dlambda = w_target on both legs
            cols = [np.full_like(lam_leg, -w_target)]
            cols += [np.zeros_like(lam_leg) for _ in range(n_elements)]
            cols += [lam_leg]
            np.savetxt(
                os.path.join(folder, "%s_%d.dat" % (name, it)),
                np.column_stack(cols),
                header="synthetic",
            )


def test_thermodynamic_integration_adds_ideal_mixing(make_calc, recorded_job):
    calc = binary_calc(make_calc)
    job, _ = recorded_job(Solid, calc)
    job.vol = 18.0**3 * (1024 / 500)
    job.k = [2.0, 3.0]
    _write_switching_files(job.simfolder, calc.n_elements, calc.n_iterations, 0.25)
    job.thermodynamic_integration()
    assert job.fmix == pytest.approx(get_ideal_mixing_fe(calc))
    assert job.fe == pytest.approx(job.fref + job.fmix + job.w + job.pv)
    assert job.fmix < 0


def test_thermodynamic_integration_without_swaps_has_no_mixing_term(make_calc, recorded_job):
    calc = binary_calc(make_calc, monte_carlo={"n_swaps": 0})
    job, _ = recorded_job(Solid, calc)
    job.vol = 18.0**3 * (1024 / 500)
    job.k = [2.0, 3.0]
    _write_switching_files(job.simfolder, calc.n_elements, calc.n_iterations, 0.25)
    job.thermodynamic_integration()
    assert job.fmix == 0
    assert job.fe == pytest.approx(job.fref + job.w + job.pv)


def test_report_has_ideal_mixing(make_calc, recorded_job):
    calc = binary_calc(make_calc)
    job, _ = recorded_job(Solid, calc)
    job.vol = 18.0**3 * (1024 / 500)
    job.k = [2.0, 3.0]
    _write_switching_files(job.simfolder, calc.n_elements, calc.n_iterations, 0.25)
    job.thermodynamic_integration()
    job.volatom = job.vol / calc._natoms  # set by run_averaging in a real run
    job.submit_report()
    assert job.report["results"]["ideal_mixing"] == pytest.approx(job.fmix)
