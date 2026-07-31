"""
Tm uncertainty from the solid/liquid free-energy crossing.

``find_tm`` rejects only ``arg == 0`` and ``arg == len - 1``, so the local
slope stencil around the crossing has to survive a crossing that lands close
to either end of the sweep: past the top it would raise IndexError in
postprocessing (after all the MD has been paid for), and below the bottom a
negative index is silently wrapped by numpy to the far end, turning the local
slope into a chord across the whole temperature range.
"""
import logging

import numpy as np
import pytest

from calphy.routines import MeltingTemp


def _job():
    """A MeltingTemp with just enough state to exercise _crossing_error."""
    job = object.__new__(MeltingTemp)
    job.logger = logging.getLogger("test_find_tm_crossing_error")
    job.logger.addHandler(logging.NullHandler())
    return job


def _curves(n, t0=850.0, t1=1250.0):
    """Curved F(T) for both phases -- S grows with T, so the curves are not
    straight and a wrapped chord is measurably different from a local slope."""
    T = np.linspace(t0, t1, n)
    f_sol = -2.10 - 2.0e-4 * T - 9.0e-5 * T * np.log(T)
    f_lqd = -2.06 - 2.4e-4 * T - 9.6e-5 * T * np.log(T)
    err = np.full(n, 1.0e-4)
    return (T, f_sol, err), (T, f_lqd, err)


def _reference_tmerr(job, arg, suberr, half):
    """Local, correctly-clamped slope difference computed independently."""
    n = len(job.solres[1])
    lo, hi = max(arg - half, 0), min(arg + half, n - 1)
    ss = (job.solres[1][hi] - job.solres[1][lo]) / (job.solres[0][hi] - job.solres[0][lo])
    ls = (job.lqdres[1][hi] - job.lqdres[1][lo]) / (job.lqdres[0][hi] - job.lqdres[0][lo])
    return suberr / (ss - ls)


@pytest.mark.parametrize("arg", [1, 5, 25, 49])
def test_low_edge_crossing_does_not_wrap(arg):
    """A crossing inside the stencil width must still use a *local* slope."""
    job = _job()
    n = 50000
    job.solres, job.lqdres = _curves(n)
    suberr = 1.0e-4

    got = job._crossing_error(arg, suberr)
    expected = _reference_tmerr(job, arg, suberr, MeltingTemp.CROSSING_STENCIL)
    assert got == pytest.approx(expected, rel=1e-12)

    # and it must differ from what the wrapped (negative-index) read gives
    lo_wrapped = arg - MeltingTemp.CROSSING_STENCIL
    hi = arg + MeltingTemp.CROSSING_STENCIL
    ss = (job.solres[1][hi] - job.solres[1][lo_wrapped]) / (
        job.solres[0][hi] - job.solres[0][lo_wrapped]
    )
    ls = (job.lqdres[1][hi] - job.lqdres[1][lo_wrapped]) / (
        job.lqdres[0][hi] - job.lqdres[0][lo_wrapped]
    )
    wrapped = suberr / (ss - ls)
    assert got != pytest.approx(wrapped, rel=1e-6)


@pytest.mark.parametrize("offset", [2, 10, 50])
def test_high_edge_crossing_does_not_raise(offset):
    """A crossing near the top of the sweep must not IndexError."""
    job = _job()
    n = 50000
    job.solres, job.lqdres = _curves(n)
    arg = n - offset
    tmerr = job._crossing_error(arg, 1.0e-4)
    assert np.isfinite(tmerr)


@pytest.mark.parametrize("n", [12, 40, 101, 500])
def test_short_sweep_is_safe_at_every_interior_index(n):
    """
    The stencil scales with the sweep, so no interior crossing can index out
    of bounds however short the sweep is.
    """
    job = _job()
    job.solres, job.lqdres = _curves(n)
    for arg in range(1, n - 1):
        tmerr = job._crossing_error(arg, 1.0e-4)
        assert np.isfinite(tmerr), f"n={n} arg={arg}"


def test_parallel_curves_report_no_error_instead_of_dividing_by_zero():
    """Identical slopes give no crossing scale; report nan, do not blow up."""
    job = _job()
    n = 1000
    T = np.linspace(850.0, 1250.0, n)
    f = -2.10 - 2.0e-4 * T
    err = np.full(n, 1.0e-4)
    job.solres = (T, f, err)
    job.lqdres = (T, f + 0.01, err)  # parallel, offset
    assert np.isnan(job._crossing_error(n // 2, 1.0e-4))


def test_degenerate_sweep_reports_no_error():
    """Too few samples to form any stencil at all."""
    job = _job()
    T = np.array([900.0, 901.0])
    job.solres = (T, np.array([-3.0, -3.001]), np.array([1e-4, 1e-4]))
    job.lqdres = (T, np.array([-3.0, -3.002]), np.array([1e-4, 1e-4]))
    # n // 20 == 0 -> half clamps to 1, lo=0, hi=1 -> still usable
    assert np.isfinite(job._crossing_error(1, 1.0e-4))
