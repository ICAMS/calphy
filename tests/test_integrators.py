import pytest
from calphy.integrators import *

def test_ideal_gas():
	a = get_ideal_gas_fe(1000, 0.07, 1000, [26], [1])
	assert np.abs(a+0.8900504315410337) < 1E-5

def test_uf():
	a = get_uhlenbeck_ford_fe(1000, 0.07, 50, 2)
	assert np.abs(a-5.37158083028874) < 1E-5


def test_integrate_rs_reports_the_worst_iteration(tmp_path):
    # Two sweeps of a constant potential energy U; the second has a backward
    # offset of 0.01 eV/atom, so max |wf - wb| / (2 lambda) = 0.005 at lambda
    # = 0.5.  The dissipation flag must see the bad sweep, not the clean one.
    import numpy as np

    lam = np.linspace(1.0, 0.5, 101)
    U = -3.0
    for i, offset in [(1, 0.0), (2, 0.01)]:
        fwd = np.column_stack((lam * U, 0 * lam, 1000 + 0 * lam, lam))
        bwd = np.column_stack((lam * (U + offset), 0 * lam, 1000 + 0 * lam, lam))[::-1]
        np.savetxt(tmp_path / ("ts.forward_%d.dat" % i), fwd)
        np.savetxt(tmp_path / ("ts.backward_%d.dat" % i), bwd)

    _, ediss = integrate_rs(str(tmp_path), -3.5, 400.0, 100, p=0, nsims=2,
                            scale_energy=True)
    assert ediss == pytest.approx(0.005)
