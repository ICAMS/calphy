import pytest
from calphy.integrators import *

def test_ideal_gas():
	a = get_ideal_gas_fe(1000, 0.07, 1000, [26], [1])
	assert np.abs(a+0.8900504315410337) < 1E-5

def test_uf():
	a = get_uhlenbeck_ford_fe(1000, 0.07, 50, 2)
	assert np.abs(a-5.37158083028874) < 1E-5

def test_solid_ref():
	a = get_einstein_crystal_fe(1000, 1000, [26], 200, [1.2], [1])
	assert np.abs(a+0.4727067261423942) < 1E-5
