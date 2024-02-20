import pytest
from calphy.integrators import *

def test_ideal_gas():
	a = get_ideal_gas_fe(1000, 0.07, 1000, [26], [1])
	assert np.abs(a+0.8900504315410337) < 1E-5

def test_uf():
	a = get_uhlenbeck_ford_fe(1000, 0.07, 50, 2)
	assert np.abs(a-5.37158083028874) < 1E-5
