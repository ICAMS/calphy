import pytest
from calphy.input import read_inputfile

def test_options():
	options = read_inputfile("examples/Cu_EAM/input.yaml")
	assert options[0]._temperature == 1300