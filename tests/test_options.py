import pytest
from calphy.input import read_yamlfile

def test_options():
	options = read_yamlfile("examples/Cu_EAM/input.yaml")
	assert options["calculations"][0]["temperature"] == 1300