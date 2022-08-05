import pytest
import calphy.helpers as ch
import numpy as np

def test_nones():
	a = [None, 1, 2]
	assert ch.check_if_any_is_none(a) == True
	assert ch.check_if_any_is_not_none(a) == True

	b = [None, None, None]
	assert ch.check_if_any_is_none(b) == True
	assert ch.check_if_any_is_not_none(b) == False

	c = None
	assert ch.check_if_any_is_none(c) == True
	assert ch.check_if_any_is_not_none(c) == False

	d = 1
	assert ch.check_if_any_is_none(d) == False
	assert ch.check_if_any_is_not_none(d) == True

	d = [1, 2, 3]
	assert ch.check_if_any_is_none(d) == False
	assert ch.check_if_any_is_not_none(d) == True

def test_replace_nones():
	a = [None, 1, 2]
	b = [3, 5, 6]
	c = ch.replace_nones(a, b)
	assert c[0] == 3

def test_validate_spring_constants():

	d = [1, 2, 4]
	e = ch.validate_spring_constants(d)
	assert e[0] == 1

	d = [1, np.NaN, 4]
	e = ch.validate_spring_constants(d)
	assert e[1] == 1
