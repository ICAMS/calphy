import pytest
import calphy.helpers as ch

def test_nones():
	a = [None, 1, 2]
	assert ch.check_if_none(a) == False

	b = [None, None, None]
	assert ch.check_if_none(b) == True

	c = None
	assert ch.check_if_none(c) == True

	d = 1
	assert ch.check_if_none(d) == False

def test_replace_nones():
	a = [None, 1, 2]
	b = [3, 5, 6]
	c = ch.replace_nones(a, b)
	assert c[0] == 3
