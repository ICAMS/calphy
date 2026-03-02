import pytest
import tempfile
import os
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

	d = [1, np.nan, 4]
	e = ch.validate_spring_constants(d)
	assert e[1] == 1

def test_prepare_log_no_handler_accumulation():
    """
    Calling prepare_log twice for the same file should not accumulate handlers.
    Each call should reset to exactly one file handler.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        logfile = os.path.join(tmpdir, "test.log")

        logger1 = ch.prepare_log(logfile)
        assert len(logger1.handlers) == 1, "Should have exactly 1 handler after first call"

        logger2 = ch.prepare_log(logfile)
        assert len(logger2.handlers) == 1, "Should still have exactly 1 handler after second call"
        assert logger1 is logger2, "Should return the same logger object"

def test_prepare_log_no_cross_contamination():
    """
    Two independent calculations logging to different files must not
    write into each other's log files.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        log1 = os.path.join(tmpdir, "calc1.log")
        log2 = os.path.join(tmpdir, "calc2.log")

        logger1 = ch.prepare_log(log1)
        logger1.info("message from calc1")

        logger2 = ch.prepare_log(log2)
        logger2.info("message from calc2")

        # Flush all handlers
        for h in logger1.handlers:
            h.flush()
        for h in logger2.handlers:
            h.flush()

        content1 = open(log1).read()
        content2 = open(log2).read()

        assert "message from calc1" in content1, "calc1.log should contain calc1 message"
        assert "message from calc2" not in content1, "calc1.log must NOT contain calc2 message"
        assert "message from calc2" in content2, "calc2.log should contain calc2 message"
        assert "message from calc1" not in content2, "calc2.log must NOT contain calc1 message"
