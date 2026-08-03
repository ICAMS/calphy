"""
Irreversibility check on switching paths.

Dissipation enters the free energy directly, and a path whose structure changed
partway through it dissipates orders of magnitude more than one that merely ran
fast.  The numbers below are the real ones from a Ag ``melting_temperature``
run: a clean reversible-scaling sweep sat at 1.7e-4 eV/atom, while the same
sweep with the solid melting near the top sat at 2.3e-2 -- so the threshold has
to separate those two without crying wolf on the clean one.
"""
import logging
import os

import numpy as np
import pytest

from calphy.routines import MeltingTemp
from calphy.solid import Solid

CLEAN = 1.7e-4      # measured: reversible sweep
MELTED = 2.3e-2     # measured: solid melted during the sweep


def _log_text(tmp_path, job):
    for handler in job.logger.handlers:
        handler.flush()
    with open(os.path.join(str(tmp_path), "calphy.log")) as fh:
        return fh.read()


def test_clean_sweep_is_not_flagged(make_calc, recorded_job, tmp_path):
    """A reversible sweep must pass silently -- otherwise the warning is noise."""
    calc = make_calc("B1")
    job, _ = recorded_job(Solid, calc)
    assert job.check_dissipation(CLEAN, "sweep") is False
    assert "WARNING" not in _log_text(tmp_path, job).split("dissipation:")[-1]


def test_melted_sweep_is_flagged(make_calc, recorded_job, tmp_path):
    """The dissipation of a sweep that melted must be called out."""
    calc = make_calc("B1")
    job, _ = recorded_job(Solid, calc)
    assert job.check_dissipation(MELTED, "sweep") is True
    log = _log_text(tmp_path, job)
    assert "WARNING" in log
    assert "far from reversible" in log


def test_threshold_separates_the_two_measured_cases(make_calc, recorded_job):
    """The shipped default sits between the clean and the melted measurement."""
    calc = make_calc("B1")
    default = calc.tolerance.dissipation
    assert CLEAN < default < MELTED, (
        "default tolerance.dissipation=%r does not separate a clean sweep "
        "(%.1e) from a melted one (%.1e)" % (default, CLEAN, MELTED)
    )


def test_check_can_be_disabled(make_calc, recorded_job):
    calc = make_calc("B1", tolerance={"dissipation": 0})
    job, _ = recorded_job(Solid, calc)
    assert job.check_dissipation(MELTED, "sweep") is False


def test_sign_is_ignored(make_calc, recorded_job):
    """Dissipation is a magnitude; a negative value is just as irreversible."""
    calc = make_calc("B1")
    job, _ = recorded_job(Solid, calc)
    assert job.check_dissipation(-MELTED, "sweep") is True


# --------------------------------------------------------------------------- #
# MeltingTemp surfacing
# --------------------------------------------------------------------------- #
class _Job:
    def __init__(self, phase, ediss, high):
        self.calc = type("C", (), {"reference_phase": phase})()
        self.ediss = ediss
        self.ediss_high = high


def _melting_temp(sol, lqd, threshold=1e-3):
    job = object.__new__(MeltingTemp)
    job.logger = logging.getLogger("test_dissipation_check")
    job.logger.handlers = []
    job.logger.propagate = True
    job.soljob, job.lqdjob = sol, lqd
    job.calc = type("C", (), {"tolerance": type("T", (), {"dissipation": threshold})()})()
    return job


def test_melting_temp_repeats_the_warning(caplog):
    """
    The sub-calculations log into their own folders, which nobody reads when
    the point of the mode is one number at the end.
    """
    job = _melting_temp(_Job("solid", MELTED, True), _Job("liquid", 3e-4, False))
    with caplog.at_level(logging.WARNING):
        job._report_sweep_quality()
    text = " ".join(r.message for r in caplog.records)
    assert "unreliable" in text
    assert "solid" in text
    assert "liquid" not in text.split("tolerance.dissipation")[0]


def test_melting_temp_silent_when_both_sweeps_are_clean(caplog):
    job = _melting_temp(_Job("solid", CLEAN, False), _Job("liquid", 3e-4, False))
    with caplog.at_level(logging.WARNING):
        job._report_sweep_quality()
    assert not [r for r in caplog.records if r.levelname == "WARNING"]


def test_melting_temp_survives_jobs_without_the_attribute(caplog):
    """A job that never ran a sweep has no verdict; that is not an error."""
    sol = _Job("solid", 0, False)
    del sol.ediss_high
    job = _melting_temp(sol, _Job("liquid", 3e-4, False))
    job._report_sweep_quality()  # must not raise


# --------------------------------------------------------------------------- #
# harvest
# --------------------------------------------------------------------------- #
def test_gather_results_surfaces_the_verdict(tmp_path):
    """A harvested frame can be filtered on the flag without re-deriving it."""
    import yaml
    from calphy.postprocessing import gather_results

    root = str(tmp_path)
    for name, high in (("run_ok", False), ("run_bad", True)):
        d = os.path.join(root, name)
        os.makedirs(d)
        with open(os.path.join(d, "input_file.yaml"), "w") as fh:
            yaml.safe_dump({"calculations": [{
                "mode": "ts", "temperature": 900, "pressure": 0,
                "reference_phase": "solid", "phase_name": "",
                "element": ["Ag"], "reference_composition": 0.0,
            }]}, fh)
        with open(os.path.join(d, "report.yaml"), "w") as fh:
            yaml.safe_dump({
                "results": {
                    "free_energy": -3.4,
                    "ts_dissipation": MELTED if high else CLEAN,
                    "ts_dissipation_high": high,
                },
                "input": {"element": "Ag", "concentration": "1.0"},
            }, fh)
        t = np.linspace(850, 1250, 20)
        np.savetxt(os.path.join(d, "temperature_sweep.dat"),
                   np.column_stack((t, np.full(20, -3.4), np.full(20, 1e-5))))

    df = gather_results(root).set_index("calculation")
    # pandas may hand these back as numpy bools; compare by value
    assert bool(df.loc["run_bad"]["ts_dissipation_high"]) is True
    assert bool(df.loc["run_ok"]["ts_dissipation_high"]) is False
    assert df.loc["run_bad"]["ts_dissipation"] == pytest.approx(MELTED)
    # the flag is usable as a filter, which is the point of surfacing it
    flagged = df[df["ts_dissipation_high"] == True]  # noqa: E712
    assert list(flagged.index) == ["run_bad"]
