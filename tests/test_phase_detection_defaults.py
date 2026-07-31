"""
Structural phase-stability detection defaults.

The melt / solidification checks stay OFF by default -- every mode except
``melting_temperature`` merely warns -- while ``melting_temperature``, whose
temperature bracket is advanced *only* by MeltedError / SolidifiedError, turns
them on for its own sub-calculations.
"""

import os
import tempfile

import numpy as np
import pytest
import yaml

from calphy.input import read_inputfile
from calphy.routines import MeltingTemp
from calphy.solid import Solid


# --------------------------------------------------------------------------- #
# detection defaults
# --------------------------------------------------------------------------- #
STRUCTURE = """LAMMPS data file

1 atoms
1 atom types

0.0 4.05 xlo xhi
0.0 4.05 ylo yhi
0.0 4.05 zlo zhi

Masses

1 26.98

Atoms

1 1 0.0 0.0 0.0
"""


def _melting_temp_job(tmpdir, tolerance=None):
    """Build a MeltingTemp whose sub-calculations are ready to inspect."""
    structure_file = os.path.join(tmpdir, "conf.data")
    with open(structure_file, "w") as fh:
        fh.write(STRUCTURE)

    calc = {
        "mode": "melting_temperature",
        "lattice": structure_file,
        "reference_phase": "solid",
        "temperature": 1000,
        "pressure": 0,
        "element": ["Al"],
        "mass": [26.98],
        "pair_style": ["eam/alloy"],
        "pair_coeff": ["* * Al.eam.alloy Al"],
    }
    if tolerance is not None:
        calc["tolerance"] = tolerance

    inputfile = os.path.join(tmpdir, "input.yaml")
    with open(inputfile, "w") as fh:
        yaml.safe_dump({"calculations": [calc]}, fh)

    calculations = read_inputfile(inputfile)
    prev = os.getcwd()
    os.chdir(tmpdir)
    try:
        job = MeltingTemp(calculation=calculations[0], simfolder=tmpdir)
        job.prepare_calcs()
    finally:
        os.chdir(prev)
    return job


def test_detection_off_by_default():
    """Every other mode keeps the checks unreachable -- the shipped default."""
    from calphy.input import Tolerance

    tol = Tolerance()
    assert tol.solid_fraction == 0.0, "melt check must stay off by default"
    assert tol.liquid_fraction == 1.0, "solidify check must stay off by default"


def test_melting_temperature_enables_detection():
    """
    melting_temperature cannot bracket Tm without the checks: run_jobs signals
    'melted'/'froze' only via the exceptions, so they must be live for both
    sub-calculations.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        job = _melting_temp_job(tmpdir)

        assert len(job.calculations) == 2
        for sub in job.calculations:
            assert sub.tolerance.solid_fraction > 0, (
                "melt detection unreachable for the %s sub-calculation"
                % sub.reference_phase
            )
            assert sub.tolerance.liquid_fraction < 1, (
                "solidify detection unreachable for the %s sub-calculation"
                % sub.reference_phase
            )


def test_melting_temperature_respects_explicit_tolerance():
    """An explicit user setting is never overridden, including a deliberate off."""
    with tempfile.TemporaryDirectory() as tmpdir:
        job = _melting_temp_job(
            tmpdir, tolerance={"solid_fraction": 0.42, "liquid_fraction": 0.11}
        )
        for sub in job.calculations:
            assert sub.tolerance.solid_fraction == pytest.approx(0.42)
            assert sub.tolerance.liquid_fraction == pytest.approx(0.11)


@pytest.mark.parametrize(
    "phase, tolerance, expect_warning",
    [
        ("solid", {"solid_fraction": 0.0}, True),
        ("solid", {"solid_fraction": 0.7}, False),
        ("liquid", {"liquid_fraction": 1.0}, True),
        ("liquid", {"liquid_fraction": 0.05}, False),
    ],
)
def test_disabled_detection_is_warned(make_calc, recorded_job, tmp_path,
                                      phase, tolerance, expect_warning):
    """A disabled check must be visible in the log after the fact.

    calphy's loggers write to their own file handler and do not propagate, so
    assert against the log file the job actually produces.
    """
    from calphy.liquid import Liquid

    scenario = "B1" if phase == "solid" else "B4"
    calc = make_calc(scenario, tolerance=tolerance)
    job, _ = recorded_job(Solid if phase == "solid" else Liquid, calc)

    for handler in job.logger.handlers:
        handler.flush()
    with open(os.path.join(str(tmp_path), "calphy.log")) as fh:
        log = fh.read()

    warned = "WARNING" in log and "DISABLED" in log
    assert warned is expect_warning, log
