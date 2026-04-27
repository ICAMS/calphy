import pytest
from calphy.input import read_inputfile
import os


def test_pressure_coupling_default_iso():
    """Scalar pressure uses 'iso' barostat by default."""
    calculations = read_inputfile(os.path.join(os.getcwd(), "tests/input.yaml"))
    calc = calculations[0]
    assert calc._pressure_coupling == "iso"


def test_pressure_coupling_explicit_tri():
    """Explicit pressure_coupling: tri is accepted and stored."""
    calculations = read_inputfile(os.path.join(os.getcwd(), "tests/inp_tri.yaml"))
    calc = calculations[0]
    assert calc._pressure_coupling == "tri"
    assert calc.pressure_coupling == "tri"


def test_pressure_coupling_invalid_raises():
    """An unrecognised pressure_coupling value raises a validation error."""
    from calphy.input import Calculation

    with pytest.raises(Exception):
        Calculation(
            element=["Al"],
            mass=[26.98],
            reference_phase="solid",
            pair_style=["eam/alloy"],
            pair_coeff=["* * tests/Al99.eam.alloy Al"],
            lattice="FCC",
            lattice_constant=4.05,
            temperature=800.0,
            pressure_coupling="bad_value",
        )
