"""
Schema-level tests for the Dammak quantum thermal bath surface
(mode: fe-qtb). These do not run MD; they just verify that calphy's
input layer parses QTB inputs correctly, sets the internal _qtb
flag, and rejects unsupported combinations.
"""

import pytest
import yaml
from calphy.input import read_inputfile


def _write_yaml(tmp_path, payload):
    fn = tmp_path / "input.yaml"
    fn.write_text(yaml.safe_dump(payload))
    return str(fn)


def _base_calc():
    return {
        "element": "Cu",
        "lattice": "FCC",
        "lattice_constant": 3.615,
        "mass": 63.546,
        "md": {"timestep": 0.001},
        "n_equilibration_steps": 1000,
        "n_iterations": 1,
        "n_switching_steps": 1500,
        "pair_coeff": "* * tests/Cu01.eam.alloy Cu",
        "pair_style": "eam/alloy",
        "pressure": 0.0,
        "queue": {"scheduler": "local", "cores": 1},
        "reference_phase": "solid",
        "repeat": [4, 4, 4],
        "temperature": 300.0,
    }


def test_mode_fe_qtb_sets_qtb_flag(tmp_path):
    """mode: fe-qtb sets _qtb=True and normalises mode back to 'fe'."""
    calc = _base_calc()
    calc["mode"] = "fe-qtb"
    fn = _write_yaml(tmp_path, {"calculations": [calc]})
    [opts] = read_inputfile(fn)
    assert opts._qtb is True
    assert opts.mode == "fe"  # normalised so downstream dispatch works
    # QTB defaults populated
    qtb = opts.quantum_thermal_bath
    assert qtb.thermostat_damping == pytest.approx(0.1)
    assert qtb.f_max == pytest.approx(200.0)
    assert qtb.n_f == 100


def test_mode_fe_does_not_set_qtb_flag(tmp_path):
    calc = _base_calc()
    calc["mode"] = "fe"
    fn = _write_yaml(tmp_path, {"calculations": [calc]})
    [opts] = read_inputfile(fn)
    assert opts._qtb is False
    assert opts.mode == "fe"


def test_qtb_user_parameters_round_trip(tmp_path):
    """User-set QTB knobs survive the schema and are visible downstream."""
    calc = _base_calc()
    calc["mode"] = "fe-qtb"
    calc["quantum_thermal_bath"] = {
        "thermostat_damping": 0.05,
        "barostat_damping": 0.2,
        "f_max": 30.0,
        "n_f": 50,
    }
    fn = _write_yaml(tmp_path, {"calculations": [calc]})
    [opts] = read_inputfile(fn)
    qtb = opts.quantum_thermal_bath
    assert qtb.thermostat_damping == pytest.approx(0.05)
    assert qtb.barostat_damping == pytest.approx(0.2)
    assert qtb.f_max == pytest.approx(30.0)
    assert qtb.n_f == 50


def test_fe_qtb_rejects_liquid_reference(tmp_path):
    """The Uhlenbeck-Ford liquid reference is classical; fe-qtb forbids it."""
    calc = _base_calc()
    calc["mode"] = "fe-qtb"
    calc["reference_phase"] = "liquid"
    fn = _write_yaml(tmp_path, {"calculations": [calc]})
    with pytest.raises(Exception) as exc_info:
        read_inputfile(fn)
    # Should mention fe-qtb and liquid in the message
    msg = str(exc_info.value).lower()
    assert "fe-qtb" in msg or "qtb" in msg
    assert "liquid" in msg
