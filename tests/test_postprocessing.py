"""Tests for calphy.postprocessing.gather_results.

Builds minimal calculation folders by hand (input_file.yaml + report.yaml,
plus temperature_sweep.dat / ts.forward_i.dat / ts.backward_i.dat for the
ts mode) rather than running real calphy jobs, mirroring the on-disk shape
Calculation.dump() produces -- see calphy/input.py for the phase_name
(default "") and reference_composition (default 0.0) field defaults.
"""
import os

import numpy as np
import yaml

from calphy.postprocessing import gather_results

INPUT_TS = {
    "calculations": [
        {
            "mode": "ts",
            "temperature": [500, 600],
            "pressure": 0,
            "reference_phase": "solid",
            "phase_name": "",
            "reference_composition": 0.0,
        }
    ]
}

INPUT_FE = {
    "calculations": [
        {
            "mode": "fe",
            "temperature": 500,
            "pressure": 0,
            "reference_phase": "solid",
            "phase_name": "",
            "reference_composition": 0.0,
        }
    ]
}


def _write_yaml(path, data):
    with open(path, "w") as fh:
        yaml.dump(data, fh)


def _make_ts_calc(folder, with_replicas=True):
    # "replicas" = independent forward/backward switching runs (nsims in
    # calphy/integrators.py), each written as a ts.forward_i.dat/ts.backward_i.dat
    # pair; with_replicas=False omits them to test the no-file fallback.
    os.makedirs(folder, exist_ok=True)
    _write_yaml(os.path.join(folder, "input_file.yaml"), INPUT_TS)
    _write_yaml(
        os.path.join(folder, "report.yaml"),
        {
            "input": {"element": "Al Cu", "concentration": "0.5 0.5"},
            "results": {"free_energy": -3.5, "ts_dissipation": 0.002},
        },
    )
    t = np.array([500.0, 600.0])
    f = np.array([-3.5, -3.6])
    ferr = np.array([0.001, 0.0012])
    np.savetxt(os.path.join(folder, "temperature_sweep.dat"), np.column_stack((t, f, ferr)))

    if with_replicas:
        lam = np.linspace(1, 0.5, 5)
        dx = np.arange(5, dtype=float) + 1.0
        p = np.zeros(5)
        vol = np.ones(5)
        np.savetxt(
            os.path.join(folder, "ts.forward_1.dat"),
            np.column_stack((dx, p, vol, lam)),
            header="x",
        )
        np.savetxt(
            os.path.join(folder, "ts.backward_1.dat"),
            np.column_stack((dx, p, vol, lam)),
            header="x",
        )
    return t, f, ferr, lam if with_replicas else None


def _make_fe_calc(folder, with_report=True):
    os.makedirs(folder, exist_ok=True)
    _write_yaml(os.path.join(folder, "input_file.yaml"), INPUT_FE)
    if with_report:
        _write_yaml(
            os.path.join(folder, "report.yaml"),
            {
                "input": {"element": "Al Cu", "concentration": "0.5 0.5"},
                "results": {"free_energy": -3.4, "dissipation": 0.0005},
            },
        )


def test_gather_results_ts_mode_with_replicas(tmp_path):
    """ts mode picks up free_energy_error, ts_dissipation, and -- when the
    switching data is asked for -- per-replica forward/backward energy-diff +
    lambda arrays."""
    mainfolder = tmp_path / "calcs"
    t, f, ferr, lam = _make_ts_calc(str(mainfolder / "run1"))

    df = gather_results(str(mainfolder), include_sweep_data=True)
    row = df.loc[df.calculation == "run1"].iloc[0]

    assert row.status == "True"
    assert row.calculation_mode == "ts"
    np.testing.assert_allclose(row.free_energy_error, ferr)
    assert np.isnan(row.dissipation)
    assert row.ts_dissipation == 0.002

    assert row.forward_energy_diff is not None
    assert row.backward_energy_diff is not None
    assert len(row.forward_energy_diff) == 1
    dx = np.arange(5, dtype=float) + 1.0
    np.testing.assert_allclose(row.forward_energy_diff[0], dx / lam)
    np.testing.assert_allclose(row.forward_lambda[0], lam)
    np.testing.assert_allclose(row.backward_lambda[0], lam)


def test_gather_results_ts_mode_without_replicas(tmp_path):
    """Without ts.forward_i.dat/ts.backward_i.dat files, the diff/lambda
    columns fall back to None but free_energy_error still parses."""
    mainfolder = tmp_path / "calcs"
    t, f, ferr, _ = _make_ts_calc(str(mainfolder / "run1"), with_replicas=False)

    df = gather_results(str(mainfolder))
    row = df.loc[df.calculation == "run1"].iloc[0]

    assert row.status == "True"
    np.testing.assert_allclose(row.free_energy_error, ferr)
    assert row.forward_energy_diff is None
    assert row.backward_energy_diff is None
    assert row.forward_lambda is None
    assert row.backward_lambda is None


def test_gather_results_fe_mode(tmp_path):
    """fe mode has no temperature sweep: free_energy_error is fixed at 0.0
    and dissipation is read straight from report.yaml; no diff arrays."""
    mainfolder = tmp_path / "calcs"
    _make_fe_calc(str(mainfolder / "run2"))

    df = gather_results(str(mainfolder))
    row = df.loc[df.calculation == "run2"].iloc[0]

    assert row.status == "True"
    assert row.calculation_mode == "fe"
    assert row.free_energy == -3.4
    assert row.free_energy_error == 0.0
    assert row.dissipation == 0.0005
    assert np.isnan(row.ts_dissipation)
    assert row.forward_energy_diff is None
    assert row.backward_energy_diff is None


def test_gather_results_missing_report(tmp_path):
    """A calculation folder with no report.yaml (e.g. still running or
    failed) reports status False and leaves every new column at its
    NaN/None default -- no regression vs. the pre-existing columns."""
    mainfolder = tmp_path / "calcs"
    _make_fe_calc(str(mainfolder / "run3"), with_report=False)

    df = gather_results(str(mainfolder))
    row = df.loc[df.calculation == "run3"].iloc[0]

    assert row.status == "False"
    assert np.isnan(row.free_energy)
    assert np.isnan(row.free_energy_error)
    assert np.isnan(row.dissipation)
    assert np.isnan(row.ts_dissipation)
    assert row.forward_energy_diff is None
    assert row.backward_energy_diff is None
