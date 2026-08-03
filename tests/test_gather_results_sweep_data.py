"""
gather_results and the raw reversible-scaling switching data.

``ts.forward_*.dat`` / ``ts.backward_*.dat`` are written by ``fix print 1`` --
one row per MD step -- so loading them for every calculation dominates both the
runtime and the size of the returned frame.  They are therefore opt-in, and a
calculation whose sweep files are damaged (a job killed mid-sweep leaves a
ragged final row) must not take the whole gather down with it.
"""
import os

import numpy as np
import pytest
import yaml

from calphy.postprocessing import gather_results

NSTEPS = 200
SWEEP_COLS = [
    "forward_energy_diff",
    "backward_energy_diff",
    "forward_lambda",
    "backward_lambda",
]


def _make_folders(root, n_folders=3, n_replicas=1, nsteps=NSTEPS):
    """Write n_folders finished ts calculations under root."""
    lam = np.linspace(1.0, 0.68, nsteps)
    dU = -3.48 + 0.01 * np.sin(np.linspace(0, 20, nsteps))
    press = 600 * np.cos(np.linspace(0, 30, nsteps))
    vol = 6015 + np.linspace(0, 40, nsteps)
    hdr = "dU[eV/atom] press[bar] vol[A^3] lambda"

    names = []
    for k in range(n_folders):
        name = f"ts-fcc-{900 + k}-0"
        d = os.path.join(root, name)
        os.makedirs(d)
        names.append(name)
        with open(os.path.join(d, "input_file.yaml"), "w") as fh:
            yaml.safe_dump({"calculations": [{
                "mode": "ts", "temperature": 900 + k, "pressure": 0,
                "reference_phase": "solid", "phase_name": "",
                "element": ["Ag"], "reference_composition": 0.0,
            }]}, fh)
        with open(os.path.join(d, "report.yaml"), "w") as fh:
            yaml.safe_dump({
                "results": {"free_energy": -3.4, "dissipation": 1e-4},
                "input": {"element": "Ag", "concentration": "1.0"},
            }, fh)
        t = np.linspace(850, 1250, nsteps)
        np.savetxt(os.path.join(d, "temperature_sweep.dat"),
                   np.column_stack((t, dU, np.full(nsteps, 1e-5))))
        for r in range(1, n_replicas + 1):
            np.savetxt(os.path.join(d, f"ts.forward_{r}.dat"),
                       np.column_stack((dU, press, vol, lam)), header=hdr)
            np.savetxt(os.path.join(d, f"ts.backward_{r}.dat"),
                       np.column_stack((dU[::-1], press, vol, lam[::-1])),
                       header=hdr)
    return names


def _corrupt(root, name):
    """Truncate a replica mid-row, as a job killed mid-sweep would."""
    path = os.path.join(root, name, "ts.forward_1.dat")
    text = open(path).read()
    with open(path, "w") as fh:
        fh.write(text[: int(len(text) * 0.5)] + "\n-3.48 12.0\n")


def test_sweep_data_is_off_by_default(tmp_path):
    """The expensive columns exist but stay empty unless asked for."""
    root = str(tmp_path)
    _make_folders(root)
    df = gather_results(root)
    assert len(df) == 3
    for col in SWEEP_COLS:
        assert col in df.columns, f"{col} disappeared from the schema"
        assert df[col].isna().all() or all(v is None for v in df[col])


def test_default_gather_ignores_damaged_sweep_files(tmp_path):
    """
    With sweep data off, a damaged switching file is irrelevant -- it is never
    opened, so the gather is unaffected.
    """
    root = str(tmp_path)
    names = _make_folders(root)
    _corrupt(root, names[1])
    df = gather_results(root)
    assert len(df) == 3
    assert df["error_code"].isna().all() or all(
        e is None for e in df["error_code"]
    )


def test_sweep_data_loads_when_requested(tmp_path):
    root = str(tmp_path)
    _make_folders(root, n_folders=2, n_replicas=2)
    df = gather_results(root, include_sweep_data=True)
    for col in SWEEP_COLS:
        for cell in df[col]:
            assert cell is not None
            assert len(cell) == 2, "one array per replica expected"
            assert all(len(a) == NSTEPS for a in cell)
    # energy differential is dU/lambda
    fd = df["forward_energy_diff"].iloc[0][0]
    fl = df["forward_lambda"].iloc[0][0]
    assert np.all(np.isfinite(fd))
    assert fl[0] == pytest.approx(1.0)


@pytest.mark.parametrize("stride", [2, 10, 50])
def test_sweep_data_stride_reduces_samples(tmp_path, stride):
    root = str(tmp_path)
    _make_folders(root, n_folders=1)
    df = gather_results(root, include_sweep_data=True, sweep_data_stride=stride)
    expected = len(range(0, NSTEPS, stride))
    for col in SWEEP_COLS:
        assert len(df[col].iloc[0][0]) == expected


def test_damaged_sweep_file_does_not_abort_the_gather(tmp_path):
    """
    One unreadable calculation must cost only itself: every other folder still
    comes back, and the failure is recorded rather than raised.
    """
    root = str(tmp_path)
    names = _make_folders(root, n_folders=3)
    _corrupt(root, names[1])

    with pytest.warns(RuntimeWarning, match="could not read switching data"):
        df = gather_results(root, include_sweep_data=True)

    assert len(df) == 3, "a damaged folder took the whole gather down"
    df = df.set_index("calculation")
    bad = df.loc[names[1]]
    assert bad["forward_energy_diff"] is None
    assert "unreadable switching data" in str(bad["error_code"])
    # the intact folders are untouched
    for good in (names[0], names[2]):
        assert df.loc[good]["forward_energy_diff"] is not None
        assert df.loc[good]["error_code"] is None
    # free energies survive everywhere, including the damaged calculation
    assert all(len(f) == NSTEPS for f in df["free_energy"])


def test_invalid_stride_rejected(tmp_path):
    with pytest.raises(ValueError, match="sweep_data_stride"):
        gather_results(str(tmp_path), sweep_data_stride=0)
