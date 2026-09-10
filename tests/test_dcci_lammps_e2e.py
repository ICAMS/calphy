"""End-to-end ``mode: dcci`` on the real LAMMPS binary (DCCI_PLAN Part 3).

Cu (Mishin Cu01), 500 atoms, from a rough coexistence point at 0 bar to 50 kbar
in six blocks of 500 steps, out and back.  Deliberately tiny: it checks that
the whole mode runs, writes its files, produces a line with the right sign and
a slope consistent with dH/(T dV), and closes within a loose hysteresis bound.
The physics validation with converged settings lives in the example notebook.

Marked ``lammps`` and ``slow``; skipped when no ``lmp`` binary resolves.
"""
import os
import glob
import sys
import subprocess

import numpy as np
import pytest
import yaml

from calphy.runner import resolve_lammps_executable

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

try:
    BINARY = resolve_lammps_executable(None)
except ValueError:
    BINARY = None

pytestmark = [
    pytest.mark.lammps,
    pytest.mark.slow,
    pytest.mark.skipif(BINARY is None, reason="no lmp binary resolves"),
]

INPUT = {
    "calculations": [
        {
            "mode": "dcci",
            "element": "Cu",
            "mass": 63.546,
            "lattice": "FCC",
            "lattice_constant": 3.615,
            "repeat": [5, 5, 5],
            "temperature": 1340.0,
            "pressure": [0.0, 50000.0],
            "pair_style": "eam/alloy",
            "pair_coeff": "* * %s Cu" % os.path.join(HERE, "Cu01.eam.alloy"),
            "md": {"timestep": 0.001, "n_small_steps": 1000},
            "n_equilibration_steps": 2000,
            "n_switching_steps": 3000,
            "n_iterations": 1,
            "dcci": {"n_block_steps": 500},
            "queue": {"scheduler": "local", "cores": 1},
        }
    ]
}


@pytest.mark.parametrize("parallel", [False, True], ids=["sequential", "parallel_cells"])
def test_dcci_cu_out_and_back(tmp_path, parallel):
    import copy

    inp = copy.deepcopy(INPUT)
    if parallel:
        # two cells at once, one core each: exercises the threaded block loop
        inp["calculations"][0]["dcci"]["parallel_cells"] = True
        inp["calculations"][0]["queue"]["cores"] = 2
    with open(tmp_path / "input.yaml", "w") as fh:
        yaml.safe_dump(inp, fh, sort_keys=False)
    env = dict(os.environ, PYTHONPATH=REPO)
    proc = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.argv=['calphy_kernel','-i','input.yaml','-k','0']; "
         "from calphy.queuekernel import main; main()"],
        cwd=str(tmp_path), capture_output=True, text=True, timeout=1500, env=env,
    )
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-4000:]

    (folder,) = [d for d in glob.glob(str(tmp_path / "dcci-*")) if os.path.isdir(d)]
    for name in ("coexistence_line.dat", "report.yaml", "metadata.yaml",
                 "dcci.forward_1.dat", "dcci.backward_1.dat",
                 "solid/conf.dcci.forward_1.data", "liquid/conf.dcci.backward_1.data",
                 "solid/dcci_forward_1.log.lammps", "liquid/dcci_backward_1.log.lammps"):
        assert os.path.exists(os.path.join(folder, name)), name

    fwd = np.loadtxt(os.path.join(folder, "dcci.forward_1.dat"))
    assert fwd.shape[0] >= 3 and np.all(np.isfinite(fwd))
    temp, p_real, dpdt = fwd[:, 3], fwd[:, 5], fwd[:, 10]
    # copper melts at higher temperature under pressure
    assert np.all(np.diff(temp) > 0) and np.all(np.diff(p_real) > 0)
    assert p_real[-1] >= 50000.0
    # the produced line has the slope the block averages say it should
    fd = np.gradient(p_real, temp)
    assert np.allclose(dpdt[1:-1], fd[1:-1], rtol=0.3)
    # of the right order for Cu (about 40 K/GPa, i.e. ~250 bar/K)
    assert 150 < dpdt.mean() < 400

    line = np.loadtxt(os.path.join(folder, "coexistence_line.dat"))
    assert line.shape[0] == fwd.shape[0] + 1
    # the first point is the starting pressure; its temperature averages the
    # forward start (T0 exactly) with where the backward sweep came back to
    assert line[0, 0] == 0.0 and line[0, 3] == 1340.0
    assert abs(line[0, 1] - 1340.0) < 100.0
    with open(os.path.join(folder, "report.yaml")) as fh:
        report = yaml.safe_load(fh)
    # loose: six blocks of 500 steps is far from converged, but it must close
    assert abs(report["results"]["hysteresis"]) < 100.0
    assert report["results"]["n_blocks_used"] == fwd.shape[0]
