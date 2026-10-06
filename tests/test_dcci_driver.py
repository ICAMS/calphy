"""The d-CCI block loop driven against a fake LAMMPS runner.

The runner records every command and answers ``read_timeseries`` with synthetic
block averages for a system whose solid and liquid differ by constant du and dv,
so the integrated line can be checked against the analytic solution and the
emitted command stream against the design (one ramping fix npt per sweep,
``run n start 0 stop N`` blocks, lambda redefined per block).
"""
import os

import numpy as np
import pytest
import yaml

import calphy.helpers as ph
from calphy.dcci import DynamicCCI, cce_slope
from calphy.input import Calculation
from calphy.integrators import EV_A3_TO_BAR
from calphy.liquid import Liquid
from calphy.solid import Solid

HERE = os.path.dirname(os.path.abspath(__file__))

T0 = 1300.0
P0, P1 = 0.0, 100000.0
N_BLOCK, N_SWEEP = 500, 2000
U_S, U_L, V_S, V_L = -3.0, -2.9, 12.0, 12.6

BASE = dict(
    element=["Cu"], mass=[63.546], mode="dcci", temperature=T0,
    pressure=[P0, P1], lattice="FCC", lattice_constant=3.615, repeat=[2, 2, 2],
    pair_style="eam/alloy", pair_coeff="* * %s Cu" % os.path.join(HERE, "Cu01.eam.alloy"),
    n_equilibration_steps=100, n_switching_steps=N_SWEEP,
    dcci={"n_block_steps": N_BLOCK, "stop_at_target_pressure": False},
)


class FakeRunner:
    """Records commands; synthesises ave/time rows from the lambda ramps it saw."""

    def __init__(self, directory):
        self.directory = directory
        self.commands = []
        self.blocks = []          # (lam_start, lam_end) per block run
        self._pending = None
        self.closed = False

    def command(self, s):
        cmd = " ".join(s.split())
        self.commands.append(cmd)
        t = cmd.split()
        if t[:3] == ["variable", "lam", "equal"] and "step" in t[3]:
            a, rest = t[3].split("+(", 1)
            b = rest.split(")")[0]
            self._pending = (float(a), float(a) + float(b))
        elif t[0] == "run" and self._pending is not None and "start" in t:
            self.blocks.append(self._pending)

    def sync(self):
        pass

    def close(self):
        self.closed = True

    def rotate_logs(self, name):
        pass

    def read_timeseries(self, name, usecols=None):
        phase = os.path.basename(self.directory)
        u, v = (U_S, V_S) if phase == "solid" else (U_L, V_L)
        rows = []
        for k, (l0, l1) in enumerate(self.blocks):
            rows.append([(k + 1) * N_BLOCK, u, v, 0.5 * (l0 + l1)])
        return np.array(rows)


@pytest.fixture
def driver(tmp_path, monkeypatch, request):
    monkeypatch.chdir(tmp_path)
    runners = {}
    overrides = getattr(request, "param", {})

    def fake_create_object(calc, directory):
        r = FakeRunner(directory)
        runners.setdefault(os.path.basename(directory), []).append(r)
        return r

    def fake_averaging(self):
        self.lx = self.ly = self.lz = 14.46
        self.volatom = self.lx**3 / self.natoms
        open(os.path.join(self.simfolder, "conf.equilibration.data"), "w").close()

    def fake_solid_fraction(path):
        # the solid stays solid, the liquid stays liquid
        return 10**6 if "/solid/" in path else 0

    monkeypatch.setattr(ph, "create_object", fake_create_object)
    monkeypatch.setattr(Solid, "run_averaging", fake_averaging)
    monkeypatch.setattr(Liquid, "run_averaging", fake_averaging)
    monkeypatch.setattr(ph, "find_solid_fraction", fake_solid_fraction)

    data = {**BASE, **overrides}
    if "dcci" in overrides:
        data["dcci"] = {**BASE["dcci"], **overrides["dcci"]}
    calc = Calculation(**data)
    simfolder = calc.create_folders()
    job = DynamicCCI(calculation=calc, simfolder=simfolder)
    job.runners = runners
    return job


def test_forward_sweep_integrates_the_analytic_line(driver):
    job = driver
    job.prepare_cells()
    job.equilibrate_cells()
    state = job.run_sweep("forward", iteration=1)

    n_blocks = N_SWEEP // N_BLOCK
    assert state.n_blocks == n_blocks
    assert state.p_rs == pytest.approx(P1)

    data = np.loadtxt(os.path.join(job.simfolder, "dcci.forward_1.dat"))
    assert data.shape == (n_blocks, 11)
    lam, temp, p_rs, p_real = data[:, 2], data[:, 3], data[:, 4], data[:, 5]
    # constant du, dv: lambda is linear in P_RS with the analytic slope
    f = cce_slope(U_S, U_L, V_S, V_L)
    h = (P1 - P0) / n_blocks / EV_A3_TO_BAR
    assert lam == pytest.approx(1.0 + f * h * np.arange(1, n_blocks + 1))
    assert temp == pytest.approx(T0 / lam)
    assert p_real == pytest.approx(p_rs / lam)
    assert np.all(np.diff(temp) > 0) and np.all(np.diff(p_real) > 0)
    # the diagnostic slope of the produced line matches the finite-difference one
    dpdt = data[:, 10]
    fd = np.gradient(p_real, temp)
    assert np.allclose(dpdt[1:-1], fd[1:-1], rtol=0.05)
    assert state.lam == pytest.approx(lam[-1])

    for phase in ("solid", "liquid"):
        (run,) = job.runners[phase]
        cmds = run.commands
        npt = [c for c in cmds if " npt " in c]
        assert len(npt) == 2, npt
        # equilibration at constant P_RS, then a single ramp for the whole sweep
        assert npt[0].split()[-3:-1] == ["%f" % P0, "%f" % P0]
        assert npt[1].split()[-3:-1] == ["%f" % P0, "%f" % P1]
        blocks = [c for c in cmds if c.startswith("run ") and "start 0 stop %d" % N_SWEEP in c]
        assert len(blocks) == n_blocks
        assert cmds.count("reset_timestep 0") == 1
        assert any(c.startswith("pair_style hybrid/scaled v_lam") for c in cmds)
        assert any(c.startswith("fix fav all ave/time 1 %d %d v_u_real v_vol_atom v_lam" % (N_BLOCK, N_BLOCK)) for c in cmds)
        assert any(c.startswith("write_data") and "conf.dcci.forward_1.data" in c for c in cmds)
        assert run.closed
        # the lambda ramp of each block starts where the previous corrected value ended
        starts = [b[0] for b in run.blocks]
        assert starts == pytest.approx([1.0] + list(lam[:-1]))


def test_forward_sweep_stops_at_target_pressure(driver):
    job = driver
    job.calc.dcci.stop_at_target_pressure = True
    job.prepare_cells()
    job.equilibrate_cells()
    # the ramp is applied to the scaled pressure; the real pressure P_RS/lambda
    # runs ahead of it and reaches pressure[1] before the ramp does
    state = job.run_sweep("forward", iteration=1)
    data = np.atleast_2d(np.loadtxt(os.path.join(job.simfolder, "dcci.forward_1.dat")))
    assert state.n_blocks == data.shape[0] == 3 < N_SWEEP // N_BLOCK
    assert data[-1, 5] >= P1 > data[-2, 5]
    assert state.p_rs == pytest.approx(data[-1, 4]) and state.p_rs < P1
    (run,) = job.runners["solid"]
    assert len([c for c in run.commands if c.startswith("run ") and "start 0" in c]) == 3


def test_trapezoid_corrector(driver):
    job = driver
    job.calc.dcci.integrator = "trapezoid"
    job.prepare_cells()
    job.equilibrate_cells()
    job.run_sweep("forward", iteration=1)
    data = np.loadtxt(os.path.join(job.simfolder, "dcci.forward_1.dat"))
    # constant slope: trapezoid and euler coincide
    f = cce_slope(U_S, U_L, V_S, V_L)
    h = (P1 - P0) / (N_SWEEP // N_BLOCK) / EV_A3_TO_BAR
    assert data[:, 2] == pytest.approx(1.0 + f * h * np.arange(1, data.shape[0] + 1))


def test_cells_are_built_at_the_starting_point(driver):
    job = driver
    job.prepare_cells()
    for phase in ("solid", "liquid"):
        c = job.jobs[phase].calc
        assert c.reference_phase == phase
        assert (c._pressure, c._pressure_stop) == (P0, P0)
        assert c._temperature == T0
        assert c.tolerance.solid_fraction == 0.7 and c.tolerance.liquid_fraction == 0.05
        assert os.path.isdir(os.path.join(job.simfolder, phase))
    with open(os.path.join(job.simfolder, "input_file.yaml")) as fh:
        assert yaml.safe_load(fh)["calculations"][0]["mode"] == "dcci"


def test_full_cycle_closes_and_reports(driver):
    job = driver
    job.calculate_coexistence_line()

    n_blocks = N_SWEEP // N_BLOCK
    fwd = np.loadtxt(os.path.join(job.simfolder, "dcci.forward_1.dat"))
    bwd = np.loadtxt(os.path.join(job.simfolder, "dcci.backward_1.dat"))
    assert fwd.shape == bwd.shape == (n_blocks, 11)
    # constant du, dv: the backward sweep retraces the forward one exactly
    assert bwd[-1, 2] == pytest.approx(1.0)
    assert bwd[-1, 4] == pytest.approx(P0) and bwd[-1, 5] == pytest.approx(P0)
    assert bwd[:, 4][::-1] == pytest.approx(np.concatenate(([P0], fwd[:-1, 4])))

    line = np.loadtxt(os.path.join(job.simfolder, "coexistence_line.dat"))
    assert line.shape == (n_blocks + 1, 5)
    assert line[0, :2] == pytest.approx([P0, T0])
    assert line[:, 1] == pytest.approx(np.concatenate(([T0], fwd[:, 3])))
    assert np.all(line[:, 2] < 1e-6)
    assert job.hysteresis == pytest.approx(0.0, abs=1e-6) and not job.hysteresis_high

    with open(os.path.join(job.simfolder, "report.yaml")) as fh:
        report = yaml.safe_load(fh)
    assert report["results"]["hysteresis_high"] is False
    assert report["results"]["pressure_reached"] == pytest.approx(line[-1, 0])
    assert report["results"]["temperature_at_pressure_reached"] == pytest.approx(line[-1, 1])
    assert report["results"]["n_blocks_used"] == n_blocks
    assert report["input"]["pressure_stop"] == P1
    assert os.path.exists(os.path.join(job.simfolder, "metadata.yaml"))
    assert os.path.exists(os.path.join(job.simfolder, "input_configuration.data"))

    # backward sessions start from the forward end configuration at its lambda
    for phase in ("solid", "liquid"):
        fwd_run, bwd_run = job.runners[phase]
        assert any("conf.dcci.forward_1.data" in c and c.startswith("read_data") for c in bwd_run.commands)
        first_lam = [c for c in bwd_run.commands if c.startswith("variable lam equal")][0]
        assert float(first_lam.split()[3]) == pytest.approx(fwd[-1, 2])
        npt = [c for c in bwd_run.commands if " npt " in c]
        assert npt[1].split()[-3:-1] == ["%f" % P1, "%f" % P0]
        assert not any("remap" in c for c in bwd_run.commands)
        assert any("velocity" in c for c in fwd_run.commands)
        assert not any("velocity" in c for c in bwd_run.commands)


def test_integrate_dcci_hysteresis_and_error(tmp_path):
    from calphy.dcci import integrate_dcci

    t0, p0 = 1000.0, 0.0
    cols = 11
    p_f = np.array([10000.0, 20000.0, 30000.0])
    t_f = np.array([1010.0, 1020.0, 1030.0])
    offset = 4.0
    for i in (1, 2):
        fwd = np.zeros((3, cols)); fwd[:, 3] = t_f; fwd[:, 5] = p_f
        # the backward sweep runs 2 K (iteration 1) / 6 K (iteration 2) hot
        d = offset - 2.0 if i == 1 else offset + 2.0
        # backward rows at 20, 10, 0 kbar: the forward line there, shifted by d
        bwd = np.zeros((3, cols)); bwd[:, 5] = p_f[::-1] - 10000.0
        bwd[:, 3] = np.interp(bwd[:, 5], [p0, *p_f], [t0, *t_f]) + d
        np.savetxt(tmp_path / ("dcci.forward_%d.dat" % i), fwd)
        np.savetxt(tmp_path / ("dcci.backward_%d.dat" % i), bwd)

    hyst, (p, t, err, tf, tb) = integrate_dcci(str(tmp_path), t0, p0, nsims=2, return_values=True)
    assert hyst == pytest.approx(offset)
    assert p == pytest.approx([p0, *p_f])
    assert tf == pytest.approx([t0, *t_f])
    # the backward sweep starts exactly where the forward one ended
    shift = np.array([offset, offset, offset, 0.0])
    assert tb == pytest.approx(tf + shift)
    assert t == pytest.approx(tf + shift / 2)
    # half the hysteresis, plus the spread of the two iterations' means
    stderr = np.array([np.std([s / 2 - 1.0, s / 2 + 1.0], ddof=1) / np.sqrt(2) if s else 0.0 for s in shift])
    assert err == pytest.approx(np.sqrt((shift / 2) ** 2 + stderr**2))
    assert os.path.exists(tmp_path / "coexistence_line.dat")


@pytest.mark.parametrize(
    "driver", [{"dcci": {"parallel_cells": True}, "queue": {"cores": 4}}], indirect=True
)
def test_parallel_cells_give_the_same_line(driver):
    """The cells only talk at block boundaries, so running them concurrently
    must reproduce the sequential result exactly, with half the cores each."""
    job = driver
    assert job.parallel and job.cell_cores == 2
    job.calculate_coexistence_line()
    for phase in ("solid", "liquid"):
        assert job.jobs[phase].calc.queue.cores == 2
        # both sweeps happened in each cell, run as one block per exchange
        fwd_run, bwd_run = job.runners[phase]
        assert len([c for c in fwd_run.commands if c.startswith("run ") and "start 0" in c]) == N_SWEEP // N_BLOCK

    n_blocks = N_SWEEP // N_BLOCK
    fwd = np.loadtxt(os.path.join(job.simfolder, "dcci.forward_1.dat"))
    f = cce_slope(U_S, U_L, V_S, V_L)
    h = (P1 - P0) / n_blocks / EV_A3_TO_BAR
    assert fwd[:, 2] == pytest.approx(1.0 + f * h * np.arange(1, n_blocks + 1))
    assert job.hysteresis == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize(
    "driver", [{"dcci": {"parallel_cells": True}, "queue": {"cores": 1}}], indirect=True
)
def test_parallel_cells_falls_back_to_sequential_on_one_core(driver):
    job = driver
    assert not job.parallel and job.cell_cores == 1
    job.prepare_cells()
    assert all(job.jobs[p].calc.queue.cores == 1 for p in ("solid", "liquid"))
    with open(os.path.join(job.simfolder, "calphy.log")) as fh:
        assert "parallel_cells needs queue.cores >= 2" in fh.read()
