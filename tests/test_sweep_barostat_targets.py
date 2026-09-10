"""Barostat and thermostat targets of the temperature-sweep fixes at finite pressure.

Reversible scaling (mode ``ts``) samples the real system at (T0/lambda, P) only
when the barostat imposes the *scaled* pressure lambda*P (de Koning, Antonelli,
Yip, J. Chem. Phys. 115, 11025 (2001), Eq. 11) -- the sweep fix has to ramp
P -> (T0/Tf)*P and the middle equilibration of the backward sweep has to sit at
(T0/Tf)*P.  ``tscale`` and the pre-flight range scan heat the *real* system and
have to hold the real pressure instead.  None of this is visible at pressure 0,
where the golden command streams are recorded, so the targets are pinned here at
50 kbar with the same RecordingRunner setup.
"""
import pytest

from calphy.solid import Solid

P0 = 50000.0  # bar
T0, TF = 400.0, 700.0  # K, from the B6/B8/B12 baseline inputs
N_SWEEP = 5000  # n_switching_steps of those inputs


def _npt_fixes(commands):
    """(index, id, Tstart, Tstop, Pstart, Pstop) of every ``fix ... npt`` command."""
    out = []
    for i, cmd in enumerate(commands):
        t = cmd.split()
        if len(t) > 3 and t[0] == "fix" and t[3] == "npt":
            k = t.index("iso")
            out.append((i, t[1], float(t[5]), float(t[6]), float(t[k + 1]), float(t[k + 2])))
    return out


def _sweep_runs(commands):
    return [i for i, cmd in enumerate(commands) if cmd == "run %d" % N_SWEEP]


def _last_fix_before(fixes, index):
    return [f for f in fixes if f[0] < index][-1]


def _job(make_calc, recorded_job, scenario, **overrides):
    calc = make_calc(scenario, pressure=P0, **overrides)
    job, rec = recorded_job(Solid, calc)
    job.lx = job.ly = job.lz = 18.075
    return job, rec


def test_ts_forward_sweep_ramps_scaled_pressure(make_calc, recorded_job):
    job, rec = _job(make_calc, recorded_job, "B6")
    job._reversible_scaling_forward(iteration=1)
    cmds = rec.commands

    (run_idx,) = _sweep_runs(cmds)
    fixes = _npt_fixes(cmds)
    sweep = _last_fix_before(fixes, run_idx)

    # Warm start and COM-constrained equilibration hold the real pressure.
    for f in fixes:
        if f[0] < sweep[0]:
            assert f[4] == f[5] == P0
    # The sweep ramps the scaled pressure lambda*P from P to (T0/Tf)*P.
    assert sweep[1] == "f1"
    assert sweep[2] == sweep[3] == T0
    assert sweep[4] == P0
    assert sweep[5] == pytest.approx(P0 * T0 / TF)
    # The COM-corrected temperature is rebound to the replaced fix.
    assert "fix_modify f1 temp tcm" in cmds[sweep[0] + 1 : run_idx]


def test_ts_backward_sweep_ramps_scaled_pressure_back(make_calc, recorded_job):
    job, rec = _job(make_calc, recorded_job, "B6")
    job._reversible_scaling_backward(iteration=1)
    cmds = rec.commands

    (run_idx,) = _sweep_runs(cmds)
    fixes = _npt_fixes(cmds)
    pf = P0 * T0 / TF

    # Middle equilibration at lambda = T0/Tf sits at the scaled pressure pf.
    middle = fixes[0]
    assert middle[4] == pytest.approx(pf)
    assert middle[5] == pytest.approx(pf)
    # The sweep ramps pf -> P.
    sweep = _last_fix_before(fixes, run_idx)
    assert sweep[4] == pytest.approx(pf)
    assert sweep[5] == P0
    assert "fix_modify f1 temp tcm" in cmds[sweep[0] + 1 : run_idx]


def test_ts_nvt_sweep_has_no_barostat(make_calc, recorded_job):
    job, rec = _job(make_calc, recorded_job, "B6", npt=False)
    job._reversible_scaling_forward(iteration=1)
    cmds = rec.commands
    (run_idx,) = _sweep_runs(cmds)
    assert _npt_fixes(cmds) == []
    # Only the warm-start thermostat is replaced; the COM-constrained nvt fix
    # stays in place for the sweep (nothing to ramp without a barostat).
    assert cmds[:run_idx].count("unfix f1") == 1
    nvt = [c for c in cmds[:run_idx] if c.startswith("fix f1 all nvt")]
    assert len(nvt) == 2 and "fixedpoint" in nvt[-1]


def test_tscale_holds_real_pressure_and_reverses_temperature(make_calc, recorded_job):
    job, rec = _job(make_calc, recorded_job, "B8")
    job.temperature_scaling(iteration=1)
    cmds = rec.commands

    fixes = _npt_fixes(cmds)
    assert fixes, "no npt fix recorded"
    for f in fixes:
        assert f[4] == f[5] == P0

    fwd_run, bwd_run = _sweep_runs(cmds)
    fwd = _last_fix_before(fixes, fwd_run)
    bwd = _last_fix_before(fixes, bwd_run)
    assert (fwd[2], fwd[3]) == (T0, TF)
    assert (bwd[2], bwd[3]) == (TF, T0)


def test_prescan_holds_real_pressure(make_calc, recorded_job):
    job, rec = _job(make_calc, recorded_job, "B12")
    job.scan_temperature_range()
    fixes = _npt_fixes(rec.commands)
    assert fixes, "no npt fix recorded"
    for f in fixes:
        assert f[4] == f[5] == P0
    ramp = [f for f in fixes if f[1] == "f2"]
    assert len(ramp) == 1 and (ramp[0][2], ramp[0][3]) == (T0, TF)
