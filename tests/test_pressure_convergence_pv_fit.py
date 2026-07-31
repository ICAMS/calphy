"""
The P-V fit / box-rescale path of ``Phase.run_pressure_convergence``.

Samples recorded before a ``change_box`` must never leak into the mean
pressure, the stored box, or the P-V fit, and the scale factor must be referred
to the volume the box actually has when ``change_box`` runs.
"""
import numpy as np
import pytest

from calphy.phase import Phase
from calphy.solid import Solid


def test_fit_volume_scale_uses_current_volume():
    """
    change_box scales the box as it stands now, so the factor must be referred
    to the current volume rather than the windowed mean that fed the fit.
    """
    # exactly linear P(V): P = 0 at V = 16.0
    v0 = 16.0
    pv_history = [(v0 - 0.4, 400.0), (v0 - 0.3, 300.0), (v0 - 0.2, 200.0)]

    # referred to the mean-volume history entry (the old behaviour)
    default = Phase._fit_volume_scale(pv_history, 0.0)
    assert default == pytest.approx((v0 / (v0 - 0.2)) ** (1 / 3), rel=1e-9)

    # referred to where the box really is now
    current = 15.5
    explicit = Phase._fit_volume_scale(pv_history, 0.0, current)
    assert explicit == pytest.approx((v0 / current) ** (1 / 3), rel=1e-9)
    assert explicit != pytest.approx(default, rel=1e-6)


class _PVStub:
    """
    Minimal runner stub with an exactly-linear equation of state.

    ``run`` appends one cycle of samples at the current box and relaxes the box
    slightly toward equilibrium (so the P-V history has distinct volumes, as a
    real barostat would); ``change_box`` scales it. No noise, so the fit is
    exact and any contamination from pre-rescale samples shows up as a clean
    numerical discrepancy.
    """

    # ``relax`` is deliberately far slower than the sampling window: the box
    # must still be well away from equilibrium when the P-V fit first fires,
    # so that the rescale -- not the barostat -- is what closes the gap. That
    # is the regime in which stale pre-rescale samples actually corrupt the
    # averages.
    def __init__(self, natoms, ncount, v0=16.0, bulk=1.0e6, start_frac=0.97,
                 relax=0.002):
        self.natoms = natoms
        self.ncount = ncount
        self.v0 = v0
        self.bulk = bulk
        self.relax = relax
        self.volatom = v0 * start_frac
        self.rows = []
        self.commands = []
        self.n_change_box = 0

    def _pressure(self):
        return -self.bulk * (self.volatom - self.v0) / self.v0

    def command(self, cmd):
        self.commands.append(" ".join(cmd.split()))
        tokens = cmd.split()
        if tokens[0] == "run":
            for _ in range(self.ncount):
                # barostat nudges the box toward equilibrium each sample
                self.volatom += self.relax * (self.v0 - self.volatom)
                length = (self.volatom * self.natoms) ** (1.0 / 3.0)
                self.rows.append([len(self.rows), length, length, length,
                                  self._pressure()])
        elif tokens[0] == "change_box":
            scale = float(tokens[tokens.index("scale") + 1])
            self.volatom *= scale ** 3
            self.n_change_box += 1

    def sync(self):
        pass

    def close(self):
        pass

    def rotate_logs(self, stage):
        pass

    def read_timeseries(self, name, usecols=None):
        arr = np.array(self.rows, dtype=float)
        if usecols is not None:
            arr = arr[:, list(usecols)]
        return arr


def test_pressure_convergence_discards_pre_rescale_samples(make_calc, recorded_job):
    """
    After a change_box the stored box must describe the *new* cell only.

    The stub's equation of state puts zero pressure exactly at ``v0``, so a
    correct run converges with the stored vol/atom on ``v0``. If samples taken
    before the rescale still counted, the stored box would sit between the old
    and new volumes instead.
    """
    calc = make_calc("B1")
    # tight tolerance so the loop has to reach equilibrium rather than stop on
    # the first cycle, which is what brings the P-V fit into play at all
    calc.tolerance.pressure = 1.0
    job, _ = recorded_job(Solid, calc)

    ncount = int(calc.md.n_small_steps) // int(
        calc.md.n_every_steps * calc.md.n_repeat_steps
    )
    stub = _PVStub(natoms=job.natoms, ncount=ncount)
    job.run_pressure_convergence(stub)

    assert stub.n_change_box > 0, "P-V fit path was never exercised"
    assert job.volatom == pytest.approx(stub.v0, rel=2e-4), (
        "stored vol/atom %.6f is not the equilibrium volume %.6f -- "
        "pre-rescale samples leaked into the average"
        % (job.volatom, stub.v0)
    )
    # and the stored box is self-consistent with the stored volume
    assert job.lx * job.ly * job.lz / job.natoms == pytest.approx(
        job.volatom, rel=2e-3
    )


def test_pressure_convergence_skips_transient_cycle_after_rescale(
    make_calc, recorded_job
):
    """The cycle straight after a rescale is pure transient and must not
    declare convergence."""
    calc = make_calc("B1")
    calc.tolerance.pressure = 1.0
    job, _ = recorded_job(Solid, calc)

    ncount = int(calc.md.n_small_steps) // int(
        calc.md.n_every_steps * calc.md.n_repeat_steps
    )
    stub = _PVStub(natoms=job.natoms, ncount=ncount)
    job.run_pressure_convergence(stub)

    runs = [c for c in stub.commands if c.startswith("run")]
    changes = [i for i, c in enumerate(stub.commands) if c.startswith("change_box")]
    # every change_box is followed by at least one more run before the loop ends
    for idx in changes:
        assert any(c.startswith("run") for c in stub.commands[idx + 1:]), (
            "a rescale was the last thing the loop did; the new box was never "
            "sampled before the stored value was taken"
        )
    assert len(runs) > len(changes)
