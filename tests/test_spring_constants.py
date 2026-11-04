"""Spring-constant estimate from the MSD in Solid.

A small MSD (low T, heavy or stiff elements) used to be rounded to 0 Å² and
replaced by 1 Å², giving a k that was orders of magnitude too soft. These tests
feed a synthetic msd.dat through the estimate and the convergence loop.
"""
import logging
from types import SimpleNamespace

import numpy as np
import pytest

from calphy.integrators import kb
from calphy.solid import Solid


class FakeLmp:
    """Serves msd.dat as `rows` more samples per run command."""

    def __init__(self, columns, rows):
        self.columns = [np.asarray(c, dtype=float) for c in columns]
        self.rows = rows
        self.nruns = 0

    def command(self, cmd):
        if cmd.startswith("run"):
            self.nruns += 1

    def sync(self):
        pass

    def read_timeseries(self, name):
        assert name == "msd.dat"
        n = self.nruns * self.rows
        steps = np.arange(n, dtype=float)
        return np.column_stack([steps] + [c[:n] for c in self.columns])


def make_solid(elements, temperature, spring_constants=None, n_cycles=20,
               tolerance=0.1):
    job = Solid.__new__(Solid)
    job.logger = logging.getLogger("test_spring_constants")
    job.calc = SimpleNamespace(
        element=elements,
        n_elements=len(elements),
        _temperature=temperature,
        _qtb=False,
        spring_constants=spring_constants,
        n_equilibration_steps=25000,
        md=SimpleNamespace(
            n_every_steps=10,
            n_repeat_steps=10,
            n_small_steps=10000,
            n_cycles=n_cycles,
            thermostat_damping=[0.1, 0.1],
        ),
        tolerance=SimpleNamespace(spring_constant=tolerance),
    )
    return job


def noisy_msd(mean, n, rel_noise=0.05, seed=0):
    rng = np.random.default_rng(seed)
    return mean * (1 + rel_noise * rng.standard_normal(n))


def test_small_msd_gives_stiff_spring():
    # bcc W at 100 K: MSD ~0.004 Å² must give k = 3kT/0.004 ~ 6.46, not 3kT/1
    job = make_solid(["W"], 100.0)
    lmp = FakeLmp([np.full(300, 0.004)], rows=100)
    lmp.nruns = 3
    k_mean, _ = job.analyse_spring_constants(lmp)
    assert k_mean[0] == pytest.approx(3 * kb * 100.0 / 0.004)
    assert k_mean[0] == pytest.approx(6.46, abs=0.01)


def test_msd_is_not_rounded():
    # 0.0149 and 0.0051 both rounded to 0.01 before
    job = make_solid(["A", "B"], 300.0)
    lmp = FakeLmp([np.full(300, 0.0149), np.full(300, 0.0051)], rows=100)
    lmp.nruns = 3
    k_mean, _ = job.analyse_spring_constants(lmp)
    assert k_mean[0] == pytest.approx(3 * kb * 300.0 / 0.0149)
    assert k_mean[1] == pytest.approx(3 * kb * 300.0 / 0.0051)


def test_zero_msd_single_element_raises():
    job = make_solid(["W"], 100.0)
    lmp = FakeLmp([np.zeros(300)], rows=100)
    with pytest.raises(ValueError):
        job.run_spring_constant_convergence(lmp)
    # no usable estimate at all: stop on the first cycle instead of looping
    assert lmp.nruns == 1


def test_zero_msd_empty_element_takes_other_k():
    # an element with no atoms has MSD 0; it gets the other element's k
    job = make_solid(["Cu", "Zr"], 100.0)
    lmp = FakeLmp([noisy_msd(0.004, 2000), np.zeros(2000)], rows=100)
    job.run_spring_constant_convergence(lmp)
    assert np.isfinite(job.k).all()
    assert job.k[0] == pytest.approx(6.46, rel=0.05)
    assert job.k[1] == job.k[0]


def test_user_spring_constant_overrides_estimate():
    job = make_solid(["Cu", "Zr"], 100.0, spring_constants=[None, 2.0])
    lmp = FakeLmp([noisy_msd(0.004, 2000), np.zeros(2000)], rows=100)
    job.run_spring_constant_convergence(lmp)
    assert job.k[0] == pytest.approx(6.46, rel=0.05)
    assert job.k[1] == 2.0


def test_convergence_stops_after_two_cycles():
    # a stationary MSD converges as soon as two estimates can be compared
    job = make_solid(["W"], 100.0)
    lmp = FakeLmp([noisy_msd(0.004, 2000)], rows=100)
    job.run_spring_constant_convergence(lmp)
    assert lmp.nruns == 2
    assert job.k[0] == pytest.approx(6.46, rel=0.05)


def test_convergence_is_bounded_by_n_cycles():
    # an MSD that keeps drifting never meets the tolerance: use the last
    # estimate after n_cycles instead of iterating further or failing
    job = make_solid(["W"], 100.0, n_cycles=5)
    drifting = np.repeat(0.004 * 1.5 ** np.arange(5), 100)
    lmp = FakeLmp([drifting], rows=100)
    job.run_spring_constant_convergence(lmp)
    assert lmp.nruns == 5
    assert np.isfinite(job.k[0])
