"""md.seed: one master seed controls every stochastic step of a job.

Phase.__init__ seeds np.random from md.seed (drawing a fresh seed when the
key is unset) and backfills the value into the simfolder copy of the input
file, so every run -- seeded or not -- is reproducible after the fact.
"""
import os

import pytest
import yaml
from pydantic import ValidationError

from calphy.solid import Solid


@pytest.fixture
def solid_job(make_calc, recorded_job):
    def _build(seed="pin"):
        calc = make_calc(
            "B1", tolerance={"pressure": 1e12, "spring_constant": 1e12}
        )
        job, rec = recorded_job(Solid, calc, fixtures=["avg.dat", "msd.dat"])
        if seed == "pin":
            pass                      # recorded_job pins md.seed = 42
        else:
            # rebuild the job with a different/absent seed: recorded_job set
            # 42 on calc before Phase copied it, so make a fresh job
            calc.md.seed = seed
            job = Solid(calculation=calc, simfolder=job.simfolder)
        return job, rec

    return _build


def test_same_seed_gives_identical_command_stream(solid_job):
    job1, rec1 = solid_job()
    job1.run_averaging()
    stream1 = list(rec1.commands)

    job2, rec2 = solid_job()
    job2.run_averaging()
    assert stream1 == list(rec2.commands)


def test_different_seed_changes_the_stream(solid_job):
    job1, rec1 = solid_job()
    job1.run_averaging()
    stream1 = list(rec1.commands)

    job2, rec2 = solid_job(seed=99)
    job2.run_averaging()
    assert stream1 != list(rec2.commands)


def test_unset_seed_is_backfilled_and_recorded(solid_job):
    job, _ = solid_job(seed=None)
    assert isinstance(job.calc.md.seed, int) and job.calc.md.seed > 0
    with open(os.path.join(job.simfolder, "input_file.yaml")) as fh:
        dumped = yaml.safe_load(fh)
    assert dumped["calculations"][0]["md"]["seed"] == job.calc.md.seed


def test_unset_seed_differs_between_jobs(solid_job):
    job1, _ = solid_job(seed=None)
    job2, _ = solid_job(seed=None)
    assert job1.calc.md.seed != job2.calc.md.seed


def test_seed_does_not_mutate_the_callers_calculation(make_calc, recorded_job):
    calc = make_calc("B1", tolerance={"pressure": 1e12, "spring_constant": 1e12})
    job, _ = recorded_job(Solid, calc, fixtures=["avg.dat", "msd.dat"])
    calc.md.seed = None
    Solid(calculation=calc, simfolder=job.simfolder)
    assert calc.md.seed is None       # Phase works on a deepcopy


def test_qtb_seed_key_removed_with_hint(make_calc):
    with pytest.raises(ValidationError, match="use md.seed"):
        make_calc("B1", quantum_thermal_bath={"seed": 1234})


def test_seed_must_be_positive(make_calc):
    with pytest.raises(ValidationError):
        make_calc("B1", md={"timestep": 0.001, "seed": -3})
