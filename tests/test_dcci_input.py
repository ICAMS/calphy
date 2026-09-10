"""Input validation and dispatch of ``mode: dcci``."""
import os

import pytest
from pydantic import ValidationError

from calphy.input import Calculation
from calphy.runner import required_styles

HERE = os.path.dirname(os.path.abspath(__file__))

BASE = dict(
    element=["Cu"], mass=[63.546], mode="dcci", temperature=1350.0,
    pressure=[0.0, 100000.0], lattice="FCC", lattice_constant=3.615,
    pair_style="eam/alloy", pair_coeff="* * %s Cu" % os.path.join(HERE, "Cu01.eam.alloy"),
)


@pytest.fixture
def build(tmp_path, monkeypatch):
    # structure generation during validation writes a data file into cwd
    monkeypatch.chdir(tmp_path)

    def _build(**extra):
        return Calculation(**{**BASE, **extra})

    return _build


def error_of(build, **extra):
    with pytest.raises(ValidationError) as ei:
        build(**extra)
    return str(ei.value)


def test_valid_dcci_input(build):
    calc = build()
    assert calc.mode == "dcci"
    assert calc._temperature == 1350.0 and calc._temperature_stop == 1350.0
    assert (calc._pressure, calc._pressure_stop) == (0.0, 100000.0)
    # the mode always builds a solid and a liquid cell; the value only names the folder
    assert calc.reference_phase == "solid"
    assert calc.create_identifier().startswith("dcci-")
    assert calc.create_identifier().endswith("-solid-1350-0")


def test_dcci_block_defaults(build):
    d = build().dcci
    assert d.n_block_steps == 1000
    assert d.integrator == "euler"
    assert d.stop_at_target_pressure is True
    assert d.n_check_blocks == 0
    assert d.hysteresis_tolerance == 5.0
    assert d.parallel_cells is False


def test_dcci_block_overrides_and_hints(build):
    calc = build(dcci={"n_block_steps": 200, "integrator": "trapezoid"})
    assert calc.dcci.n_block_steps == 200 and calc.dcci.integrator == "trapezoid"
    assert "integrator" in error_of(build, dcci={"integrator": "rk4"})
    assert "did you mean 'n_block_steps'" in error_of(build, dcci={"n_block_step": 5})
    assert "n_block_steps" in error_of(build, dcci={"n_block_steps": 0})


def test_dcci_parallel_cells_needs_two_cores(build):
    from calphy.dcci import DynamicCCI

    calc = build(dcci={"parallel_cells": True}, queue={"cores": 4})
    job = DynamicCCI(calculation=calc, simfolder=calc.create_folders())
    assert job.parallel and job.cell_cores == 2

    calc = build(dcci={"parallel_cells": True}, queue={"cores": 1}, folder_prefix="one")
    job = DynamicCCI(calculation=calc, simfolder=calc.create_folders())
    assert not job.parallel and job.cell_cores == 1

    calc = build(queue={"cores": 4}, folder_prefix="seq")
    job = DynamicCCI(calculation=calc, simfolder=calc.create_folders())
    assert not job.parallel and job.cell_cores == 4


def test_dcci_needs_pressure_range(build):
    assert "pressure: [P_start, P_stop]" in error_of(build, pressure=0.0)
    assert "pressure: [P_start, P_stop]" in error_of(build, pressure=[100000.0])


def test_dcci_needs_explicit_coexistence_temperature(build):
    assert "temperature: <T_coex>" in error_of(build, temperature=[1300.0, 1400.0])
    # the mendeleev melting-point guess is not a coexistence point of the potential
    assert "temperature: <T_coex>" in error_of(build, temperature=0)


def test_dcci_rejects_swaps_and_nvt(build):
    assert "monte_carlo swaps" in error_of(build, monte_carlo={"n_swaps": 10})
    assert "npt: True" in error_of(build, npt=False)


def test_dcci_requires_hybrid_scaled(build):
    styles = required_styles(build())
    assert "hybrid/scaled" in styles["pair"]
    assert {"npt", "ave/time"} <= styles["fix"]


def test_dcci_dispatch(build, tmp_path):
    from calphy.dcci import DynamicCCI
    from calphy.queuekernel import setup_calculation, run_calculation

    calc = build()
    job = setup_calculation(calc)
    assert isinstance(job, DynamicCCI)
    assert os.path.isdir(job.simfolder) and job.simfolder.startswith(str(tmp_path))
    assert os.path.exists(os.path.join(job.simfolder, "input_file.yaml"))
    assert (job.t0, job.p_start, job.p_stop) == (1350.0, 0.0, 100000.0)

    called = []
    job.calculate_coexistence_line = lambda: called.append(True)
    run_calculation(job)
    assert called == [True]
