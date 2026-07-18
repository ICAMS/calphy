"""Strict input validation: unknown keys are hard errors with hints.

A typo in an input key must never silently fall back to a default -- that is
a silent scientific error. Every input model forbids extra keys and the error
carries the best available hint (did-you-mean, wrong-block, removed-in-v2).
"""
import os

import pytest
import yaml
from pydantic import ValidationError

from calphy.input import Calculation, read_inputfile

BASE = dict(
    element=["Cu"], mass=[63.5], mode="fe", temperature=500, pressure=0,
    lattice="FCC", lattice_constant=3.6, reference_phase="solid",
    pair_style="eam/alloy", pair_coeff="* * Cu.eam.alloy Cu",
)


def build(**extra):
    return Calculation(**{**BASE, **extra})


def error_of(**extra):
    with pytest.raises(ValidationError) as ei:
        build(**extra)
    return str(ei.value)


# --------------------------------------------------------------------------- #
# unknown keys raise, with hints
# --------------------------------------------------------------------------- #
def test_typo_top_level_did_you_mean():
    msg = error_of(n_iteration=99)
    assert "unknown input key 'n_iteration'" in msg
    assert "did you mean 'n_iterations'?" in msg


def test_typo_in_nested_block_did_you_mean():
    msg = error_of(md={"timesteps": 0.005})
    assert "unknown input key 'timesteps' in the 'md:' block" in msg
    assert "did you mean 'timestep'?" in msg


def test_wrong_block_hint_top_level_key():
    # a real md key given at the calculation level points at the md block
    msg = error_of(timestep=0.005)
    assert "'timestep' belongs in the 'md:' block" in msg


def test_wrong_block_hint_nested_key():
    # a real top-level key inside a block points back at the calculation block
    msg = error_of(md={"element": "Cu"})
    assert "'element' belongs in the calculation block" in msg


def test_removed_key_migration_message():
    msg = error_of(savefile=True)
    assert "job-state pickling was removed in calphy v2" in msg


def test_unknown_key_without_any_hint():
    msg = error_of(totally_bogus=1)
    assert "unknown input key 'totally_bogus' in the calculation block" in msg


def test_multiple_unknown_keys_all_reported():
    msg = error_of(n_iteration=1, md={"timesteps": 0.1})
    assert "n_iteration" in msg
    # nested block errors surface separately from top-level ones
    msg2 = error_of(md={"timesteps": 0.1, "n_cycle": 5})
    assert "timesteps" in msg2 and "n_cycle" in msg2


def test_valid_input_still_passes():
    calc = build(n_iterations=3, md={"timestep": 0.002})
    assert calc.n_iterations == 3
    assert calc.md.timestep == 0.002


# --------------------------------------------------------------------------- #
# file-level checks
# --------------------------------------------------------------------------- #
def _write(tmp_path, data, name="input.yaml"):
    f = tmp_path / name
    with open(f, "w") as fh:
        yaml.safe_dump(data, fh)
    return str(f)


def _calc_dict(**extra):
    d = {k: v for k, v in BASE.items()}
    d.update(extra)
    return d


def test_read_inputfile_rejects_unknown_key(tmp_path):
    f = _write(tmp_path, {"calculations": [_calc_dict(n_iteration=5)]})
    with pytest.raises(ValidationError, match="did you mean 'n_iterations'"):
        read_inputfile(f)


def test_read_inputfile_rejects_stray_top_level_key(tmp_path):
    f = _write(tmp_path, {"calculations": [_calc_dict()], "element": ["Cu"]})
    with pytest.raises(ValueError, match="unknown top-level key"):
        read_inputfile(f)


def test_read_inputfile_legacy_format_clear_error(tmp_path):
    # pre-v1 format: top-level keys, no calculations list
    f = _write(tmp_path, {"element": ["Cu"], "mass": [63.5], "calculations_old": []})
    with pytest.raises(ValueError, match="legacy calphy input format"):
        read_inputfile(f)


def test_read_inputfile_no_calculations_key(tmp_path):
    f = _write(tmp_path, {"something": 1})
    with pytest.raises(ValueError, match="no 'calculations:' list"):
        read_inputfile(f)


def test_fast_mode_skips_validation(tmp_path):
    # validate=False (many many input calcs) must stay cheap: no pydantic
    # validation, so even a typo key passes through un-checked here. Full
    # validation always happens when the calculation actually runs.
    f = _write(tmp_path, {"calculations": [_calc_dict(n_iteration=5)]})
    calcs = read_inputfile(f, validate=False)
    assert len(calcs) == 1


def test_repo_example_inputs_are_key_clean():
    # every shipped example must parse without unknown-key errors (other
    # errors, like missing structure files or API keys, are fine here)
    import glob

    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    files = sorted(
        glob.glob(os.path.join(root, "examples", "**", "*.yaml"), recursive=True)
        + glob.glob(os.path.join(root, "docs", "source", "examples", "*", "*.yaml"))
        + glob.glob(os.path.join(here, "baselines", "inputs", "*.yaml"))
    )
    assert files, "no example inputs found"
    bad = []
    for f in files:
        try:
            read_inputfile(f, validate=True)
        except Exception as exc:
            if "unknown input key" in str(exc) or "removed in calphy" in str(exc):
                bad.append((f, str(exc)))
    assert not bad, bad
