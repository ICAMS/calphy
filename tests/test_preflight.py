"""Unit tests for the preflight capability check (Part 4).  No LAMMPS binary."""
import logging
import types

import pytest

import calphy.runner as R

FAKE_BIN = "/fake/lmp"

# Canned `lmp -h` style listings -----------------------------------------------
FULL_H = """\
Large-scale Atomic/Molecular Massively Parallel Simulator - dummy

* Pair styles:

eam        eam/alloy  eam/fs     ufm        hybrid/scaled
lj/cut     zero

* Bond styles:

zero

* Fix styles:

nve        nvt        npt        langevin   ave/time
print      momentum   ti/spring  atom/swap  qtb

* Compute styles:

msd        temp/com   pair       pe
"""

MINIMAL_H = """\
Large-scale Atomic/Molecular Massively Parallel Simulator - dummy

* Pair styles:

eam        eam/alloy  eam/fs     lj/cut

* Fix styles:

nve        nvt        npt        langevin   ave/time   print   momentum

* Compute styles:

msd        temp/com   pair
"""


@pytest.fixture(autouse=True)
def _clear_cache():
    R._STYLE_CACHE.clear()
    yield
    R._STYLE_CACHE.clear()


def install_fake_h(monkeypatch, stdout, counter=None):
    def _run(argv, capture_output=None, text=None, timeout=None):
        if counter is not None:
            counter.append(argv)
        return types.SimpleNamespace(returncode=0, stdout=stdout, stderr="")
    monkeypatch.setattr(R.subprocess, "run", _run)


# --------------------------------------------------------------------------- #
# Parser
# --------------------------------------------------------------------------- #
def test_parse_styles_extracts_categories():
    styles = R._parse_styles(FULL_H)
    assert {"eam/alloy", "ufm", "hybrid/scaled"} <= styles["pair"]
    assert {"nve", "ti/spring", "atom/swap", "qtb"} <= styles["fix"]
    assert {"msd", "temp/com", "pair"} <= styles["compute"]
    # tokens from other sections (bond) must not leak into our categories
    assert "zero" not in styles["fix"] and "zero" not in styles["compute"]


def test_parse_styles_empty_on_garbage():
    assert R._parse_styles("no styles here\n") == {"pair": set(), "fix": set(), "compute": set()}


# --------------------------------------------------------------------------- #
# preflight: pass / fail
# --------------------------------------------------------------------------- #
def test_solid_fe_passes_on_full_build(make_calc, monkeypatch):
    install_fake_h(monkeypatch, FULL_H)
    R.preflight(make_calc("B1"), FAKE_BIN)      # no raise


def test_solid_fe_missing_ti_spring(make_calc, monkeypatch):
    install_fake_h(monkeypatch, MINIMAL_H)
    with pytest.raises(ValueError) as ei:
        R.preflight(make_calc("B1"), FAKE_BIN)
    msg = str(ei.value)
    assert "ti/spring" in msg and "EXTRA-FIX" in msg


def test_liquid_missing_ufm(make_calc, monkeypatch):
    install_fake_h(monkeypatch, MINIMAL_H)
    with pytest.raises(ValueError) as ei:
        R.preflight(make_calc("B4"), FAKE_BIN)
    msg = str(ei.value)
    assert "ufm" in msg and "EXTRA-PAIR" in msg
    assert "hybrid/scaled" in msg


def test_qtb_missing_qtb_style(make_calc, monkeypatch):
    install_fake_h(monkeypatch, MINIMAL_H)
    calc = make_calc("B11")
    assert calc._qtb is True
    with pytest.raises(ValueError) as ei:
        R.preflight(calc, FAKE_BIN)
    assert "qtb" in str(ei.value) and "QTB" in str(ei.value)


def test_mc_swaps_needs_atom_swap(make_calc, monkeypatch):
    install_fake_h(monkeypatch, MINIMAL_H)
    calc = make_calc("B6", monte_carlo={"n_swaps": 5, "forward_swap_types": [1, 2]})
    with pytest.raises(ValueError) as ei:
        R.preflight(calc, FAKE_BIN)
    assert "atom/swap" in str(ei.value) and "MC" in str(ei.value)


def test_required_styles_ts_needs_hybrid_scaled(make_calc):
    req = R.required_styles(make_calc("B6"))      # ts / solid
    assert "hybrid/scaled" in req["pair"]
    assert "ti/spring" in req["fix"]


# --------------------------------------------------------------------------- #
# caching, skip-hatch, unreadable binary
# --------------------------------------------------------------------------- #
def test_cache_avoids_second_subprocess(monkeypatch):
    calls = []
    install_fake_h(monkeypatch, FULL_H, counter=calls)
    R.get_binary_styles(FAKE_BIN)
    R.get_binary_styles(FAKE_BIN)
    assert len(calls) == 1                          # second call served from cache


def test_skip_preflight_env(make_calc, monkeypatch):
    calls = []
    install_fake_h(monkeypatch, MINIMAL_H, counter=calls)
    monkeypatch.setenv("CALPHY_SKIP_PREFLIGHT", "1")
    R.preflight(make_calc("B1"), FAKE_BIN)          # would fail if it ran
    assert calls == []                              # short-circuited before subprocess


def test_unreadable_h_output_skips_with_warning(make_calc, monkeypatch, caplog):
    install_fake_h(monkeypatch, "")                 # e.g. a binary that crashes on -h
    with caplog.at_level(logging.WARNING, logger="calphy.runner"):
        R.preflight(make_calc("B1"), FAKE_BIN)      # needs ti/spring, but skips
    assert any("skipping" in r.message.lower() for r in caplog.records)


def test_query_styles_handles_missing_binary():
    # OSError from a non-existent binary -> empty style sets, no crash
    styles = R._query_styles("/definitely/not/a/binary")
    assert styles == {"pair": set(), "fix": set(), "compute": set()}
