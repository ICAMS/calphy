"""A disabled phase-stability check must not run the structure analysis.

`tolerance.solid_fraction: 0` (or `liquid_fraction: 1`) makes the comparison in
check_if_melted / check_if_solidfied impossible to trigger, but the solid fraction
was still computed with pyscal3 -- which segfaults in find.neighbors on some LAMMPS
triclinic dumps (Mg2Si antifluorite, 60-degree primitive cell), killing the run.
"""
from types import SimpleNamespace

import pytest

import calphy.phase as phase_mod
from calphy.phase import Phase


def _fake(solid_fraction=0.0, liquid_fraction=1.0):
    tol = SimpleNamespace(solid_fraction=solid_fraction, liquid_fraction=liquid_fraction)
    return SimpleNamespace(calc=SimpleNamespace(tolerance=tol), simfolder="/nonexistent", natoms=10)


def _explode(_):
    raise AssertionError("find_solid_fraction must not be called when the check is disabled")


def test_disabled_melt_check_skips_structure_analysis(monkeypatch):
    monkeypatch.setattr(phase_mod.ph, "find_solid_fraction", _explode)
    Phase.check_if_melted(_fake(solid_fraction=0.0), lmp=None, filename="traj.temp.dat")


def test_disabled_solidification_check_skips_structure_analysis(monkeypatch):
    monkeypatch.setattr(phase_mod.ph, "find_solid_fraction", _explode)
    Phase.check_if_solidfied(_fake(liquid_fraction=1.0), lmp=None, filename="traj.temp.dat")


def test_enabled_melt_check_still_runs(monkeypatch):
    called = []
    monkeypatch.setattr(phase_mod.ph, "find_solid_fraction", lambda f: called.append(f) or 10)
    Phase.check_if_melted(_fake(solid_fraction=0.7), lmp=None, filename="traj.temp.dat")
    assert called == ["/nonexistent/traj.temp.dat"]


def test_enabled_solidification_check_still_runs(monkeypatch):
    called = []
    monkeypatch.setattr(phase_mod.ph, "find_solid_fraction", lambda f: called.append(f) or 0)
    Phase.check_if_solidfied(_fake(liquid_fraction=0.05), lmp=None, filename="traj.temp.dat")
    assert called == ["/nonexistent/traj.temp.dat"]
