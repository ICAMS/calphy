"""set_composition_scaling_potential: the pair lists of a composition_scaling
calculation become the initial and final potential of the alchemical switch,
component by component."""
import os

import pytest

from calphy.input import Calculation
from calphy.composition_transformation import CompositionTransformation
from calphy.routines import set_composition_scaling_potential

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ZRCU = os.path.join(REPO, "examples", "example_10", "ZrCu.data")
POT = os.path.join(REPO, "examples", "potentials", "ZrCu.eam.fs")

BASE = dict(
    element=["Zr", "Cu"], mass=[91.224, 63.546], lattice=ZRCU, lattice_constant=0.0,
    mode="composition_scaling",
    composition_scaling={"output_chemical_composition": {"Cu": 532, "Zr": 492}},
    temperature=800, pressure=0, reference_phase="solid",
)


@pytest.fixture
def in_tmp(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)


def _prepare(**extra):
    calc = Calculation(**{**BASE, **extra})
    comp = CompositionTransformation(calc)
    n = set_composition_scaling_potential(calc, comp)
    return calc, comp, n


def test_plain_potential_becomes_two_copies(in_tmp):
    calc, comp, n = _prepare(pair_style="eam/fs", pair_coeff="* * %s Zr Cu" % POT)
    assert n == 1
    assert calc.pair_style == ["eam/fs", "eam/fs"]
    assert calc._pair_style_names == ["eam/fs", "eam/fs"]
    old, new = comp.update_pair_coeff("* * %s Zr Cu" % POT)
    assert calc.pair_coeff == [old, new]
    # the elements follow the per-type mappings (Cu-Cu, Zr-Cu, Zr-Zr): the
    # initial mapping keeps the 20 transformed atoms Zr, the final names them Cu
    assert old.split()[3:] == comp.pair_list_old == ["Cu", "Zr", "Zr"]
    assert new.split()[3:] == comp.pair_list_new == ["Cu", "Cu", "Zr"]


def test_overlay_components_are_rewritten_pairwise(in_tmp):
    calc, comp, n = _prepare(
        pair_mode="overlay",
        pair_style=["eam/fs", "zero 5.0"],
        pair_coeff=["* * eam/fs %s Zr Cu" % POT, "* * zero"],
    )
    assert n == 2
    assert calc.pair_style == ["eam/fs", "zero 5.0", "eam/fs", "zero 5.0"]
    assert calc._pair_style_with_options == calc.pair_style
    assert calc._pair_style_names == ["eam/fs", "zero", "eam/fs", "zero"]
    assert calc.pair_coeff == [
        "* * eam/fs %s Cu Zr Zr" % POT,
        "* * zero",
        "* * eam/fs %s Cu Cu Zr" % POT,
        "* * zero",
    ]


def test_extra_entries_of_a_plain_potential_are_dropped(in_tmp):
    """Without pair_mode overlay only the first pair_style/pair_coeff describe
    the potential; a stray second entry must not leak into the end states."""
    calc, comp, n = _prepare(
        pair_style=["eam/fs", "eam/alloy"],
        pair_coeff=["* * %s Zr Cu" % POT, "* * other.eam.alloy Zr Cu"],
    )
    assert n == 1
    assert calc.pair_style == ["eam/fs", "eam/fs"]
    assert len(calc.pair_coeff) == 2
