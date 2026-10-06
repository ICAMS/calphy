"""Tests for the ASE-backed built-in lattices in calphy.structures."""

import numpy as np
import pytest

from calphy.structures import (
    BUILTIN_LATTICES,
    IDEAL_C_OVER_A,
    canonical_lattice,
    default_lattice,
    make_lattice,
)


@pytest.mark.parametrize(
    "name, atoms_per_cell",
    [("simple_cubic", 1), ("bcc", 2), ("fcc", 4), ("diamond", 8), ("a15", 8)],
)
def test_cubic_lattices(name, atoms_per_cell):
    atoms = make_lattice(name, "Cu", 3.6, [2, 3, 4])
    assert len(atoms) == atoms_per_cell * 24
    assert np.allclose(atoms.cell, np.diag([7.2, 10.8, 14.4]))
    assert set(atoms.get_chemical_symbols()) == {"Cu"}
    assert all(atoms.pbc)


def test_fcc_nearest_neighbor_distance():
    atoms = make_lattice("fcc", "Cu", 3.6, [3, 3, 3])
    d = atoms.get_distances(0, range(1, len(atoms)), mic=True)
    assert np.isclose(d.min(), 3.6 / np.sqrt(2))
    assert np.sum(np.isclose(d, d.min())) == 12


def test_hcp_orthogonal_cell_ideal_c_over_a():
    atoms = make_lattice("hcp", "Mg", 3.2, [3, 3, 3])
    assert len(atoms) == 4 * 27
    assert np.allclose(atoms.cell.angles(), 90)
    assert np.allclose(
        atoms.cell.lengths(), [9.6, 9.6 * np.sqrt(3), 9.6 * IDEAL_C_OVER_A]
    )
    d = atoms.get_distances(0, range(1, len(atoms)), mic=True)
    # ideal hcp: 12 nearest neighbours at a
    assert np.sum(np.isclose(d, 3.2)) == 12


def test_hcp_c_over_a():
    atoms = make_lattice("hcp", "Zn", 2.66, [1, 1, 1], c_over_a=1.856)
    assert np.isclose(atoms.cell[2, 2], 2.66 * 1.856)


def test_a15_neighbors():
    # beta-W: 2a sites have 12 neighbours at a*sqrt(5)/4, 6c sites form
    # chains with 2 neighbours at a/2
    atoms = make_lattice("a15", "W", 5.05, [2, 2, 2])
    d = np.sort(atoms.get_all_distances(mic=True), axis=1)[:, 1]
    assert np.isclose(d.min(), 5.05 / 2)
    assert np.isclose(d.max(), 5.05 * np.sqrt(5) / 4)
    assert np.sum(np.isclose(d, 5.05 / 2)) == 6 * 8


def test_c_over_a_rejected_for_cubic():
    with pytest.raises(ValueError, match="hcp"):
        make_lattice("fcc", "Cu", 3.6, [1, 1, 1], c_over_a=1.6)


def test_names():
    assert canonical_lattice("FCC") == "fcc"
    assert canonical_lattice("sc") == "simple_cubic"
    assert canonical_lattice("dhcp") is None
    assert set(BUILTIN_LATTICES) == {"fcc", "bcc", "hcp", "diamond", "simple_cubic", "a15"}
    with pytest.raises(ValueError, match="Unknown lattice"):
        make_lattice("b2", "Cu", 3.0, [1, 1, 1])


@pytest.mark.parametrize(
    "element, lattice, a",
    [("Cu", "fcc", 3.61), ("Fe", "bcc", 2.87), ("Mg", "hcp", 3.21),
     ("Si", "diamond", 5.43), ("Po", "simple_cubic", 3.35), ("W", "bcc", 3.16)],
)
def test_default_lattice(element, lattice, a):
    assert default_lattice(element) == (lattice, a)


@pytest.mark.parametrize("element", ["He", "Ga", "Xx"])
def test_default_lattice_unknown(element):
    # He: no crystal reference state; Ga: orthorhombic; Xx: not an element
    assert default_lattice(element) is None


# --- through the input layer -------------------------------------------------

from ase.io import read

from calphy.input import Calculation


def _calc(**kw):
    base = dict(mass=1.0, mode="fe", temperature=1000, pressure=0,
                reference_phase="solid", pair_style="eam", pair_coeff="* * x")
    return Calculation(**{**base, **kw})


def test_input_default_lattice(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calc = _calc(element="W")
    assert calc._original_lattice == "bcc"
    assert calc.lattice_constant == 3.16
    assert calc.repeat == [5, 5, 5]
    assert calc._natoms == 250
    assert len(read(calc.lattice, format="lammps-data", atom_style="atomic")) == 250


def test_input_sc_alias(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calc = _calc(element="Cu", lattice="SC", lattice_constant=2.5, repeat=[2, 2, 2])
    assert calc._original_lattice == "simple_cubic"
    assert calc._natoms == 8


def test_input_lattice_constant_from_element(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calc = _calc(element="Fe", lattice="bcc", repeat=[4, 4, 4])
    assert calc.lattice_constant == 2.87
    assert calc._natoms == 128


@pytest.mark.parametrize("lattice", ["dhcp", "b2", "no_such_file.data"])
def test_input_unknown_lattice(tmp_path, monkeypatch, lattice):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="neither a built-in lattice"):
        _calc(element="Cu", lattice=lattice, lattice_constant=3.0)


def test_input_c_over_a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calc = _calc(element="Zn", lattice="hcp", lattice_constant=2.66,
                 c_over_a=1.856, repeat=[2, 2, 2])
    atoms = read(calc.lattice, format="lammps-data", atom_style="atomic")
    assert np.isclose(atoms.cell[2, 2], 2 * 2.66 * 1.856)


def test_input_c_over_a_with_default_lattice(tmp_path, monkeypatch):
    # Mg defaults to hcp, so c_over_a applies
    monkeypatch.chdir(tmp_path)
    calc = _calc(element="Mg", c_over_a=1.624, repeat=[1, 1, 1])
    atoms = read(calc.lattice, format="lammps-data", atom_style="atomic")
    assert np.isclose(atoms.cell[2, 2] / atoms.cell[0, 0], 1.624)


def test_input_default_hcp_is_ideal(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calc = _calc(element="Mg", lattice="hcp", lattice_constant=3.2, repeat=[1, 1, 1])
    atoms = read(calc.lattice, format="lammps-data", atom_style="atomic")
    assert np.isclose(atoms.cell[2, 2], 3.2 * IDEAL_C_OVER_A)


@pytest.mark.parametrize("kw", [dict(element="Cu", lattice="fcc", lattice_constant=3.6),
                                dict(element="Cu")])
def test_input_c_over_a_rejected_for_non_hcp(tmp_path, monkeypatch, kw):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="c_over_a only applies"):
        _calc(c_over_a=1.6, **kw)


def test_input_c_over_a_rejected_for_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    make_lattice("hcp", "Mg", 3.2, [2, 2, 2]).write("mg.data", format="lammps-data")
    with pytest.raises(ValueError, match="built-in hcp"):
        _calc(element="Mg", lattice="mg.data", c_over_a=1.6)


def test_input_c_over_a_positive(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="greater than 0"):
        _calc(element="Mg", lattice="hcp", lattice_constant=3.2, c_over_a=0)
