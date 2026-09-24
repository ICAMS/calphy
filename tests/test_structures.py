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
    [("simple_cubic", 1), ("bcc", 2), ("fcc", 4), ("diamond", 8)],
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


def test_c_over_a_rejected_for_cubic():
    with pytest.raises(ValueError, match="hcp"):
        make_lattice("fcc", "Cu", 3.6, [1, 1, 1], c_over_a=1.6)


def test_names():
    assert canonical_lattice("FCC") == "fcc"
    assert canonical_lattice("sc") == "simple_cubic"
    assert canonical_lattice("dhcp") is None
    assert set(BUILTIN_LATTICES) == {"fcc", "bcc", "hcp", "diamond", "simple_cubic"}
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
