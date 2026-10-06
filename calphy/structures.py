"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021-2026 (c) Sarath Menon, Yury Lysogorskiy, Ralf Drautz
Interdisciplinary Centre for Advanced Materials Simulation (ICAMS),
Ruhr University Bochum, 44801 Bochum, Germany

calphy is published and distributed under the Academic Software Licence v1.0 (ASL).
calphy is distributed in the hope that it will be useful for non-commercial academic
research, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the LICENSE file for details.

The ASL permits academic non-commercial use only. Contact
sarath.menon@ruhr-uni-bochum.de to enquire about commercial use rights.

More information about the program can be found in:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz.
"Automated Free Energy Calculation from Atomistic Simulations." Physical Review Materials 5(10), 2021
DOI: 10.1103/PhysRevMaterials.5.103801

For more information contact:
sarath.menon@ruhr-uni-bochum.de

Built-in single-element crystal lattices, built with ASE.
"""

import numpy as np
from ase.build import bulk
from ase.spacegroup import crystal
from ase.data import atomic_numbers, reference_states

IDEAL_C_OVER_A = np.sqrt(8.0 / 3.0)

# accepted lattice names -> canonical name
_ALIASES = {
    "fcc": "fcc",
    "bcc": "bcc",
    "hcp": "hcp",
    "diamond": "diamond",
    "simple_cubic": "simple_cubic",
    "sc": "simple_cubic",
    "a15": "a15",
}

# canonical name -> ase.build.bulk crystalstructure
_ASE_NAMES = {
    "fcc": "fcc",
    "bcc": "bcc",
    "hcp": "hcp",
    "diamond": "diamond",
    "simple_cubic": "sc",
}

BUILTIN_LATTICES = tuple(_ASE_NAMES) + ("a15",)


def canonical_lattice(name):
    """
    Return the canonical name of a built-in lattice, or None if ``name``
    is not one.
    """
    return _ALIASES.get(name.lower())


def default_lattice(element):
    """
    Ground-state lattice and lattice constant of an element.

    Taken from ``ase.data.reference_states``; only elements whose reference
    state is one of the built-in lattices are supported.

    Returns
    -------
    (str, float) or None
        Canonical lattice name and lattice constant a in Angstrom, or None
        if no built-in lattice is known for the element.
    """
    z = atomic_numbers.get(element)
    if z is None:
        return None
    state = reference_states[z]
    if not state:
        return None
    name = canonical_lattice(state["symmetry"])
    if name is None:
        return None
    return name, state["a"]


def make_lattice(name, element, lattice_constant, repeat, c_over_a=None):
    """
    Build a single-element crystal in an orthogonal cell.

    Parameters
    ----------
    name : str
        Built-in lattice name, see ``BUILTIN_LATTICES`` (``sc`` is accepted
        for ``simple_cubic``). ``a15`` is the single-element A15 (beta-W)
        structure.
    element : str
        Chemical symbol.
    lattice_constant : float
        Lattice constant a in Angstrom.
    repeat : sequence of 3 ints
        Repetitions of the unit cell along x, y and z.
    c_over_a : float, optional
        c/a ratio for hcp; the ideal value sqrt(8/3) if not given.

    Returns
    -------
    ase.Atoms
        Cubic cells for the cubic lattices; for hcp the 4-atom orthorhombic
        cell with edges a, sqrt(3) a, c.
    """
    canonical = canonical_lattice(name)
    if canonical is None:
        raise ValueError(
            f"Unknown lattice '{name}', available: {', '.join(BUILTIN_LATTICES)}"
        )
    if canonical == "hcp":
        if c_over_a is None:
            c_over_a = IDEAL_C_OVER_A
        atoms = bulk(
            element,
            "hcp",
            a=lattice_constant,
            c=c_over_a * lattice_constant,
            orthorhombic=True,
        )
    else:
        if c_over_a is not None:
            raise ValueError("c_over_a can only be used with the hcp lattice")
        if canonical == "a15":
            # Pm-3n, Wyckoff sites 2a and 6d, 8 atoms per cubic cell
            atoms = crystal(
                [element, element],
                basis=[(0, 0, 0), (0.25, 0.5, 0)],
                spacegroup=223,
                cellpar=[lattice_constant] * 3 + [90] * 3,
            )
        else:
            atoms = bulk(
                element, _ASE_NAMES[canonical], a=lattice_constant, cubic=True
            )
    return atoms.repeat(tuple(repeat))
