"""
Test for composition transformation bug fix.

This test verifies that composition_transformation.py correctly handles
LAMMPS data files where all atoms are of a single type that needs to be
partially transformed to another element.

Bug: When reading LAMMPS data files, the code wasn't providing Z_of_type
mapping, causing ASE to incorrectly identify atom types. Additionally, the
code was using pyscal's numeric types instead of element species names.
"""

import pytest
import numpy as np
import os
import tempfile
from calphy.composition_transformation import CompositionTransformation
from calphy.phase_diagram import SimpleCalculation


def verify_lammps_data_file(
    filepath, expected_natoms, expected_ntypes, expected_type_counts
):
    """
    Verify contents of a LAMMPS data file.

    Parameters
    ----------
    filepath : str
        Path to the LAMMPS data file
    expected_natoms : int
        Expected total number of atoms
    expected_ntypes : int
        Expected number of atom types declared in header
    expected_type_counts : dict
        Dictionary mapping type number to expected count
        e.g., {1: 1000, 2: 600, 3: 400}
    """
    assert os.path.exists(filepath), f"File {filepath} does not exist"

    with open(filepath, "r") as f:
        lines = f.readlines()

    # Check header
    header_natoms = None
    header_ntypes = None
    atoms_section_start = None

    for i, line in enumerate(lines):
        if "atoms" in line and header_natoms is None:
            header_natoms = int(line.split()[0])
        if "atom types" in line:
            header_ntypes = int(line.split()[0])
        if "Atoms" in line:
            atoms_section_start = i + 1
            break

    assert (
        header_natoms == expected_natoms
    ), f"Expected {expected_natoms} atoms in header, got {header_natoms}"
    assert (
        header_ntypes == expected_ntypes
    ), f"Expected {expected_ntypes} atom types in header, got {header_ntypes}"

    # Parse atoms section and count types
    type_counts = {}
    atom_count = 0

    for line in lines[atoms_section_start:]:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 5:  # atom_id type x y z
            atom_type = int(parts[1])
            type_counts[atom_type] = type_counts.get(atom_type, 0) + 1
            atom_count += 1

    assert (
        atom_count == expected_natoms
    ), f"Expected {expected_natoms} atoms in data, got {atom_count}"

    # Verify type counts
    for type_num, expected_count in expected_type_counts.items():
        actual_count = type_counts.get(type_num, 0)
        assert (
            actual_count == expected_count
        ), f"Expected {expected_count} atoms of type {type_num}, got {actual_count}"

    # Verify all types are within declared range
    for type_num in type_counts.keys():
        assert (
            1 <= type_num <= expected_ntypes
        ), f"Type {type_num} out of range [1, {expected_ntypes}]"

    return type_counts


@pytest.fixture
def mg_structure_file():
    """Create a temporary LAMMPS data file with only Mg atoms (all type 2)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        filepath = os.path.join(tmpdir, "Mg.lammps.data")

        # Simple HCP Mg structure - all atoms are type 2
        natoms = 2048
        content = f"""(written by ASE)

{natoms} atoms
2 atom types

0.0                   25.68  xlo xhi
0.0      44.479064738368763  ylo yhi
0.0                41.93544  zlo zhi

Atoms # atomic

"""

        # Generate atom positions for a simple HCP lattice
        nx, ny, nz = 8, 8, 16  # Should give 8*8*16*2 = 2048 atoms
        lx, ly, lz = 25.68, 44.479064738368763, 41.93544

        atom_id = 1
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    for sublat in range(2):
                        x = (ix + sublat * 0.5) * lx / nx
                        y = (iy + sublat * 0.5) * ly / ny
                        z = (iz + sublat * 0.5) * lz / nz
                        content += (
                            f"  {atom_id:4d}   2  {x:20.15f}  {y:20.15f}  {z:20.15f}\n"
                        )
                        atom_id += 1

        with open(filepath, "w") as f:
            f.write(content)

        yield filepath


def test_composition_transformation_single_type_file(mg_structure_file):
    """
    Test composition transformation when starting structure has all atoms
    of a single type (type 2 = Mg) and we want to add Al atoms.

    Input: 2048 Mg atoms (all type 2 in LAMMPS file)
    Output: 102 Al + 1946 Mg atoms

    This test reproduces the bug where the transformation fails with
    "ValueError: high <= 0" due to incorrect type/species handling.
    """
    # Create a simple calculation object
    calc = SimpleCalculation(
        lattice=mg_structure_file,
        element=["Al", "Mg"],
        input_chemical_composition={"Al": 0, "Mg": 2048},
        output_chemical_composition={"Al": 102, "Mg": 1946},
    )

    # Create composition transformation
    comp = CompositionTransformation(calc)

    # Verify basic properties
    assert comp.natoms == 2048
    assert comp.typedict == {"Al": 1, "Mg": 2}

    # Verify input/output compositions
    assert comp.input_chemical_composition == {"Al": 0, "Mg": 2048}
    assert comp.output_chemical_composition == {"Al": 102, "Mg": 1946}

    # Verify transformation requirements
    assert comp.to_remove == {"Mg": 102}
    assert comp.to_add == {"Al": 102}

    # Check that we have exactly one transformation: Mg -> Al
    assert len(comp.transformation_list) == 1
    assert comp.transformation_list[0]["primary_element"] == "Mg"
    assert comp.transformation_list[0]["secondary_element"] == "Al"
    assert comp.transformation_list[0]["count"] == 102

    # Verify pair lists
    assert comp.pair_list_old == ["Mg", "Mg"]
    assert comp.pair_list_new == ["Al", "Mg"]

    # Test pair_coeff update
    pair_coeff_input = "pair_coeff * * /path/to/potential.yace Al Mg"
    pc_old, pc_new = comp.update_pair_coeff(pair_coeff_input)

    # Verify pair coeffs are different and properly formatted
    # pc_old should have Mg because all atoms start as Mg
    # pc_new should have both Al and Mg after transformation
    assert "Mg" in pc_old
    assert "Al" in pc_new and "Mg" in pc_new

    # Write structure and verify
    output_structure = os.path.join(
        os.path.dirname(mg_structure_file), "transformed.lammps.data"
    )
    comp.write_structure(output_structure)

    assert os.path.exists(output_structure)

    # Read and verify the written structure
    import re
    with open(output_structure, "r") as f:
        content = f.read()
        # Check for "2048 atoms" with flexible whitespace (space or tab)
        assert re.search(r"2048\s+atoms", content), "Expected '2048 atoms' in output file"

    # Verify entropy contribution calculation
    entropy = comp.entropy_contribution
    assert isinstance(entropy, float)
    assert entropy < 0  # Should be negative for this transformation


def test_composition_transformation_with_pace_potential(mg_structure_file):
    """
    Test with PACE potential format (as in the original bug report).
    Verifies that pair_coeff commands are correctly generated for PACE potentials.
    """
    # Create a simple calculation object
    calc = SimpleCalculation(
        lattice=mg_structure_file,
        element=["Al", "Mg"],
        input_chemical_composition={"Al": 0, "Mg": 2048},
        output_chemical_composition={"Al": 102, "Mg": 1946},
    )

    # Create composition transformation
    comp = CompositionTransformation(calc)

    # Test PACE potential format
    pair_coeff_input = "pair_coeff * * /home/user/AlMgZn.yace Al Mg"
    pc_old, pc_new = comp.update_pair_coeff(pair_coeff_input)

    # Both should contain the potential file path
    assert "/home/user/AlMgZn.yace" in pc_old
    assert "/home/user/AlMgZn.yace" in pc_new

    # Verify element specifications
    assert "Mg" in pc_old
    assert "Al" in pc_new and "Mg" in pc_new

    # Write and verify structure
    output_structure = os.path.join(
        os.path.dirname(mg_structure_file), "transformed_pace.lammps.data"
    )
    comp.write_structure(output_structure)

    # Verify file content: 2048 atoms, 2 types (Al and Mg)
    # Type 1: 102 Al (Mg→Al), Type 2: 1946 Mg (Mg→Mg)
    type_counts = verify_lammps_data_file(
        output_structure,
        expected_natoms=2048,
        expected_ntypes=2,
        expected_type_counts={1: 102, 2: 1946},
    )


def test_composition_transformation_swap_types(mg_structure_file):
    """
    Test get_swap_types method returns correct type mappings.
    """
    # Create a simple calculation object
    calc = SimpleCalculation(
        lattice=mg_structure_file,
        element=["Al", "Mg"],
        input_chemical_composition={"Al": 0, "Mg": 2048},
        output_chemical_composition={"Al": 102, "Mg": 1946},
    )

    # Create composition transformation
    comp = CompositionTransformation(calc)

    # Get swap types
    swap_types = comp.get_swap_types()

    assert isinstance(swap_types, list)
    assert len(swap_types) == 2

    # Write and verify structure
    output_structure = os.path.join(
        os.path.dirname(mg_structure_file), "transformed_swap.lammps.data"
    )
    comp.write_structure(output_structure)

    # Verify file content: 2048 atoms, 2 types
    # Type 1: 102 Al (Mg→Al), Type 2: 1946 Mg (Mg→Mg)
    type_counts = verify_lammps_data_file(
        output_structure,
        expected_natoms=2048,
        expected_ntypes=2,
        expected_type_counts={1: 102, 2: 1946},
    )


@pytest.fixture
def al_structure_file():
    """Create a temporary LAMMPS data file with only Al atoms (all type 1)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        filepath = os.path.join(tmpdir, "Al.lammps.data")

        natoms = 2000
        content = f"""(written by ASE)

{natoms} atoms
2 atom types

0.0                   32.40  xlo xhi
0.0                   32.40  ylo yhi
0.0                   32.40  zlo zhi

Atoms # atomic

"""

        # Generate atom positions for a simple FCC lattice
        nx, ny, nz = 10, 10, 10
        lx, ly, lz = 32.40, 32.40, 32.40

        atom_id = 1
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    for basis in [
                        (0, 0, 0),
                        (0.5, 0.5, 0),
                        (0.5, 0, 0.5),
                        (0, 0.5, 0.5),
                    ]:
                        if atom_id > natoms:
                            break
                        x = (ix + basis[0]) * lx / nx
                        y = (iy + basis[1]) * ly / ny
                        z = (iz + basis[2]) * lz / nz
                        content += (
                            f"  {atom_id:4d}   1  {x:20.15f}  {y:20.15f}  {z:20.15f}\n"
                        )
                        atom_id += 1
                    if atom_id > natoms:
                        break
                if atom_id > natoms:
                    break
            if atom_id > natoms:
                break

        with open(filepath, "w") as f:
            f.write(content)

        yield filepath


@pytest.fixture
def almg_structure_file():
    """Create a temporary LAMMPS data file with 50% Al and 50% Mg."""
    with tempfile.TemporaryDirectory() as tmpdir:
        filepath = os.path.join(tmpdir, "AlMg.lammps.data")

        natoms = 2000
        n_al = natoms // 2
        content = f"""(written by ASE)

{natoms} atoms
2 atom types

0.0                   32.40  xlo xhi
0.0                   32.40  ylo yhi
0.0                   32.40  zlo zhi

Atoms # atomic

"""

        nx, ny, nz = 10, 10, 10
        lx, ly, lz = 32.40, 32.40, 32.40

        atom_id = 1
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    for basis in [
                        (0, 0, 0),
                        (0.5, 0.5, 0),
                        (0.5, 0, 0.5),
                        (0, 0.5, 0.5),
                    ]:
                        if atom_id > natoms:
                            break
                        x = (ix + basis[0]) * lx / nx
                        y = (iy + basis[1]) * ly / ny
                        z = (iz + basis[2]) * lz / nz
                        atom_type = 1 if atom_id <= n_al else 2
                        content += f"  {atom_id:4d}   {atom_type}  {x:20.15f}  {y:20.15f}  {z:20.15f}\n"
                        atom_id += 1
                    if atom_id > natoms:
                        break
                if atom_id > natoms:
                    break
            if atom_id > natoms:
                break

        with open(filepath, "w") as f:
            f.write(content)

        yield filepath


def test_100_percent_transformation(al_structure_file):
    """
    Test 100% composition transformation: Pure Al → Pure Mg.

    Edge case where all atoms transform from one element to another.
    """
    natoms = 2000

    calc = SimpleCalculation(
        lattice=al_structure_file,
        element=["Al", "Mg"],
        input_chemical_composition={"Al": natoms, "Mg": 0},
        output_chemical_composition={"Al": 0, "Mg": natoms},
    )

    comp = CompositionTransformation(calc)

    # Verify transformation
    assert comp.to_remove == {"Al": natoms}
    assert comp.to_add == {"Mg": natoms}

    # Should have only one mapping: Al→Mg
    assert len(comp.unique_mappings) == 1
    assert comp.unique_mappings[0] == "Al-Mg"

    # Pair lists - matches actual structure types (1 type: pure Al)
    # Ensures consistency: data file has 1 type, pair_coeff maps 1 element
    assert comp.pair_list_old == ["Al"]
    assert comp.pair_list_new == ["Mg"]

    # Test pair_coeff
    pair_coeff_input = "pair_coeff * * AlMg.yace Al Mg"
    pc_old, pc_new = comp.update_pair_coeff(pair_coeff_input)

    assert "Al" in pc_old
    assert "Mg" in pc_new

    # Swap types - for 100% transformation, only the transforming type
    # Uses the target element's type from typedict
    swap_types = comp.get_swap_types()
    assert len(swap_types) == 1
    # Type number depends on how mappingdict assigns it

    # Write and verify structure
    output_structure = os.path.join(
        os.path.dirname(al_structure_file), "transformed_100pct.lammps.data"
    )
    comp.write_structure(output_structure)

    # Verify file content: 2000 atoms, 1 type (all Al→Mg)
    type_counts = verify_lammps_data_file(
        output_structure,
        expected_natoms=2000,
        expected_ntypes=1,
        expected_type_counts={1: 2000},
    )


def test_partial_transformation_al_to_almg(al_structure_file):
    """
    Test partial transformation: 100% Al → 80% Al + 20% Mg.

    This tests the case where we add a new element to a pure structure.
    """
    natoms = 2000
    n_al_final = int(0.80 * natoms)
    n_mg_final = int(0.20 * natoms)

    calc = SimpleCalculation(
        lattice=al_structure_file,
        element=["Al", "Mg"],
        input_chemical_composition={"Al": natoms, "Mg": 0},
        output_chemical_composition={"Al": n_al_final, "Mg": n_mg_final},
    )

    comp = CompositionTransformation(calc)

    # Verify transformation
    assert comp.to_remove == {"Al": natoms - n_al_final}
    assert comp.to_add == {"Mg": n_mg_final}

    # Should have two mappings: Al→Al and Al→Mg
    assert len(comp.unique_mappings) == 2
    assert "Al-Al" in comp.unique_mappings
    assert "Al-Mg" in comp.unique_mappings

    # Pair lists
    assert comp.pair_list_old == ["Al", "Al"]
    assert comp.pair_list_new == ["Al", "Mg"]

    # Swap types should include both types (same initial element)
    # Order: [transforming_type, conserved_type] = [Al→Mg, Al→Al]
    swap_types = comp.get_swap_types()
    assert len(swap_types) == 2
    assert swap_types == [2, 1]

    # Write and verify structure
    output_structure = os.path.join(
        os.path.dirname(al_structure_file), "transformed_partial.lammps.data"
    )
    comp.write_structure(output_structure)

    # Verify file content: 2000 atoms, 2 types
    # Type 1: 1600 Al (Al→Al), Type 2: 400 Mg (Al→Mg)
    type_counts = verify_lammps_data_file(
        output_structure,
        expected_natoms=2000,
        expected_ntypes=2,
        expected_type_counts={1: 1600, 2: 400},
    )


def test_enrichment_transformation(almg_structure_file):
    """
    Test enrichment: 50% Al + 50% Mg → 80% Al + 20% Mg.

    This tests transforming existing Mg atoms to Al.
    """
    natoms = 2000
    n_al_initial = natoms // 2
    n_mg_initial = natoms // 2
    n_al_final = int(0.80 * natoms)
    n_mg_final = int(0.20 * natoms)

    calc = SimpleCalculation(
        lattice=almg_structure_file,
        element=["Al", "Mg"],
        input_chemical_composition={"Al": n_al_initial, "Mg": n_mg_initial},
        output_chemical_composition={"Al": n_al_final, "Mg": n_mg_final},
    )

    comp = CompositionTransformation(calc)

    # Verify transformation
    assert comp.to_remove == {"Mg": n_mg_initial - n_mg_final}
    assert comp.to_add == {"Al": n_al_final - n_al_initial}

    # Should have three mappings: Al→Al, Mg→Al, Mg→Mg
    assert len(comp.unique_mappings) == 3
    assert "Al-Al" in comp.unique_mappings
    assert "Mg-Al" in comp.unique_mappings
    assert "Mg-Mg" in comp.unique_mappings

    # Pair lists
    assert comp.pair_list_old == ["Al", "Mg", "Mg"]
    assert comp.pair_list_new == ["Al", "Al", "Mg"]

    # Swap types should be the two outcomes for Mg atoms (Mg→Al and Mg→Mg)
    # Order: [transforming_type, conserved_type] = [Mg→Al, Mg→Mg]
    swap_types = comp.get_swap_types()
    assert len(swap_types) == 2
    # Ordered as [type 2 (Mg→Al), type 3 (Mg→Mg)]
    assert swap_types == [2, 3]

    # Write and verify structure
    output_structure = os.path.join(
        os.path.dirname(almg_structure_file), "transformed_enrichment.lammps.data"
    )
    comp.write_structure(output_structure)

    # Verify file content: 2000 atoms, 3 types
    # Type 1: 1000 Al (Al→Al), Type 2: 600 Al (Mg→Al), Type 3: 400 Mg (Mg→Mg)
    type_counts = verify_lammps_data_file(
        output_structure,
        expected_natoms=2000,
        expected_ntypes=3,
        expected_type_counts={1: 1000, 2: 600, 3: 400},
    )


def test_depletion_transformation(almg_structure_file):
    """
    Test depletion: 50% Al + 50% Mg → 20% Al + 80% Mg.

    This tests transforming existing Al atoms to Mg.
    """
    natoms = 2000
    n_al_initial = natoms // 2
    n_mg_initial = natoms // 2
    n_al_final = int(0.20 * natoms)
    n_mg_final = int(0.80 * natoms)

    calc = SimpleCalculation(
        lattice=almg_structure_file,
        element=["Al", "Mg"],
        input_chemical_composition={"Al": n_al_initial, "Mg": n_mg_initial},
        output_chemical_composition={"Al": n_al_final, "Mg": n_mg_final},
    )

    comp = CompositionTransformation(calc)

    # Verify transformation
    assert comp.to_remove == {"Al": n_al_initial - n_al_final}
    assert comp.to_add == {"Mg": n_mg_final - n_mg_initial}

    # Should have three mappings: Al→Al, Al→Mg, Mg→Mg
    assert len(comp.unique_mappings) == 3
    assert "Al-Al" in comp.unique_mappings
    assert "Al-Mg" in comp.unique_mappings
    assert "Mg-Mg" in comp.unique_mappings

    # Pair lists
    assert comp.pair_list_old == ["Al", "Al", "Mg"]
    assert comp.pair_list_new == ["Al", "Mg", "Mg"]

    # Swap types should be the two Al types
    # Order: [transforming_type, conserved_type] = [Al→Mg, Al→Al]
    swap_types = comp.get_swap_types()
    assert len(swap_types) == 2
    assert swap_types == [2, 1]

    # Write and verify structure
    output_structure = os.path.join(
        os.path.dirname(almg_structure_file), "transformed_depletion.lammps.data"
    )
    comp.write_structure(output_structure)

    # Verify file content: 2000 atoms, 3 types
    # Type 1: 400 Al (Al→Al), Type 2: 600 Mg (Al→Mg), Type 3: 1000 Mg (Mg→Mg)
    type_counts = verify_lammps_data_file(
        output_structure,
        expected_natoms=2000,
        expected_ntypes=3,
        expected_type_counts={1: 400, 2: 600, 3: 1000},
    )


def test_small_transformation(mg_structure_file):
    """
    Test very small transformation: only 1 atom changes.

    Edge case with minimal transformation.
    """
    natoms = 2048

    calc = SimpleCalculation(
        lattice=mg_structure_file,
        element=["Al", "Mg"],
        input_chemical_composition={"Al": 0, "Mg": natoms},
        output_chemical_composition={"Al": 1, "Mg": natoms - 1},
    )

    comp = CompositionTransformation(calc)

    # Verify transformation
    assert comp.to_remove == {"Mg": 1}
    assert comp.to_add == {"Al": 1}

    # Should have two mappings
    assert len(comp.unique_mappings) == 2
    assert "Mg-Al" in comp.unique_mappings
    assert "Mg-Mg" in comp.unique_mappings

    # Transformation list should have 1 atom
    assert len(comp.transformation_list) == 1
    assert comp.transformation_list[0]["count"] == 1

    # Write and verify structure
    output_structure = os.path.join(
        os.path.dirname(mg_structure_file), "transformed_small.lammps.data"
    )
    comp.write_structure(output_structure)

    # Verify file content: 2048 atoms, 2 types
    # Type 1: 1 Al (Mg→Al), Type 2: 2047 Mg (Mg→Mg)
    type_counts = verify_lammps_data_file(
        output_structure,
        expected_natoms=2048,
        expected_ntypes=2,
        expected_type_counts={1: 1, 2: 2047},
    )


def test_three_element_transformation(almg_structure_file):
    """
    Test transformation involving three elements: Al-Mg → Al-Mg-Cu.

    This tests adding a third element to an existing binary alloy.
    """
    natoms = 2000
    n_al_initial = natoms // 2
    n_mg_initial = natoms // 2

    # Final: 40% Al, 40% Mg, 20% Cu
    n_al_final = int(0.40 * natoms)
    n_mg_final = int(0.40 * natoms)
    n_cu_final = int(0.20 * natoms)

    calc = SimpleCalculation(
        lattice=almg_structure_file,
        element=["Al", "Mg", "Cu"],
        input_chemical_composition={"Al": n_al_initial, "Mg": n_mg_initial, "Cu": 0},
        output_chemical_composition={
            "Al": n_al_final,
            "Mg": n_mg_final,
            "Cu": n_cu_final,
        },
    )

    comp = CompositionTransformation(calc)

    # Verify transformations
    assert comp.to_remove["Al"] == n_al_initial - n_al_final
    assert comp.to_remove["Mg"] == n_mg_initial - n_mg_final
    assert comp.to_add == {"Cu": n_cu_final}

    # Should have multiple mappings including Cu
    assert any("Cu" in mapping for mapping in comp.unique_mappings)

    # Verify structure can be written
    output_structure = os.path.join(
        os.path.dirname(almg_structure_file), "AlMgCu.lammps.data"
    )
    comp.write_structure(output_structure)
    assert os.path.exists(output_structure)

    # Verify file content: 2000 atoms total with correct type distribution
    # Multiple types for Al→Al, Al→Cu, Mg→Mg, Mg→Cu transformations
    # Expected: 800 Al, 800 Mg, 400 Cu total
    type_counts = verify_lammps_data_file(
        output_structure,
        expected_natoms=2000,
        expected_ntypes=len(comp.unique_mappings),
        expected_type_counts={},  # Complex mapping, just verify total
    )
    # Verify total count is correct
    assert sum(type_counts.values()) == 2000
