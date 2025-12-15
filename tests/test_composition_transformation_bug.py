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
                        content += f"  {atom_id:4d}   2  {x:20.15f}  {y:20.15f}  {z:20.15f}\n"
                        atom_id += 1
        
        with open(filepath, 'w') as f:
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
        output_chemical_composition={"Al": 102, "Mg": 1946}
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
    assert comp.transformation_list[0]['primary_element'] == 'Mg'
    assert comp.transformation_list[0]['secondary_element'] == 'Al'
    assert comp.transformation_list[0]['count'] == 102
    
    # Verify pair lists
    assert comp.pair_list_old == ['Mg', 'Mg']
    assert comp.pair_list_new == ['Al', 'Mg']
    
    # Test pair_coeff update
    pair_coeff_input = "pair_coeff * * /path/to/potential.yace Al Mg"
    pc_old, pc_new = comp.update_pair_coeff(pair_coeff_input)
    
    # Verify pair coeffs are different and properly formatted
    # pc_old should have Mg because all atoms start as Mg
    # pc_new should have both Al and Mg after transformation
    assert "Mg" in pc_old
    assert "Al" in pc_new and "Mg" in pc_new
    
    # Write structure and verify
    output_structure = os.path.join(os.path.dirname(mg_structure_file), "transformed.lammps.data")
    comp.write_structure(output_structure)
    
    assert os.path.exists(output_structure)
    
    # Read and verify the written structure
    with open(output_structure, 'r') as f:
        content = f.read()
        assert "2048 atoms" in content
    
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
        output_chemical_composition={"Al": 102, "Mg": 1946}
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


def test_composition_transformation_swap_types(mg_structure_file):
    """
    Test get_swap_types method returns correct type mappings.
    """
    # Create a simple calculation object
    calc = SimpleCalculation(
        lattice=mg_structure_file,
        element=["Al", "Mg"],
        input_chemical_composition={"Al": 0, "Mg": 2048},
        output_chemical_composition={"Al": 102, "Mg": 1946}
    )
    
    # Create composition transformation
    comp = CompositionTransformation(calc)
    
    # Get swap types
    swap_types = comp.get_swap_types()
    
    assert isinstance(swap_types, list)
    assert len(swap_types) == 2
