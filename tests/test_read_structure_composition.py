import pytest
import os
from calphy.phase_diagram import read_structure_composition


# Get the directory where test files are located
TEST_DIR = os.path.dirname(os.path.abspath(__file__))


def test_read_structure_composition_with_masses():
    """
    Test read_structure_composition with a LAMMPS data file that has a Masses section.
    The file has 1 Al atom (type 1) and 1 Ni atom (type 2).
    """
    lattice_file = os.path.join(TEST_DIR, "test_data_with_masses.data")
    element_list = ["Al", "Ni"]
    
    composition = read_structure_composition(lattice_file, element_list)
    
    # Should have 1 atom of each type
    assert composition["Al"] == 1
    assert composition["Ni"] == 1
    assert sum(composition.values()) == 2


def test_read_structure_composition_without_masses():
    """
    Test read_structure_composition with a LAMMPS data file without a Masses section.
    The file has 1 Al atom (type 1) and 1 Ni atom (type 2).
    """
    lattice_file = os.path.join(TEST_DIR, "test_data_without_masses.data")
    element_list = ["Al", "Ni"]
    
    composition = read_structure_composition(lattice_file, element_list)
    
    # Should have 1 atom of each type
    assert composition["Al"] == 1
    assert composition["Ni"] == 1
    assert sum(composition.values()) == 2


def test_read_structure_composition_element_ordering():
    """
    Test that element ordering matters: element[0] corresponds to LAMMPS type 1,
    element[1] corresponds to LAMMPS type 2.
    """
    lattice_file = os.path.join(TEST_DIR, "test_data_with_masses.data")
    
    # Correct ordering: Al is type 1, Ni is type 2
    composition1 = read_structure_composition(lattice_file, ["Al", "Ni"])
    assert composition1["Al"] == 1
    assert composition1["Ni"] == 1
    
    # Reversed ordering: Ni is type 1, Al is type 2
    # This should give reversed counts
    composition2 = read_structure_composition(lattice_file, ["Ni", "Al"])
    assert composition2["Ni"] == 1  # Now Ni is assigned to type 1
    assert composition2["Al"] == 1  # Now Al is assigned to type 2


def test_read_structure_composition_missing_element():
    """
    Test that elements not present in the structure get count 0.
    """
    lattice_file = os.path.join(TEST_DIR, "test_data_with_masses.data")
    element_list = ["Al", "Ni", "Cu"]  # Cu is not in the file
    
    composition = read_structure_composition(lattice_file, element_list)
    
    assert composition["Al"] == 1
    assert composition["Ni"] == 1
    assert composition["Cu"] == 0
    assert sum(composition.values()) == 2


def test_read_structure_composition_consistency():
    """
    Test that both files (with and without Masses section) give the same result.
    """
    file_with_masses = os.path.join(TEST_DIR, "test_data_with_masses.data")
    file_without_masses = os.path.join(TEST_DIR, "test_data_without_masses.data")
    element_list = ["Al", "Ni"]
    
    comp1 = read_structure_composition(file_with_masses, element_list)
    comp2 = read_structure_composition(file_without_masses, element_list)
    
    # Both should give identical results
    assert comp1 == comp2
    assert comp1["Al"] == comp2["Al"] == 1
    assert comp1["Ni"] == comp2["Ni"] == 1
