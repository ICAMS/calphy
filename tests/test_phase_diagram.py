import pytest
import numpy as np
from calphy.phase_diagram import read_structure_composition


def test_read_structure_composition_single_element():
    """Test reading a structure with only one element type"""
    # Use existing test structure file
    structure_file = 'tests/conf1.data'
    elements = ['Cu', 'Zr']
    
    comp = read_structure_composition(structure_file, elements)
    
    # Should have Cu atoms and 0 Zr atoms
    assert isinstance(comp, dict)
    assert 'Cu' in comp
    assert 'Zr' in comp
    assert comp['Cu'] > 0
    assert comp['Zr'] == 0
    assert sum(comp.values()) > 0


def test_read_structure_composition_element_not_in_structure():
    """Test that elements not present in structure get count 0"""
    structure_file = 'tests/conf1.data'
    elements = ['Cu', 'Zr', 'Al', 'Ni']
    
    comp = read_structure_composition(structure_file, elements)
    
    # Should have all elements with appropriate counts
    assert len(comp) == 4
    assert comp['Cu'] > 0
    assert comp['Zr'] == 0
    assert comp['Al'] == 0
    assert comp['Ni'] == 0


def test_read_structure_composition_element_order_matters():
    """Test that element list order determines type mapping"""
    structure_file = 'tests/conf1.data'
    
    # Different element orders should give different results
    # because element[0] maps to LAMMPS type 1, element[1] to type 2, etc.
    elements1 = ['Cu', 'Zr']
    elements2 = ['Zr', 'Cu']
    
    comp1 = read_structure_composition(structure_file, elements1)
    comp2 = read_structure_composition(structure_file, elements2)
    
    # With Cu first, Cu should have the atoms
    assert comp1['Cu'] > 0
    assert comp1['Zr'] == 0
    
    # With Zr first, Zr should have the atoms (because type 1 maps to first element)
    assert comp2['Zr'] > 0
    assert comp2['Cu'] == 0


def test_read_structure_composition_total_atoms():
    """Test that total atom count is preserved regardless of element order"""
    structure_file = 'tests/conf1.data'
    
    elements1 = ['Cu', 'Zr']
    elements2 = ['Zr', 'Cu']
    elements3 = ['Cu', 'Zr', 'Al']
    
    comp1 = read_structure_composition(structure_file, elements1)
    comp2 = read_structure_composition(structure_file, elements2)
    comp3 = read_structure_composition(structure_file, elements3)
    
    # Total should be the same regardless of element list
    total1 = sum(comp1.values())
    total2 = sum(comp2.values())
    total3 = sum(comp3.values())
    
    assert total1 == total2 == total3
    assert total1 > 0


def test_read_structure_composition_invalid_file():
    """Test that invalid file path raises appropriate error"""
    with pytest.raises(Exception):
        read_structure_composition('nonexistent_file.data', ['Cu', 'Zr'])


def test_read_structure_composition_empty_element_list():
    """Test behavior with empty element list"""
    structure_file = 'tests/conf1.data'
    elements = []
    
    comp = read_structure_composition(structure_file, elements)
    
    # Should return empty dict
    assert isinstance(comp, dict)
    assert len(comp) == 0


def test_read_structure_composition_single_element_list():
    """Test with single element in list"""
    structure_file = 'tests/conf1.data'
    elements = ['Cu']
    
    comp = read_structure_composition(structure_file, elements)
    
    assert len(comp) == 1
    assert 'Cu' in comp
    assert comp['Cu'] > 0


def test_composition_transformation_100_percent():
    """Test that 100% composition transformation is allowed (pure phase endpoint)"""
    from calphy.composition_transformation import CompositionTransformation
    
    # Create a test calculation that transforms all atoms from one element to another
    class TestCalc:
        def __init__(self):
            self.lattice = 'tests/conf1.data'
            self.element = ['Cu', 'Al']
            self.composition_scaling = type('obj', (object,), {
                '_input_chemical_composition': {'Cu': 500, 'Al': 0},
                '_output_chemical_composition': {'Cu': 0, 'Al': 500},
                'input_chemical_composition': property(lambda s: s._input_chemical_composition),
                'output_chemical_composition': property(lambda s: s._output_chemical_composition),
                'restrictions': []
            })()
    
    calc = TestCalc()
    
    # Should not raise an error
    comp = CompositionTransformation(calc)
    
    # Verify the transformation is set up correctly
    assert 'Cu' in comp.to_remove
    assert 'Al' in comp.to_add
    assert comp.to_remove['Cu'] == 500
    assert comp.to_add['Al'] == 500
    assert len(comp.transformation_list) == 1
    assert comp.transformation_list[0]['primary_element'] == 'Cu'
    assert comp.transformation_list[0]['secondary_element'] == 'Al'
    assert comp.transformation_list[0]['count'] == 500
