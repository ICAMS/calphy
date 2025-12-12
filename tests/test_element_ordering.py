import pytest
from calphy.input import Calculation, read_inputfile


def test_correct_element_ordering():
    """Test that correct element ordering is accepted"""
    calc_dict = {
        'element': ['Cu', 'Zr'],
        'mass': [63.546, 91.224],
        'pair_coeff': ['* * potential.eam.fs Cu Zr'],
        'pair_style': ['eam/fs'],
        'mode': 'fe',
        'temperature': 1000,
        'pressure': 0,
        'lattice': 'tests/conf1.data',  # Use existing data file
    }
    
    # Should not raise any exception
    calc = Calculation(**calc_dict)
    assert calc.element == ['Cu', 'Zr']
    assert calc.mass == [63.546, 91.224]


def test_wrong_element_ordering():
    """Test that mismatched element ordering is rejected"""
    calc_dict = {
        'element': ['Cu', 'Zr'],
        'mass': [63.546, 91.224],
        'pair_coeff': ['* * potential.eam.fs Zr Cu'],  # Wrong order!
        'pair_style': ['eam/fs'],
        'mode': 'fe',
        'temperature': 1000,
        'pressure': 0,
        'lattice': 'fcc',
        'lattice_constant': 3.61,
    }
    
    with pytest.raises(ValueError) as exc_info:
        calc = Calculation(**calc_dict)
    
    assert 'ordering mismatch' in str(exc_info.value).lower()


def test_single_element_no_ordering_issue():
    """Test that single element systems don't trigger ordering validation"""
    calc_dict = {
        'element': ['Cu'],
        'mass': [63.546],
        'pair_coeff': ['* * potential.eam Cu'],
        'pair_style': ['eam'],
        'mode': 'fe',
        'temperature': 1000,
        'pressure': 0,
        'lattice': 'fcc',
        'lattice_constant': 3.61,
    }
    
    # Should not raise any exception for single element
    calc = Calculation(**calc_dict)
    assert calc.element == ['Cu']


def test_no_elements_in_pair_coeff():
    """Test that pair_coeff without element names skips validation"""
    calc_dict = {
        'element': ['Cu', 'Zr'],
        'mass': [63.546, 91.224],
        'pair_coeff': ['* * potential.eam'],  # No elements specified
        'pair_style': ['eam'],
        'mode': 'fe',
        'temperature': 1000,
        'pressure': 0,
        'lattice': 'tests/conf1.data',  # Use existing data file
    }
    
    # Should not raise exception when pair_coeff has no elements
    calc = Calculation(**calc_dict)
    assert calc.element == ['Cu', 'Zr']


def test_element_mismatch():
    """Test that completely different elements in pair_coeff are rejected"""
    calc_dict = {
        'element': ['Cu', 'Zr'],
        'mass': [63.546, 91.224],
        'pair_coeff': ['* * potential.eam.fs Al Ni'],  # Different elements!
        'pair_style': ['eam/fs'],
        'mode': 'fe',
        'temperature': 1000,
        'pressure': 0,
        'lattice': 'fcc',
        'lattice_constant': 3.61,
    }
    
    with pytest.raises(ValueError) as exc_info:
        calc = Calculation(**calc_dict)
    
    assert 'element mismatch' in str(exc_info.value).lower()


def test_three_element_ordering():
    """Test ordering validation works for 3+ elements"""
    # Correct ordering
    calc_dict = {
        'element': ['Cu', 'Zr', 'Al'],
        'mass': [63.546, 91.224, 26.982],
        'pair_coeff': ['* * potential.eam.fs Cu Zr Al'],
        'pair_style': ['eam/fs'],
        'mode': 'fe',
        'temperature': 1000,
        'pressure': 0,
        'lattice': 'tests/conf1.data',  # Use existing data file
    }
    
    calc = Calculation(**calc_dict)
    assert calc.element == ['Cu', 'Zr', 'Al']
    
    # Wrong ordering
    calc_dict_wrong = {
        'element': ['Cu', 'Zr', 'Al'],
        'mass': [63.546, 91.224, 26.982],
        'pair_coeff': ['* * potential.eam.fs Al Cu Zr'],  # Different order
        'pair_style': ['eam/fs'],
        'mode': 'fe',
        'temperature': 1000,
        'pressure': 0,
        'lattice': 'tests/conf1.data',  # Use existing data file
    }
    
    with pytest.raises(ValueError) as exc_info:
        calc = Calculation(**calc_dict_wrong)
    
    assert 'ordering mismatch' in str(exc_info.value).lower()


def test_existing_example_files():
    """Test that existing example files still load correctly"""
    # Single element examples should work
    calcs = read_inputfile('tests/input.yaml')
    assert len(calcs) > 0
    assert calcs[0].element is not None
