"""
Tests for phase diagram preparation validation
"""

import pytest
import os
import tempfile
import yaml
from calphy.phase_diagram import prepare_inputs_for_phase_diagram


def test_binary_system_validation_too_many_elements():
    """Test that phase diagram preparation rejects systems with more than 2 elements"""
    
    phase_data = {
        'phases': [{
            'phase_name': 'test_phase',
            'element': ['Cu', 'Zr', 'Al'],  # 3 elements - should fail
            'mass': [63.546, 91.224, 26.982],
            'reference_phase': 'solid',
            'lattice': 'test.data',
            'pair_coeff': ['* * Cu-Zr-Al.eam.alloy Cu Zr Al'],
            'composition': {
                'reference_element': 'Cu',
                'range': [0.3, 0.7],
                'interval': 0.1,
                'reference': 0.5
            },
            'temperature': {
                'range': [300, 500],
                'interval': 100
            }
        }]
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.safe_dump(phase_data, f)
        temp_file = f.name
    
    try:
        with pytest.raises(ValueError, match="currently supports only binary systems.*Found 3 elements"):
            prepare_inputs_for_phase_diagram(temp_file)
    finally:
        os.unlink(temp_file)


def test_binary_system_validation_single_element():
    """Test that phase diagram preparation rejects single element systems"""
    
    phase_data = {
        'phases': [{
            'phase_name': 'test_phase',
            'element': ['Cu'],  # 1 element - should fail
            'mass': [63.546],
            'reference_phase': 'solid',
            'lattice': 'test.data',
            'pair_coeff': ['* * Cu.eam'],
            'composition': {
                'reference_element': 'Cu',
                'range': [1.0],
                'interval': 0.1,
                'reference': 1.0
            },
            'temperature': {
                'range': [300],
                'interval': 100
            }
        }]
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.safe_dump(phase_data, f)
        temp_file = f.name
    
    try:
        with pytest.raises(ValueError, match="currently supports only binary systems.*Found 1 element"):
            prepare_inputs_for_phase_diagram(temp_file)
    finally:
        os.unlink(temp_file)


def test_element_ordering_validation_in_phase_diagram():
    """Test that phase diagram preparation validates element ordering matches pair_coeff"""
    
    phase_data = {
        'phases': [{
            'phase_name': 'test_phase',
            'element': ['Cu', 'Zr'],
            'mass': [63.546, 91.224],
            'reference_phase': 'solid',
            'lattice': 'test.data',
            'pair_coeff': ['* * CuZr.eam.alloy Zr Cu'],  # Wrong order - should fail
            'composition': {
                'reference_element': 'Cu',
                'range': [0.3, 0.7],
                'interval': 0.1,
                'reference': 0.5
            },
            'temperature': {
                'range': [300, 500],
                'interval': 100
            }
        }]
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.safe_dump(phase_data, f)
        temp_file = f.name
    
    try:
        with pytest.raises(ValueError, match="Element ordering mismatch"):
            prepare_inputs_for_phase_diagram(temp_file)
    finally:
        os.unlink(temp_file)


def test_element_ordering_validation_correct():
    """Test that correct element ordering passes validation"""
    
    # Create a minimal LAMMPS data file
    data_content = """Test structure

2 atoms
2 atom types

0 10 xlo xhi
0 10 ylo yhi
0 10 zlo zhi

Masses

1 63.546
2 91.224

Atoms

1 1 0 0 0
2 2 5 5 5
"""
    
    temp_dir = tempfile.mkdtemp()
    data_file = os.path.join(temp_dir, 'test.data')
    yaml_file = os.path.join(temp_dir, 'test.yaml')
    
    with open(data_file, 'w') as f:
        f.write(data_content)
    
    phase_data = {
        'phases': [{
            'phase_name': 'test_phase',
            'element': ['Cu', 'Zr'],
            'mass': [63.546, 91.224],
            'reference_phase': 'solid',
            'lattice': data_file,
            'pair_coeff': ['* * CuZr.eam.alloy Cu Zr'],  # Correct order
            'composition': {
                'reference_element': 'Cu',
                'range': [0.5],
                'interval': 0.1,
                'reference': 0.5
            },
            'temperature': {
                'range': [300],
                'interval': 100
            }
        }]
    }
    
    with open(yaml_file, 'w') as f:
        yaml.safe_dump(phase_data, f)
    
    original_dir = os.getcwd()
    try:
        os.chdir(temp_dir)
        # Should not raise - correct ordering
        prepare_inputs_for_phase_diagram(yaml_file)
        # If it succeeds, output file should exist
        output_file = 'test_phase_test.yaml'
        assert os.path.exists(output_file)
    finally:
        os.chdir(original_dir)
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_phase_no_pair_coeff_skips_validation():
    """Test that phases without pair_coeff skip element ordering validation"""
    
    # Create a minimal LAMMPS data file
    data_content = """Test structure

2 atoms
2 atom types

0 10 xlo xhi
0 10 ylo yhi
0 10 zlo zhi

Masses

1 63.546
2 91.224

Atoms

1 1 0 0 0
2 2 5 5 5
"""
    
    temp_dir = tempfile.mkdtemp()
    data_file = os.path.join(temp_dir, 'test.data')
    yaml_file = os.path.join(temp_dir, 'test.yaml')
    
    with open(data_file, 'w') as f:
        f.write(data_content)
    
    phase_data = {
        'phases': [{
            'phase_name': 'test_phase',
            'element': ['Cu', 'Zr'],
            'mass': [63.546, 91.224],
            'reference_phase': 'solid',
            'lattice': data_file,
            # No pair_coeff - validation should be skipped
            'composition': {
                'reference_element': 'Cu',
                'range': [0.5],
                'interval': 0.1,
                'reference': 0.5
            },
            'temperature': {
                'range': [300],
                'interval': 100
            }
        }]
    }
    
    with open(yaml_file, 'w') as f:
        yaml.safe_dump(phase_data, f)
    
    original_dir = os.getcwd()
    try:
        os.chdir(temp_dir)
        # Should not raise even without pair_coeff
        prepare_inputs_for_phase_diagram(yaml_file)
        output_file = 'test_phase_test.yaml'
        assert os.path.exists(output_file)
    finally:
        os.chdir(original_dir)
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)
