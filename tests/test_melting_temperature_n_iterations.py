"""
Test that n_iterations is properly preserved in melting_temperature mode.

This test verifies the fix for the bug where melting_temperature calculations
with n_iterations > 1 would fail because n_iterations wasn't passed to the
internal ts calculations.
"""

import pytest
import tempfile
import os
import yaml
from calphy.routines import MeltingTemp
from calphy.input import read_inputfile


def test_melting_temperature_preserves_n_iterations():
    """
    Test that n_iterations is preserved when melting_temperature creates
    internal ts calculations for solid and liquid phases.
    """
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a dummy structure file
        structure_file = os.path.join(tmpdir, "conf.data")
        with open(structure_file, 'w') as f:
            f.write("""LAMMPS data file

1 atoms
1 atom types

0.0 4.05 xlo xhi
0.0 4.05 ylo yhi
0.0 4.05 zlo zhi

Masses

1 26.98

Atoms

1 1 0.0 0.0 0.0
""")
        
        # Create a minimal input file for melting_temperature with n_iterations > 1
        input_data = {
            "calculations": [
                {
                    "mode": "melting_temperature",
                    "lattice": structure_file,
                    "state": "solid",
                    "temperature": 1000,
                    "pressure": 0,
                    "n_iterations": 3,  # This should be preserved
                    "element": ["Al"],
                    "mass": [26.98],
                    "pair_style": ["eam/alloy"],
                    "pair_coeff": ["* * Al.eam.alloy Al"],
                }
            ]
        }
        
        input_file = os.path.join(tmpdir, "input.yaml")
        with open(input_file, "w") as f:
            yaml.safe_dump(input_data, f)
        
        # Read the input file
        calculations = read_inputfile(input_file)
        calc = calculations[0]
        
        # Create MeltingTemp object
        melt = MeltingTemp(calculation=calc, simfolder=tmpdir)
        
        # Prepare the internal calculations
        melt.prepare_calcs()
        
        # Verify that both solid and liquid calculations have n_iterations = 3
        assert len(melt.calculations) == 2, "Should create 2 calculations (solid and liquid)"
        
        solid_calc = melt.calculations[0]
        liquid_calc = melt.calculations[1]
        
        assert solid_calc.mode == "ts", "First calculation should be ts mode"
        assert liquid_calc.mode == "ts", "Second calculation should be ts mode"
        
        assert solid_calc.reference_phase == "solid", "First calculation should be solid"
        assert liquid_calc.reference_phase == "liquid", "Second calculation should be liquid"
        
        # Critical assertion: n_iterations should be preserved
        assert solid_calc.n_iterations == 3, f"Solid calc should have n_iterations=3, got {solid_calc.n_iterations}"
        assert liquid_calc.n_iterations == 3, f"Liquid calc should have n_iterations=3, got {liquid_calc.n_iterations}"


def test_melting_temperature_default_n_iterations():
    """
    Test that when n_iterations is not specified, it defaults to 1.
    """
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a dummy structure file
        structure_file = os.path.join(tmpdir, "conf.data")
        with open(structure_file, 'w') as f:
            f.write("""LAMMPS data file

1 atoms
1 atom types

0.0 4.05 xlo xhi
0.0 4.05 ylo yhi
0.0 4.05 zlo zhi

Masses

1 26.98

Atoms

1 1 0.0 0.0 0.0
""")
        
        # Create input without n_iterations specified
        input_data = {
            "calculations": [
                {
                    "mode": "melting_temperature",
                    "lattice": structure_file,
                    "state": "solid",
                    "temperature": 1000,
                    "pressure": 0,
                    # n_iterations not specified - should default to 1
                    "element": ["Al"],
                    "mass": [26.98],
                    "pair_style": ["eam/alloy"],
                    "pair_coeff": ["* * Al.eam.alloy Al"],
                }
            ]
        }
        
        input_file = os.path.join(tmpdir, "input.yaml")
        with open(input_file, "w") as f:
            yaml.safe_dump(input_data, f)
        
        calculations = read_inputfile(input_file)
        calc = calculations[0]
        
        melt = MeltingTemp(calculation=calc, simfolder=tmpdir)
        melt.prepare_calcs()
        
        # Should default to 1
        assert melt.calculations[0].n_iterations == 1
        assert melt.calculations[1].n_iterations == 1


def test_melting_temperature_n_iterations_values():
    """
    Test various values of n_iterations to ensure they're all preserved correctly.
    """
    
    for n_iter in [1, 2, 5, 10]:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a dummy structure file
            structure_file = os.path.join(tmpdir, "conf.data")
            with open(structure_file, 'w') as f:
                f.write("""LAMMPS data file

1 atoms
1 atom types

0.0 4.05 xlo xhi
0.0 4.05 ylo yhi
0.0 4.05 zlo zhi

Masses

1 26.98

Atoms

1 1 0.0 0.0 0.0
""")
            
            input_data = {
                "calculations": [
                    {
                        "mode": "melting_temperature",
                        "lattice": structure_file,
                        "state": "solid",
                        "temperature": 1000,
                        "pressure": 0,
                        "n_iterations": n_iter,
                        "element": ["Al"],
                        "mass": [26.98],
                        "pair_style": ["eam/alloy"],
                        "pair_coeff": ["* * Al.eam.alloy Al"],
                    }
                ]
            }
            
            input_file = os.path.join(tmpdir, "input.yaml")
            with open(input_file, "w") as f:
                yaml.safe_dump(input_data, f)
            
            calculations = read_inputfile(input_file)
            calc = calculations[0]
            
            melt = MeltingTemp(calculation=calc, simfolder=tmpdir)
            melt.prepare_calcs()
            
            assert melt.calculations[0].n_iterations == n_iter, \
                f"Solid calc should have n_iterations={n_iter}"
            assert melt.calculations[1].n_iterations == n_iter, \
                f"Liquid calc should have n_iterations={n_iter}"
