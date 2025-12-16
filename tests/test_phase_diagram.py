import pytest
import numpy as np
from calphy.phase_diagram import read_structure_composition


def test_read_structure_composition_single_element():
    """Test reading a structure with only one element type"""
    # Use existing test structure file
    structure_file = "tests/conf1.data"
    elements = ["Cu", "Zr"]

    comp = read_structure_composition(structure_file, elements)

    # Should have Cu atoms and 0 Zr atoms
    assert isinstance(comp, dict)
    assert "Cu" in comp
    assert "Zr" in comp
    assert comp["Cu"] > 0
    assert comp["Zr"] == 0
    assert sum(comp.values()) > 0


def test_read_structure_composition_element_not_in_structure():
    """Test that elements not present in structure get count 0"""
    structure_file = "tests/conf1.data"
    elements = ["Cu", "Zr", "Al", "Ni"]

    comp = read_structure_composition(structure_file, elements)

    # Should have all elements with appropriate counts
    assert len(comp) == 4
    assert comp["Cu"] > 0
    assert comp["Zr"] == 0
    assert comp["Al"] == 0
    assert comp["Ni"] == 0


def test_read_structure_composition_element_order_matters():
    """Test that element list order determines type mapping"""
    structure_file = "tests/conf1.data"

    # Different element orders should give different results
    # because element[0] maps to LAMMPS type 1, element[1] to type 2, etc.
    elements1 = ["Cu", "Zr"]
    elements2 = ["Zr", "Cu"]

    comp1 = read_structure_composition(structure_file, elements1)
    comp2 = read_structure_composition(structure_file, elements2)

    # With Cu first, Cu should have the atoms
    assert comp1["Cu"] > 0
    assert comp1["Zr"] == 0

    # With Zr first, Zr should have the atoms (because type 1 maps to first element)
    assert comp2["Zr"] > 0
    assert comp2["Cu"] == 0


def test_read_structure_composition_total_atoms():
    """Test that total atom count is preserved regardless of element order"""
    structure_file = "tests/conf1.data"

    elements1 = ["Cu", "Zr"]
    elements2 = ["Zr", "Cu"]
    elements3 = ["Cu", "Zr", "Al"]

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
        read_structure_composition("nonexistent_file.data", ["Cu", "Zr"])


def test_read_structure_composition_empty_element_list():
    """Test behavior with empty element list"""
    structure_file = "tests/conf1.data"
    elements = []

    comp = read_structure_composition(structure_file, elements)

    # Should return empty dict
    assert isinstance(comp, dict)
    assert len(comp) == 0


def test_read_structure_composition_single_element_list():
    """Test with single element in list"""
    structure_file = "tests/conf1.data"
    elements = ["Cu"]

    comp = read_structure_composition(structure_file, elements)

    assert len(comp) == 1
    assert "Cu" in comp
    assert comp["Cu"] > 0


def test_composition_transformation_100_percent():
    """Test that 100% composition transformation is allowed (pure phase endpoint)"""
    from calphy.composition_transformation import CompositionTransformation

    # Create a test calculation that transforms all atoms from one element to another
    class TestCalc:
        def __init__(self):
            self.lattice = "tests/conf1.data"
            self.element = ["Cu", "Al"]
            self.composition_scaling = type(
                "obj",
                (object,),
                {
                    "_input_chemical_composition": {"Cu": 500, "Al": 0},
                    "_output_chemical_composition": {"Cu": 0, "Al": 500},
                    "input_chemical_composition": property(
                        lambda s: s._input_chemical_composition
                    ),
                    "output_chemical_composition": property(
                        lambda s: s._output_chemical_composition
                    ),
                    "restrictions": [],
                },
            )()

    calc = TestCalc()

    # Should not raise an error
    comp = CompositionTransformation(calc)

    # Verify the transformation is set up correctly
    assert "Cu" in comp.to_remove
    assert "Al" in comp.to_add
    assert comp.to_remove["Cu"] == 500
    assert comp.to_add["Al"] == 500
    assert len(comp.transformation_list) == 1
    assert comp.transformation_list[0]["primary_element"] == "Cu"
    assert comp.transformation_list[0]["secondary_element"] == "Al"
    assert comp.transformation_list[0]["count"] == 500


def test_composition_transformation_single_atom_type_pair_coeff():
    """Test that pair_coeff matches actual atom types in structure.

    Pure Cu structure (1 atom type) should generate pair_coeff with 1 element mapping,
    regardless of calc.element count. This ensures LAMMPS consistency: the number of
    elements in pair_coeff must match the number of atom types in the data file header.
    """
    from calphy.composition_transformation import CompositionTransformation

    # Create a test calculation - pure Cu structure transforming to Al
    class TestCalc:
        def __init__(self):
            self.lattice = "tests/conf1.data"
            self.element = ["Cu", "Al"]  # 2 elements in config
            self.composition_scaling = type(
                "obj",
                (object,),
                {
                    "_input_chemical_composition": {"Cu": 500, "Al": 0},
                    "_output_chemical_composition": {"Cu": 0, "Al": 500},
                    "input_chemical_composition": property(
                        lambda s: s._input_chemical_composition
                    ),
                    "output_chemical_composition": property(
                        lambda s: s._output_chemical_composition
                    ),
                    "restrictions": [],
                },
            )()

    calc = TestCalc()
    comp = CompositionTransformation(calc)

    # Structure has 1 actual atom type (pure Cu), so pair_coeff has 1 element
    # This matches the LAMMPS data file which declares 1 atom type
    assert (
        len(comp.pair_list_old) == 1
    ), f"Expected 1 element in pair_list_old, got {len(comp.pair_list_old)}: {comp.pair_list_old}"
    assert (
        len(comp.pair_list_new) == 1
    ), f"Expected 1 element in pair_list_new, got {len(comp.pair_list_new)}: {comp.pair_list_new}"
    assert comp.pair_list_old == ["Cu"]
    assert comp.pair_list_new == ["Al"]


def test_create_composition_array_single_non_reference():
    """Test that single composition value different from reference is not marked as reference"""
    from calphy.phase_diagram import _create_composition_array
    
    # Test: reference=0.0, single comp=1.0 (should NOT be marked as reference)
    comp_arr, is_reference = _create_composition_array(
        comp_range=1.0,
        interval=0.1,
        reference=0.0
    )
    
    assert len(comp_arr) == 1
    assert comp_arr[0] == 1.0
    assert is_reference[0] == False, "Single composition at 1.0 should not be marked as reference when reference=0.0"


def test_create_composition_array_single_is_reference():
    """Test that single composition value equal to reference IS marked as reference"""
    from calphy.phase_diagram import _create_composition_array
    
    # Test: reference=0.0, single comp=0.0 (SHOULD be marked as reference)
    comp_arr, is_reference = _create_composition_array(
        comp_range=0.0,
        interval=0.1,
        reference=0.0
    )
    
    assert len(comp_arr) == 1
    assert comp_arr[0] == 0.0
    assert is_reference[0] == True, "Single composition at 0.0 should be marked as reference when reference=0.0"


def test_create_composition_array_range():
    """Test composition array creation with range"""
    from calphy.phase_diagram import _create_composition_array
    
    # Test: reference=0.0, range [0.0, 1.0], interval=0.25
    comp_arr, is_reference = _create_composition_array(
        comp_range=[0.0, 1.0],
        interval=0.25,
        reference=0.0
    )
    
    assert len(comp_arr) == 5  # 0.0, 0.25, 0.5, 0.75, 1.0
    assert is_reference[0] == True  # 0.0 is reference
    assert all(not ref for ref in is_reference[1:])  # Others are not reference

