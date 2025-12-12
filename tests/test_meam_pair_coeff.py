"""
Tests for MEAM pair_coeff handling in composition transformation
"""

import pytest
from calphy.composition_transformation import CompositionTransformation


def test_meam_two_file_format():
    """Test that MEAM potentials with two files (library + parameter) are handled correctly"""

    # We'll test update_pair_coeff directly by creating a minimal mock
    # that has the necessary attributes but doesn't require file I/O

    class MinimalMock:
        def __init__(self):
            self.element = ["Si", "Al"]
            # For transformation from 50-50 to 30-70 Si-Al
            # Old: Si Al (equal amounts, one of each)
            # New: Si Si Al Al Al (30% Si, 70% Al - roughly)
            self.pair_list_old = ["Si", "Al"]
            self.pair_list_new = ["Si", "Al", "Al"]  # More Al in output

    comp = MinimalMock()
    comp.__class__.iselement = CompositionTransformation.iselement
    comp.__class__.update_pair_coeff = CompositionTransformation.update_pair_coeff

    # MEAM format: pair_coeff * * library.meam El1 El2 parameter.meam El1 El2
    pair_coeff = "* * ../potentials/library.meam Si Al ../potentials/si.meam Si Al"

    pc_old, pc_new = comp.update_pair_coeff(pair_coeff)

    print(f"\nOld: {pc_old}")
    print(f"New: {pc_new}")

    # Both element lists should be updated
    # The files should be preserved
    assert "../potentials/library.meam" in pc_old
    assert "../potentials/si.meam" in pc_old
    assert "../potentials/library.meam" in pc_new
    assert "../potentials/si.meam" in pc_new

    # Extract and verify element lists
    old_parts = pc_old.split()
    new_parts = pc_new.split()

    # Find file positions
    lib_idx_old = old_parts.index("../potentials/library.meam")
    param_idx_old = old_parts.index("../potentials/si.meam")
    lib_idx_new = new_parts.index("../potentials/library.meam")
    param_idx_new = new_parts.index("../potentials/si.meam")

    # Elements between library file and parameter file
    elements_after_lib_old = old_parts[lib_idx_old + 1 : param_idx_old]
    elements_after_lib_new = new_parts[lib_idx_new + 1 : param_idx_new]

    # Elements after parameter file (should be at end)
    elements_after_param_old = old_parts[param_idx_old + 1 :]
    elements_after_param_new = new_parts[param_idx_new + 1 :]

    # Both element lists in each command should match
    assert (
        elements_after_lib_old == elements_after_param_old
    ), f"Old command has mismatched element lists: {elements_after_lib_old} vs {elements_after_param_old}"
    assert (
        elements_after_lib_new == elements_after_param_new
    ), f"New command has mismatched element lists: {elements_after_lib_new} vs {elements_after_param_new}"


def test_regular_pair_coeff_still_works():
    """Ensure regular (single-file) pair_coeff still works after MEAM fix"""

    class MinimalMock:
        def __init__(self):
            self.element = ["Cu", "Zr"]
            self.pair_list_old = ["Cu", "Zr"]
            self.pair_list_new = ["Cu", "Zr", "Zr"]  # More Zr

    comp = MinimalMock()
    comp.__class__.iselement = CompositionTransformation.iselement
    comp.__class__.update_pair_coeff = CompositionTransformation.update_pair_coeff

    # Regular EAM format: pair_coeff * * potential.eam.alloy El1 El2
    pair_coeff = "* * CuZr.eam.alloy Cu Zr"

    pc_old, pc_new = comp.update_pair_coeff(pair_coeff)

    print(f"\nOld: {pc_old}")
    print(f"New: {pc_new}")

    assert "CuZr.eam.alloy" in pc_old
    assert "CuZr.eam.alloy" in pc_new

    # Should have element specification
    assert "Cu" in pc_old and "Zr" in pc_old
    assert "Cu" in pc_new and "Zr" in pc_new


def test_meam_three_elements():
    """Test MEAM with three elements"""

    class MinimalMock:
        def __init__(self):
            self.element = ["Al", "Si", "Mg"]
            self.pair_list_old = ["Al", "Si", "Mg"]
            self.pair_list_new = ["Al", "Al", "Si", "Mg"]  # More Al

    comp = MinimalMock()
    comp.__class__.iselement = CompositionTransformation.iselement
    comp.__class__.update_pair_coeff = CompositionTransformation.update_pair_coeff

    pair_coeff = "* * library.meam Al Si Mg param.meam Al Si Mg"

    pc_old, pc_new = comp.update_pair_coeff(pair_coeff)

    print(f"\nOld: {pc_old}")
    print(f"New: {pc_new}")

    # Verify both files are present
    assert "library.meam" in pc_old and "library.meam" in pc_new
    assert "param.meam" in pc_old and "param.meam" in pc_new

    # Verify elements appear in both positions
    old_parts = pc_old.split()
    new_parts = pc_new.split()

    lib_idx_old = old_parts.index("library.meam")
    param_idx_old = old_parts.index("param.meam")

    elements_after_lib = old_parts[lib_idx_old + 1 : param_idx_old]
    elements_after_param = old_parts[param_idx_old + 1 :]

    # Both element lists should match in the old command
    assert elements_after_lib == elements_after_param
    assert len(elements_after_lib) == 3
