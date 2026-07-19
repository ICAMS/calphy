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


def test_create_composition_array_direct_values():
    """Test that a direct 'values' list bypasses range/interval"""
    from calphy.phase_diagram import _create_composition_array
    
    direct = [0.0, 0.1, 0.3, 0.6, 1.0]  # non-equidistant
    comp_arr, is_reference = _create_composition_array(
        comp_range=None,
        interval=None,
        reference=0.0,
        values=direct
    )
    
    assert len(comp_arr) == 5
    assert list(comp_arr) == direct
    assert is_reference[0] == True   # 0.0 is reference
    assert all(not ref for ref in is_reference[1:])  # others are not reference


def test_create_composition_array_direct_values_reference_in_middle():
    """Test direct values array where reference is not the first entry"""
    from calphy.phase_diagram import _create_composition_array
    
    direct = [0.2, 0.5, 1.0]
    comp_arr, is_reference = _create_composition_array(
        comp_range=None,
        interval=None,
        reference=0.5,
        values=direct
    )
    
    assert len(comp_arr) == 3
    assert is_reference[0] == False
    assert is_reference[1] == True   # 0.5 matches reference
    assert is_reference[2] == False


def test_create_composition_array_values_overrides_range():
    """values kwarg takes priority over range/interval when both are supplied"""
    from calphy.phase_diagram import _create_composition_array
    
    direct = [0.0, 0.5, 1.0]
    comp_arr, is_reference = _create_composition_array(
        comp_range=[0.0, 1.0],
        interval=0.1,          # would give 11 points if used
        reference=0.0,
        values=direct
    )
    
    # Should use direct values (3 points), not range/interval (11 points)
    assert len(comp_arr) == 3
    assert list(comp_arr) == direct


def test_create_temperature_array_direct_values():
    """Test that a direct 'values' list bypasses range/interval"""
    from calphy.phase_diagram import _create_temperature_array

    direct = [1000, 1200, 1500, 1800]  # non-equidistant
    temp_arr = _create_temperature_array(temp_range=None, interval=None, values=direct)

    assert len(temp_arr) == 4
    assert list(temp_arr) == [float(t) for t in direct]


def test_create_temperature_array_values_overrides_range():
    """values kwarg takes priority over range/interval when both are supplied"""
    from calphy.phase_diagram import _create_temperature_array

    direct = [1300, 1600, 1800]
    temp_arr = _create_temperature_array(
        temp_range=[1300, 1800],
        interval=100,          # would give 6 points if used
        values=direct
    )

    # Should use direct values (3 points), not range/interval (6 points)
    assert len(temp_arr) == 3
    assert list(temp_arr) == [float(t) for t in direct]


def test_create_temperature_array_range_unchanged():
    """Existing range/interval behaviour is unaffected when values is not given"""
    from calphy.phase_diagram import _create_temperature_array

    temp_arr = _create_temperature_array(temp_range=[1300, 1500], interval=100)

    assert len(temp_arr) == 3  # 1300, 1400, 1500
    assert temp_arr[0] == 1300.0
    assert temp_arr[-1] == 1500.0


# ---------------------------------------------------------------------------
# Helpers for PhaseDiagram to_pickle / from_pickle tests
# ---------------------------------------------------------------------------

def _make_phase_diagram(calculated=True):
    """
    Build a PhaseDiagram instance without touching the filesystem by
    bypassing __init__ and setting attributes directly.

    Parameters
    ----------
    calculated : bool
        If True, populate the post-calculate() attributes with dummy data.
        If False, leave them as None (pre-calculate state).
    """
    import pandas as _pd
    from calphy.phase_diagram import PhaseDiagram

    obj = object.__new__(PhaseDiagram)
    obj.reference_element = "Ag"
    obj.phases = ["cufcc", "agfcc", "lqd"]
    obj.composition_intervals = {
        "cufcc": (0.0, 0.5),
        "agfcc": (0.5, 1.0),
        "lqd":   (0.0, 1.0),
    }
    obj.df = _pd.DataFrame({
        "phase": ["cufcc", "agfcc", "lqd"],
        "composition": [0.2, 0.8, 0.5],
        "temperature": [1000, 1000, 1000],
        "free_energy": [-1.0, -1.2, -0.9],
        "status": ["True", "True", "True"],
    })

    if calculated:
        obj.tangents = [[[0.1, 0.4]], [[0.2, 0.6]]]
        obj.temperatures = [800, 900]
        obj.tangent_types = [["cufcc-lqd"], ["agfcc-lqd"]]
        obj._calc_kwargs = {
            "fit_order": 4,
            "method": "polynomial",
            "boundary_trim": 0.1,
            "remove_self_tangents_for": [],
            "ideal_configurational_entropy": True,
            "end_weight": 3,
            "end_indices": 4,
        }
    else:
        obj.tangents = None
        obj.temperatures = None
        obj.tangent_types = None
        obj._calc_kwargs = {}

    return obj


# ---------------------------------------------------------------------------
# PhaseDiagram.to_pickle / PhaseDiagram.from_pickle tests
# ---------------------------------------------------------------------------

def test_phase_diagram_save_creates_file(tmp_path):
    """to_pickle() should create a file at the given path."""
    from calphy.phase_diagram import PhaseDiagram

    pd_obj = _make_phase_diagram(calculated=True)
    out = tmp_path / "pd.pkl"
    pd_obj.to_pickle(str(out))
    assert out.exists()
    assert out.stat().st_size > 0


def test_phase_diagram_load_restores_calculated(tmp_path):
    """from_pickle() should restore a fully-calculated PhaseDiagram exactly."""
    from calphy.phase_diagram import PhaseDiagram

    pd_obj = _make_phase_diagram(calculated=True)
    out = str(tmp_path / "pd.pkl")
    pd_obj.to_pickle(out)

    loaded = PhaseDiagram.from_pickle(out)

    assert isinstance(loaded, PhaseDiagram)
    assert loaded.reference_element == pd_obj.reference_element
    assert loaded.phases == pd_obj.phases
    assert loaded.composition_intervals == pd_obj.composition_intervals
    assert loaded.tangents == pd_obj.tangents
    assert loaded.temperatures == pd_obj.temperatures
    assert loaded.tangent_types == pd_obj.tangent_types
    assert loaded._calc_kwargs == pd_obj._calc_kwargs
    assert loaded.df.shape == pd_obj.df.shape


def test_phase_diagram_save_load_pre_calculate(tmp_path):
    """to_pickle/from_pickle round-trip works even when calculate() has not been called."""
    from calphy.phase_diagram import PhaseDiagram

    pd_obj = _make_phase_diagram(calculated=False)
    out = str(tmp_path / "pd_pre.pkl")
    pd_obj.to_pickle(out)

    loaded = PhaseDiagram.from_pickle(out)

    assert loaded.tangents is None
    assert loaded.temperatures is None
    assert loaded.tangent_types is None
    assert loaded.phases == pd_obj.phases
    assert loaded.reference_element == pd_obj.reference_element


def test_phase_diagram_load_missing_file(tmp_path):
    """from_pickle() should raise FileNotFoundError for a non-existent path."""
    from calphy.phase_diagram import PhaseDiagram

    with pytest.raises(FileNotFoundError):
        PhaseDiagram.from_pickle(str(tmp_path / "does_not_exist.pkl"))


def test_phase_diagram_repr_after_load(tmp_path):
    """__repr__ should reflect the correct state after loading."""
    from calphy.phase_diagram import PhaseDiagram

    pd_calc = _make_phase_diagram(calculated=True)
    pd_uncalc = _make_phase_diagram(calculated=False)

    out_calc = str(tmp_path / "calc.pkl")
    out_uncalc = str(tmp_path / "uncalc.pkl")
    pd_calc.to_pickle(out_calc)
    pd_uncalc.to_pickle(out_uncalc)

    loaded_calc = PhaseDiagram.from_pickle(out_calc)
    loaded_uncalc = PhaseDiagram.from_pickle(out_uncalc)

    assert "calculated" in repr(loaded_calc)
    assert "not calculated" in repr(loaded_uncalc)


# ---------------------------------------------------------------------------
# Helpers for to_tdb / from_tdb tests
# ---------------------------------------------------------------------------

def _make_tdb_phase_diagram():
    """
    Build a minimal PhaseDiagram with enough synthetic Au-Cu data for
    :meth:`PhaseDiagram.build_calphad_surface` (and hence
    :meth:`to_tdb`) to fit without touching real MD output.

    Two full-range phases (fcc and lqd) plus one limited-range ordered
    phase (aucu in the window x_Cu ∈ [0.45, 0.55]), each populated with
    multi-T G(x, T) arrays at a handful of compositions including the
    pure endpoints required by ``build_calphad_surface``.
    """
    import numpy as np
    import pandas as _pd
    from calphy.phase_diagram import PhaseDiagram

    obj = object.__new__(PhaseDiagram)
    obj.reference_element = "Cu"
    obj.phases = ["fcc", "lqd", "aucu"]
    obj.composition_intervals = {
        "fcc": (0.0, 1.0),
        "lqd": (0.0, 1.0),
        "aucu": (0.45, 0.55),
    }

    T = np.linspace(400.0, 1400.0, 11)
    rows = []

    # Full-range phases: pure endpoints + several interior compositions.
    # G(x, T) chosen to give:
    #   - fcc more stable at low T
    #   - lqd more stable at high T (crosses fcc near ~1200 K at endpoints)
    #   - a mild attractive RK excess so the fit produces non-trivial L_k.
    for ph, base, slope, l0 in (
        ("fcc", -4.5, -3.0e-4, -0.05),
        ("lqd", -4.4, -4.0e-4, -0.20),
    ):
        for x in [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]:
            G = base + slope * T + l0 * x * (1.0 - x)
            rows.append({
                "phase": ph,
                "composition": float(x),
                "temperature": T.copy(),
                "free_energy": G.copy(),
                "status": "True",
                "is_reference": x in (0.0, 1.0),
            })

    # Limited-range ordered phase: data only inside its composition window.
    # Slightly more stable than fcc at x ≈ 0.5 to give the bounded fit
    # something to anchor onto.
    for x in [0.45, 0.475, 0.5, 0.525, 0.55]:
        G = -4.5 - 3.0e-4 * T - 0.07 * x * (1.0 - x)
        rows.append({
            "phase": "aucu",
            "composition": float(x),
            "temperature": T.copy(),
            "free_energy": G.copy(),
            "status": "True",
            "is_reference": False,
        })

    obj.df = _pd.DataFrame(rows)
    obj.tangents = None
    obj.temperatures = None
    obj.tangent_types = None
    obj._calc_kwargs = {}
    obj._calphad_surfaces = {}
    return obj


def test_phase_diagram_to_tdb_writes_sgte_structure(tmp_path):
    """to_tdb emits the expected SGTE header, FUNCTION and PHASE sections."""
    from calphy.phase_diagram import PhaseDiagram

    pd_obj = _make_tdb_phase_diagram()
    out = tmp_path / "synthetic.tdb"
    pd_obj.to_tdb(str(out), elements=("Au", "Cu"))

    assert out.exists() and out.stat().st_size > 0
    text = out.read_text()

    # Header conventions
    assert "TEMP_LIM" in text
    assert "DEFINE_SYSTEM_DEFAULT" in text
    assert "TYPE_DEFINITION % SEQ" in text
    assert "$ CALPHY_TDB_METADATA" in text

    # GHSER FUNCTION blocks emitted once per element
    assert "FUNCTION GHSERAU" in text
    assert "FUNCTION GHSERCU" in text

    # Phase declarations: solid solutions use (M):VA, liquid stays 1-sub
    assert "PHASE FCC % 2 1.0 1.0" in text
    assert "PHASE LQD % 1 1.0" in text
    assert "PHASE AUCU % 2 1.0 1.0" in text

    # Pure-element G of the host phase references GHSER
    assert "+GHSERAU#" in text
    assert "+GHSERCU#" in text


def test_phase_diagram_from_tdb_roundtrips_metadata(tmp_path):
    """from_tdb recovers phases, composition intervals and surfaces."""
    from calphy.phase_diagram import PhaseDiagram

    pd_obj = _make_tdb_phase_diagram()
    out = tmp_path / "synthetic.tdb"
    pd_obj.to_tdb(str(out), elements=("Au", "Cu"))

    loaded = PhaseDiagram.from_tdb(str(out))

    assert loaded.reference_element == "Cu"
    assert sorted(loaded.phases) == ["aucu", "fcc", "lqd"]
    assert loaded.composition_intervals["fcc"] == (0.0, 1.0)
    assert loaded.composition_intervals["lqd"] == (0.0, 1.0)
    assert loaded.composition_intervals["aucu"] == (0.45, 0.55)
    # Surfaces are present for solution + limited-range phases
    assert "fcc" in loaded._calphad_surfaces
    assert "lqd" in loaded._calphad_surfaces
    assert "aucu" in loaded._calphad_surfaces
    # df is empty by design — raw F(x,T) data is not stored in the TDB
    assert loaded.df.shape[0] == 0


def test_phase_diagram_from_tdb_missing_metadata(tmp_path):
    """from_tdb rejects a TDB that lacks the $ CALPHY_TDB_METADATA header."""
    from calphy.phase_diagram import PhaseDiagram

    bare = tmp_path / "bare.tdb"
    bare.write_text(
        "TYPE_DEFINITION % SEQ *!\n"
        "ELEMENT AU FCC_A1 0.0 0.0 0.0 !\n"
        "PHASE FCC % 1 1.0 !\n"
    )
    with pytest.raises(ValueError, match="CALPHY_TDB_METADATA"):
        PhaseDiagram.from_tdb(str(bare))


def test_phase_diagram_from_tdb_missing_file(tmp_path):
    """from_tdb raises FileNotFoundError for a non-existent path."""
    from calphy.phase_diagram import PhaseDiagram

    with pytest.raises(FileNotFoundError):
        PhaseDiagram.from_tdb(str(tmp_path / "does_not_exist.tdb"))


def test_to_tdb_poly6_L_parameters_and_metadata(tmp_path):
    """Default poly6 mode writes full T-dependence in full-range L params."""
    pd_obj = _make_tdb_phase_diagram()
    out = tmp_path / "poly6.tdb"
    pd_obj.to_tdb(str(out), elements=("Au", "Cu"))
    text = out.read_text()

    # metadata records the L form
    import json
    from calphy.phase_diagram import _TDB_METADATA_PREFIX
    meta_line = next(
        l for l in text.splitlines() if l.startswith(_TDB_METADATA_PREFIX)
    )
    meta = json.loads(meta_line[len(_TDB_METADATA_PREFIX):])
    assert meta["L_temperature_form"] == "poly6"
    assert meta["limited_L_n_terms"] == 3
    assert "aucu" in meta["limited_rk_order"]

    # at least one FCC L parameter carries a T*LN(T) term (poly6, not linear)
    fcc_L_lines = [
        l for l in text.splitlines() if l.startswith("PARAMETER L(FCC")
    ]
    assert fcc_L_lines
    assert any("T*LN(T)" in l for l in fcc_L_lines)


def test_to_tdb_linear_mode_and_validation(tmp_path):
    """L_temperature_form='linear' refits a+bT; invalid values raise."""
    pd_obj = _make_tdb_phase_diagram()
    out = tmp_path / "linear.tdb"
    pd_obj.to_tdb(str(out), elements=("Au", "Cu"), L_temperature_form="linear")
    text = out.read_text()
    fcc_L_lines = [
        l for l in text.splitlines() if l.startswith("PARAMETER L(FCC")
    ]
    assert fcc_L_lines
    assert not any("LN(T)" in l for l in fcc_L_lines)

    with pytest.raises(ValueError, match="L_temperature_form"):
        pd_obj.to_tdb(str(out), elements=("Au", "Cu"), L_temperature_form="cubic")


def test_from_tdb_full_range_surfaces_roundtrip_exactly(tmp_path):
    """poly6 mode: from_tdb recovers full-range G(x,T) to write precision."""
    from calphy.phase_diagram import PhaseDiagram, _eval_calphad_surface_at

    pd_obj = _make_tdb_phase_diagram()
    pd_obj.build_calphad_surface(rk_order=3)
    out = tmp_path / "exact.tdb"
    pd_obj.to_tdb(str(out), elements=("Au", "Cu"))

    loaded = PhaseDiagram.from_tdb(str(out))
    xg = np.linspace(0.01, 0.99, 51)
    for ph in ("fcc", "lqd"):
        for T in (450.0, 900.0, 1350.0):
            g0 = _eval_calphad_surface_at(pd_obj._calphad_surfaces[ph], xg, T)
            g1 = _eval_calphad_surface_at(loaded._calphad_surfaces[ph], xg, T)
            assert np.max(np.abs(g1 - g0)) < 1e-5  # eV/atom


def test_limited_range_bounded_fit_convex_well():
    """The bounded fit keeps the ordered-phase well convex (no W bottom)."""
    from calphy.phase_diagram import (
        _phase_data_as_points,
        _fit_limited_range_surface_bounded,
        _eval_calphad_surface_at,
    )

    pd_obj = _make_tdb_phase_diagram()
    pd_obj.build_calphad_surface(rk_order=3)
    host = pd_obj._calphad_surfaces["fcc"]
    points = _phase_data_as_points(pd_obj.df, "aucu")

    fit = _fit_limited_range_surface_bounded(
        points, host, rk_order=4, n_T_terms=3,
        constraint_T_range=(298.15, 1400.0),
    )
    # L_coeffs comes back zero-padded to the 6-term CALPHAD layout
    assert fit["L_coeffs"].shape == (4, 6)
    assert np.allclose(fit["L_coeffs"][:, 3:], 0.0)
    assert fit["constraint_violation_max"] < 1e-4

    xg = np.linspace(0.451, 0.549, 101)
    for T in (500.0, 700.0, 900.0):
        g_ph = _eval_calphad_surface_at(fit, xg, T)
        g_host = _eval_calphad_surface_at(host, xg, T)
        diff = g_ph - g_host
        # single local minimum: sign of the derivative changes at most once
        sign_changes = np.sum(np.abs(np.diff(np.sign(np.diff(diff)))) > 0)
        assert sign_changes <= 1

    with pytest.raises(ValueError, match="n_T_terms"):
        _fit_limited_range_surface_bounded(points, host, n_T_terms=4)


def test_to_tdb_pycalphad_reproduces_surface(tmp_path):
    """pycalphad evaluates the TDB to the same G(x,T) as calphy's surface."""
    pycalphad = pytest.importorskip("pycalphad")
    from pycalphad import Database, calculate
    from calphy.phase_diagram import _eval_calphad_surface_at, EV_TO_J_MOL

    pd_obj = _make_tdb_phase_diagram()
    pd_obj.build_calphad_surface(rk_order=3)
    out = tmp_path / "pycalphad.tdb"
    pd_obj.to_tdb(str(out), elements=("Au", "Cu"))

    db = Database(str(out))
    assert {"FCC", "LQD", "AUCU"} <= set(db.phases)

    xg = np.linspace(0.05, 0.95, 19)
    for ph, tdb_name in (("fcc", "FCC"), ("lqd", "LQD")):
        n_subl = len(db.phases[tdb_name].constituents)
        cols = [1.0 - xg, xg]
        if n_subl == 2:
            cols.append(np.ones_like(xg))
        pts = np.column_stack(cols)
        for T in (500.0, 1000.0, 1350.0):
            res = calculate(
                db, ["AU", "CU", "VA"], tdb_name, T=T, P=101325, N=1,
                points=pts, output="GM",
            )
            g_tdb = res.GM.values.flatten() / EV_TO_J_MOL
            g_surf = _eval_calphad_surface_at(
                pd_obj._calphad_surfaces[ph], xg, T
            )
            assert np.max(np.abs(g_tdb - g_surf)) < 1e-4  # eV/atom

