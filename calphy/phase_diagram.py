import numpy as np
from tqdm.notebook import trange
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import itertools
from itertools import combinations
import math
import copy
import os
import pickle
from calphy.composition_transformation import CompositionTransformation
import yaml
import matplotlib.patches as mpatches

from calphy.integrators import kb

from scipy.spatial import ConvexHull
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter


colors = [
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
]

matcolors = {
    "amber": {
        50: "#fff8e1",
        100: "#ffecb3",
        200: "#ffe082",
        300: "#ffd54f",
        400: "#ffca28",
        500: "#ffc107",
        600: "#ffb300",
        700: "#ffa000",
        800: "#ff8f00",
        900: "#ff6f00",
    },
    "blue_grey": {
        50: "#ECEFF1",
        100: "#CFD8DC",
        200: "#B0BEC5",
        300: "#90A4AE",
        400: "#78909C",
        500: "#607D8B",
        600: "#546E7A",
        700: "#455A64",
        800: "#37474F",
        900: "#263238",
    },
    "blue": {
        50: "#E3F2FD",
        100: "#BBDEFB",
        200: "#90CAF9",
        300: "#64B5F6",
        400: "#42A5F5",
        500: "#2196F3",
        600: "#1E88E5",
        700: "#1976D2",
        800: "#1565C0",
        900: "#0D47A1",
    },
    "brown": {
        50: "#EFEBE9",
        100: "#D7CCC8",
        200: "#BCAAA4",
        300: "#A1887F",
        400: "#8D6E63",
        500: "#795548",
        600: "#6D4C41",
        700: "#5D4037",
        800: "#4E342E",
        900: "#3E2723",
    },
    "cyan": {
        50: "#E0F7FA",
        100: "#B2EBF2",
        200: "#80DEEA",
        300: "#4DD0E1",
        400: "#26C6DA",
        500: "#00BCD4",
        600: "#00ACC1",
        700: "#0097A7",
        800: "#00838F",
        900: "#006064",
    },
    "deep_orange": {
        50: "#FBE9E7",
        100: "#FFCCBC",
        200: "#FFAB91",
        300: "#FF8A65",
        400: "#FF7043",
        500: "#FF5722",
        600: "#F4511E",
        700: "#E64A19",
        800: "#D84315",
        900: "#BF360C",
    },
    "deep_purple": {
        50: "#EDE7F6",
        100: "#D1C4E9",
        200: "#B39DDB",
        300: "#9575CD",
        400: "#7E57C2",
        500: "#673AB7",
        600: "#5E35B1",
        700: "#512DA8",
        800: "#4527A0",
        900: "#311B92",
    },
    "green": {
        50: "#E8F5E9",
        100: "#C8E6C9",
        200: "#A5D6A7",
        300: "#81C784",
        400: "#66BB6A",
        500: "#4CAF50",
        600: "#43A047",
        700: "#388E3C",
        800: "#2E7D32",
        900: "#1B5E20",
    },
    "grey": {
        50: "#FAFAFA",
        100: "#F5F5F5",
        200: "#EEEEEE",
        300: "#E0E0E0",
        400: "#BDBDBD",
        500: "#9E9E9E",
        600: "#757575",
        700: "#616161",
        800: "#424242",
        900: "#212121",
    },
    "indigo": {
        50: "#E8EAF6",
        100: "#C5CAE9",
        200: "#9FA8DA",
        300: "#7986CB",
        400: "#5C6BC0",
        500: "#3F51B5",
        600: "#3949AB",
        700: "#303F9F",
        800: "#283593",
        900: "#1A237E",
    },
    "light_blue": {
        50: "#E1F5FE",
        100: "#B3E5FC",
        200: "#81D4FA",
        300: "#4FC3F7",
        400: "#29B6F6",
        500: "#03A9F4",
        600: "#039BE5",
        700: "#0288D1",
        800: "#0277BD",
        900: "#01579B",
    },
    "light_green": {
        50: "#F1F8E9",
        100: "#DCEDC8",
        200: "#C5E1A5",
        300: "#AED581",
        400: "#9CCC65",
        500: "#8BC34A",
        600: "#7CB342",
        700: "#689F38",
        800: "#558B2F",
        900: "#33691E",
    },
    "lime": {
        50: "#F9FBE7",
        100: "#F0F4C3",
        200: "#E6EE9C",
        300: "#DCE775",
        400: "#D4E157",
        500: "#CDDC39",
        600: "#C0CA33",
        700: "#AFB42B",
        800: "#9E9D24",
        900: "#827717",
    },
    "orange": {
        50: "#FFF3E0",
        100: "#FFE0B2",
        200: "#FFCC80",
        300: "#FFB74D",
        400: "#FFA726",
        500: "#FF9800",
        600: "#FB8C00",
        700: "#F57C00",
        800: "#EF6C00",
        900: "#E65100",
    },
    "pink": {
        50: "#FCE4EC",
        100: "#F8BBD0",
        200: "#F48FB1",
        300: "#F06292",
        400: "#EC407A",
        500: "#E91E63",
        600: "#D81B60",
        700: "#C2185B",
        800: "#AD1457",
        900: "#880E4F",
    },
    "purple": {
        50: "#F3E5F5",
        100: "#E1BEE7",
        200: "#CE93D8",
        300: "#BA68C8",
        400: "#AB47BC",
        500: "#9C27B0",
        600: "#8E24AA",
        700: "#7B1FA2",
        800: "#6A1B9A",
        900: "#4A148C",
    },
    "red": {
        50: "#FFEBEE",
        100: "#FFCDD2",
        200: "#EF9A9A",
        300: "#E57373",
        500: "#F44336",
        600: "#E53935",
        700: "#D32F2F",
        800: "#C62828",
        900: "#B71C1C",
    },
    "teal": {
        50: "#E0F2F1",
        100: "#B2DFDB",
        200: "#80CBC4",
        300: "#4DB6AC",
        400: "#26A69A",
        500: "#009688",
        600: "#00897B",
        700: "#00796B",
        800: "#00695C",
        900: "#004D40",
    },
    "yellow": {
        50: "#FFFDE7",
        100: "#FFF9C4",
        200: "#FFF59D",
        300: "#FFF176",
        400: "#FFEE58",
        500: "#FFEB3B",
        600: "#FDD835",
        700: "#FBC02D",
        800: "#F9A825",
        900: "#F57F17",
    },
}


def read_structure_composition(lattice_file, element_list):
    """
    Read a LAMMPS data file and determine the input chemical composition.

    Parameters
    ----------
    lattice_file : str
        Path to the LAMMPS data file
    element_list : list
        List of element symbols in order (element[0] = type 1, element[1] = type 2, etc.)

    Returns
    -------
    dict
        Dictionary mapping element symbols to atom counts
        Elements not present in the structure will have count 0
    """
    from ase.io import read
    from collections import Counter

    # Read the structure file
    structure = read(lattice_file, format="lammps-data", style="atomic")

    # Get the species/types from the structure
    # ASE reads LAMMPS types as species strings ('1', '2', etc.)
    if "species" in structure.arrays:
        types_in_structure = structure.arrays["species"]
    else:
        # Fallback: get atomic numbers and convert to strings
        types_in_structure = [str(x) for x in structure.get_atomic_numbers()]

    # Count atoms by type
    type_counts = Counter(types_in_structure)

    # Build composition mapping element names to counts
    # element[0] corresponds to LAMMPS type '1', element[1] to type '2', etc.
    input_chemical_composition = {}
    for idx, element in enumerate(element_list):
        lammps_type = str(idx + 1)  # LAMMPS types are 1-indexed
        input_chemical_composition[element] = type_counts.get(lammps_type, 0)

    return input_chemical_composition


# Constants for phase diagram preparation
COMPOSITION_TOLERANCE = 1e-5


def _create_composition_array(comp_range, interval, reference, values=None):
    """
    Create composition array from range specification or a direct list of values.

    Parameters
    ----------
    comp_range : list or scalar
        Composition range [min, max] or single value. Ignored when *values* is given.
    interval : float
        Composition interval. Ignored when *values* is given.
    reference : float
        Reference composition value
    values : list, optional
        Explicit list of compositions to use instead of deriving them from
        *comp_range* and *interval*. Allows non-equidistant spacing.

    Returns
    -------
    tuple
        (comp_arr, is_reference) - composition array and boolean array marking reference compositions
    """
    # Direct values array takes priority over range/interval
    if values is not None:
        comp_arr = np.asarray(values, dtype=float)
        is_reference = np.abs(comp_arr - reference) < COMPOSITION_TOLERANCE
        return comp_arr, is_reference

    # Convert to list if scalar
    if not isinstance(comp_range, list):
        comp_range = [comp_range]

    if len(comp_range) == 2:
        comp_arr = np.arange(comp_range[0], comp_range[-1], interval)
        last_val = comp_range[-1]
        if last_val not in comp_arr:
            comp_arr = np.append(comp_arr, last_val)
        is_reference = np.abs(comp_arr - reference) < COMPOSITION_TOLERANCE
    elif len(comp_range) == 1:
        comp_arr = [comp_range[0]]
        # Check if this single composition equals the reference
        is_reference = [np.abs(comp_range[0] - reference) < COMPOSITION_TOLERANCE]
    else:
        raise ValueError("Composition range should be scalar or list of two values!")

    return comp_arr, is_reference


def _create_temperature_array(temp_range, interval, values=None):
    """
    Create temperature array from range specification or a direct list of values.

    Parameters
    ----------
    temp_range : list or scalar
        Temperature range [min, max] or single value. Ignored when *values* is given.
    interval : float
        Temperature interval. Ignored when *values* is given.
    values : list, optional
        Explicit list of temperatures to use instead of deriving them from
        *temp_range* and *interval*. Allows non-equidistant spacing.

    Returns
    -------
    ndarray
        Temperature array
    """
    # Direct values array takes priority over range/interval
    if values is not None:
        return np.asarray(values, dtype=float)

    # Convert to list if scalar
    if not isinstance(temp_range, list):
        temp_range = [temp_range]

    if len(temp_range) == 2:
        ntemps = int((temp_range[-1] - temp_range[0]) / interval) + 1
        temp_arr = np.linspace(temp_range[0], temp_range[-1], ntemps, endpoint=True)
    elif len(temp_range) == 1:
        temp_arr = [temp_range[0]]
    else:
        raise ValueError("Temperature range should be scalar or list of two values!")

    return temp_arr


def _add_temperature_calculations(calc_dict, temp_arr, all_calculations):
    """
    Helper to add calculations for each temperature point.

    Parameters
    ----------
    calc_dict : dict
        Base calculation dictionary
    temp_arr : array
        Array of temperatures
    all_calculations : list
        List to append calculations to
    """
    for temp in temp_arr:
        calc_for_temp = copy.deepcopy(calc_dict)
        calc_for_temp["temperature"] = int(temp)
        all_calculations.append(calc_for_temp)


def fix_data_file(datafile, nelements):
    """
    Change the atom types keyword in the structure file
    """
    lines = []
    with open(datafile, "r") as fin:
        for line in fin:
            if "atom types" in line:
                lines.append(f"{nelements} atom types\n")
            else:
                lines.append(line)
    outfile = datafile + "mod.data"
    with open(outfile, "w") as fout:
        for line in lines:
            fout.write(line)
    return outfile


class CScale:
    def __init__(self):
        self._input_chemical_composition = None
        self._output_chemical_composition = None
        self.restrictions = []

    @property
    def input_chemical_composition(self):
        return self._input_chemical_composition

    @property
    def output_chemical_composition(self):
        return self._output_chemical_composition


class SimpleCalculation:
    """
    Simple calc class
    """

    def __init__(
        self, lattice, element, input_chemical_composition, output_chemical_composition
    ):
        self.lattice = lattice
        self.element = element
        self.composition_scaling = CScale()
        self.composition_scaling._input_chemical_composition = (
            input_chemical_composition
        )
        self.composition_scaling._output_chemical_composition = (
            output_chemical_composition
        )


def prepare_inputs_for_phase_diagram(inputyamlfile, calculation_base_name=None):
    with open(inputyamlfile, "r") as fin:
        data = yaml.safe_load(fin)

    if calculation_base_name is None:
        calculation_base_name = inputyamlfile

    for phase in data["phases"]:
        # Validate binary system assumption
        n_elements = len(phase["element"])
        if n_elements != 2:
            raise ValueError(
                f"Phase diagram preparation currently supports only binary systems. "
                f"Found {n_elements} elements: {phase['element']}"
            )

        # Validate element ordering consistency with pair_coeff
        # This ensures element[0] -> type 1, element[1] -> type 2
        if "pair_coeff" in phase:
            from calphy.input import _extract_elements_from_pair_coeff

            # pair_coeff can be a list or a string - handle both
            pair_coeff = phase["pair_coeff"]
            if isinstance(pair_coeff, list):
                pair_coeff = pair_coeff[0] if pair_coeff else None
            pair_coeff_elements = _extract_elements_from_pair_coeff(pair_coeff)
            if pair_coeff_elements != phase["element"]:
                raise ValueError(
                    f"Element ordering mismatch for phase '{phase.get('phase_name', 'unnamed')}'!\n"
                    f"Elements in 'element' field: {phase['element']}\n"
                    f"Elements from pair_coeff: {pair_coeff_elements}\n"
                    f"These must match exactly in order (element[0] -> LAMMPS type 1, element[1] -> type 2)."
                )

        phase_reference_state = phase["reference_phase"]
        phase_name = phase["phase_name"]

        comps = phase["composition"]
        reference_element = comps["reference_element"]
        if "use_composition_scaling" in comps.keys():
            use_composition_scaling = bool(comps["use_composition_scaling"])
        else:
            use_composition_scaling = True
        if str(phase_reference_state) == "liquid":
            use_composition_scaling = False

        other_element_list = copy.deepcopy(phase["element"])
        other_element_list.remove(reference_element)
        other_element = other_element_list[0]

        # Auto-calculate reference composition from the lattice file
        # if the user did not provide it explicitly.
        if "reference" not in comps:
            lattice_file = phase["lattice"]
            struct_comp = read_structure_composition(lattice_file, phase["element"])
            total_atoms = sum(struct_comp.values())
            if total_atoms > 0:
                ref_comp = struct_comp.get(reference_element, 0) / total_atoms
            else:
                ref_comp = 0.0
            comps["reference"] = ref_comp

        # Create composition array using helper function
        # A direct 'values' list in the input takes priority over range/interval
        comp_arr, is_reference = _create_composition_array(
            comps.get("range"),
            comps.get("interval"),
            comps["reference"],
            values=comps.get("values"),
        )
        ncomps = len(comp_arr)

        # Create temperature array using helper function
        temps = phase["temperature"]
        # A direct 'values' list in the input takes priority over range/interval
        temp_arr = _create_temperature_array(
            temps.get("range"), temps.get("interval"), values=temps.get("values")
        )
        ntemps = len(temp_arr)

        all_calculations = []

        for count, comp in enumerate(comp_arr):
            # check if ref comp equals given comp
            if is_reference[count]:
                # copy the dict
                calc = copy.deepcopy(phase)

                # pop extra keys which are not needed
                # we dont kick out phase_name
                extra_keys = ["composition", "monte_carlo"]
                for key in extra_keys:
                    _ = calc.pop(key, None)

                # update file if needed
                outfile = fix_data_file(calc["lattice"], len(calc["element"]))

                # add ref phase, needed
                calc["reference_phase"] = phase_reference_state
                calc["reference_composition"] = comps["reference"]
                calc["mode"] = "fe"
                calc["folder_prefix"] = f"{phase_name}-{comp:.2f}"
                calc["lattice"] = outfile

                # Add calculations for each temperature
                _add_temperature_calculations(calc, temp_arr, all_calculations)
            else:
                # off stoichiometric
                # copy the dict
                calc = copy.deepcopy(phase)

                # read the structure file to determine input composition automatically
                input_chemical_composition = read_structure_composition(
                    calc["lattice"], calc["element"]
                )

                # calculate total number of atoms from structure
                n_atoms = sum(input_chemical_composition.values())

                if n_atoms == 0:
                    raise ValueError(
                        f"No atoms found in structure file {calc['lattice']}"
                    )

                # find number of atoms of second species based on target composition
                # we follow the convention that composition is always given with the reference element
                output_chemical_composition = {}
                n_species_b = int(np.round(comp * n_atoms, decimals=0))
                output_chemical_composition[reference_element] = n_species_b

                n_species_a = int(n_atoms - n_species_b)
                output_chemical_composition[other_element] = n_species_a

                # Note: Pure phases (n_species_a == 0 or n_species_b == 0) are allowed
                # Composition transformation can handle 100% replacement

                # good, now we need to write such a structure out; likely better to use working directory for that
                folder_prefix = f"{phase_name}-{comp:.2f}"
                calc["reference_composition"] = comps["reference"]
                # if solid, its very easy; kinda
                # if calc['reference_phase'] == 'solid':
                if use_composition_scaling:
                    # this is solid , and comp scale is turned on
                    # pop extra keys which are not needed
                    # we dont kick out phase_name
                    extra_keys = ["composition", "reference_phase"]
                    for key in extra_keys:
                        _ = calc.pop(key, None)

                    # just submit comp scales
                    # add ref phase, needed
                    calc["mode"] = "composition_scaling"
                    calc["folder_prefix"] = folder_prefix
                    calc["composition_scaling"] = {}
                    calc["composition_scaling"][
                        "output_chemical_composition"
                    ] = output_chemical_composition

                else:
                    # manually create a mixed structure - not that the pair style is always ok :)

                    outfile = os.path.join(
                        os.getcwd(),
                        os.path.basename(calc["lattice"]) + folder_prefix + ".comp.mod",
                    )
                    # print(f'finding comp trf from {input_chemical_composition} to {output_chemical_composition}')
                    # write_structure(calc['lattice'], input_chemical_composition, output_chemical_composition, outfile)

                    simplecalc = SimpleCalculation(
                        calc["lattice"],
                        calc["element"],
                        input_chemical_composition,
                        output_chemical_composition,
                    )
                    compsc = CompositionTransformation(simplecalc)
                    compsc.write_structure(outfile, for_fe_mode=True)

                    # pop extra keys which are not needed
                    # we dont kick out phase name
                    extra_keys = ["composition"]
                    for key in extra_keys:
                        _ = calc.pop(key, None)

                    # add ref phase, needed
                    calc["mode"] = "fe"
                    calc["folder_prefix"] = folder_prefix
                    calc["lattice"] = outfile

                # Add calculations for each temperature
                _add_temperature_calculations(calc, temp_arr, all_calculations)

        # finish and write up the file
        output_data = {"calculations": all_calculations}
        base_name = os.path.basename(calculation_base_name)
        for rep in [".yml", ".yaml"]:
            base_name = base_name.replace(rep, "")

        outfile_phase = phase_name + "_" + base_name + ".yaml"
        with open(outfile_phase, "w") as fout:
            yaml.safe_dump(output_data, fout)
        print(
            f"Total {len(all_calculations)} calculations found for phase {phase_name}, written to {outfile_phase}"
        )


def _get_temp_arg(tarr, temp, threshold=1e-1):
    if tarr is None:
        return None
    arg = np.argsort(np.abs(tarr - temp))[0]

    th = np.abs(tarr - temp)[arg]
    if th > threshold:
        arg = None
    return arg


def _is_val_ok(val):
    if val is None:
        return False
    elif math.isnan(val):
        return False
    else:
        return True


def _get_fe_at_args(arr, args):
    fes = []
    for count, x in enumerate(args):
        if _is_val_ok(x):
            fes.append(arr[count][int(x)])
        else:
            fes.append(None)
    return fes


def _calculate_configurational_entropy(x, correction=0):
    if correction == 0:
        s = np.array(
            [(c * np.log(c) + (1 - c) * np.log(1 - c)) if 1 > c > 0 else 0 for c in x]
        )
    else:
        arg = np.argsort(np.abs(x - correction))[0]
        left_side = x[: arg + 1]
        right_side = x[arg:]

        if len(left_side) > 0:
            left_side = left_side / left_side[-1]
            s_left = np.array(
                [
                    (c * np.log(c) + (1 - c) * np.log(1 - c)) if 1 > c > 0 else 0
                    for c in left_side
                ]
            )

        if len(right_side) > 0:
            right_side = right_side - right_side[0]
            right_side = right_side / right_side[-1]
            s_right = np.array(
                [
                    (c * np.log(c) + (1 - c) * np.log(1 - c)) if 1 > c > 0 else 0
                    for c in right_side
                ]
            )

        if len(left_side) == 0:
            return s_right
        elif len(right_side) == 0:
            return s_left
        else:
            return np.concatenate((s_left, s_right[1:]))
    return -s


def _get_free_energy_fit(
    composition,
    free_energy,
    fit_order=5,
    end_weight=3,
    end_indices=4,
    method="polynomial",
):
    """
    Fit free energy as a function of composition.

    Parameters
    ----------
    composition : array-like
        Composition values (between 0 and 1).
    free_energy : array-like
        Free energy values at each composition.
    fit_order : int
        Number of Redlich-Kister coefficients (method="redlich-kister")
        or polynomial order (method="polynomial").
    end_weight : float
        Extra weight given to endpoint data points.
    end_indices : int
        Number of points at each end to receive extra weight.
    method : str
        "redlich-kister" or "polynomial" (default).

    Returns
    -------
    If method="polynomial": polynomial coefficient array (np.polyfit style).
    If method="redlich-kister": dict with keys "L" (RK coefficients),
        "F0" and "F1" (endpoint free energies), "x0" and "x1" (endpoint
        compositions).
    """
    weights = np.ones_like(free_energy)
    weights[0:end_indices] = end_weight
    weights[-end_indices:] = end_weight

    if method == "polynomial":
        fit = np.polyfit(composition, free_energy, fit_order, w=weights)
        return fit

    # --- Redlich-Kister ---
    x = np.asarray(composition, dtype=float)
    F = np.asarray(free_energy, dtype=float)

    # Endpoint values (use actual data at boundaries)
    x0, x1 = x[0], x[-1]
    F0, F1 = F[0], F[-1]

    # Determine which endpoints are actual pure components (x=0 or x=1).
    # Only at pure-component boundaries should the excess be forced to zero.
    tol = 1e-3
    left_is_pure = (x0 < tol) or (x0 > 1.0 - tol)
    right_is_pure = (x1 < tol) or (x1 > 1.0 - tol)

    # Linear reference between endpoints
    F_lin = F0 + (F1 - F0) * (x - x0) / (x1 - x0)
    F_excess = F - F_lin

    # Shifted coordinate so that x0->0, x1->1
    xi = (x - x0) / (x1 - x0)

    if left_is_pure and right_is_pure:
        # Standard RK: excess vanishes at both endpoints via xi*(1-xi) prefactor
        n_coeffs = fit_order
        prefactor = xi * (1 - xi)
        basis = np.column_stack(
            [prefactor * (1 - 2 * xi) ** k for k in range(n_coeffs)]
        )
    elif right_is_pure:
        # Only pin excess=0 at right endpoint (x1 is pure component)
        # Use (1-xi) prefactor: vanishes at xi=1 but free at xi=0
        n_coeffs = fit_order
        prefactor = 1 - xi
        basis = np.column_stack([prefactor * xi**k for k in range(n_coeffs)])
    elif left_is_pure:
        # Only pin excess=0 at left endpoint (x0 is pure component)
        # Use xi prefactor: vanishes at xi=0 but free at xi=1
        n_coeffs = fit_order
        prefactor = xi
        basis = np.column_stack([prefactor * (1 - xi) ** k for k in range(n_coeffs)])
    else:
        # Neither endpoint is pure — use unconstrained polynomial for excess
        n_coeffs = fit_order
        basis = np.column_stack([xi**k for k in range(n_coeffs)])

    # Weighted least-squares
    W = np.diag(weights)
    L_coeffs, _, _, _ = np.linalg.lstsq(W @ basis, W @ F_excess, rcond=None)

    return {
        "L": L_coeffs,
        "F0": F0,
        "F1": F1,
        "x0": x0,
        "x1": x1,
        "left_is_pure": left_is_pure,
        "right_is_pure": right_is_pure,
    }


def _eval_free_energy_fit(fit, composition):
    """
    Evaluate a free energy fit at the given compositions.

    Parameters
    ----------
    fit : array or dict
        Output of ``_get_free_energy_fit``.
        If array → polynomial (np.polyval).
        If dict  → Redlich-Kister evaluation.
    composition : array-like
        Composition values to evaluate at.

    Returns
    -------
    F : ndarray
        Free energy values.
    """
    if isinstance(fit, dict):
        x = np.asarray(composition, dtype=float)
        x0, x1 = fit["x0"], fit["x1"]
        F0, F1 = fit["F0"], fit["F1"]
        L = fit["L"]
        left_is_pure = fit.get("left_is_pure", True)
        right_is_pure = fit.get("right_is_pure", True)

        F_lin = F0 + (F1 - F0) * (x - x0) / (x1 - x0)
        xi = (x - x0) / (x1 - x0)

        if left_is_pure and right_is_pure:
            prefactor = xi * (1 - xi)
            F_excess = sum(L[k] * prefactor * (1 - 2 * xi) ** k for k in range(len(L)))
        elif right_is_pure:
            prefactor = 1 - xi
            F_excess = sum(L[k] * prefactor * xi**k for k in range(len(L)))
        elif left_is_pure:
            prefactor = xi
            F_excess = sum(L[k] * prefactor * (1 - xi) ** k for k in range(len(L)))
        else:
            F_excess = sum(L[k] * xi**k for k in range(len(L)))

        return F_lin + F_excess
    else:
        return np.polyval(fit, composition)


def get_phase_free_energy(
    df,
    phase,
    temp,
    composition_interval=(0, 1),
    ideal_configurational_entropy=False,
    entropy_correction=0.0,
    fit_order=5,
    composition_grid=10000,
    composition_cutoff=None,
    reset_value=1,
    plot=False,
    end_weight=3,
    end_indices=4,
    method="polynomial",
):
    """
    Get the free energy of a phase as a function of composition.

    Parameters
    ----------
    df: Pandas dataframe
        Dataframe consisting of values from simulation. Should contain at least columns composition, phase, `free_energy` and `temperature`.
        `energy_free` and `temperature` should be arrays of equal length, generally an output from reversible scaling calculation.

    phase: str
        phase for which calculation is to be done. Should be present in `df`.

    temp: float
        temperature at which the free energy curves are to be calculated.

    composition_interval: tuple, optional
        If provided, this composition interval is considered. Default (0, 1)

    ideal_configuration_entropy: bool, optional\
        If True, add the ideal configurational entropy. See Notes. Default False.

    entropy_correction: float, optional.
        The composition of the ordered phase. See Notes. Default None.

    fit_order: int, optional
        Order of the polynomial fit used for fitting free energy as a function of composition. Default 5.

    composition_grid: int, optional
        Number of composition points to be used for fitting. Default 10000.

    composition_cutoff: float, optional
        term for correcting incomplete data. If two consecutive composition values are separated by more than `composition_cutoff`,
        it is reset to `reset_value`. Default None.

    reset_value: float, optional
        see above. Default 1.

    plot: bool, optional
        If True, plot the calculated free energy curves.

    method: str, optional
        Fitting method for F(x). "polynomial" (default) uses np.polyfit;
        "redlich-kister" uses the Redlich-Kister expansion
        F_excess = x(1-x) * sum_k L_k*(1-2x)^k.

    Returns
    -------
    result_dict: dict
        contains keys: "phase", "temperature", "composition", "free_energy", and "entropy".

    Notes
    -----
    To be added
    """
    df_phase = df.loc[df["phase"] == phase]
    # drop Nones
    df_phase = df_phase.sort_values(by="composition")
    df_phase = df_phase[
        (df_phase["composition"] >= composition_interval[0])
        & (df_phase["composition"] <= composition_interval[1])
    ]

    composition = df_phase["composition"].values
    args = df_phase["temperature"].apply(_get_temp_arg, args=(temp,))
    fes = _get_fe_at_args(df_phase["free_energy"].values, args)

    # Track which rows already had ideal entropy applied by fix_composition_scaling
    # and which rows are composition_scaling (the only ones we should touch).
    if "_entropy_corrected" in df_phase.columns:
        ec_flags_all = df_phase["_entropy_corrected"].values
    else:
        ec_flags_all = np.zeros(len(df_phase), dtype=bool)

    if "calculation_mode" in df_phase.columns:
        is_comp_scaling_all = (
            df_phase["calculation_mode"] == "composition_scaling"
        ).values
    else:
        is_comp_scaling_all = np.zeros(len(df_phase), dtype=bool)

    # filter out None values
    composition = np.array(
        [composition[count] for count, x in enumerate(fes) if x is not None]
    )
    ec_flags = np.array(
        [ec_flags_all[count] for count, x in enumerate(fes) if x is not None]
    )
    is_comp_scaling = np.array(
        [is_comp_scaling_all[count] for count, x in enumerate(fes) if x is not None]
    )
    fes = np.array([x for x in fes if x is not None])

    if (len(fes) == 0) or (fes is None):
        warnings.warn("Some temperatures could not be found!")
    elif len(fes) <= fit_order:
        warnings.warn(
            f"Not enough data points ({len(fes)}) for fit order {fit_order} for phase '{phase}'. "
            f"Need at least {fit_order+1} points. Returning None."
        )
        return None
    else:
        entropy_term_arr = (
            kb
            * temp
            * _calculate_configurational_entropy(
                composition, correction=entropy_correction
            )
        )
        # Only modify composition_scaling rows; fe-mode rows (liquid or
        # solid) already contain whatever configurational entropy the
        # simulation sampled, so we never add/remove the analytical term.
        if ideal_configurational_entropy:
            # Add entropy to comp_scaling rows that don't have it yet
            mask = is_comp_scaling & ~ec_flags
            fes = fes.copy()
            fes[mask] -= entropy_term_arr[mask]
            entropy_term = entropy_term_arr
        else:
            # Remove entropy from comp_scaling rows that already have it
            mask = is_comp_scaling & ec_flags
            fes = fes.copy()
            fes[mask] += entropy_term_arr[mask]
            entropy_term = []

        fe_fit = _get_free_energy_fit(
            composition,
            fes,
            fit_order=fit_order,
            end_weight=end_weight,
            end_indices=end_indices,
            method=method,
        )
        # Use the requested composition_interval for the evaluation grid
        # so the range is consistent even when some endpoint data is
        # missing at certain temperatures.
        comp_lo = (
            composition_interval[0]
            if composition_interval is not None
            else np.min(composition)
        )
        comp_hi = (
            composition_interval[1]
            if composition_interval is not None
            else np.max(composition)
        )
        compfine = np.linspace(comp_lo, comp_hi, composition_grid)

        # now fit on the comp grid again
        fe = _eval_free_energy_fit(fe_fit, compfine)

        if composition_cutoff is not None:
            distances = [np.min(np.abs(c - composition)) for c in compfine]
            filters = [
                x for x in range(len(distances)) if distances[x] > composition_cutoff
            ]
            fe[filters] = reset_value

        if plot:
            plt.scatter(composition, fes, s=4, label=f"{phase}-calc.", color="#e57373")
            plt.plot(compfine, fe, label=f"{phase}-fit", color="#b71c1c")
            plt.xlabel("x")
            plt.ylabel("F (eV/atom)")
            plt.legend()

        return {
            "phase": phase,
            "temperature": temp,
            "composition": compfine,
            "free_energy": fe,
            "entropy": entropy_term,
            "raw_composition": composition,
            "raw_free_energy": fes,
        }
    return None


def get_free_energy_mixing(dict_list, threshold=1e-3, boundary_trim=0.1):
    """
    Input is a list of dictionaries

    Get free energy of mixing by subtracting end member values.
    End members are chosen automatically.

    Parameters
    ----------
    dict_list : list of dict
        Phase free-energy dictionaries (output of ``get_phase_free_energy``).
    threshold : float
        Tolerance for matching end-member compositions (default 1e-3).
    boundary_trim : float
        Composition width to trim from the boundaries of partial-range
        phases.  A partial-range phase is one whose composition range
        does not reach the global minimum or maximum.  Trimming removes
        the edge region where the global linear reference can produce
        artefactual dips in F_mix.  Set to ``0`` to disable.
    """
    dict_list = np.atleast_1d(dict_list)

    dict_list = np.array([dct for dct in dict_list if dct is not None])

    # we have to get min_comp from all possible values
    min_comp = np.min([np.min(d["composition"]) for d in dict_list])
    max_comp = np.max([np.max(d["composition"]) for d in dict_list])

    # now left ref will be min fe value from all dicts, corresponds to min_comp
    min_fe = []
    max_fe = []
    for d in dict_list:
        diff = np.abs(d["composition"] - min_comp)
        arg = np.argsort(diff)[0]
        if diff[arg] < threshold:
            min_fe.append(d["free_energy"][arg])
        diff = np.abs(d["composition"] - max_comp)
        arg = np.argsort(diff)[0]
        if diff[arg] < threshold:
            max_fe.append(d["free_energy"][arg])

    # lists are grabbed, now get the references
    left_ref = np.min(min_fe)
    right_ref = np.min(max_fe)

    # print(left_ref, right_ref)
    # now once again, loop through, and add the diff
    for d in dict_list:
        # adjust ref based on composition demands
        scaled_comp = d["composition"] / max_comp
        right_ref_scaled = right_ref * scaled_comp
        left_ref_scaled = left_ref * (1 - scaled_comp)

        # print(d["free_energy"][-1])
        # print((right_ref_scaled + left_ref_scaled)[-1])
        ref = d["free_energy"] - (right_ref_scaled + left_ref_scaled)
        d["free_energy_mix"] = ref

    # Compute mixing energy for raw data points (if present)
    for d in dict_list:
        if "raw_composition" in d and "raw_free_energy" in d:
            rc = d["raw_composition"]
            scaled_rc = rc / max_comp
            raw_ref = right_ref * scaled_rc + left_ref * (1 - scaled_rc)
            d["raw_free_energy_mix"] = d["raw_free_energy"] - raw_ref

    # Trim boundary points from partial-range phases
    if boundary_trim > 0:
        for d in dict_list:
            comp = d["composition"]
            c_min, c_max = np.min(comp), np.max(comp)
            left_partial = c_min > min_comp + threshold
            right_partial = c_max < max_comp - threshold

            if left_partial or right_partial:
                # Cap trim so it never exceeds 1/4 of the phase's own range,
                # preventing narrow phases (e.g. ordered compounds) from being
                # trimmed to an empty array.
                phase_range = c_max - c_min
                trim = min(float(boundary_trim), phase_range / 4.0)

                mask = np.ones(len(comp), dtype=bool)
                if left_partial:
                    mask &= comp >= c_min + trim
                if right_partial:
                    mask &= comp <= c_max - trim
                if mask.any():
                    d["composition"] = comp[mask]
                    d["free_energy"] = d["free_energy"][mask]
                    d["free_energy_mix"] = d["free_energy_mix"][mask]

    return dict_list


TABLEAU10 = [
    "#4E79A7",
    "#F28E2B",
    "#E15759",
    "#76B7B2",
    "#59A14F",
    "#EDC948",
    "#B07AA1",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AC",
]


def create_color_list(phases):
    combinations_list = ["-".join(pair) for pair in combinations(phases, 2)]
    same_element_pairs = ["-".join([item, item]) for item in phases]
    final_combinations = same_element_pairs + combinations_list

    color_dict = {}

    for count, combination in enumerate(final_combinations):
        color_hex = TABLEAU10[count % len(TABLEAU10)]
        color_dict[combination] = color_hex
        raw = combination.split("-")
        if raw[0] != raw[1]:
            reversecombo = f"{raw[1]}-{raw[0]}"
            color_dict[reversecombo] = color_hex
    return color_dict


def get_tangent_type(dict_list, tangent, energy):
    left_c = tangent[0]
    right_c = tangent[1]

    left_e = energy[0]
    right_e = energy[1]

    left_phase = None
    right_phase = None

    left_values = []
    left_phases = []
    right_values = []
    right_phases = []

    for d in dict_list:
        if len(d["composition"]) == 0:
            continue
        diff = np.abs(left_c - d["composition"])
        arg = np.argsort(diff)[0]
        if diff[arg] < 1e-5:
            a = np.abs(left_e - d["free_energy_mix"][arg])
            left_values.append(a)
            left_phases.append(d["phase"])
        diff = np.abs(right_c - d["composition"])
        arg = np.argsort(diff)[0]
        if diff[arg] < 1e-5:
            a = np.abs(right_e - d["free_energy_mix"][arg])
            right_values.append(a)
            right_phases.append(d["phase"])

    # now check min values
    left_min_arg = np.argmin(left_values)
    if left_values[left_min_arg] < 1e-5:
        # this is ok
        left_phase = left_phases[left_min_arg]

    right_min_arg = np.argmin(right_values)
    if right_values[right_min_arg] < 1e-5:
        # this is ok
        right_phase = right_phases[right_min_arg]

    phase_str = f"{left_phase}-{right_phase}"
    return phase_str


def get_common_tangents(
    dict_list, peak_cutoff=0.01, plot=False, remove_self_tangents_for=[]
):
    """
    Get common tangent constructions using convex hull method
    """
    points = np.vstack(
        [np.column_stack((d["composition"], d["free_energy_mix"])) for d in dict_list]
    )

    # if color_dict is None:
    #    color_dict = create_color_list(dict_list)

    # make common tangent constructions
    # term checks if two different phases are stable at the end points, then common tangent is needed
    hull = ConvexHull(points)
    convex_points = []
    convex_x = []
    for simplex, equation in zip(hull.simplices, hull.equations):
        # Select lower convex hull facets: outward normal has negative y-component
        if equation[1] < 0:
            convex_points.extend(points[simplex, 1])
            convex_x.extend(points[simplex, 0])

    dist = np.diff(np.sort(convex_x))
    dist = np.where(dist > peak_cutoff)[0]
    sargs = np.argsort(convex_x)
    convex_x = np.array(convex_x)
    convex_points = np.array(convex_points)

    tangents = []
    energies = []
    tangent_types = []
    phases = []

    for d in dist:
        t = [convex_x[sargs][d], convex_x[sargs][d + 1]]
        e = [convex_points[sargs][d], convex_points[sargs][d + 1]]
        phase_str = get_tangent_type(dict_list, t, e)

        remove = False
        ps = phase_str.split("-")
        if ps[0] == ps[1]:
            if ps[0] in remove_self_tangents_for:
                remove = True

        if not remove:
            tangents.append(t)
            energies.append(e)
            tangent_types.append(phase_str)
            phases.append(phase_str.split("-"))

    if plot:
        for d in dict_list:
            plt.plot(
                d["composition"],
                d["free_energy_mix"],
                color=colors[np.random.randint(len(colors))],
            )
        for t, e in zip(tangents, energies):
            plt.plot(t, e, color="black", ls="dashed")
        plt.ylim(top=max(0.0, np.max(convex_points) * 1.1))

    return (
        np.array(tangents),
        np.array(energies),
        np.array(tangent_types),
        np.array(phases),
    )


def _get_single_phase_boundaries(tangents, temperatures, tangent_types):
    """
    Derive single-phase region composition intervals from common-tangent data.

    At each temperature the tangent construction identifies coexistence windows
    ``(x_left, x_right)``.  The single-phase regions are the complementary
    composition intervals that lie *outside* every coexistence window.

    Returns
    -------
    dict[str, list[tuple]]
        Mapping phase_name -> list of ``(T, x_lo, x_hi)`` sorted by temperature.
    """
    from collections import defaultdict

    phase_intervals = defaultdict(list)

    for count, T in enumerate(temperatures):
        T_tangents = tangents[count]
        T_types = np.atleast_1d(tangent_types[count])
        if len(T_tangents) == 0:
            continue

        T_tangents = np.atleast_2d(T_tangents)
        order = np.argsort(T_tangents[:, 0])
        sorted_tangs = T_tangents[order]
        sorted_types = T_types[order]

        # Build a flat boundary list: alternating (x, phase) for left then right
        # of each coexistence window, sorted by x.
        boundaries = []
        for tang, ttype in zip(sorted_tangs, sorted_types):
            x_left, x_right = float(tang[0]), float(tang[1])
            parts = str(ttype).split("-")
            left_p = parts[0] if parts[0] != "None" else None
            right_p = parts[1] if parts[1] != "None" else None
            boundaries.append((x_left, left_p))
            boundaries.append((x_right, right_p))

        # Single-phase region before the first coexistence window
        x0, p0 = boundaries[0]
        if p0 is not None:
            phase_intervals[p0].append((T, 0.0, x0))

        # Single-phase regions *between* consecutive coexistence windows
        for i in range(1, len(boundaries) - 1, 2):
            x_start, p_start = boundaries[i]  # right edge of current window
            x_end, p_end = boundaries[i + 1]  # left edge of next window
            if p_start is not None and x_end > x_start:
                phase_intervals[p_start].append((T, x_start, x_end))

        # Single-phase region after the last coexistence window
        x_last, p_last = boundaries[-1]
        if p_last is not None:
            phase_intervals[p_last].append((T, x_last, 1.0))

    for p in phase_intervals:
        phase_intervals[p].sort(key=lambda r: r[0])

    return dict(phase_intervals)


def plot_phase_diagram(
    tangents,
    temperature,
    tangent_types,
    phases,
    edgecolor="#37474f",
    linewidth=1,
    linestyle="-",
    fill=True,
    alpha=0.35,
    border_lw=2,
    smooth_boundary=0,
    color_phases=False,
    figsize=None,
    ax=None,
):
    """
    Plot a binary phase diagram.

    Parameters
    ----------
    tangents : list of arrays
        Tangent composition pairs at each temperature, output of the
        phase-diagram loop.
    temperature : list
        Temperature value for each entry in *tangents*.
    tangent_types : list of arrays
        Phase-pair labels (e.g. ``"cufcc-lqd"``) for every tangent.
    phases : list of str
        Ordered phase names used to build the colour palette.
    edgecolor : str
        Colour for polygon borders and the figure frame.
    linewidth : float
        Line width when *fill* is False (legacy horizontal-line mode).
    linestyle : str
        Line style when *fill* is False.
    fill : bool
        If True (default), render two-phase regions as filled polygons
        with coloured borders.  If False, fall back to horizontal lines.
    alpha : float
        Fill opacity for polygons (0–1).
    border_lw : float
        Line width of the polygon borders.
    smooth_boundary : int
        Savitzky-Golay window size (odd integer) for smoothing polygon
        boundaries.  Set to 0 (default) to disable.  A value of 11
        is a good starting point.
    color_phases : bool
        If True, fill *single-phase* regions with per-phase colours instead
        of filling two-phase coexistence regions.  Two-phase regions are
        left uncoloured (white background), and the legend lists each phase
        individually.  Default False (coexistence-region colouring).
    figsize : tuple or None
        Figure size.  Defaults to (7, 5).
    ax : matplotlib Axes or None
        If given, draw on this axes instead of creating a new figure.

    Returns
    -------
    fig : matplotlib Figure
    ax : matplotlib Axes
    """
    if figsize is None:
        figsize = (7, 5)

    color_dict = create_color_list(phases)
    # Per-phase colour dict (used when color_phases=True)
    phase_color_dict = {p: TABLEAU10[i % len(TABLEAU10)] for i, p in enumerate(phases)}

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    if color_phases:
        # ---- single-phase region colouring ----
        phase_intervals = _get_single_phase_boundaries(
            tangents, temperature, tangent_types
        )

        for phase, intervals in phase_intervals.items():
            color = phase_color_dict.get(phase, TABLEAU10[0])
            intervals.sort(key=lambda r: r[0])
            Ts = [r[0] for r in intervals]
            x_lo = [r[1] for r in intervals]
            x_hi = [r[2] for r in intervals]

            if smooth_boundary > 0 and len(x_lo) > smooth_boundary:
                x_lo = list(savgol_filter(x_lo, smooth_boundary, 3))
                x_hi = list(savgol_filter(x_hi, smooth_boundary, 3))

            poly_x = x_lo + x_hi[::-1]
            poly_T = Ts + Ts[::-1]
            ax.fill(poly_x, poly_T, color=color, alpha=alpha)
            ax.plot(x_lo, Ts, color=edgecolor, lw=border_lw, solid_capstyle="round")
            ax.plot(x_hi, Ts, color=edgecolor, lw=border_lw, solid_capstyle="round")

        legend_patches = [
            mpatches.Patch(color=phase_color_dict[p], label=p)
            for p in phases
            if p in phase_color_dict
        ]
        ax.legend(handles=legend_patches, loc="center left", bbox_to_anchor=(1, 0.5))
        ax.set_xlabel("Composition")
        ax.set_ylabel("T (K)")
        fig.tight_layout()
        return fig, ax

    if not fill:
        # ---- legacy: horizontal lines ----
        for count, x in enumerate(tangents):
            for c, a in enumerate(x):
                ax.plot(
                    np.array(a),
                    [temperature[count], temperature[count]],
                    linestyle,
                    lw=linewidth,
                    c=color_dict[tangent_types[count][c]],
                )
    else:
        # ---- filled polygons ----
        # Collect boundary curves per region type.
        # Each region type accumulates (T, x_left, x_right) triples.
        from collections import defaultdict

        region_data = defaultdict(list)

        for count, x in enumerate(tangents):
            T = temperature[count]
            for c, a in enumerate(x):
                label = tangent_types[count][c]
                region_data[label].append((T, a[0], a[1]))

        for label, rows in region_data.items():
            rows.sort(key=lambda r: r[0])
            Ts = [r[0] for r in rows]
            x_left = [r[1] for r in rows]
            x_right = [r[2] for r in rows]

            # Optionally smooth the boundary curves
            if smooth_boundary > 0 and len(x_left) > smooth_boundary:
                x_left = list(savgol_filter(x_left, smooth_boundary, 3))
                x_right = list(savgol_filter(x_right, smooth_boundary, 3))

            # Build a closed polygon: go up along the left boundary,
            # then back down along the right boundary.
            poly_x = x_left + x_right[::-1]
            poly_T = Ts + Ts[::-1]

            color = color_dict[label]
            ax.fill(poly_x, poly_T, color=color, alpha=alpha)
            # Draw the left and right boundary lines
            ax.plot(x_left, Ts, color=edgecolor, lw=border_lw, solid_capstyle="round")
            ax.plot(x_right, Ts, color=edgecolor, lw=border_lw, solid_capstyle="round")

    # Build legend only for regions that actually appear in the data
    seen_labels = set()
    for tt in tangent_types:
        for label in np.atleast_1d(tt):
            seen_labels.add(label)

    legend_patches = []
    seen_colors = set()
    for label in seen_labels:
        color = color_dict.get(label, TABLEAU10[0])
        # Deduplicate reversed pairs (e.g. lqd-agfcc == agfcc-lqd)
        canonical = "-".join(sorted(label.split("-")))
        if canonical not in seen_colors:
            legend_patches.append(mpatches.Patch(color=color, label=label))
            seen_colors.add(canonical)
    ax.legend(handles=legend_patches, loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_xlabel("Composition")
    ax.set_ylabel("T (K)")
    fig.tight_layout()
    return fig, ax


# ---------------------------------------------------------------------------
# CALPHAD surface helpers
# ---------------------------------------------------------------------------

def _fit_calphad_poly6(T, G):
    """
    Fit the standard six-term CALPHAD polynomial to G(T) data.

        G(T) = a + b·T + c·T·ln T + d·T² + e·T³ + f·T⁻¹

    Parameters
    ----------
    T : array-like  Temperature in K.
    G : array-like  Free energy in eV/atom.

    Returns
    -------
    coeffs : ndarray, shape (6,)  [a, b, c, d, e, f]
    """
    T = np.asarray(T, dtype=float)
    G = np.asarray(G, dtype=float)
    A = np.column_stack(
        [np.ones_like(T), T, T * np.log(T), T ** 2, T ** 3, 1.0 / T]
    )
    coeffs, _, _, _ = np.linalg.lstsq(A, G, rcond=None)
    return coeffs


def _eval_calphad_poly6(coeffs, T):
    """
    Evaluate the six-term CALPHAD polynomial.

    Parameters
    ----------
    coeffs : array-like  [a, b, c, d, e, f] from :func:`_fit_calphad_poly6`.
    T : float or array-like  Temperature in K.

    Returns
    -------
    G : ndarray  Free energy in eV/atom.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    return (
        coeffs[0]
        + coeffs[1] * T
        + coeffs[2] * T * np.log(T)
        + coeffs[3] * T ** 2
        + coeffs[4] * T ** 3
        + coeffs[5] / T
    )


def _eval_calphad_surface_at(surface, x, T):
    """
    Evaluate the CALPHAD-decomposed G(x, T) surface for a single phase at
    temperature *T* over a composition array *x*.

    The model is::

        G(x, T) = (1-x)·G_A(T) + x·G_B(T)
                  + k_B·T·[x·ln x + (1-x)·ln(1-x)]
                  + x·(1-x)·Σ_k L_k·(1-2x)^k

    Parameters
    ----------
    surface : dict  Output entry from :meth:`PhaseDiagram.build_calphad_surface`.
    x : array-like  Composition grid (0 to 1).
    T : float       Temperature in K.

    Returns
    -------
    G : ndarray  Free energy in eV/atom.
    """
    x = np.asarray(x, dtype=float)
    T = float(T)
    coeffs_A = surface["coeffs_A"]
    coeffs_B = surface["coeffs_B"]
    L = surface["L_coeffs"]

    G_A = float(_eval_calphad_poly6(coeffs_A, T))
    G_B = float(_eval_calphad_poly6(coeffs_B, T))
    G_lin = (1.0 - x) * G_A + x * G_B

    # Ideal configurational entropy  (→ 0 at x=0 and x=1)
    log_x = np.where(x > 0, np.log(np.maximum(x, 1e-300)), 0.0)
    log_1mx = np.where(x < 1, np.log(np.maximum(1.0 - x, 1e-300)), 0.0)
    G_ideal = kb * T * (x * log_x + (1.0 - x) * log_1mx)

    # Redlich-Kister excess
    pf = x * (1.0 - x)
    G_xs = sum(L[k] * pf * (1.0 - 2.0 * x) ** k for k in range(len(L)))

    return G_lin + G_ideal + G_xs


class PhaseDiagram:
    """
    High-level class for computing and plotting a binary phase diagram.

    Wraps the full workflow — data gathering, cleaning, free-energy
    fitting, common-tangent construction, and plotting — into a single
    object.

    Parameters
    ----------
    folders : dict
        Mapping of phase name → folder path (or list of folder paths).
        Pass a list to merge multiple simulation folders into a single
        phase — the raw data is combined *before* any fitting, so there
        are no cross-folder inconsistencies::

            # single folder per phase
            {'cufcc': 'cufcc', 'agfcc': 'agfcc', 'lqd': 'lqd'}

            # aufcc and cufcc merged into one 'fcc' phase
            {'fcc': ['aufcc_folder', 'cufcc_folder'], 'lqd': 'lqd'}

    reference_element : str
        The element whose fraction is used as the composition axis
        (e.g. ``'Ag'``).
    composition_intervals : dict, optional
        Per-phase composition bounds, e.g.
        ``{'fcc': (0, 1), 'lqd': (0, 1)}``.
        Phases not listed are auto-detected from the data: the
        interval is set to the ``(min, max)`` of available
        compositions for that phase.
    smooth : bool
        If True (default), smooth F(T) data with thermodynamic basis
        during the ``clean_df`` step.

    Examples
    --------
    >>> pd = PhaseDiagram(
    ...     folders={'fcc': ['aufcc', 'cufcc'], 'lqd': 'lqd'},
    ...     reference_element='Au',
    ... )
    >>> pd.calculate(T_range=(400, 1400), T_step=5, fit_order=4,
    ...              method='redlich-kister')
    >>> fig, ax = pd.plot()
    """

    def __init__(
        self, folders, reference_element, composition_intervals=None, smooth=True
    ):
        from calphy.postprocessing import (
            gather_results,
            clean_df,
            fix_composition_scaling,
        )

        self.reference_element = reference_element
        self.composition_intervals = composition_intervals or {}
        self.phases = list(folders.keys())

        # ---- gather & clean ----
        # Each value in `folders` can be a single path or a list of paths.
        # IMPORTANT: fix_composition_scaling anchors each composition_scaling
        # row to its phase's reference (is_reference=True) row.  A single
        # folder may contain calculations with *different* phase_name values
        # (e.g. lqdAg + lqdCu in the same folder), each having its own
        # reference point.  We therefore assign a unique temp name per
        # (folder_index, original_phase_name) group so that clean_df and
        # fix_composition_scaling always see exactly one reference per group.
        # Only after both steps are all groups merged under the user's key.
        temp_to_user = {}  # temporary_phase_name -> user key
        dfs = []
        for phase, folder_spec in folders.items():
            folder_list = (
                [folder_spec] if isinstance(folder_spec, str) else list(folder_spec)
            )
            for idx, folder in enumerate(folder_list):
                df = gather_results(
                    folder, reduce_composition=True, extract_phase_prefix=True
                )

                # Assign one temp name per distinct original phase_name so
                # each group (with its own reference row) is processed
                # independently through fix_composition_scaling.
                def _make_temp(orig, phase=phase, idx=idx):
                    return f"__calphy_{phase}_{idx}_{orig}__"

                for orig_name in df["phase_name"].unique():
                    temp_to_user[_make_temp(orig_name)] = phase
                df["phase_name"] = df["phase_name"].apply(_make_temp)
                dfs.append(df)

        combined = pd.concat(dfs, ignore_index=True)
        combined = combined.loc[combined.status == "True"]

        dfc = clean_df(
            combined, reference_element, combine_direct_calculations=True, smooth=smooth
        )
        dfc = fix_composition_scaling(dfc, add_ideal_entropy=True)

        # Merge temp-named groups back to user keys and label the phase column.
        user_dfc = {}
        for temp_name, phase_df in dfc.items():
            user_key = temp_to_user[temp_name]
            if user_key in user_dfc:
                user_dfc[user_key] = pd.concat(
                    [user_dfc[user_key], phase_df], ignore_index=True
                )
            else:
                user_dfc[user_key] = phase_df

        for key, val in user_dfc.items():
            val["phase"] = key

        self.df = pd.concat([val for val in user_dfc.values()])

        # ---- auto-detect composition intervals from data ----
        for phase in self.phases:
            if phase not in self.composition_intervals:
                df_p = self.df.loc[self.df["phase"] == phase, "composition"]
                if len(df_p) > 0:
                    self.composition_intervals[phase] = (
                        float(df_p.min()),
                        float(df_p.max()),
                    )

        # ---- state filled by calculate() ----
        self.tangents = None
        self.temperatures = None
        self.tangent_types = None
        self._calc_kwargs = {}
        self._calphad_surfaces = {}

    # ------------------------------------------------------------------
    # CALPHAD surface decomposition
    # ------------------------------------------------------------------

    def build_calphad_surface(self, rk_order=3, endpoint_tol=0.05):
        """
        Fit a CALPHAD-decomposed Gibbs energy surface G(x, T) for each phase.

        The model mirrors the pycalphad/TDB representation::

            G(x, T) = (1-x)·G_A(T) + x·G_B(T)
                      + k_B·T·[x·ln x + (1-x)·ln(1-x)]
                      + x·(1-x)·Σ_k L_k·(1-2x)^k

        *G_A(T)* and *G_B(T)* are six-term CALPHAD polynomials fitted to the
        pure-endpoint calphy data.  *L_k* are Redlich-Kister interaction
        parameters fitted to temperature-averaged excess free energies at
        intermediate compositions (T-independent RK approximation).

        Because G(x, T) is a smooth analytic surface in both x and T, phase
        boundaries computed from it via :meth:`calculate` (with
        ``calphad_surface=True``) are substantially smoother and more
        physically consistent than those from the per-temperature polynomial
        fits in the default mode.  This is the root reason why the
        pycalphad/TDB route yields better phase diagrams: it decomposes G into
        physically motivated components rather than fitting a raw polynomial
        slice-by-slice in composition at each T.

        Call this method *before* ``calculate(calphad_surface=True)``, or pass
        ``calphad_surface=True`` directly to ``calculate()`` (which will call
        this automatically).

        Parameters
        ----------
        rk_order : int
            Number of Redlich-Kister parameters to fit (default 3).
        endpoint_tol : float
            Composition window for locating pure endpoints: rows with
            ``composition <= endpoint_tol`` are candidates for the A endpoint
            and rows with ``composition >= 1 - endpoint_tol`` for the B
            endpoint (default 0.05).

        Returns
        -------
        dict
            Keyed by phase name.  Each value is either ``None`` (surface could
            not be built) or a dict with keys ``'coeffs_A'``, ``'coeffs_B'``,
            ``'L_coeffs'``, ``'x_data'``, ``'G_xs_data'``, ``'rk_order'``.
        """
        surfaces = {}

        for phase in self.phases:
            df_p = self.df.loc[self.df["phase"] == phase].copy()
            df_p = df_p.sort_values("composition")

            # --- Locate pure endpoints ---
            left_rows = df_p.loc[df_p["composition"] <= endpoint_tol]
            right_rows = df_p.loc[df_p["composition"] >= 1.0 - endpoint_tol]

            if len(left_rows) == 0 or len(right_rows) == 0:
                warnings.warn(
                    f"Phase '{phase}': no pure endpoints found within "
                    f"composition tolerance {endpoint_tol}. "
                    "Skipping CALPHAD surface for this phase."
                )
                surfaces[phase] = None
                continue

            # Prefer is_reference (fe-mode) rows; fall back to closest-to-endpoint.
            def _pick_ep(rows, target):
                if "is_reference" in rows.columns:
                    ref = rows.loc[rows["is_reference"] == True]
                    if len(ref) > 0:
                        rows = ref
                return rows.iloc[
                    (rows["composition"] - target).abs().argsort().iloc[0]
                ]

            left_row = _pick_ep(left_rows, 0.0)
            right_row = _pick_ep(right_rows, 1.0)

            T_A = np.asarray(left_row["temperature"], dtype=float)
            G_A_data = np.asarray(left_row["free_energy"], dtype=float)
            T_B = np.asarray(right_row["temperature"], dtype=float)
            G_B_data = np.asarray(right_row["free_energy"], dtype=float)

            # --- Fit CALPHAD 6-term polynomials to pure endpoints ---
            coeffs_A = _fit_calphad_poly6(T_A, G_A_data)
            coeffs_B = _fit_calphad_poly6(T_B, G_B_data)

            # --- Extract excess free energy at intermediate compositions ---
            x_data = []
            G_xs_data = []
            intermediate = df_p.loc[
                (df_p["composition"] > endpoint_tol)
                & (df_p["composition"] < 1.0 - endpoint_tol)
            ]

            for _, row in intermediate.iterrows():
                x = float(row["composition"])
                T_row = np.asarray(row["temperature"], dtype=float)
                G_row = np.asarray(row["free_energy"], dtype=float)

                # Linear (mechanical-mixture) reference
                G_A_at_T = _eval_calphad_poly6(coeffs_A, T_row)
                G_B_at_T = _eval_calphad_poly6(coeffs_B, T_row)
                G_lin = (1.0 - x) * G_A_at_T + x * G_B_at_T

                # Ideal configurational entropy
                G_ideal = kb * T_row * (
                    x * np.log(x) + (1.0 - x) * np.log(1.0 - x)
                )

                # Excess = total − linear reference − ideal entropy.
                # Average over temperature (T-independent RK approximation).
                G_xs_arr = G_row - G_lin - G_ideal
                x_data.append(x)
                G_xs_data.append(float(np.mean(G_xs_arr)))

            if len(x_data) < rk_order:
                warnings.warn(
                    f"Phase '{phase}': only {len(x_data)} intermediate "
                    f"compositions available, need ≥ {rk_order} for the RK "
                    "fit. Skipping CALPHAD surface for this phase."
                )
                surfaces[phase] = None
                continue

            x_arr = np.array(x_data)
            Gxs_arr = np.array(G_xs_data)

            # --- Fit Redlich-Kister parameters ---
            # G_xs(x) = x·(1-x)·Σ_k L_k·(1-2x)^k  enforces G_xs→0 at x=0,1
            prefactor = x_arr * (1.0 - x_arr)
            basis = np.column_stack(
                [prefactor * (1.0 - 2.0 * x_arr) ** k for k in range(rk_order)]
            )
            L_coeffs, _, _, _ = np.linalg.lstsq(basis, Gxs_arr, rcond=None)

            L_str = "  ".join(
                f"L{k}={L_coeffs[k] * 96485:.1f} J/mol" for k in range(rk_order)
            )
            rms = float(
                np.sqrt(np.mean((Gxs_arr - basis @ L_coeffs) ** 2)) * 96485
            )
            print(f"[{phase}]  {L_str}   RMS_xs={rms:.1f} J/mol")

            surfaces[phase] = {
                "coeffs_A": coeffs_A,
                "coeffs_B": coeffs_B,
                "L_coeffs": L_coeffs,
                "x_data": x_arr,
                "G_xs_data": Gxs_arr,
                "rk_order": rk_order,
            }

        self._calphad_surfaces = surfaces
        return surfaces

    # ------------------------------------------------------------------
    # Core computation
    # ------------------------------------------------------------------

    def calculate(
        self,
        T_range=(400, 1400),
        T_step=5,
        fit_order=4,
        method="polynomial",
        boundary_trim=0.1,
        remove_self_tangents_for=None,
        ideal_configurational_entropy=True,
        end_weight=3,
        end_indices=4,
        calphad_surface=False,
        rk_order=3,
        composition_grid=10000,
    ):
        """
        Compute common-tangent constructions across a temperature range.

        Parameters
        ----------
        T_range : tuple
            (T_min, T_max) in Kelvin.
        T_step : float
            Temperature increment.
        fit_order : int
            Polynomial / Redlich-Kister order for F(x) fits.
        method : str
            ``'polynomial'`` or ``'redlich-kister'``.
        boundary_trim : float
            Amount to trim from partial-range phase boundaries.
            ``'auto'`` (default) computes 2× the average composition
            spacing per phase.
        remove_self_tangents_for : list of str, optional
            Phase names for which same-phase tangent constructions
            should be discarded.
        ideal_configurational_entropy : bool
            Include ideal configurational entropy in the free energies.
            Defaults to True.  Set to False to remove the ideal
            configurational-entropy contribution (e.g. for ordered phases
            or when no swap moves were performed).
        end_weight : int
            Weight for endpoints in the fit.
        end_indices : int
            Number of endpoint indices to weight.
        calphad_surface : bool
            If True, use a CALPHAD-decomposed G(x, T) surface for the phase
            boundary calculation instead of the default per-temperature
            polynomial fits.  The surface is::

                G(x,T) = (1-x)·G_A(T) + x·G_B(T)
                          + k_B·T·[x ln x + (1-x) ln(1-x)]
                          + x(1-x)·Σ_k L_k·(1-2x)^k

            This mirrors the pycalphad/TDB approach and yields smoother,
            more physically consistent phase boundaries because G is smooth
            in both x and T simultaneously.  If :meth:`build_calphad_surface`
            has not been called yet it is invoked automatically using
            *rk_order*.  Default False (legacy per-T polynomial mode).
        rk_order : int
            Number of Redlich-Kister parameters for the CALPHAD surface.
            Only used when *calphad_surface* is True and the surface has not
            been pre-built.  Default 3.
        composition_grid : int
            Number of composition points on the evaluation grid when
            *calphad_surface* is True.  Default 10000.
        """
        if remove_self_tangents_for is None:
            remove_self_tangents_for = []

        # Build the CALPHAD surface on demand.
        if calphad_surface and not self._calphad_surfaces:
            self.build_calphad_surface(rk_order=rk_order)

        self._calc_kwargs = dict(
            fit_order=fit_order,
            method=method,
            boundary_trim=boundary_trim,
            remove_self_tangents_for=remove_self_tangents_for,
            ideal_configurational_entropy=ideal_configurational_entropy,
            end_weight=end_weight,
            end_indices=end_indices,
            calphad_surface=calphad_surface,
        )

        temps_arr = np.arange(T_range[0], T_range[1], T_step)
        tangents = []
        temps = []
        tangent_types = []

        for t in temps_arr:
            dict_list = []
            for phase in self.phases:
                ci = self.composition_intervals.get(phase, (0, 1))

                if calphad_surface and self._calphad_surfaces.get(phase) is not None:
                    # Use the smooth analytical CALPHAD surface.
                    compfine = np.linspace(ci[0], ci[1], composition_grid)
                    G_arr = _eval_calphad_surface_at(
                        self._calphad_surfaces[phase], compfine, t
                    )
                    d = {
                        "phase": phase,
                        "temperature": t,
                        "composition": compfine,
                        "free_energy": G_arr,
                    }
                else:
                    d = get_phase_free_energy(
                        self.df,
                        phase,
                        t,
                        ideal_configurational_entropy=ideal_configurational_entropy,
                        composition_interval=ci,
                        fit_order=fit_order,
                        method=method,
                        end_weight=end_weight,
                        end_indices=end_indices,
                    )

                if d is not None:
                    dict_list.append(d)

            if len(dict_list) > 0:
                dc = get_free_energy_mixing(dict_list, boundary_trim=boundary_trim)
                tn, _, cn, _ = get_common_tangents(
                    dc, remove_self_tangents_for=remove_self_tangents_for
                )
                tangents.append(tn)
                temps.append(t)
                tangent_types.append(cn)

        self.tangents = tangents
        self.temperatures = temps
        self.tangent_types = tangent_types


    # ------------------------------------------------------------------
    # Phase diagram plot
    # ------------------------------------------------------------------

    def plot(
        self,
        fill=True,
        alpha=0.2,
        border_lw=2,
        smooth_boundary=11,
        color_phases=False,
        figsize=None,
        ax=None,
        **kwargs,
    ):
        """
        Plot the full phase diagram.

        Parameters
        ----------
        color_phases : bool
            If True, fill single-phase regions with per-phase colours instead
            of filling two-phase coexistence regions.  Default False.

        All other keyword arguments are forwarded to
        :func:`plot_phase_diagram`.

        Returns
        -------
        fig : matplotlib Figure
        ax : matplotlib Axes
        """
        self._require_calculated()
        return plot_phase_diagram(
            self.tangents,
            self.temperatures,
            self.tangent_types,
            self.phases,
            fill=fill,
            alpha=alpha,
            border_lw=border_lw,
            smooth_boundary=smooth_boundary,
            color_phases=color_phases,
            figsize=figsize,
            ax=ax,
            **kwargs,
        )

    # ------------------------------------------------------------------
    # Free-energy curves at a single temperature
    # ------------------------------------------------------------------

    def plot_free_energy(self, T, show_data=False, figsize=None, ax=None):
        """
        Plot free-energy curves F(x) for all phases at temperature *T*.

        Parameters
        ----------
        T : float
            Temperature in Kelvin.
        show_data : bool
            If True, overlay the raw data points on each phase curve.
        figsize : tuple, optional
        ax : matplotlib Axes, optional

        Returns
        -------
        fig : matplotlib Figure
        ax : matplotlib Axes
        """
        kw = self._calc_kwargs
        if figsize is None:
            figsize = (7, 5)
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        color_dict = create_color_list(self.phases)

        for phase in self.phases:
            ci = self.composition_intervals.get(phase, (0, 1))
            d = get_phase_free_energy(
                self.df,
                phase,
                T,
                ideal_configurational_entropy=kw.get(
                    "ideal_configurational_entropy", True
                ),
                composition_interval=ci,
                fit_order=kw.get("fit_order", 4),
                method=kw.get("method", "polynomial"),
                end_weight=kw.get("end_weight", 3),
                end_indices=kw.get("end_indices", 4),
            )
            if d is not None:
                color = color_dict.get(f"{phase}-{phase}", TABLEAU10[0])
                ax.plot(
                    d["composition"], d["free_energy"], label=phase, color=color, lw=2
                )
                if show_data and "raw_composition" in d:
                    ax.scatter(
                        d["raw_composition"],
                        d["raw_free_energy"],
                        s=40,
                        zorder=5,
                        color=color,
                        edgecolors="black",
                        lw=0.8,
                    )

        ax.set_xlabel("Composition")
        ax.set_ylabel("F (eV/atom)")
        ax.set_title(f"T = {T} K")
        ax.legend()
        fig.tight_layout()
        return fig, ax

    # ------------------------------------------------------------------
    # Free-energy of mixing + tangent lines at one temperature
    # ------------------------------------------------------------------

    def plot_free_energy_mixing(self, T, show_data=False, figsize=None, ax=None):
        """
        Plot free energy of mixing F_mix(x) with common-tangent lines
        at temperature *T*.

        Parameters
        ----------
        T : float
            Temperature in Kelvin.
        show_data : bool
            If True, overlay the raw data points on each phase curve.

        Returns
        -------
        fig : matplotlib Figure
        ax : matplotlib Axes
        """
        kw = self._calc_kwargs
        if figsize is None:
            figsize = (7, 5)
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        color_dict = create_color_list(self.phases)

        dict_list = []
        for phase in self.phases:
            ci = self.composition_intervals.get(phase, (0, 1))
            d = get_phase_free_energy(
                self.df,
                phase,
                T,
                ideal_configurational_entropy=kw.get(
                    "ideal_configurational_entropy", True
                ),
                composition_interval=ci,
                fit_order=kw.get("fit_order", 4),
                method=kw.get("method", "polynomial"),
                end_weight=kw.get("end_weight", 3),
                end_indices=kw.get("end_indices", 4),
            )
            if d is not None:
                dict_list.append(d)

        if len(dict_list) == 0:
            warnings.warn(f"No valid data at T={T}")
            return fig, ax

        dc = get_free_energy_mixing(
            dict_list, boundary_trim=kw.get("boundary_trim", 0.1)
        )

        for d in dc:
            phase = d["phase"]
            color = color_dict.get(f"{phase}-{phase}", TABLEAU10[0])
            ax.plot(
                d["composition"], d["free_energy_mix"], label=phase, color=color, lw=2
            )
            if show_data and "raw_free_energy_mix" in d:
                ax.scatter(
                    d["raw_composition"],
                    d["raw_free_energy_mix"],
                    s=40,
                    zorder=5,
                    color=color,
                    edgecolors="black",
                    lw=0.8,
                )

        tn, en, _, _ = get_common_tangents(
            dc, remove_self_tangents_for=kw.get("remove_self_tangents_for", [])
        )

        for t, e in zip(tn, en):
            ax.plot(t, e, color="black", ls="dashed", lw=1.5)

        ax.set_xlabel("Composition")
        ax.set_ylabel(r"$\Delta F_\mathrm{mix}$ (eV/atom)")
        ax.set_title(f"T = {T} K")
        ax.legend()
        fig.tight_layout()
        return fig, ax

    # ------------------------------------------------------------------
    # Raw data vs fit at one temperature
    # ------------------------------------------------------------------

    def plot_data_vs_fit(self, phase, T, figsize=None, ax=None):
        """
        Compare raw free-energy data points with the fitted curve
        for a single phase at temperature *T*.

        Returns
        -------
        fig : matplotlib Figure
        ax : matplotlib Axes
        """
        kw = self._calc_kwargs
        if figsize is None:
            figsize = (7, 5)
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        ci = self.composition_intervals.get(phase, (0, 1))
        df_phase = self.df.loc[self.df["phase"] == phase].copy()
        df_phase = df_phase.sort_values(by="composition")
        df_phase = df_phase[
            (df_phase["composition"] >= ci[0]) & (df_phase["composition"] <= ci[1])
        ]

        composition = df_phase["composition"].values
        args = df_phase["temperature"].apply(_get_temp_arg, args=(T,))
        fes = _get_fe_at_args(df_phase["free_energy"].values, args)

        comp_raw = np.array(
            [composition[i] for i, x in enumerate(fes) if x is not None]
        )
        fe_raw = np.array([x for x in fes if x is not None])

        if len(fe_raw) == 0:
            warnings.warn(f"No data for phase '{phase}' at T={T}")
            return fig, ax

        # Track which rows already had ideal entropy applied and which
        # are composition_scaling (the only ones we should touch).
        if "_entropy_corrected" in df_phase.columns:
            ec_flags_all = df_phase["_entropy_corrected"].values
        else:
            ec_flags_all = np.zeros(len(df_phase), dtype=bool)
        if "calculation_mode" in df_phase.columns:
            is_cs_all = (df_phase["calculation_mode"] == "composition_scaling").values
        else:
            is_cs_all = np.zeros(len(df_phase), dtype=bool)
        ec_flags = np.array(
            [ec_flags_all[i] for i, x in enumerate(fes) if x is not None]
        )
        is_cs = np.array([is_cs_all[i] for i, x in enumerate(fes) if x is not None])

        ice = kw.get("ideal_configurational_entropy", True)
        entropy_term = kb * T * _calculate_configurational_entropy(comp_raw)
        if ice:
            mask = is_cs & ~ec_flags
            fe_raw = fe_raw.copy()
            fe_raw[mask] -= entropy_term[mask]
        else:
            mask = is_cs & ec_flags
            fe_raw = fe_raw.copy()
            fe_raw[mask] += entropy_term[mask]

        # Plot raw data
        ax.scatter(
            comp_raw,
            fe_raw,
            s=40,
            zorder=5,
            label="data",
            edgecolors="black",
            facecolors="white",
            lw=1.5,
        )

        # Plot fit
        d = get_phase_free_energy(
            self.df,
            phase,
            T,
            ideal_configurational_entropy=ice,
            composition_interval=ci,
            fit_order=kw.get("fit_order", 4),
            method=kw.get("method", "polynomial"),
            end_weight=kw.get("end_weight", 3),
            end_indices=kw.get("end_indices", 4),
        )
        if d is not None:
            ax.plot(
                d["composition"], d["free_energy"], color="#E15759", lw=2, label="fit"
            )

        ax.set_xlabel("Composition")
        ax.set_ylabel("F (eV/atom)")
        ax.set_title(f"{phase} — T = {T} K")
        ax.legend()
        fig.tight_layout()
        return fig, ax

    # ------------------------------------------------------------------
    # CALPHAD surface diagnostics
    # ------------------------------------------------------------------

    def plot_calphad_surface_fit(self, figsize=None):
        """
        Diagnostic plot for the CALPHAD surface fit.

        Shows two subplots per phase:

        * **Left**: fitted CALPHAD G(T) polynomials for the pure endpoints
          (A and B) overlaid on the raw calphy data.
        * **Right**: Redlich-Kister excess G_xs(x) fitted curve overlaid on
          the temperature-averaged excess data points.

        Requires :meth:`build_calphad_surface` to have been called first.

        Parameters
        ----------
        figsize : tuple or None

        Returns
        -------
        fig : matplotlib Figure
        axes : ndarray of Axes
        """
        if not self._calphad_surfaces:
            raise RuntimeError(
                "No CALPHAD surface data. Call build_calphad_surface() first."
            )

        n_phases = sum(
            1 for v in self._calphad_surfaces.values() if v is not None
        )
        if n_phases == 0:
            raise RuntimeError("No valid CALPHAD surfaces were built.")

        if figsize is None:
            figsize = (6 * n_phases, 4)
        fig, axes = plt.subplots(1, n_phases * 2, figsize=figsize)
        if n_phases * 2 == 1:
            axes = [axes]
        axes = np.atleast_1d(axes)

        col = 0
        for phase in self.phases:
            surf = self._calphad_surfaces.get(phase)
            if surf is None:
                continue

            df_p = self.df.loc[self.df["phase"] == phase]
            endpoint_tol = 0.05

            # --- Left panel: G(T) endpoint fits ---
            ax_gt = axes[col]
            col += 1

            for label, comp_lo, comp_hi, coeffs, color in [
                ("A (x≈0)", 0.0, endpoint_tol, surf["coeffs_A"], "#4E79A7"),
                ("B (x≈1)", 1.0 - endpoint_tol, 1.0, surf["coeffs_B"], "#E15759"),
            ]:
                ep_rows = df_p.loc[
                    (df_p["composition"] >= comp_lo)
                    & (df_p["composition"] <= comp_hi)
                ]
                for _, row in ep_rows.iterrows():
                    T_raw = np.asarray(row["temperature"], dtype=float)
                    G_raw = np.asarray(row["free_energy"], dtype=float)
                    ax_gt.scatter(T_raw, G_raw, s=12, alpha=0.6, color=color)
                    T_fit = np.linspace(T_raw.min(), T_raw.max(), 300)
                    ax_gt.plot(
                        T_fit,
                        _eval_calphad_poly6(coeffs, T_fit),
                        lw=2,
                        color=color,
                        label=label,
                    )

            ax_gt.set_xlabel("T (K)")
            ax_gt.set_ylabel("G (eV/atom)")
            ax_gt.set_title(f"{phase} — endpoint G(T) fits")
            ax_gt.legend(fontsize=8)

            # --- Right panel: RK excess G_xs(x) ---
            ax_rk = axes[col]
            col += 1

            x_data = surf["x_data"]
            Gxs_data = surf["G_xs_data"]
            L = surf["L_coeffs"]
            rk_order = surf["rk_order"]

            ax_rk.scatter(
                x_data,
                Gxs_data * 96485,
                s=40,
                zorder=5,
                label="data (avg over T)",
                edgecolors="black",
                facecolors="white",
                lw=1.5,
            )

            x_fine = np.linspace(0, 1, 500)
            pf = x_fine * (1.0 - x_fine)
            G_xs_fit = sum(L[k] * pf * (1.0 - 2.0 * x_fine) ** k for k in range(rk_order))
            ax_rk.plot(
                x_fine,
                G_xs_fit * 96485,
                lw=2,
                color="#59A14F",
                label=f"RK fit (order {rk_order})",
            )

            ax_rk.axhline(0, color="black", lw=0.8, ls="--")
            ax_rk.set_xlabel("Composition")
            ax_rk.set_ylabel(r"$G_{xs}$ (J/mol)")
            ax_rk.set_title(f"{phase} — RK excess")
            ax_rk.legend(fontsize=8)

        fig.tight_layout()
        return fig, axes

    # ------------------------------------------------------------------
    # Data convergence overview
    # ------------------------------------------------------------------

    def plot_convergence(self, figsize=None):
        """
        Show which (phase, composition, temperature) calculations
        succeeded.  Each phase gets a subplot with temperature on
        the y-axis and composition on the x-axis.  Successful runs
        are shown as filled circles; missing data as open circles.

        Returns
        -------
        fig : matplotlib Figure
        axes : array of matplotlib Axes
        """
        n = len(self.phases)
        if figsize is None:
            figsize = (4 * n, 5)
        fig, axes = plt.subplots(1, n, figsize=figsize, sharey=True)
        if n == 1:
            axes = [axes]

        for ax, phase in zip(axes, self.phases):
            df_p = self.df.loc[self.df["phase"] == phase]
            for _, row in df_p.iterrows():
                comp = row["composition"]
                tarr = np.atleast_1d(row["temperature"])
                farr = np.atleast_1d(row["free_energy"])
                ok = ~np.isnan(farr.astype(float))
                ax.scatter(
                    np.full(ok.sum(), comp), tarr[ok], c="#4E79A7", s=15, zorder=3
                )
                if (~ok).any():
                    ax.scatter(
                        np.full((~ok).sum(), comp),
                        tarr[~ok],
                        facecolors="none",
                        edgecolors="#E15759",
                        s=15,
                        zorder=3,
                    )
            ax.set_title(phase)
            ax.set_xlabel("Composition")

        axes[0].set_ylabel("T (K)")
        fig.tight_layout()
        return fig, axes

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _require_calculated(self):
        if self.tangents is None:
            raise RuntimeError("No tangent data. Call .calculate() first.")

    def __repr__(self):
        n = len(self.phases)
        status = "calculated" if self.tangents is not None else "not calculated"
        return (
            f"PhaseDiagram(phases={self.phases}, "
            f"ref='{self.reference_element}', {status})"
        )

    def save(self, filename):
        """
        Save the PhaseDiagram object to a file using pickle.

        Parameters
        ----------
        filename : str
            Path to the output file (e.g. ``'phase_diagram.pkl'``).
        """
        with open(filename, "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def load(filename):
        """
        Load a PhaseDiagram object previously saved with :meth:`save`.

        Parameters
        ----------
        filename : str
            Path to the pickle file.

        Returns
        -------
        PhaseDiagram
        """
        with open(filename, "rb") as f:
            return pickle.load(f)

    def dump_df(self, filename):
        """
        Save the merged, reference-corrected DataFrame to a pickle file.

        This captures the data *after* gathering, cleaning, and reference
        correction but *before* any F(x) fitting, making it the ideal
        checkpoint for inspecting or re-using raw data across sessions.

        The DataFrame can be reloaded with ``pandas.read_pickle(filename)``.

        Parameters
        ----------
        filename : str
            Path to the output file (e.g. ``'raw_data.pkl'``).
        """
        self.df.to_pickle(filename)
