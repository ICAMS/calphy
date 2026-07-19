import numpy as np
from tqdm.notebook import trange
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import itertools
from itertools import combinations
from collections import defaultdict
import math
import copy
import os
import pickle
from calphy.composition_transformation import CompositionTransformation
import yaml
import matplotlib.patches as mpatches
import re
import json
from fractions import Fraction
from datetime import datetime

from calphy.integrators import kb

from scipy.spatial import ConvexHull
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter


EV_TO_J_MOL = 96485.33212331002
R_GAS_CONSTANT = 8.31446261815324

_SER_REF = {
    "AG": "FCC_A1",
    "AL": "FCC_A1",
    "AU": "FCC_A1",
    "CO": "HCP_A3",
    "CR": "BCC_A2",
    "CU": "FCC_A1",
    "FE": "BCC_A2",
    "MG": "HCP_A3",
    "MN": "BCC_A2",
    "MO": "BCC_A2",
    "NI": "FCC_A1",
    "PB": "FCC_A1",
    "PD": "FCC_A1",
    "PT": "FCC_A1",
    "SI": "DIA_A4",
    "SN": "BCT_A5",
    "TI": "HCP_A3",
    "V": "BCC_A2",
    "W": "BCC_A2",
    "ZN": "HCP_A3",
    "ZR": "HCP_A3",
}

_ELEMENT_SYMBOLS = [
    "AC", "AG", "AL", "AM", "AR", "AS", "AT", "AU", "B", "BA", "BE", "BH",
    "BI", "BK", "BR", "C", "CA", "CD", "CE", "CF", "CL", "CM", "CN", "CO",
    "CR", "CS", "CU", "DB", "DS", "DY", "ER", "ES", "EU", "F", "FE", "FL",
    "FM", "FR", "GA", "GD", "GE", "H", "HE", "HF", "HG", "HO", "HS", "I",
    "IN", "IR", "K", "KR", "LA", "LI", "LR", "LU", "LV", "MC", "MD", "MG",
    "MN", "MO", "MT", "N", "NA", "NB", "ND", "NE", "NH", "NI", "NO", "NP",
    "O", "OG", "OS", "P", "PA", "PB", "PD", "PM", "PO", "PR", "PT", "PU",
    "RA", "RB", "RE", "RF", "RG", "RH", "RN", "RU", "S", "SB", "SC", "SE",
    "SG", "SI", "SM", "SN", "SR", "TA", "TB", "TC", "TE", "TH", "TI", "TL",
    "TM", "TS", "U", "V", "W", "XE", "Y", "YB", "ZN", "ZR",
]

_TDB_METADATA_PREFIX = "$ CALPHY_TDB_METADATA "


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
    dict_list, peak_cutoff=0.003, plot=False, remove_self_tangents_for=[]
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
    A = np.column_stack([np.ones_like(T), T, T * np.log(T), T**2, T**3, 1.0 / T])
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
    T = np.asarray(T, dtype=float)
    return (
        coeffs[0]
        + coeffs[1] * T
        + coeffs[2] * T * np.log(T)
        + coeffs[3] * T**2
        + coeffs[4] * T**3
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

    # Redlich-Kister excess — L can be 1-D (T-independent) or 2-D (T-dependent)
    # Shape (rk_order,)   → L_k constant
    # Shape (rk_order, 2) → L_k(T) = a_k + b_k·T  (legacy)
    # Shape (rk_order, 6) → L_k(T) = CALPHAD poly6 (preferred)
    pf = x * (1.0 - x)
    L = np.asarray(L)
    if L.ndim == 2 and L.shape[1] == 6:
        G_xs = sum(_eval_calphad_poly6(L[k], T) * pf * (1.0 - 2.0 * x) ** k for k in range(len(L)))
    elif L.ndim == 2:
        G_xs = sum((L[k, 0] + L[k, 1] * T) * pf * (1.0 - 2.0 * x) ** k for k in range(len(L)))
    else:
        G_xs = sum(L[k] * pf * (1.0 - 2.0 * x) ** k for k in range(len(L)))

    return G_lin + G_ideal + G_xs


def _phase_is_liquid(phase_name):
    """Heuristic: does this user-supplied phase name refer to a liquid?

    Used by :meth:`PhaseDiagram.to_tdb` to decide which full-range phases
    skip the SGTE-conventional ``(M):(VA)`` 2-sublattice form (vacancy
    interstitial sublattice), which is universal for crystalline solid
    solutions but never used for liquids.
    """
    n = str(phase_name).lower().strip()
    return n in {"l", "liq", "lqd", "liquid"} or n.startswith("liq")


def _tdb_phase_name(name):
    phase_name = re.sub(r"[^0-9A-Za-z_]+", "_", str(name).upper()).strip("_")
    if not phase_name:
        raise ValueError("Phase names must not be empty")
    if phase_name[0].isdigit():
        phase_name = f"P_{phase_name}"
    return phase_name


def _format_tdb_expr_poly6(coeffs):
    a, b, c, d, e, f = np.asarray(coeffs, dtype=float)
    terms = []
    if abs(a) > 1e-10:
        terms.append(f"{a:+.8g}")
    if abs(b) > 1e-10:
        terms.append(f"{b:+.8g}*T")
    if abs(c) > 1e-10:
        terms.append(f"{c:+.8g}*T*LN(T)")
    if abs(d) > 1e-14:
        terms.append(f"{d:+.8g}*T**2")
    if abs(e) > 1e-18:
        terms.append(f"{e:+.8g}*T**3")
    if abs(f) > 1e-6:
        terms.append(f"{f:+.8g}*T**(-1)")
    return "".join(terms) if terms else "0"


def _format_tdb_expr_linear(a, b):
    terms = []
    if abs(a) > 1e-10:
        terms.append(f"{a:+.8g}")
    if abs(b) > 1e-10:
        terms.append(f"{b:+.8g}*T")
    return "".join(terms) if terms else "0"


_TDB_TERM_RE = re.compile(
    r"([+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)"  # signed number, optional exponent
    r"((?:\*T\*\*\(?-?\d+\)?|\*T\*(?:LN|ln)\(T\)|\*T)?)"  # optional factor
)


def _parse_tdb_poly6_expr(expr):
    """Parse a TDB polynomial expression into six CALPHAD-poly coefficients.

    Recognised terms: constant, ``*T``, ``*T*LN(T)``, ``*T**2``, ``*T**3``,
    ``*T**(-1)``.  Whitespace is ignored.  Unrecognised factors are silently
    dropped — sufficient for round-tripping :func:`_format_tdb_expr_poly6`
    output, not a general TDB parser.
    """
    text = expr.replace(" ", "").replace("\t", "").replace("\n", "")
    coeffs = {"a": 0.0, "b": 0.0, "c": 0.0, "d": 0.0, "e": 0.0, "f": 0.0}
    pos = 0
    while pos < len(text):
        m = _TDB_TERM_RE.match(text, pos)
        if not m or m.end() == pos:
            pos += 1
            continue
        num_str = m.group(1)
        factor = m.group(2).upper()
        try:
            val = float(num_str)
        except ValueError:
            pos = m.end()
            continue
        if factor == "":
            coeffs["a"] += val
        elif factor == "*T":
            coeffs["b"] += val
        elif factor.endswith("*LN(T)"):
            coeffs["c"] += val
        elif factor in ("*T**2", "*T**(2)"):
            coeffs["d"] += val
        elif factor in ("*T**3", "*T**(3)"):
            coeffs["e"] += val
        elif factor in ("*T**(-1)", "*T**-1"):
            coeffs["f"] += val
        pos = m.end()
    return np.array(
        [coeffs["a"], coeffs["b"], coeffs["c"], coeffs["d"], coeffs["e"], coeffs["f"]]
    )


def _phase_data_as_points(df, phase):
    rows = []
    df_phase = df.loc[df["phase"] == phase]
    for _, row in df_phase.iterrows():
        temperatures = np.atleast_1d(np.asarray(row["temperature"], dtype=float))
        free_energies = np.atleast_1d(np.asarray(row["free_energy"], dtype=float))
        if len(temperatures) != len(free_energies):
            continue
        composition = float(row["composition"])
        for temperature, free_energy in zip(temperatures, free_energies):
            if np.isfinite(temperature) and np.isfinite(free_energy):
                rows.append((composition, float(temperature), float(free_energy)))
    return pd.DataFrame(rows, columns=["composition", "temperature", "free_energy"])


def _infer_binary_elements(phases, reference_element):
    ref = reference_element.upper()
    found = set()
    element_set = set(_ELEMENT_SYMBOLS)
    for phase in phases:
        token = re.sub(r"[^A-Za-z]", "", str(phase)).upper()
        if ref not in token:
            continue
        pos = 0
        while pos < len(token):
            two = token[pos : pos + 2]
            one = token[pos : pos + 1]
            if two in element_set:
                found.add(two)
                pos += 2
            elif one in element_set:
                found.add(one)
                pos += 1
            else:
                pos += 1
    found.add(ref)
    others = sorted(found - {ref})
    if len(others) == 1:
        return [others[0], ref]
    raise ValueError(
        "Could not infer both binary elements from phase names. "
        "Pass elements=('A', reference_element) explicitly."
    )


def _stoichiometry_from_reference_composition(x_reference, max_denominator=8):
    frac = Fraction(float(x_reference)).limit_denominator(max_denominator)
    n_ref = frac.numerator
    total = frac.denominator
    m_other = total - n_ref
    if m_other <= 0 or n_ref <= 0:
        raise ValueError("Limited-range compound phases require an interior stoichiometry")
    return int(m_other), int(n_ref)


def _stoichiometry_from_phase_formula(phase, el_a, el_b):
    token = re.sub(r"[^0-9A-Za-z]", "", str(phase)).upper()
    pattern = re.compile(f"({re.escape(el_a)}|{re.escape(el_b)})([0-9]*)")
    counts = {el_a: 0, el_b: 0}
    matched = ""
    for element, number in pattern.findall(token):
        counts[element] += int(number) if number else 1
        matched += element + number
    if matched != token or counts[el_a] <= 0 or counts[el_b] <= 0:
        return None
    return int(counts[el_a]), int(counts[el_b])


def _infer_compound_stoichiometry(phase, points, el_a, el_b):
    from_formula = _stoichiometry_from_phase_formula(phase, el_a, el_b)
    if from_formula is not None:
        return from_formula
    x_ref = float(np.median(points["composition"]))
    return _stoichiometry_from_reference_composition(x_ref)


def _safe_xlogx(y):
    y = np.asarray(y, dtype=float)
    return np.where(y > 0.0, y * np.log(np.maximum(y, 1e-300)), 0.0)


def _site_ratio_weights(site_ratios):
    site_ratios = np.asarray(site_ratios, dtype=float)
    if site_ratios.ndim != 1 or len(site_ratios) < 2:
        raise ValueError("A pseudo-sublattice compound model needs at least two sublattices")
    if np.any(site_ratios <= 0):
        raise ValueError("Sublattice site ratios must be positive")
    return site_ratios / np.sum(site_ratios)


def _normalise_compound_sublattice_model(model):
    if isinstance(model, dict):
        model = model.get("site_ratios")
    if model is None:
        raise ValueError("compound_sublattice_models entries must define site_ratios")
    site_ratios = tuple(int(x) if float(x).is_integer() else float(x) for x in model)
    _site_ratio_weights(site_ratios)
    return site_ratios


def _infer_compound_sublattice_model(m_other, n_ref, mode="compact"):
    if mode == "compact":
        return (int(m_other), int(n_ref))
    if mode == "expanded":
        return tuple([1] * int(m_other + n_ref))
    raise ValueError("compound_sublattice_mode must be 'compact' or 'expanded'")


def _sublattice_endmember_bits(n_sublattices):
    return np.asarray(list(itertools.product([0, 1], repeat=n_sublattices)), dtype=int)


def _sublattice_grid(n_sublattices, grid_points):
    grid = np.linspace(0.0, 1.0, int(grid_points))
    if n_sublattices == 2:
        return grid[:, None]
    n_combinations = int(grid_points) ** (n_sublattices - 1)
    if n_combinations > 250000:
        raise ValueError(
            "Pseudo-sublattice grid is too large. Reduce pseudo_grid_points or use fewer sublattices."
        )
    return np.asarray(list(itertools.product(grid, repeat=n_sublattices - 1)), dtype=float)


def _pseudo_sublattice_minimum(theta, x_reference, temperature, site_ratios, grid_points):
    x_reference = np.asarray(x_reference, dtype=float)
    temperature = np.asarray(temperature, dtype=float)
    ratios = _site_ratio_weights(site_ratios)
    n_sublattices = len(ratios)
    endmember_bits = _sublattice_endmember_bits(n_sublattices)
    theta = np.asarray(theta, dtype=float).reshape(len(endmember_bits), 2)

    known_grid = _sublattice_grid(n_sublattices, grid_points)
    y_known = np.broadcast_to(known_grid[None, :, :], (len(x_reference), len(known_grid), n_sublattices - 1))
    y_last = (
        x_reference[:, None]
        - np.sum(y_known * ratios[:-1][None, None, :], axis=2)
    ) / ratios[-1]
    y_ref = np.concatenate([y_known, y_last[:, :, None]], axis=2)
    valid = (y_ref >= -1e-12) & (y_ref <= 1.0 + 1e-12)
    valid = np.all(valid, axis=2)
    y_ref = np.clip(y_ref, 0.0, 1.0)
    y_other = 1.0 - y_ref

    g_ref = np.zeros((len(x_reference), len(known_grid)), dtype=float)
    for idx, bits in enumerate(endmember_bits):
        probability = np.prod(np.where(bits[None, None, :] == 1, y_ref, y_other), axis=2)
        g_ref += probability * (theta[idx, 0] + theta[idx, 1] * temperature[:, None])
    g_id = R_GAS_CONSTANT * temperature[:, None] * (
        np.sum(ratios[None, None, :] * (_safe_xlogx(y_ref) + _safe_xlogx(y_other)), axis=2)
    )
    values = g_ref + g_id
    values[~valid] = np.inf
    result = np.min(values, axis=1)
    return result


def _fit_linear_tdb_parameter(temperature, free_energy_j_mol):
    design = np.column_stack([np.ones_like(temperature), temperature])
    coeffs, _, _, _ = np.linalg.lstsq(design, free_energy_j_mol, rcond=None)
    return coeffs


def _fit_pseudo_sublattice_phase(
    points,
    site_ratios,
    target_x_reference=None,
    max_points=600,
    grid_points=81,
    wrong_endmember_penalty=50000.0,
    regularization=1e-8,
):
    from scipy.optimize import least_squares

    if len(points) < 4:
        raise ValueError("At least four finite points are needed to fit a pseudo-sublattice phase")

    fit_points = points.sort_values(["composition", "temperature"]).copy()
    if len(fit_points) > max_points:
        keep = np.linspace(0, len(fit_points) - 1, max_points).round().astype(int)
        fit_points = fit_points.iloc[np.unique(keep)]

    x = fit_points["composition"].to_numpy(dtype=float)
    temperature = fit_points["temperature"].to_numpy(dtype=float)
    target = fit_points["free_energy"].to_numpy(dtype=float) * EV_TO_J_MOL
    site_ratios = _normalise_compound_sublattice_model(site_ratios)
    ratios = _site_ratio_weights(site_ratios)
    if target_x_reference is None:
        target_x_reference = float(np.median(points["composition"]))
    stoich = float(target_x_reference)
    distances = np.abs(points["composition"].to_numpy(dtype=float) - stoich)
    ordered_points = points.iloc[distances <= max(np.min(distances), 1e-12)]
    ordered_temperature = ordered_points["temperature"].to_numpy(dtype=float)
    ordered_energy = ordered_points["free_energy"].to_numpy(dtype=float) * EV_TO_J_MOL
    ordered_a, ordered_b = _fit_linear_tdb_parameter(ordered_temperature, ordered_energy)
    high_a = ordered_a + float(wrong_endmember_penalty)
    endmember_bits = _sublattice_endmember_bits(len(site_ratios))
    endmember_compositions = endmember_bits @ ratios
    ordered_mask = np.abs(endmember_compositions - stoich) <= max(
        np.min(np.abs(endmember_compositions - stoich)), 1e-12
    )
    theta0 = []
    for is_ordered in ordered_mask:
        theta0.extend([ordered_a if is_ordered else high_a, ordered_b])
    theta0 = np.asarray(theta0, dtype=float)
    scale = max(1000.0, float(np.nanstd(target)))

    def residual(theta):
        pred = _pseudo_sublattice_minimum(theta, x, temperature, site_ratios, grid_points)
        resid = (pred - target) / scale
        if regularization > 0:
            denom = np.maximum(np.abs(theta0), 1.0)
            reg = np.sqrt(regularization) * (theta - theta0) / denom
            resid = np.concatenate([resid, reg])
        return resid

    result = least_squares(residual, theta0, max_nfev=400)
    pred = _pseudo_sublattice_minimum(result.x, x, temperature, site_ratios, grid_points)
    rms = float(np.sqrt(np.mean((pred - target) ** 2)))
    return {
        "theta": result.x,
        "site_ratios": list(site_ratios),
        "rms_j_mol": rms,
        "n_fit_points": int(len(fit_points)),
        "success": bool(result.success),
        "message": result.message,
    }


def _fit_compound_two_sublattice(points, host_surface, m_first, n_second,
                                  anti_site_penalty_j_mol=80000.0,
                                  pure_sublattice_penalty_j_mol=80000.0):
    """
    Fit a 2-sublattice ``(A,B):(A,B)`` compound with site ratios ``(m, n)``
    where the stoichiometric endmember ``A:B`` corresponds to mole fraction
    ``x_A = m / (m + n)``.

    Only the stoichiometric endmember ``G(A:B; 0)`` is fit to calphy data
    (linear in T). The pure-element-on-each-sublattice endmembers ``A:A``
    and ``B:B`` are set to the host solution phase's ``G_A(T)`` and
    ``G_B(T)`` plus a large per-atom penalty ``pure_sublattice_penalty_j_mol``
    so that the CEF cannot lower G below the host through configurational
    disordering away from the stoichiometric column. The anti-site
    endmember ``G(B:A; 0)`` equals the stoichiometric one plus
    ``anti_site_penalty_j_mol``. Together these two penalties confine the
    compound to a narrow stability window around ``x_A = m/(m+n)``.

    All endmember energies are returned in J / mol-formula-unit (multiplied
    by ``m + n``) since pycalphad's TDB convention writes ``PARAMETER
    G(phase, ...; 0)`` per formula unit.
    """
    m, n = int(m_first), int(n_second)
    N = m + n
    # Composition convention: data column is mole fraction of the reference
    # element (the second-sublattice element), so the stoichiometric
    # composition in data coordinates is n/(m+n).
    x_stoich_data = n / N
    x_stoich_first = m / N  # x of the first (non-reference) element

    x = points["composition"].to_numpy(dtype=float)
    T_all = points["temperature"].to_numpy(dtype=float)
    G_all = points["free_energy"].to_numpy(dtype=float)  # eV / atom

    distances = np.abs(x - x_stoich_data)
    near = distances <= max(np.min(distances), 1e-12)
    T_stoich = T_all[near]
    G_stoich = G_all[near] * EV_TO_J_MOL  # J / mol-atom

    G_per_formula = G_stoich * N
    a_AB, b_AB = _fit_linear_tdb_parameter(T_stoich, G_per_formula)

    coeffs_A_per_formula = np.asarray(host_surface["coeffs_A"], dtype=float) * EV_TO_J_MOL * N
    coeffs_B_per_formula = np.asarray(host_surface["coeffs_B"], dtype=float) * EV_TO_J_MOL * N
    pure_penalty_per_formula = float(pure_sublattice_penalty_j_mol) * N

    pred = a_AB + b_AB * T_stoich
    rms = float(np.sqrt(np.mean((pred - G_per_formula) ** 2)) / N)

    return {
        "site_ratios": (m, n),
        "stoichiometry": (m, n),
        "x_stoich": x_stoich_first,
        "theta_AB": (a_AB, b_AB),  # stoichiometric endmember, fit to data
        "theta_AA_poly6": coeffs_A_per_formula,  # pure A on both: host G_A * N
        "theta_BB_poly6": coeffs_B_per_formula,  # pure B on both: host G_B * N
        "pure_sublattice_penalty": pure_penalty_per_formula,
        "anti_site_penalty": float(anti_site_penalty_j_mol) * N,
        "rms_j_mol": rms,
        "n_fit_points": int(len(T_stoich)),
    }


def _fit_limited_range_surface(points, host_surface, rk_order=3,
                                anchor_outside_window=True,
                                anchor_margin=0.02,
                                anchor_count=8,
                                anchor_weight=1.0):
    """
    Fit a 1-sublattice (A,B) CALPHAD surface to a limited-composition-range
    phase, sharing the pure-element reference functions ``G_A(T)`` and
    ``G_B(T)`` with a host solution phase.

    Only the temperature-dependent Redlich-Kister excess parameters
    ``L_k(T) = a_k + b_k T`` are fit. This mirrors
    :meth:`PhaseDiagram.build_calphad_surface` but for phases that do not
    span the full composition range. The shared SER guarantees a consistent
    reference state with the host phase so the compound's stability against
    the host is determined purely by the excess term, which is exactly the
    quantity calphy data resolves in the compound's window.

    To keep the polynomial from blowing up outside the data window, anchor
    points are appended at evenly spaced compositions in
    ``[0, x_min - anchor_margin]`` and ``[x_max + anchor_margin, 1]`` where
    the compound's free energy is set equal to the host's free energy at the
    same ``(x, T)``. This sets the excess to zero outside the window so the
    compound coincides with the host (and is therefore not spuriously stable
    far from its stoichiometric range).

    The model evaluated by pycalphad is::

        G(x, T) = (1-x) G_A(T) + x G_B(T)
                  + R T [x ln x + (1-x) ln(1-x)]
                  + x (1-x) Sum_k (a_k + b_k T) (1 - 2x)^k
    """
    x = points["composition"].to_numpy(dtype=float)
    T = points["temperature"].to_numpy(dtype=float)
    G_data = points["free_energy"].to_numpy(dtype=float)  # eV/atom

    coeffs_A = host_surface["coeffs_A"]
    coeffs_B = host_surface["coeffs_B"]

    def _g_lin_ideal(x_arr, T_arr):
        G_A = _eval_calphad_poly6(coeffs_A, T_arr)
        G_B = _eval_calphad_poly6(coeffs_B, T_arr)
        G_lin = (1.0 - x_arr) * G_A + x_arr * G_B
        log_x = np.where(x_arr > 0, np.log(np.maximum(x_arr, 1e-300)), 0.0)
        log_1mx = np.where(x_arr < 1, np.log(np.maximum(1.0 - x_arr, 1e-300)), 0.0)
        G_ideal = kb * T_arr * (x_arr * log_x + (1.0 - x_arr) * log_1mx)
        return G_lin, G_ideal

    G_lin, G_ideal = _g_lin_ideal(x, T)
    G_xs = G_data - G_lin - G_ideal
    weights = np.ones_like(G_xs)

    if anchor_outside_window:
        x_min = float(np.min(x))
        x_max = float(np.max(x))
        unique_T = np.unique(T)
        anchors_x = []
        if x_min - anchor_margin > 1e-3:
            anchors_x.extend(np.linspace(1e-3, x_min - anchor_margin, anchor_count))
        if 1.0 - (x_max + anchor_margin) > 1e-3:
            anchors_x.extend(np.linspace(x_max + anchor_margin, 1.0 - 1e-3, anchor_count))
        if anchors_x:
            anchor_x = np.repeat(anchors_x, len(unique_T))
            anchor_T = np.tile(unique_T, len(anchors_x))
            # Anchor: compound G == host G at (x, T) -> excess == 0.
            x = np.concatenate([x, anchor_x])
            T = np.concatenate([T, anchor_T])
            G_xs = np.concatenate([G_xs, np.zeros_like(anchor_x)])
            weights = np.concatenate([weights, anchor_weight * np.ones_like(anchor_x)])

    pf = x * (1.0 - x)
    cols = []
    for k in range(rk_order):
        pk = pf * (1.0 - 2.0 * x) ** k
        cols.append(pk)
        cols.append(pk * T)
    basis = np.column_stack(cols)
    W = np.sqrt(weights)
    L_flat, _, _, _ = np.linalg.lstsq(basis * W[:, None], G_xs * W, rcond=None)
    L_coeffs = L_flat.reshape(rk_order, 2)
    rms = float(np.sqrt(np.mean((basis @ L_flat - G_xs) ** 2)) * EV_TO_J_MOL)
    return {
        "coeffs_A": np.asarray(coeffs_A, dtype=float),
        "coeffs_B": np.asarray(coeffs_B, dtype=float),
        "L_coeffs": L_coeffs,
        "rms_j_mol": rms,
        "n_fit_points": int(len(x)),
        "rk_order": int(rk_order),
    }


def _limited_phase_pseudo_points(
    df, phase, composition_interval, fit_order=4, n_T_max=41, n_x=41,
    ideal_configurational_entropy=True,
):
    """
    Build dense pseudo-data for a limited-range phase from calphy's own
    per-temperature polynomial fits (:func:`get_phase_free_energy` — the
    same smoothing :meth:`PhaseDiagram.calculate` applies to narrow
    phases).

    Fitting the TDB Redlich-Kister model against ~11 raw composition
    points leaves the interior of a high-order RK well underdetermined:
    the outside-window inequality constraints then tilt the well minimum
    away from the data (0.02+ in x for AuCu MD data), which shows up as
    a "boot-shaped" single-phase window near the order-disorder dome
    top.  Sampling the per-T polynomial fit densely across the window
    pins the interior shape — the RK fit then tracks the well-minimum
    location and depth of calphy's own construction to ~0.003 in x.

    Parameters
    ----------
    df : DataFrame
        The PhaseDiagram dataframe.
    phase : str
        Phase name.
    composition_interval : tuple
        (x_lo, x_hi) window for the phase.
    fit_order : int
        Polynomial order of the per-T fit (default 4 — the same default
        :meth:`PhaseDiagram.calculate` uses).
    n_T_max : int
        Maximum number of temperature slices to sample (evenly spaced
        over the unique data temperatures).
    n_x : int
        Composition grid points per temperature slice.

    Returns
    -------
    DataFrame with columns (composition, temperature, free_energy), or
    ``None`` if no temperature slice could be fit (too few compositions
    for *fit_order*) — callers should fall back to the raw data points.
    """
    sub = df.loc[df["phase"] == phase]
    all_T = sorted(
        {
            float(t)
            for arr in sub["temperature"]
            for t in np.atleast_1d(np.asarray(arr, dtype=float))
            if np.isfinite(t)
        }
    )
    if not all_T:
        return None
    if len(all_T) > n_T_max:
        idx = np.linspace(0, len(all_T) - 1, int(n_T_max)).round().astype(int)
        all_T = [all_T[i] for i in np.unique(idx)]

    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for T in all_T:
            d = get_phase_free_energy(
                df,
                phase,
                T,
                composition_interval=tuple(composition_interval),
                ideal_configurational_entropy=ideal_configurational_entropy,
                fit_order=fit_order,
                composition_grid=int(n_x),
            )
            if d is None:
                continue
            for x, G in zip(d["composition"], d["free_energy"]):
                if np.isfinite(x) and np.isfinite(G):
                    rows.append((float(x), float(T), float(G)))
    if not rows:
        return None
    return pd.DataFrame(
        rows, columns=["composition", "temperature", "free_energy"]
    )


def _limited_phase_min_depth_by_T(points, host_surface, fit=None):
    """
    Per-temperature minimum of ``G_phase - G_host`` (eV/atom) across the
    phase's composition points — the depth of the ordered-phase well
    below the host.  With ``fit=None`` the depth comes from the data
    points themselves; otherwise the fit surface is evaluated at the
    same (x, T) points.

    Returns (T_sorted, depth) arrays.
    """
    x = points["composition"].to_numpy(dtype=float)
    T = points["temperature"].to_numpy(dtype=float)
    if fit is None:
        G = points["free_energy"].to_numpy(dtype=float)
    else:
        G = np.array(
            [
                _eval_calphad_surface_at(fit, np.array([xi]), ti)[0]
                for xi, ti in zip(x, T)
            ]
        )
    G_host = np.array(
        [
            _eval_calphad_surface_at(host_surface, np.array([xi]), ti)[0]
            for xi, ti in zip(x, T)
        ]
    )
    diff = G - G_host
    T_unique = np.unique(T)
    depth = np.array([np.min(diff[T == t]) for t in T_unique])
    return T_unique, depth


def _first_zero_crossing(T, depth):
    """T where depth first crosses from <=0 to >0 (linear interp), or None."""
    for i in range(len(depth) - 1):
        if depth[i] <= 0.0 < depth[i + 1]:
            f = depth[i] / (depth[i] - depth[i + 1])
            return float(T[i] + f * (T[i + 1] - T[i]))
    return None


def _fit_limited_range_surface_bounded(
    points,
    host_surface,
    rk_order=3,
    n_T_terms=2,
    n_constraint_x=32,
    n_constraint_T=24,
    constraint_T_range=None,
    constraint_margin=0.02,
    bound_slack=-0.001,
    convexity_constraint=True,
    n_convexity_x=25,
    ridge_lambda=1e-8,
):
    """
    Fit a 1-sublattice (A,B) RK surface to a limited-window phase, with the
    inequality constraint ``G_phase(x, T) >= G_host(x, T) - bound_slack``
    enforced at sample points *outside* the data window.

    Inside the window the L parameters are free to follow the data (which
    sits below the host).  Outside the window, the polynomial is forced to
    sit at or above the host, so the limited-range phase can never become
    spuriously stable far from its real composition range — addressing the
    failure mode where the RK tail drops below G_host (and pycalphad picks
    the phase) at x → 0 or x → 1.

    Mathematically: since the phase and host share G_A(T), G_B(T) and the
    ideal-entropy term, ``G_phase - G_host == G_xs_phase - G_xs_host``.  The
    constraint reduces to a linear inequality
    ``basis(x, T) @ L >= G_xs_host(x, T) - bound_slack`` for each constraint
    point.  We solve the constrained QP via SLSQP, warm-started from the
    unconstrained least-squares solution.

    Parameters
    ----------
    points : DataFrame
        Per (composition, temperature) data points for the phase.
    host_surface : dict
        Output of :meth:`PhaseDiagram.build_calphad_surface` for the host
        (full-range) phase.  Must contain ``coeffs_A``, ``coeffs_B``,
        ``L_coeffs``.
    rk_order : int
        RK order (default 3).
    n_T_terms : int
        Temperature basis terms per L_k: ``2`` → ``a + b·T`` (default),
        ``3`` → ``a + b·T + c·T·ln T``.  Three terms capture the
        curvature of the ordered-phase stability window in T (e.g. the
        steepening of the well below ~500 K) that a linear form misses.
    n_constraint_x : int
        How many composition points to enforce the inequality at, total
        across the two outside-window regions (default 20).
    constraint_margin : float
        Distance from the window edge at which to start placing constraints.
    bound_slack : float
        Slack on the inequality (eV/atom); positive values let G_phase dip
        slightly below G_host, negative values force a margin.  Default
        ``-0.001`` (a 1 meV/atom strict margin): the constraint is only
        enforced at sample points, and between samples the polynomial can
        wiggle a fraction of a meV below the host — enough for pycalphad
        to draw hairline spurious lenses near the pure endpoints, where
        the phase and host energies are degenerate by construction.  The
        margin buries those wiggles.
    convexity_constraint : bool
        If True (default), additionally require ``d²G/dx² >= 0`` at
        ``n_convexity_x`` compositions *inside* the data window (crossed
        with the constraint temperatures).  A narrow ordered phase is a
        single-phase field — its G(x) well must be convex.  Without this,
        a high-order RK fit can develop a slightly W-shaped well bottom,
        which pycalphad resolves as a spurious miscibility gap inside the
        phase (extra ``#2`` composition sets in the diagram).  The
        constraint is linear in L (G_xs is linear in L and the ideal
        term's curvature is fixed), so it slots into the same SLSQP
        problem.
    ridge_lambda : float
        Scaled Tikhonov ridge on the L coefficients (default ``1e-8``).
        Without it, high-order fits can satisfy the data and every
        sampled constraint with absurd cancelling magnitudes (L terms of
        ±thousands of eV summing to meV) whose residual wiggles
        *between* the constraint samples dip below the host and draw
        hairline phantom lenses.  Each coefficient is scaled by its
        natural magnitude (1, T_mid, T_mid·lnT_mid) so the penalty is
        dimensionally even.  Keep it weak: narrow windows legitimately
        need large high-order L (``(1-2x)^k`` is tiny there), and a
        strong ridge biases the well's T-slope, shifting dissolution
        temperatures by tens of K.

    Returns
    -------
    dict
        Same keys as :func:`_fit_limited_range_surface` plus
        ``constraint_violation_max`` (largest violation in eV/atom, ≤ 0 if
        feasible) and ``n_constraints``.  ``L_coeffs`` is returned in the
        zero-padded six-term CALPHAD layout ``(rk_order, 6)`` — columns
        ``[a, b·T, c·T·lnT, 0, 0, 0]`` — so it evaluates directly with
        :func:`_eval_calphad_surface_at`.
    """
    from scipy.optimize import minimize

    x = points["composition"].to_numpy(dtype=float)
    T = points["temperature"].to_numpy(dtype=float)
    G_data = points["free_energy"].to_numpy(dtype=float)

    coeffs_A = np.asarray(host_surface["coeffs_A"], dtype=float)
    coeffs_B = np.asarray(host_surface["coeffs_B"], dtype=float)

    def _g_lin_ideal(x_arr, T_arr):
        G_A = _eval_calphad_poly6(coeffs_A, T_arr)
        G_B = _eval_calphad_poly6(coeffs_B, T_arr)
        G_lin = (1.0 - x_arr) * G_A + x_arr * G_B
        log_x = np.where(x_arr > 0, np.log(np.maximum(x_arr, 1e-300)), 0.0)
        log_1mx = np.where(x_arr < 1, np.log(np.maximum(1.0 - x_arr, 1e-300)), 0.0)
        G_ideal = kb * T_arr * (x_arr * log_x + (1.0 - x_arr) * log_1mx)
        return G_lin, G_ideal

    if n_T_terms not in (2, 3):
        raise ValueError(f"n_T_terms must be 2 or 3, got {n_T_terms}")

    def _basis(x_arr, T_arr):
        pf = x_arr * (1.0 - x_arr)
        cols = []
        for k in range(rk_order):
            pk = pf * (1.0 - 2.0 * x_arr) ** k
            cols.append(pk)
            cols.append(pk * T_arr)
            if n_T_terms >= 3:
                cols.append(pk * T_arr * np.log(T_arr))
        return np.column_stack(cols)

    # ---- data side: linear LSQ target for G_xs at data points ----
    G_lin, G_ideal = _g_lin_ideal(x, T)
    G_xs_target = G_data - G_lin - G_ideal
    A_data = _basis(x, T)

    # ---- ridge rows: suppress huge cancelling L magnitudes ----
    if ridge_lambda and ridge_lambda > 0:
        T_mid = float(np.mean(T))
        scales = [1.0, T_mid, T_mid * np.log(T_mid)][:n_T_terms]
        ridge = np.sqrt(float(ridge_lambda)) * np.diag(
            np.tile(np.asarray(scales), rk_order)
        )
        A_data = np.vstack([A_data, ridge])
        G_xs_target = np.concatenate(
            [G_xs_target, np.zeros(ridge.shape[0])]
        )

    # ---- constraint side: sample points outside the window ----
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    half = max(1, n_constraint_x // 2)
    # Each outside region gets a linear grid plus a geometric grid
    # clustered toward the pure endpoint: G_phase - G_host tends to zero
    # there by construction (shared GHSER, x(1-x) excess prefactor), so
    # sub-meV polynomial wiggles between widely spaced samples are
    # enough to create spurious hairline lenses near x = 0 / x = 1.
    cx_list = []
    if x_min - constraint_margin > 1e-3:
        hi = x_min - constraint_margin
        cx_list.append(np.linspace(1e-3, hi, half))
        cx_list.append(np.geomspace(1e-3, hi, half))
    if 1.0 - (x_max + constraint_margin) > 1e-3:
        lo = x_max + constraint_margin
        cx_list.append(np.linspace(lo, 1.0 - 1e-3, half))
        cx_list.append(1.0 - np.geomspace(1e-3, 1.0 - lo, half))
    constraint_x = (
        np.unique(np.concatenate(cx_list)) if cx_list else np.empty(0, dtype=float)
    )

    # Sample T uniformly across the constraint range, not just at the data
    # temperatures.  pycalphad evaluates G(x, T) at temperatures the data may
    # never have seen (extrapolation beyond the highest sampled T), so the
    # constraint G_phase >= G_host must hold across the full TDB validity
    # window or the polynomial will undershoot the host where unconstrained.
    if constraint_T_range is None:
        constraint_T_range = (float(np.min(T)), float(np.max(T)))
    constraint_T = np.linspace(
        float(constraint_T_range[0]),
        float(constraint_T_range[1]),
        int(n_constraint_T),
    )
    if len(constraint_x) and len(constraint_T):
        cx, cT = np.meshgrid(constraint_x, constraint_T, indexing="ij")
        cx = cx.ravel()
        cT = cT.ravel()
        A_con = _basis(cx, cT)
        # Vectorised host G_xs at constraint points (T-element-wise).  The
        # host's L_coeffs may be 1-D (T-independent), 2-D linear (a,b), or
        # 2-D six-term — match _eval_calphad_surface_at's logic.
        L_host = np.asarray(host_surface["L_coeffs"])
        pf_c = cx * (1.0 - cx)
        G_xs_host_c = np.zeros_like(cx)
        if L_host.ndim == 2 and L_host.shape[1] == 6:
            for k in range(L_host.shape[0]):
                L_k_T = _eval_calphad_poly6(L_host[k], cT)
                G_xs_host_c += pf_c * (1.0 - 2.0 * cx) ** k * L_k_T
        elif L_host.ndim == 2:
            for k in range(L_host.shape[0]):
                L_k_T = L_host[k, 0] + L_host[k, 1] * cT
                G_xs_host_c += pf_c * (1.0 - 2.0 * cx) ** k * L_k_T
        else:
            for k in range(len(L_host)):
                G_xs_host_c += pf_c * (1.0 - 2.0 * cx) ** k * L_host[k]
        rhs_con = G_xs_host_c - float(bound_slack)
    else:
        A_con = np.zeros((0, n_T_terms * rk_order))
        rhs_con = np.zeros(0)
    n_outside = len(rhs_con)

    # ---- convexity constraints inside the data window ----
    # Require G''(x, T) >= 0 with G = G_lin + G_ideal + G_xs.  G_lin is
    # linear in x, G_ideal'' = kb·T·(1/x + 1/(1-x)) is fixed, and G_xs''
    # is linear in L, so the constraint is  B''(x,T) @ L >= -G_ideal''.
    # B'' is evaluated by central differences of the basis.
    if convexity_constraint:
        h = 1e-4
        cvx_x = np.linspace(x_min + h, x_max - h, int(n_convexity_x))
        cvx_T = np.linspace(
            float(constraint_T_range[0]),
            float(constraint_T_range[1]),
            int(n_constraint_T),
        )
        vx, vT = np.meshgrid(cvx_x, cvx_T, indexing="ij")
        vx = vx.ravel()
        vT = vT.ravel()
        B_pp = (_basis(vx + h, vT) - 2.0 * _basis(vx, vT) + _basis(vx - h, vT)) / h**2
        G_ideal_pp = kb * vT * (1.0 / vx + 1.0 / (1.0 - vx))
        A_con = np.vstack([A_con, B_pp])
        rhs_con = np.concatenate([rhs_con, -G_ideal_pp])

    # ---- unconstrained warm start ----
    L0, _, _, _ = np.linalg.lstsq(A_data, G_xs_target, rcond=None)

    # ---- constrained QP via SLSQP ----
    def obj(L):
        r = A_data @ L - G_xs_target
        return 0.5 * float(r @ r)

    def jac(L):
        r = A_data @ L - G_xs_target
        return A_data.T @ r

    if len(rhs_con):
        cons = [{
            "type": "ineq",
            "fun": lambda L, A=A_con, b=rhs_con: A @ L - b,
            "jac": lambda L, A=A_con: A,
        }]
        result = minimize(
            obj, L0, jac=jac, method="SLSQP", constraints=cons,
            options={"maxiter": 300, "ftol": 1e-10},
        )
        L_flat = result.x
        success = bool(result.success)
        # Report the G_phase >= G_host violation (eV/atom) separately from
        # the convexity rows, whose residuals live in eV/atom per x² units.
        slack_all = rhs_con - A_con @ L_flat
        violation_max = (
            float(np.max(slack_all[:n_outside])) if n_outside else 0.0
        )
        convexity_violation_max = (
            float(np.max(slack_all[n_outside:]))
            if len(slack_all) > n_outside
            else 0.0
        )
    else:
        L_flat = L0
        success = True
        violation_max = 0.0
        convexity_violation_max = 0.0

    # Zero-padded six-term CALPHAD layout: [a, b, c, 0, 0, 0] per L_k.
    L_by_k = L_flat.reshape(rk_order, n_T_terms)
    L_coeffs = np.zeros((rk_order, 6))
    L_coeffs[:, :n_T_terms] = L_by_k
    # RMS over the data rows only (excluding any ridge rows).
    n_data = len(x)
    rms = float(
        np.sqrt(
            np.mean((A_data[:n_data] @ L_flat - G_xs_target[:n_data]) ** 2)
        )
        * EV_TO_J_MOL
    )
    return {
        "coeffs_A": coeffs_A,
        "coeffs_B": coeffs_B,
        "L_coeffs": L_coeffs,
        "rms_j_mol": rms,
        "n_fit_points": int(len(x)),
        "rk_order": int(rk_order),
        "n_T_terms": int(n_T_terms),
        "n_constraints": int(len(rhs_con)),
        "constraint_violation_max": violation_max,
        "convexity_violation_max": convexity_violation_max,
        "success": success,
    }


def _fit_solution_with_ghser_reference(
    phase_data, host_surface, rk_order=3, l_T_dependence="linear"
):
    """
    Fit a 1-sublattice (A,B) solution phase using GHSER references for the
    pure endpoints, simultaneously fitting:

      * linear pure-element offsets ``(a_A, b_A)`` and ``(a_B, b_B)`` such
        that ``G_phase(pure_A, T) = GHSER_A(T) + a_A + b_A·T``,
      * Redlich-Kister excess parameters ``L_k(T)``.

    The model is::

        G_phase(x, T) = (1-x)·[GHSER_A(T) + a_A + b_A·T]
                        + x·[GHSER_B(T) + a_B + b_B·T]
                        + k_B·T·[x·ln x + (1-x)·ln(1-x)]
                        + x·(1-x)·Σ_k L_k(T)·(1-2x)^k

    GHSER_A and GHSER_B come from ``host_surface['coeffs_A'/'coeffs_B']``
    (in eV/atom, six-term CALPHAD poly).  For the host phase itself, the
    offsets should come out essentially zero — useful as a sanity check.

    Fitting both offsets and L_k together is critical: in the SGTE
    formulation the L_k parameters depend on the choice of pure-element
    reference.  Refitting offsets only and re-using the old L_k (computed
    against the phase's own pure-G fit) produces inconsistent surfaces and
    spurious miscibility gaps in pycalphad.

    Parameters
    ----------
    phase_data : DataFrame
        Per-(x, T) rows for the phase (already filtered to one phase).
    host_surface : dict
        Output of build_calphad_surface for the host (SER) phase.
    rk_order : int
        Number of RK terms (default 3).
    l_T_dependence : str
        ``'linear'`` (default): ``L_k(T) = a + b·T``.
        ``'poly6'``: same 6-term CALPHAD polynomial as G_A/G_B.

    Returns
    -------
    dict with keys:
      ``offset_A``, ``offset_B`` : ndarray(2)  (a, b) in eV/atom
      ``L_coeffs`` : ndarray(rk_order, n_T)  L_k(T) coefficients
                     (n_T = 2 for linear, 6 for poly6)
      ``rms_j_mol``, ``n_fit_points``, ``rk_order``, ``l_T_dependence``
    """
    coeffs_A = np.asarray(host_surface["coeffs_A"], dtype=float)
    coeffs_B = np.asarray(host_surface["coeffs_B"], dtype=float)

    rows = []
    for _, r in phase_data.iterrows():
        x = float(r["composition"])
        T_row = np.atleast_1d(np.asarray(r["temperature"], dtype=float))
        G_row = np.atleast_1d(np.asarray(r["free_energy"], dtype=float))
        if len(T_row) != len(G_row):
            continue
        for ti, gi in zip(T_row, G_row):
            if np.isfinite(ti) and np.isfinite(gi):
                rows.append((x, float(ti), float(gi)))
    if len(rows) < 6:
        raise ValueError("Need at least 6 points for the GHSER fit")
    x_all = np.array([r[0] for r in rows], dtype=float)
    T_all = np.array([r[1] for r in rows], dtype=float)
    G_all = np.array([r[2] for r in rows], dtype=float)

    # Subtract the GHSER and ideal-entropy contributions from data — what's
    # left is what the offsets + L_k must explain.
    G_A_host = _eval_calphad_poly6(coeffs_A, T_all)
    G_B_host = _eval_calphad_poly6(coeffs_B, T_all)
    G_lin_host = (1.0 - x_all) * G_A_host + x_all * G_B_host
    log_x = np.where(x_all > 0, np.log(np.maximum(x_all, 1e-300)), 0.0)
    log_1mx = np.where(x_all < 1, np.log(np.maximum(1.0 - x_all, 1e-300)), 0.0)
    G_ideal = kb * T_all * (x_all * log_x + (1.0 - x_all) * log_1mx)
    target = G_all - G_lin_host - G_ideal

    # Design matrix: 4 offset columns + 2*rk_order or 6*rk_order L columns.
    cols = []
    # offset_A: contributes (1-x)*1, (1-x)*T
    cols.append(1.0 - x_all)
    cols.append((1.0 - x_all) * T_all)
    # offset_B: x*1, x*T
    cols.append(x_all)
    cols.append(x_all * T_all)
    pf = x_all * (1.0 - x_all)
    if l_T_dependence == "linear":
        for k in range(rk_order):
            pk = pf * (1.0 - 2.0 * x_all) ** k
            cols.append(pk)
            cols.append(pk * T_all)
    elif l_T_dependence == "poly6":
        log_T = np.log(T_all)
        for k in range(rk_order):
            pk = pf * (1.0 - 2.0 * x_all) ** k
            cols.append(pk)
            cols.append(pk * T_all)
            cols.append(pk * T_all * log_T)
            cols.append(pk * T_all ** 2)
            cols.append(pk * T_all ** 3)
            cols.append(pk / T_all)
    else:
        raise ValueError(f"l_T_dependence must be 'linear' or 'poly6'")

    basis = np.column_stack(cols)
    sol, _, _, _ = np.linalg.lstsq(basis, target, rcond=None)

    offset_A = sol[0:2]
    offset_B = sol[2:4]
    L_flat = sol[4:]
    n_T = 2 if l_T_dependence == "linear" else 6
    L_coeffs = L_flat.reshape(rk_order, n_T)
    rms = float(np.sqrt(np.mean((basis @ sol - target) ** 2)) * EV_TO_J_MOL)
    return {
        "coeffs_A": coeffs_A,
        "coeffs_B": coeffs_B,
        "offset_A": offset_A,
        "offset_B": offset_B,
        "L_coeffs": L_coeffs,
        "rms_j_mol": rms,
        "n_fit_points": int(len(x_all)),
        "rk_order": int(rk_order),
        "l_T_dependence": l_T_dependence,
    }


def _fit_solution_linear_L_regularized(
    phase_data, host_surface, rk_order=3, lam=1.0,
    fit_pure_element_offsets=False,
):
    """
    Refit L_k(T) = a + b·T for a full-range solution phase, using ridge
    regularisation to prevent the unstable huge-cancelling-magnitudes
    behaviour of an ordinary linear LSQ on noisy liquid data.

    When ``fit_pure_element_offsets=True`` also fits linear offsets
    ``(a_A, b_A)``, ``(a_B, b_B)`` such that the phase's pure-element G
    is ``GHSER_A + a_A + b_A·T`` and ``GHSER_B + a_B + b_B·T`` (SGTE
    convention: non-host phases reference the host's GHSER with a small
    lattice-stability offset).  GHSER_A/B come from ``host_surface``;
    offsets are fit simultaneously with L_k from the same data, both
    subject to the ridge penalty.  This is the form used by COST507,
    solders.tdb, and standard SGTE databases.

    Solves::

        min   ||A·θ - target||²  +  λ · Σ θ_i² · scale_i²

    Parameters
    ----------
    phase_data : DataFrame
        Per-(x, T) rows for the phase.
    host_surface : dict
        Provides ``coeffs_A``, ``coeffs_B`` (the GHSER polynomials).
    rk_order : int
        Number of RK terms (default 3).
    lam : float
        Ridge regularisation strength.  Larger → smaller magnitude L_k
        (and offsets).  0 reduces to ordinary LSQ.  Default 1.0.
    fit_pure_element_offsets : bool
        If True, fit pure-element linear offsets simultaneously.  Used
        for non-host phases (offsets ≠ 0).  For the host phase, leave
        False so the pure G is exactly the input host_surface — offsets
        would be ~0 anyway.

    Returns
    -------
    dict with keys ``coeffs_A``, ``coeffs_B`` (= host GHSER, unchanged),
    ``offset_A``, ``offset_B`` (each ndarray(2): (a, b) in eV/atom; zero
    if not fitted), ``L_coeffs`` (rk_order × 2), ``rms_j_mol``,
    ``n_fit_points``, ``rk_order``, ``lambda_used``.
    """
    coeffs_A = np.asarray(host_surface["coeffs_A"], dtype=float)
    coeffs_B = np.asarray(host_surface["coeffs_B"], dtype=float)

    rows = []
    for _, r in phase_data.iterrows():
        x = float(r["composition"])
        T_row = np.atleast_1d(np.asarray(r["temperature"], dtype=float))
        G_row = np.atleast_1d(np.asarray(r["free_energy"], dtype=float))
        if len(T_row) != len(G_row):
            continue
        for ti, gi in zip(T_row, G_row):
            if np.isfinite(ti) and np.isfinite(gi):
                rows.append((x, float(ti), float(gi)))
    min_pts = (4 if fit_pure_element_offsets else 0) + 2 * rk_order
    if len(rows) < min_pts:
        raise ValueError(
            f"Need at least {min_pts} points for the regularised fit."
        )
    x_all = np.array([r[0] for r in rows], dtype=float)
    T_all = np.array([r[1] for r in rows], dtype=float)
    G_all = np.array([r[2] for r in rows], dtype=float)

    G_A = _eval_calphad_poly6(coeffs_A, T_all)
    G_B = _eval_calphad_poly6(coeffs_B, T_all)
    G_lin = (1.0 - x_all) * G_A + x_all * G_B
    log_x = np.where(x_all > 0, np.log(np.maximum(x_all, 1e-300)), 0.0)
    log_1mx = np.where(x_all < 1, np.log(np.maximum(1.0 - x_all, 1e-300)), 0.0)
    G_ideal = kb * T_all * (x_all * log_x + (1.0 - x_all) * log_1mx)
    target = G_all - G_lin - G_ideal

    cols = []
    scale_list = []
    T_mid = float(np.mean(T_all))
    if fit_pure_element_offsets:
        # (1-x)·a_A, (1-x)·T·b_A, x·a_B, x·T·b_B
        cols.append(1.0 - x_all)
        scale_list.append(1.0)
        cols.append((1.0 - x_all) * T_all)
        scale_list.append(T_mid)
        cols.append(x_all)
        scale_list.append(1.0)
        cols.append(x_all * T_all)
        scale_list.append(T_mid)
    pf = x_all * (1.0 - x_all)
    for k in range(rk_order):
        pk = pf * (1.0 - 2.0 * x_all) ** k
        cols.append(pk)
        scale_list.append(1.0)
        cols.append(pk * T_all)
        scale_list.append(T_mid)
    A = np.column_stack(cols)
    scale = np.array(scale_list)

    # Solve (A.T A + λ·diag(scale²)) θ = A.T b
    AtA = A.T @ A
    Atb = A.T @ target
    reg = lam * np.diag(scale ** 2)
    theta = np.linalg.solve(AtA + reg, Atb)

    if fit_pure_element_offsets:
        offset_A = theta[0:2]
        offset_B = theta[2:4]
        L_flat = theta[4:]
    else:
        offset_A = np.zeros(2)
        offset_B = np.zeros(2)
        L_flat = theta
    L_coeffs = L_flat.reshape(rk_order, 2)
    rms = float(np.sqrt(np.mean((A @ theta - target) ** 2)) * EV_TO_J_MOL)
    return {
        "coeffs_A": coeffs_A,
        "coeffs_B": coeffs_B,
        "offset_A": offset_A,
        "offset_B": offset_B,
        "L_coeffs": L_coeffs,
        "rms_j_mol": rms,
        "n_fit_points": int(len(x_all)),
        "rk_order": int(rk_order),
        "lambda_used": float(lam),
    }


def _fit_pure_element_offset(df, phase, host_surface, element_side, endpoint_tol=0.05):
    """
    Fit a linear-in-T offset (a + b·T) that approximates
    ``G_phase(pure_element, T) - GHSER_element(T)`` from data.

    Used by :meth:`PhaseDiagram.to_tdb` when writing SGTE-style FUNCTION
    GHSER blocks: a non-host phase's pure-element parameter is written as
    ``+GHSER<EL># + a + b·T`` instead of an independent 6-term polynomial.
    The linear offset is monotone in T so the host-vs-non-host pure G
    crossing (melting line) occurs exactly once — no spurious second
    crossing at high T from independent poly6 extrapolation.

    Parameters
    ----------
    df : DataFrame
        Phase-diagram dataframe with at least 'phase', 'composition',
        'temperature', 'free_energy' columns.
    phase : str
        Phase name to fit (e.g. 'lqd').
    host_surface : dict
        Output of build_calphad_surface for the SER host phase (provides
        coeffs_A, coeffs_B = GHSER_A, GHSER_B).
    element_side : str
        'A' for non-reference element (composition ≈ 0), 'B' for the
        reference element (composition ≈ 1).
    endpoint_tol : float
        Composition window for picking pure-endpoint rows.

    Returns
    -------
    (a, b) in eV/atom, or ``None`` if no pure-endpoint data exists for
    this phase.
    """
    df_phase = df.loc[df["phase"] == phase]
    if element_side == "A":
        rows = df_phase.loc[df_phase["composition"] <= endpoint_tol]
        target = 0.0
    else:
        rows = df_phase.loc[df_phase["composition"] >= 1.0 - endpoint_tol]
        target = 1.0
    if len(rows) == 0:
        return None

    # Prefer is_reference rows; fall back to closest-to-endpoint.
    if "is_reference" in rows.columns:
        ref = rows.loc[rows["is_reference"] == True]
        if len(ref) > 0:
            rows = ref
    pick = rows.iloc[(rows["composition"] - target).abs().argsort().iloc[0]]

    T = np.atleast_1d(np.asarray(pick["temperature"], dtype=float))
    G = np.atleast_1d(np.asarray(pick["free_energy"], dtype=float))
    if len(T) != len(G) or len(T) < 2:
        return None

    host_coeffs = host_surface["coeffs_A" if element_side == "A" else "coeffs_B"]
    G_host = _eval_calphad_poly6(host_coeffs, T)
    residual = G - G_host

    A = np.column_stack([np.ones_like(T), T])
    (a, b), _, _, _ = np.linalg.lstsq(A, residual, rcond=None)
    return float(a), float(b)


def _fit_line_compound_phase(points, x_reference=None):
    if len(points) < 2:
        raise ValueError("At least two finite temperature points are needed to fit a line compound")
    if x_reference is None:
        x_reference = float(np.median(points["composition"]))
    distances = np.abs(points["composition"].to_numpy(dtype=float) - x_reference)
    line_points = points.iloc[distances <= max(np.min(distances), 1e-12)]
    coeffs = _fit_calphad_poly6(
        line_points["temperature"].to_numpy(dtype=float),
        line_points["free_energy"].to_numpy(dtype=float),
    )
    predicted = _eval_calphad_poly6(coeffs, line_points["temperature"].to_numpy(dtype=float))
    rms = float(
        np.sqrt(
            np.mean((predicted - line_points["free_energy"].to_numpy(dtype=float)) ** 2)
        )
        * EV_TO_J_MOL
    )
    return {
        "coeffs": coeffs,
        "x_reference": float(x_reference),
        "rms_j_mol": rms,
        "n_fit_points": int(len(line_points)),
    }


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
                return rows.iloc[(rows["composition"] - target).abs().argsort().iloc[0]]

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
            # Collect one row per (composition, T) point so we can fit
            # temperature-dependent RK parameters: L_k(T) = a_k + b_k·T
            x_data = []    # unique composition values (for the count check)
            xs_all = []    # composition for each data point
            Ts_all = []    # temperature for each data point
            Gxs_all = []   # G_xs for each data point
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
                G_ideal = kb * T_row * (x * np.log(x) + (1.0 - x) * np.log(1.0 - x))

                # Excess = total − linear reference − ideal entropy.
                G_xs_arr = G_row - G_lin - G_ideal

                x_data.append(x)
                xs_all.extend([x] * len(T_row))
                Ts_all.extend(T_row.tolist())
                Gxs_all.extend(G_xs_arr.tolist())

            if len(x_data) < rk_order:
                warnings.warn(
                    f"Phase '{phase}': only {len(x_data)} intermediate "
                    f"compositions available, need ≥ {rk_order} for the RK "
                    "fit. Skipping CALPHAD surface for this phase."
                )
                surfaces[phase] = None
                continue

            x_fit = np.array(xs_all)
            T_fit = np.array(Ts_all)
            Gxs_vec = np.array(Gxs_all)
            pf_fit = x_fit * (1.0 - x_fit)

            # --- Fit temperature-dependent Redlich-Kister parameters ---
            # G_xs(x,T) = x(1-x) · Σ_k L_k(T) · (1-2x)^k
            # L_k(T) uses the same 6-term CALPHAD polynomial as G_A/G_B:
            #   L_k(T) = a + b·T + c·T·lnT + d·T² + e·T³ + f·T⁻¹
            # Build 6*rk_order column design matrix.
            cols = []
            log_T_fit = np.log(T_fit)
            for k in range(rk_order):
                pk = pf_fit * (1.0 - 2.0 * x_fit) ** k
                cols.append(pk)
                cols.append(pk * T_fit)
                cols.append(pk * T_fit * log_T_fit)
                cols.append(pk * T_fit ** 2)
                cols.append(pk * T_fit ** 3)
                cols.append(pk / T_fit)
            basis = np.column_stack(cols)
            L_flat, _, _, _ = np.linalg.lstsq(basis, Gxs_vec, rcond=None)
            # Shape (rk_order, 6): each row is [a, b, c, d, e, f] for L_k(T)
            L_coeffs = L_flat.reshape(rk_order, 6)

            Gxs_pred = basis @ L_flat
            rms = float(np.sqrt(np.mean((Gxs_vec - Gxs_pred) ** 2)))
            T_ref = float(np.mean(T_fit))
            L_str = "  ".join(
                f"L{k}(T_ref={T_ref:.0f})={_eval_calphad_poly6(L_coeffs[k], T_ref):.4f} eV/atom"
                for k in range(rk_order)
            )
            print(f"[{phase}]  {L_str}   RMS_xs={rms:.4f} eV/atom")

            surfaces[phase] = {
                "coeffs_A": coeffs_A,
                "coeffs_B": coeffs_B,
                "L_coeffs": L_coeffs,
                "x_data": np.array(x_data),
                "G_xs_data": np.array([np.mean(Gxs_vec[np.array(xs_all) == x]) for x in x_data]),
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
        calphad_surface=True,
        rk_order=3,
        composition_grid=10000,
        peak_cutoff=0.003,
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
            *rk_order*.  Default True.  Set to False to revert to the
            legacy per-temperature polynomial mode.
        rk_order : int
            Number of Redlich-Kister parameters for the CALPHAD surface.
            Only used when *calphad_surface* is True and the surface has not
            been pre-built.  Default 3.
        composition_grid : int
            Number of composition points on the evaluation grid when
            *calphad_surface* is True.  Default 10000.
        peak_cutoff : float
            Minimum composition-gap width required to report a two-phase
            coexistence region from the convex-hull construction.  Reduce
            this to catch narrow coexistence regions near pure endpoints
            (e.g. the fcc-liquid window close to the melting point of a
            pure component).  Default 0.003.
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
                    dc, peak_cutoff=peak_cutoff, remove_self_tangents_for=remove_self_tangents_for
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

        n_phases = sum(1 for v in self._calphad_surfaces.values() if v is not None)
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
                    (df_p["composition"] >= comp_lo) & (df_p["composition"] <= comp_hi)
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
            # Evaluate L_k at a representative temperature for the composition plot
            L = np.asarray(L)
            T_ref_plot = float(np.mean(surf["x_data"]) * 0 + 1000.0)  # use 1000 K as reference
            if L.ndim == 2 and L.shape[1] == 6:
                Lk_at_T = [_eval_calphad_poly6(L[k], T_ref_plot) for k in range(rk_order)]
            elif L.ndim == 2:
                Lk_at_T = [L[k, 0] + L[k, 1] * T_ref_plot for k in range(rk_order)]
            else:
                Lk_at_T = [L[k] for k in range(rk_order)]
            G_xs_fit = sum(
                Lk_at_T[k] * pf * (1.0 - 2.0 * x_fine) ** k for k in range(rk_order)
            )
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

    def to_pickle(self, filename):
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
    def from_pickle(filename):
        """
        Load a PhaseDiagram object previously saved with :meth:`to_pickle`.

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

    def to_df(self, filename):
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

    def to_parquet(self, filename):
        """
        Save the processed DataFrame to a Parquet file (version-independent).

        Columns that hold numpy arrays (``temperature``, ``free_energy``) are
        flattened to one row per temperature point; a ``row_id`` column tracks
        which rows belong together.  Object-level metadata
        (``reference_element``, ``phases``, ``composition_intervals``) is
        embedded in the Parquet schema so the full object can be reconstructed
        by :meth:`from_parquet`.

        Parameters
        ----------
        filename : str
            Output path (e.g. ``'phase_diagram.parquet'``).
        """
        import json
        import pyarrow as pa
        import pyarrow.parquet as pq

        scalar_cols = [
            c for c in self.df.columns if c not in ("temperature", "free_energy")
        ]

        rows = []
        for row_id, (_, row) in enumerate(self.df.iterrows()):
            temps = np.asarray(row["temperature"], dtype=float)
            fes = np.asarray(row["free_energy"], dtype=float)
            scalars = {c: row[c] for c in scalar_cols}
            for t, fe in zip(temps, fes):
                rows.append({"row_id": row_id, "temperature": t, "free_energy": fe, **scalars})

        flat_df = pd.DataFrame(rows)
        table = pa.Table.from_pandas(flat_df, preserve_index=False)

        meta = {
            b"calphy_reference_element": self.reference_element.encode(),
            b"calphy_phases": json.dumps(self.phases).encode(),
            b"calphy_composition_intervals": json.dumps(
                {k: list(v) for k, v in self.composition_intervals.items()}
            ).encode(),
        }
        existing = table.schema.metadata or {}
        table = table.replace_schema_metadata({**existing, **meta})
        pq.write_table(table, filename)

    def to_tdb(
        self,
        filename,
        elements=None,
        line_compounds=None,
        rk_order=3,
        T_min=298.15,
        T_max=None,
        L_temperature_form="poly6",
        limited_rk_order=None,
        limited_L_n_terms=3,
        limited_fit_order=4,
        full_range_ridge_lambda=100.0,
        compound_anti_site_penalty=80000.0,
        compound_pure_sublattice_penalty=80000.0,
    ):
        """
        Write a TDB file representing this phase diagram in SGTE/CALPHAD
        convention so pycalphad (or Thermo-Calc / OpenCalphad) can read it.

        TDB structure
        -------------
        - Header: ``TEMP_LIM``, ``DEFINE_SYSTEM_DEFAULT``,
          ``DEFAULT_COMMAND DEFINE_SYS_ELEMENT VA``, ``TYPE_DEFINITION``.
        - ``FUNCTION GHSER<EL>`` blocks for each element, taken from the
          *host* phase's pure-element 6-term polynomial.
        - The host phase (full-range solid solution containing both pure
          endpoints; e.g. ``fcc``) and limited-range phases are written as
          2-sublattice ``(A,B):VA`` solid solutions referencing GHSER
          (the universal SGTE convention seen in COST507/solders.tdb).
        - Non-host full-range phases (e.g. ``lqd``) keep an independent
          six-term polynomial for their pure-element G — the linear
          ``GHSER + a + b·T`` offset is too lossy for noisy MD liquid data
          (drops melting points by hundreds of K).
        - Phases listed in ``line_compounds`` are 2-sublattice ``A:B``
          with a single 6-term polynomial G(T) per formula unit.
        - Full-range phases' Redlich-Kister L parameters are written with
          the temperature dependence chosen by ``L_temperature_form``:
          the exact six-term polynomials from
          :meth:`build_calphad_surface` (``"poly6"``, default — the TDB
          then reproduces calphy's G(x,T) surface *exactly* inside the
          validity range) or a ridge-regularised ``a + b·T`` refit
          (``"linear"``, the standard SGTE convention).  Limited-range
          phases always use ``a + b·T`` from the bounded fit.

        A JSON metadata comment (``$ CALPHY_TDB_METADATA``) is embedded
        in the header so :meth:`from_tdb` can recover phase classification
        and stoichiometry without trying to infer it from the TDB content.

        Parameters
        ----------
        filename : str
            Output TDB path.
        elements : tuple of (A, B), optional
            Binary element pair.  ``B`` should be ``reference_element``.
            If omitted, the pair is inferred from phase names.
        line_compounds : list of str, optional
            Phase names to treat as line compounds (single stoichiometry,
            single G(T) polynomial).  Default: ``[]``.
        rk_order : int
            Number of Redlich-Kister terms (default 3).
        T_min, T_max : float
            Temperature validity range written in PARAMETER blocks.
            ``T_max`` defaults to the max temperature found in ``self.df``
            — pycalphad extrapolates blindly beyond this, so leaving it
            tight to the actual MD range is the conservative choice.
        L_temperature_form : str
            Temperature form for full-range solution phases' L_k(T)
            parameters.

            ``"poly6"`` (default)
                Write the six-term CALPHAD polynomials fitted by
                :meth:`build_calphad_surface` verbatim.  The TDB free
                energy of every full-range phase is then *identical* to
                calphy's own CALPHAD surface, so a pycalphad phase
                diagram computed from the TDB matches
                ``calculate(calphad_surface=True)`` inside the validity
                range.  Do not evaluate the TDB outside
                ``[T_min, T_max]`` — the polynomial tails are unbounded.
            ``"linear"``
                Legacy SGTE-style refit of each L_k as ``a + b·T`` with
                Tikhonov ridge regularisation
                (``full_range_ridge_lambda``).  Smoother extrapolation
                but drifts up to tens of meV/atom from the calphy
                surface, visibly narrowing e.g. the solid–liquid lens.
        limited_rk_order : int, optional
            RK order for limited-range phases' bounded fit.  By default
            the order is auto-selected per phase: the lowest of
            (3, 4, 5, 6) that fits the phase's pseudo-data to better
            than ~1 meV/atom RMS (wider composition windows need more
            terms — an underfit well overstabilises the phase and
            pushes its dissolution temperature up by ~100 K, an overfit
            one adds interior wiggle).  Pass an integer to force a
            fixed order instead.
        limited_L_n_terms : int
            T-basis terms per L_k in the bounded fit: ``2`` → a + b·T,
            ``3`` → a + b·T + c·T·lnT (default).  Three terms follow
            the curvature of the stability window in T.
        limited_fit_order : int
            Polynomial order of the per-temperature smoothing fit used
            to generate dense pseudo-data for limited-range phases
            (default 4 — the same default :meth:`calculate` uses for
            narrow phases).  The RK model is fit against these smooth
            per-T curves rather than the raw composition points, which
            keeps the well-minimum location — and hence the shape of
            the single-phase window near the dome top — faithful to
            calphy's own construction.
        full_range_ridge_lambda : float
            Tikhonov ridge weight used in refitting full-range solution
            phases' L_k(T) as ``a + b·T`` when
            ``L_temperature_form="linear"``.  Plain LSQ on noisy liquid
            data gives unstable huge cancelling magnitudes that produce
            wavy liquidi and spurious miscibility gaps on extrapolation;
            the ridge bounds the magnitudes.  Default ``100`` — increase
            for smoother (and less data-faithful) L; decrease for more
            data-faithful (and possibly less stable) L.
        compound_anti_site_penalty, compound_pure_sublattice_penalty : float
            Reserved for the 2-sublattice CEF (A,B):(A,B) compound code
            path (``_fit_compound_two_sublattice``).  Not exercised under
            the current default settings — limited-range phases use the
            1-sublattice bounded fit instead.

        Limitations and caveats
        -----------------------
        The TDB is a *lossy projection* of the calphy MD data:

        1. **Full-range phases are exact only in poly6 mode.** With the
           default ``L_temperature_form="poly6"`` the full-range solution
           phases' G(x,T) in the TDB equals calphy's CALPHAD surface
           exactly (same 6-term pure-element G and L_k polynomials).
           With ``"linear"``, L_k are refit as ``a + b·T`` with Tikhonov
           ridge regularisation (``full_range_ridge_lambda``), trading
           data-fit accuracy for smoother T extrapolation: the per-point
           ``G(x,T)`` drift is a few meV/atom for well-sampled phases
           (FCC) and up to ~20-40 meV/atom for noisy ones (LQD), which
           visibly distorts the solid-liquid lens.

        2. **Limited-range phases use a bounded SLSQP fit**, not the
           original ``build_calphad_surface`` (which has no surface for
           narrow phases anyway).  The bounded fit shares the host's
           pure-element G_A/G_B and fits L_k(T) only to data inside the
           composition window, with the constraint
           ``G_phase(x,T) ≥ G_host(x,T)`` enforced at sample points
           *outside* the window to prevent the polynomial from making the
           phase spuriously stable far from its data.  This usually
           produces narrower stability windows in the pycalphad diagram
           than calphy's own per-T polynomial common-tangent gives.

        3. **Plot only inside the data T range.**  ``T_max`` defaults to
           the MD data's actual maximum; extrapolating beyond — even by
           100 K — frequently exposes wild 6-term polynomial tails
           (visible as spurious LQD/FCC lenses).  Set the plot range to
           ``T_max`` or below.

        4. **The raw F(x,T) data is not stored.**  Use :meth:`to_parquet`
           if you need a lossless round-trip.

        5. **Pycalphad's TDB parser rejects inline ``REF:`` tags** in
           PARAMETER and FUNCTION blocks (despite Thermo-Calc accepting
           them in COST507 etc.), so traceability is emitted as
           ``$ REF: calphy-MD-<date>`` comment lines instead.
        """
        line_compounds = list(line_compounds or [])
        for ph in line_compounds:
            if ph not in self.phases:
                raise ValueError(
                    f"line_compound '{ph}' is not in self.phases {self.phases}"
                )
        if L_temperature_form not in ("poly6", "linear"):
            raise ValueError(
                "L_temperature_form must be 'poly6' or 'linear', got "
                f"{L_temperature_form!r}"
            )

        # ---- Element handling ----
        if elements is None:
            el_a, el_b = _infer_binary_elements(self.phases, self.reference_element)
        else:
            el_a, el_b = [e.upper() for e in elements]
        ref_upper = self.reference_element.upper()
        if el_b != ref_upper:
            el_a, el_b = el_b, el_a
        if el_b != ref_upper:
            raise ValueError(
                f"Neither element in {elements} matches reference_element "
                f"'{self.reference_element}'"
            )

        # ---- Temperature validity range ----
        if T_max is None:
            t_max = 0.0
            for _, row in self.df.iterrows():
                T = np.atleast_1d(np.asarray(row["temperature"], dtype=float))
                if len(T):
                    t_max = max(t_max, float(np.max(T)))
            T_max = t_max if t_max > 0 else 2000.0
        T_min = float(T_min)
        T_max = float(T_max)

        # ---- Classify phases ----
        def _is_full_range(ph):
            interval = self.composition_intervals.get(ph, (0.0, 1.0))
            return interval[0] <= 1e-6 and interval[1] >= 1.0 - 1e-6

        full_range = [
            ph for ph in self.phases
            if ph not in line_compounds and _is_full_range(ph)
        ]
        limited_range = [
            ph for ph in self.phases
            if ph not in line_compounds and not _is_full_range(ph)
        ]

        if not full_range and limited_range:
            raise ValueError(
                "TDB output requires at least one full-range solution phase "
                "to anchor pure-element references for limited-range phases."
            )

        # ---- Build CALPHAD surfaces for full-range phases ----
        # Always rebuild here so the rk_order kwarg actually takes effect
        # (otherwise a pre-built surface with a different rk_order would be
        # reused).  Higher rk_order = more x-curvature flexibility = more
        # noise-absorbing wiggle in the liquidus; lower rk_order = smoother.
        self.build_calphad_surface(rk_order=rk_order)

        host_surface = None
        host_phase = None
        for ph in full_range:
            s = self._calphad_surfaces.get(ph)
            if s is not None:
                host_surface = s
                host_phase = ph
                break
        if host_surface is None and limited_range:
            raise ValueError(
                "No full-range phase produced a valid calphad surface; "
                "cannot anchor limited-range phases."
            )

        # ---- Fit each phase ----
        # Every non-line-compound phase is written as a 1-sublattice (A,B)
        # substitutional solution with RK excess — same model as fcc/lqd.
        # Limited-range phases share the host (fcc) pure-element G_A(T),
        # G_B(T) so the reference state is consistent, and the L_k(T)
        # parameters are fit only to data inside the phase's composition
        # window (no anchor points pulling them toward zero, which would
        # weaken the phase inside its window).  The RK form's x(1-x)
        # factor naturally drives G_xs to 0 at the pure endpoints, so the
        # phase smoothly coincides with the host at the boundaries.
        records = []
        for ph in self.phases:
            if ph in line_compounds:
                points = _phase_data_as_points(self.df, ph)
                fit = _fit_line_compound_phase(points)
                records.append({"name": ph, "kind": "line_compound", "fit": fit})
            elif ph in full_range:
                s = self._calphad_surfaces.get(ph)
                if s is None:
                    warnings.warn(
                        f"Phase '{ph}': calphad surface unavailable; skipping."
                    )
                    continue
                # Both host and non-host full-range phases keep
                # independent 6-term G for their pure endpoints —
                # forcing a non-host phase's pure G to be GHSER +
                # linear-offset loses too much curvature for noisy MD
                # data (e.g. pure-Au LQD melting drifted by >200 K).
                # GHSER FUNCTION blocks are still emitted and used by
                # the host phase and by limited-range phases (which by
                # construction share the host's pure G).
                if L_temperature_form == "poly6":
                    # Write the build_calphad_surface fit verbatim: the
                    # TDB G(x,T) is then identical to calphy's surface.
                    records.append({"name": ph, "kind": "solution", "fit": s})
                else:
                    # Legacy: refit L_k as a + b·T with ridge.
                    phase_data = self.df.loc[self.df["phase"] == ph]
                    linear_fit = _fit_solution_linear_L_regularized(
                        phase_data, s, rk_order=rk_order,
                        lam=full_range_ridge_lambda,
                        fit_pure_element_offsets=False,
                    )
                    records.append(
                        {"name": ph, "kind": "solution", "fit": linear_fit}
                    )
            else:
                # 1-sublattice (A,B) RK with hard inequality
                # G_phase(x, T) >= G_host(x, T) enforced outside the data
                # composition window (SLSQP-bounded fit).  This gives a
                # visible solubility dome inside the window while
                # preventing the phase from being spuriously stable far
                # from its real composition range.  The 2-sublattice CEF
                # alternative with anti-site penalties — the SGTE standard
                # — collapses to a near-line-compound appearance with
                # default penalties, losing the solubility window that
                # calphy data captures.
                # Prefer dense pseudo-data from calphy's own per-T
                # polynomial fits (same smoothing calculate() uses for
                # narrow phases): it pins the RK well's interior shape,
                # keeping the well-minimum location faithful to the data
                # near the order-disorder dome top.  Fall back to the
                # raw points if the per-T fit is not possible (too few
                # compositions).
                points = _limited_phase_pseudo_points(
                    self.df,
                    ph,
                    self.composition_intervals.get(ph, (0.0, 1.0)),
                    fit_order=limited_fit_order,
                )
                if points is None:
                    points = _phase_data_as_points(self.df, ph)
                n_x = int(points["composition"].nunique())
                if limited_rk_order is not None:
                    rk_candidates = [int(limited_rk_order)]
                else:
                    # Auto-select the RK order: wider windows need more
                    # terms (au2cu spanning Δx=0.2 needs order 5 where
                    # aucu at Δx=0.1 is fine with 4).
                    rk_candidates = [
                        rk for rk in (3, 4, 5, 6) if rk <= max(1, n_x - 1)
                    ] or [max(1, n_x - 1)]

                # The phase-diagram-relevant scalar is the dissolution
                # temperature — where the well's minimum depth below the
                # host crosses zero.  Pick the RK order that reproduces
                # the pseudo-data's crossing best (an underfit well is
                # too deep at mid-T and pushes the dissolution up by
                # tens of K); fall back to RMS if the phase never
                # dissolves inside the data range.
                T_d, depth_d = _limited_phase_min_depth_by_T(
                    points, host_surface
                )
                T_cross_data = _first_zero_crossing(T_d, depth_d)
                fit = None
                best_key = None
                for rk_ph in rk_candidates:
                    trial = _fit_limited_range_surface_bounded(
                        points,
                        host_surface,
                        rk_order=rk_ph,
                        n_T_terms=limited_L_n_terms,
                        constraint_T_range=(T_min, T_max),
                    )
                    if T_cross_data is not None:
                        T_f, depth_f = _limited_phase_min_depth_by_T(
                            points, host_surface, fit=trial
                        )
                        T_cross_fit = _first_zero_crossing(T_f, depth_f)
                        dT = (
                            abs(T_cross_fit - T_cross_data)
                            if T_cross_fit is not None
                            else float("inf")
                        )
                        key = (dT, trial["rms_j_mol"])
                    else:
                        key = (0.0, trial["rms_j_mol"])
                    if best_key is None or key < best_key:
                        best_key = key
                        fit = trial
                if fit.get("constraint_violation_max", 0.0) > 1e-4:
                    warnings.warn(
                        f"Phase '{ph}': bounded fit could not fully satisfy "
                        f"G_phase >= G_host outside the window "
                        f"(max violation {fit['constraint_violation_max']:.4f} "
                        "eV/atom). Phase may still appear in unintended "
                        "regions of the pycalphad diagram."
                    )
                records.append(
                    {"name": ph, "kind": "limited_range", "fit": fit}
                )

        # ---- Compose TDB text ----
        metadata = {
            "format": "calphy.phase_diagram.to_tdb.v1",
            "reference_element": self.reference_element,
            "elements": [el_a, el_b],
            "phases": list(self.phases),
            "phase_kinds": {r["name"]: r["kind"] for r in records},
            "composition_intervals": {
                k: list(v) for k, v in self.composition_intervals.items()
            },
            "line_compounds": line_compounds,
            "host_phase": host_phase,
            "rk_order": int(rk_order),
            "L_temperature_form": L_temperature_form,
            "limited_rk_order": {
                r["name"]: int(r["fit"]["rk_order"])
                for r in records
                if r["kind"] == "limited_range"
            },
            "limited_L_n_terms": int(limited_L_n_terms),
            "T_min": T_min,
            "T_max": T_max,
            "units": "J/mol-atoms",
        }
        # Embed compound stoichiometries in metadata so from_tdb knows
        # the (m, n) used in PHASE declarations.
        compound_stoich = {}
        for rec in records:
            if rec["kind"] == "line_compound":
                m, n = _stoichiometry_from_reference_composition(
                    rec["fit"]["x_reference"]
                )
                compound_stoich[rec["name"]] = [int(m), int(n)]
            elif rec["kind"] == "compound_2sl":
                m, n = rec["fit"]["site_ratios"]
                compound_stoich[rec["name"]] = [int(m), int(n)]
        metadata["compound_stoichiometry"] = compound_stoich
        # Keep legacy alias so older from_tdb still loads line compounds.
        metadata["line_compound_stoichiometry"] = {
            k: v for k, v in compound_stoich.items() if k in line_compounds
        }

        # Reference comment, emitted as a `$` comment line above each
        # FUNCTION / PARAMETER block.  pycalphad's TDB parser rejects
        # inline `REF:` tags in PARAMETER definitions (even though Thermo-
        # Calc / SGTE-style databases such as COST507 carry them), so we
        # keep traceability via comments instead.
        ref_date = datetime.now().strftime('%Y%m%d')
        ref_comment = f"$ REF: calphy-MD-{ref_date}"

        lines = []
        lines.append("$")
        lines.append(f"$ TDB generated from calphy PhaseDiagram: {el_a}-{el_b}")
        lines.append(f"$ Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
        lines.append(_TDB_METADATA_PREFIX + json.dumps(metadata, sort_keys=True))
        lines.append("$")
        lines.append("")
        # SGTE-style header (TEMP_LIM + system defaults + type def).
        lines.append(f"TEMP_LIM {T_min:.2f} {T_max:.2f} !")
        lines.append("")
        lines.append("DEFINE_SYSTEM_DEFAULT ELEMENT 2 !")
        lines.append("DEFAULT_COMMAND DEFINE_SYS_ELEMENT VA !")
        lines.append("")
        lines.append("TYPE_DEFINITION % SEQ *!")
        lines.append("")
        lines.append("ELEMENT /-   ELECTRON_GAS     0.0 0.0 0.0 !")
        lines.append("ELEMENT VA   VACUUM           0.0 0.0 0.0 !")
        ser_a = _SER_REF.get(el_a, "FCC_A1")
        ser_b = _SER_REF.get(el_b, "FCC_A1")
        lines.append(f"ELEMENT {el_a:<4} {ser_a:<16} 0.0 0.0 0.0 !")
        lines.append(f"ELEMENT {el_b:<4} {ser_b:<16} 0.0 0.0 0.0 !")
        lines.append("")

        # FUNCTION GHSER<EL> blocks.  Standard SGTE convention: each
        # pure-element 6-term polynomial lives in one FUNCTION; all phases
        # reference it via +GHSER<EL>#.  GHSER for each element comes from
        # the host (SER) phase's fit.
        host_fit = next(
            (r["fit"] for r in records if r["name"] == host_phase), None
        )
        if host_fit is None:
            raise RuntimeError("Internal error: host phase has no fit record.")
        ghser_a = np.asarray(host_fit["coeffs_A"], dtype=float) * EV_TO_J_MOL
        ghser_b = np.asarray(host_fit["coeffs_B"], dtype=float) * EV_TO_J_MOL
        ghser_name_a = f"GHSER{el_a}"
        ghser_name_b = f"GHSER{el_b}"
        lines.append(ref_comment)
        lines.append(
            f"FUNCTION {ghser_name_a} {T_min:.2f} "
            f"{_format_tdb_expr_poly6(ghser_a)}; {T_max:.2f} N !"
        )
        lines.append(
            f"FUNCTION {ghser_name_b} {T_min:.2f} "
            f"{_format_tdb_expr_poly6(ghser_b)}; {T_max:.2f} N !"
        )
        lines.append("")

        def _format_ghser_plus_linear(ghser_name, a_b):
            """Format `+GHSER<EL># + a + b·T` with a/b in J/mol-atoms."""
            a, b = a_b
            lin = _format_tdb_expr_linear(a, b)
            if lin == "0":
                return f"+{ghser_name}#"
            # lin starts with +/- so concat is fine
            return f"+{ghser_name}#{lin}"

        # Per-phase species suffix used in PARAMETER blocks:
        #   solid solution     → "EL:VA"  (e.g. AU:VA, A,B:VA for L)
        #   liquid             → "EL"
        #   line compound      → "A:B"   (fixed stoichiometry)
        #   CEF compound       → "A:B" etc., handled per-endmember below
        def _is_solid_solution(rec):
            return (
                rec["kind"] == "solution"
                and not _phase_is_liquid(rec["name"])
            )

        # Phase declarations
        for rec in records:
            tdb_name = _tdb_phase_name(rec["name"])
            if rec["kind"] == "solution":
                if _is_solid_solution(rec):
                    # Solid solution: standard (M):(VA) 2-sublattice form
                    # (the universal SGTE/CALPHAD convention for FCC/BCC/HCP
                    # solid solutions, even when no interstitial chemistry
                    # is being modelled — the VA second sublattice is just
                    # full).
                    lines.append(f"PHASE {tdb_name} % 2 1.0 1.0 !")
                    lines.append(
                        f"CONSTITUENT {tdb_name} : {el_a},{el_b} : VA : !"
                    )
                else:
                    # Liquid: 1-sublattice as in COST507/solders.tdb.
                    lines.append(f"PHASE {tdb_name} % 1 1.0 !")
                    lines.append(f"CONSTITUENT {tdb_name} : {el_a},{el_b} : !")
            elif rec["kind"] == "limited_range":
                # Limited-range ordered/compound phases are written as
                # solid-solution-style 1-sublattice (A,B):(VA) — same
                # convention as the host phase, just with L_k fit only
                # inside the data window (with the bounded SLSQP fit).
                lines.append(f"PHASE {tdb_name} % 2 1.0 1.0 !")
                lines.append(
                    f"CONSTITUENT {tdb_name} : {el_a},{el_b} : VA : !"
                )
            elif rec["kind"] == "compound_2sl":
                # Legacy code path kept for backwards-compatibility; no
                # currently-supported `_fit_compound_two_sublattice` flow
                # produces records of this kind under the default to_tdb
                # configuration.
                m, n = compound_stoich[rec["name"]]
                lines.append(f"PHASE {tdb_name} % 2 {m} {n} !")
                lines.append(
                    f"CONSTITUENT {tdb_name} : {el_a},{el_b} : {el_a},{el_b} : !"
                )
            else:  # line_compound
                m, n = compound_stoich[rec["name"]]
                lines.append(f"PHASE {tdb_name} % 2 {m} {n} !")
                lines.append(f"CONSTITUENT {tdb_name} : {el_a} : {el_b} : !")
        lines.append("")

        # Parameters
        lines.append(ref_comment)
        for rec in records:
            tdb_name = _tdb_phase_name(rec["name"])
            fit = rec["fit"]
            if rec["kind"] == "solution":
                # Choose the species suffix (with or without `:VA`) based on
                # whether the phase is a solid solution; the L-parameter
                # interaction string is similarly suffixed.
                if _is_solid_solution(rec):
                    spec_a = f"{el_a}:VA"
                    spec_b = f"{el_b}:VA"
                    interaction = f"{el_a},{el_b}:VA"
                else:
                    spec_a = el_a
                    spec_b = el_b
                    interaction = f"{el_a},{el_b}"

                if rec["name"] == host_phase:
                    lines.append(
                        f"PARAMETER G({tdb_name},{spec_a};0) {T_min:.2f} "
                        f"+{ghser_name_a}#; {T_max:.2f} N !"
                    )
                    lines.append(
                        f"PARAMETER G({tdb_name},{spec_b};0) {T_min:.2f} "
                        f"+{ghser_name_b}#; {T_max:.2f} N !"
                    )
                else:
                    cA = np.asarray(fit["coeffs_A"], dtype=float) * EV_TO_J_MOL
                    cB = np.asarray(fit["coeffs_B"], dtype=float) * EV_TO_J_MOL
                    lines.append(
                        f"PARAMETER G({tdb_name},{spec_a};0) {T_min:.2f} "
                        f"{_format_tdb_expr_poly6(cA)}; {T_max:.2f} N !"
                    )
                    lines.append(
                        f"PARAMETER G({tdb_name},{spec_b};0) {T_min:.2f} "
                        f"{_format_tdb_expr_poly6(cB)}; {T_max:.2f} N !"
                    )
                L = np.asarray(fit["L_coeffs"], dtype=float)
                for k in range(L.shape[0]):
                    if L.ndim == 2 and L.shape[1] == 6:
                        coeffs = L[k] * EV_TO_J_MOL
                        expr = _format_tdb_expr_poly6(coeffs)
                    else:
                        a, b = L[k, 0] * EV_TO_J_MOL, L[k, 1] * EV_TO_J_MOL
                        expr = _format_tdb_expr_linear(a, b)
                    lines.append(
                        f"PARAMETER L({tdb_name},{interaction};{k}) "
                        f"{T_min:.2f} {expr}; {T_max:.2f} N !"
                    )
            elif rec["kind"] == "limited_range":
                # 1-sublattice (A,B):VA, sharing GHSER with the host
                # (the bounded fit set coeffs_A/B equal to the host's).
                # L_k(T) comes from the SLSQP-bounded fit — either
                # linear (a + b·T) or zero-padded poly6 with a T·lnT
                # term, depending on limited_L_n_terms.
                lines.append(
                    f"PARAMETER G({tdb_name},{el_a}:VA;0) {T_min:.2f} "
                    f"+{ghser_name_a}#; {T_max:.2f} N !"
                )
                lines.append(
                    f"PARAMETER G({tdb_name},{el_b}:VA;0) {T_min:.2f} "
                    f"+{ghser_name_b}#; {T_max:.2f} N !"
                )
                L = np.asarray(fit["L_coeffs"], dtype=float)
                for k in range(L.shape[0]):
                    if L.shape[1] == 6:
                        expr = _format_tdb_expr_poly6(L[k] * EV_TO_J_MOL)
                    else:
                        a = L[k, 0] * EV_TO_J_MOL
                        b = L[k, 1] * EV_TO_J_MOL
                        expr = _format_tdb_expr_linear(a, b)
                    lines.append(
                        f"PARAMETER L({tdb_name},{el_a},{el_b}:VA;{k}) "
                        f"{T_min:.2f} {expr}; {T_max:.2f} N !"
                    )
            elif rec["kind"] == "compound_2sl":
                # SGTE-style CEF endmembers, written as GHSER-referenced
                # formation energies.  See solders.tdb DO3 (Cu3Sn) lines
                # 356-368 for the canonical pattern.
                #
                # Stoichiometric A:B endmember (m A's on sublattice 1,
                # n B's on sublattice 2):
                #   G(A:B) = m·GHSER_A# + n·GHSER_B# + (a + b·T)
                # where (a + b·T) is the formation energy per formula
                # (data fit minus the GHSER reference, refit linearly).
                #
                # Anti-site B:A (n A's + m B's per formula, swapped sub-
                # lattice occupation): same formation energy plus
                # anti-site penalty, referenced to the swapped pure mix.
                #
                # Pure A:A and B:B (host-on-both-sublattices) endmembers
                # carry (m+n)·GHSER plus a large pure-sublattice penalty
                # so the CEF can never lower G by disordering away from
                # the stoichiometric column.
                m, n = compound_stoich[rec["name"]]
                N = m + n
                # Convert the stoichiometric AB data fit (in J/mol-formula)
                # into a GHSER-referenced linear formation energy: fit
                # G_data_per_formula(T) - (m·GHSER_A(T) + n·GHSER_B(T))
                # ≈ a + b·T at a handful of sample temperatures.
                a_AB_total, b_AB_total = fit["theta_AB"]  # per formula
                T_samples = np.linspace(T_min, T_max, 64)
                G_data_per_formula = a_AB_total + b_AB_total * T_samples
                G_ref_AB = (
                    m * _eval_calphad_poly6(host_fit["coeffs_A"], T_samples)
                    + n * _eval_calphad_poly6(host_fit["coeffs_B"], T_samples)
                ) * EV_TO_J_MOL
                form_AB = G_data_per_formula - G_ref_AB
                design = np.column_stack([np.ones_like(T_samples), T_samples])
                (a_form, b_form), *_ = np.linalg.lstsq(
                    design, form_AB, rcond=None
                )

                pure_penalty = float(fit["pure_sublattice_penalty"])
                antisite_penalty = float(fit["anti_site_penalty"])

                def _coef_ghser(count, ghser_name):
                    # Returns "+m*GHSER#" / "+GHSER#" depending on count.
                    if count == 1:
                        return f"+{ghser_name}#"
                    return f"+{count}*{ghser_name}#"

                # A:A — pure A on both sublattices.
                lines.append(
                    f"PARAMETER G({tdb_name},{el_a}:{el_a};0) {T_min:.2f} "
                    f"{_coef_ghser(N, ghser_name_a)}"
                    f"{_format_tdb_expr_linear(pure_penalty, 0.0)}; "
                    f"{T_max:.2f} N !"
                )
                # A:B — stoichiometric column.
                lines.append(
                    f"PARAMETER G({tdb_name},{el_a}:{el_b};0) {T_min:.2f} "
                    f"{_coef_ghser(m, ghser_name_a)}"
                    f"{_coef_ghser(n, ghser_name_b)}"
                    f"{_format_tdb_expr_linear(float(a_form), float(b_form))}; "
                    f"{T_max:.2f} N !"
                )
                # B:A — anti-site column.
                lines.append(
                    f"PARAMETER G({tdb_name},{el_b}:{el_a};0) {T_min:.2f} "
                    f"{_coef_ghser(n, ghser_name_a)}"
                    f"{_coef_ghser(m, ghser_name_b)}"
                    f"{_format_tdb_expr_linear(float(a_form) + antisite_penalty, float(b_form))}; "
                    f"{T_max:.2f} N !"
                )
                # B:B — pure B on both sublattices.
                lines.append(
                    f"PARAMETER G({tdb_name},{el_b}:{el_b};0) {T_min:.2f} "
                    f"{_coef_ghser(N, ghser_name_b)}"
                    f"{_format_tdb_expr_linear(pure_penalty, 0.0)}; "
                    f"{T_max:.2f} N !"
                )
            else:  # line_compound
                m, n = compound_stoich[rec["name"]]
                N = m + n
                coeffs = np.asarray(fit["coeffs"], dtype=float) * EV_TO_J_MOL * N
                lines.append(
                    f"PARAMETER G({tdb_name},{el_a}:{el_b};0) {T_min:.2f} "
                    f"{_format_tdb_expr_poly6(coeffs)}; {T_max:.2f} N !"
                )
            lines.append("")

        with open(filename, "w") as f:
            f.write("\n".join(lines) + "\n")

    @staticmethod
    def from_tdb(filename):
        """
        Load a :class:`PhaseDiagram` from a TDB file written by
        :meth:`to_tdb`.

        Parses the file using the embedded ``$ CALPHY_TDB_METADATA`` JSON
        header, the ``FUNCTION`` blocks (GHSER and any other inlined
        functions), and the ``PARAMETER`` expressions.  Solution and
        limited-range phases are stored in ``_calphad_surfaces``; line
        compounds are stored in ``_line_compound_fits`` (mirrors the
        output of :func:`_fit_line_compound_phase`).  Legacy CEF
        compounds are stored in ``_compound_2sl_fits`` for backwards
        compatibility.  Units are converted from J/mol-atoms back to
        eV/atom.

        The returned object can drive :meth:`calculate` (with
        ``calphad_surface=True``) and :meth:`plot` because those code
        paths only need ``_calphad_surfaces``, not raw data.

        Parameters
        ----------
        filename : str
            Path to a TDB file produced by :meth:`to_tdb`.

        Returns
        -------
        PhaseDiagram

        Limitations and caveats
        -----------------------
        1. **Requires the ``$ CALPHY_TDB_METADATA`` header.**  Arbitrary
           TDBs from external sources (COST507, solders.tdb, etc.) are
           *not* supported — phase classification (solution / line
           compound / limited range / CEF), stoichiometries and the host
           phase identity all come from that JSON header.  Raises
           ``ValueError`` if the header is missing.

        2. **The returned ``df`` is empty.**  The raw F(x,T) MD data is
           not stored in the TDB and cannot be recovered.  Methods that
           need the dataframe will fail or produce empty output:

           - :meth:`plot_data_vs_fit`
           - :meth:`plot_convergence`
           - :meth:`plot_free_energy` with ``show_data=True``
           - :meth:`plot_free_energy_mixing` with ``show_data=True``
           - The legacy non-surface path of :meth:`calculate`
             (``calphad_surface=False``)

           Use :meth:`from_parquet` if you need the raw data.

        3. **The round-trip is lossy.**  :meth:`to_tdb` applies a
           Tikhonov-regularised linear-T refit to full-range solution
           phases' L_k coefficients, so the surfaces recovered here are
           the *regularised* versions, not the original 6-term fits
           returned by :meth:`build_calphad_surface`.  Typical per-point
           ``G(x, T)`` drift inside the data T range is O(1 mJ/mol) for
           solid solutions and up to ~40 meV/atom for noisy liquids.

        4. **Limited-range phase diagrams differ from the original.**
           When the source ``PhaseDiagram`` had raw data,
           :meth:`calculate` for narrow phases fell back to per-T
           polynomial fits of the data.  Without raw data, the loaded
           object uses the bounded-fit surface only — producing somewhat
           narrower / differently-shaped solubility windows in the
           pycalphad or calphy-rendered diagram than the original.

        5. **Composition intervals are preserved verbatim** from the
           metadata header, even though the polynomial fit may not
           perfectly respect them outside the constrained sample points.
        """
        with open(filename) as f:
            text = f.read()

        metadata = None
        for line in text.splitlines():
            if line.startswith(_TDB_METADATA_PREFIX):
                metadata = json.loads(line[len(_TDB_METADATA_PREFIX):])
                break
        if metadata is None:
            raise ValueError(
                f"'{filename}' does not contain a $ CALPHY_TDB_METADATA "
                "header — only TDBs written by PhaseDiagram.to_tdb can be "
                "round-tripped."
            )

        reference_element = metadata["reference_element"]
        el_a, el_b = metadata["elements"]
        phases = list(metadata["phases"])
        phase_kinds = metadata.get("phase_kinds", {})
        composition_intervals = {
            k: tuple(v) for k, v in metadata["composition_intervals"].items()
        }
        line_compounds = list(metadata.get("line_compounds", []))
        rk_order = int(metadata.get("rk_order", 3))

        phase_to_tdb = {p: _tdb_phase_name(p) for p in phases}
        tdb_to_phase = {v: k for k, v in phase_to_tdb.items()}

        # Parse FUNCTION blocks first so +NAME# references in PARAMETER
        # expressions can be resolved by inlining the function's polynomial
        # into the caller's poly6 coefficient sum.  Tail (e.g. `REF: ...`)
        # is allowed after the `N` terminator.
        func_re = re.compile(
            r"FUNCTION\s+(\w+)\s+([\d.+eE\-]+)\s+(.*?);\s*"
            r"([\d.+eE\-]+)\s+N\b[^!]*!",
            re.DOTALL | re.IGNORECASE,
        )
        functions = {}  # name (upper) -> ndarray(6) poly6 coeffs in J/mol
        for fm in func_re.finditer(text):
            name = fm.group(1)
            expr = fm.group(3)
            functions[name.upper()] = _parse_tdb_poly6_expr(expr)

        def _resolve_function_refs(expr):
            """Inline +NAME# references; sum into a single poly6 coeff vector.

            Supports ``+FUNC#`` and ``-FUNC#`` (with optional whitespace).
            Any literal poly6 terms in the rest of the expression are
            parsed and added.  Returns ndarray(6) in J/mol.
            """
            text_clean = expr.replace(" ", "")
            total = np.zeros(6, dtype=float)
            ref_re = re.compile(r"([+-])([A-Za-z_]\w*)#")
            pos = 0
            stripped = ""
            for rm in ref_re.finditer(text_clean):
                sign = 1.0 if rm.group(1) == "+" else -1.0
                fname = rm.group(2).upper()
                if fname in functions:
                    total += sign * functions[fname]
                else:
                    # unknown function — leave literal in expr so
                    # _parse_tdb_poly6_expr can attempt to handle it
                    stripped += text_clean[pos:rm.end()]
                    pos = rm.end()
                    continue
                stripped += text_clean[pos:rm.start()]
                pos = rm.end()
            stripped += text_clean[pos:]
            total = total + _parse_tdb_poly6_expr(stripped)
            return total

        # PARAMETER regex tolerates an arbitrary tail before the final '!'
        # (e.g. `REF: calphy-MD-20260516`).
        param_re = re.compile(
            r"PARAMETER\s+([GL])\s*\(\s*([^)]+?)\s*\)\s+"
            r"([\d.+eE\-]+)\s+(.*?);\s*([\d.+eE\-]+)\s+N\b[^!]*!",
            re.DOTALL | re.IGNORECASE,
        )

        # Aggregate parameters per phase.
        # solutions[phase] = {'A_poly6': ndarray(6), 'B_poly6': ndarray(6),
        #                     'L': {k: ndarray(6)}}
        # line_compound_params[phase] = ndarray(6) (per formula unit)
        # compound_2sl_params[phase] = {'AA': c6, 'AB': c6, 'BA': c6, 'BB': c6}
        solutions = {}
        line_compound_params = {}
        compound_2sl_params = {}

        for m in param_re.finditer(text):
            ptype, args, _tmin, expr, _tmax = m.groups()
            ptype = ptype.upper()
            head, _, ord_str = args.partition(";")
            order = int(ord_str.strip()) if ord_str.strip() else 0
            phase_tdb, _, species_str = head.partition(",")
            phase_tdb = phase_tdb.strip()
            species_str = species_str.strip()
            ph = tdb_to_phase.get(phase_tdb)
            if ph is None:
                continue
            coeffs = _resolve_function_refs(expr)
            kind = phase_kinds.get(ph)
            if kind is None:
                kind = (
                    "line_compound" if ph in line_compounds
                    else "solution"
                )

            if kind == "line_compound":
                if ptype == "G" and ":" in species_str:
                    line_compound_params[ph] = coeffs
                continue

            if kind == "compound_2sl":
                # Legacy 2-sublattice CEF format — parse for backwards
                # compatibility with older TDBs written by to_tdb.
                if ptype == "G" and ":" in species_str:
                    s1, s2 = [s.strip().upper() for s in species_str.split(":")]
                    if s1 == el_a and s2 == el_a:
                        key = "AA"
                    elif s1 == el_a and s2 == el_b:
                        key = "AB"
                    elif s1 == el_b and s2 == el_a:
                        key = "BA"
                    elif s1 == el_b and s2 == el_b:
                        key = "BB"
                    else:
                        continue
                    compound_2sl_params.setdefault(ph, {})[key] = coeffs
                continue

            # solution and limited_range both store as 1-sublattice (A,B);
            # strip the optional ":VA" suffix that solid solutions carry.
            entry = solutions.setdefault(ph, {"L": {}})
            spec = species_str.strip()
            if spec.upper().endswith(":VA"):
                spec = spec[: -len(":VA")]
            if ptype == "G":
                target = spec.strip()
                if target.upper() == el_a:
                    entry["A_poly6"] = coeffs
                elif target.upper() == el_b:
                    entry["B_poly6"] = coeffs
            elif ptype == "L":
                entry["L"][order] = coeffs

        # ---- Reconstruct surface dicts (convert J/mol-atoms -> eV/atom) ----
        surfaces = {}
        for ph, entry in solutions.items():
            cA = entry.get("A_poly6")
            cB = entry.get("B_poly6")
            if cA is None or cB is None:
                surfaces[ph] = None
                continue
            L_dict = entry["L"]
            max_k = max(L_dict.keys()) if L_dict else -1
            L_rows = [L_dict.get(k, np.zeros(6)) for k in range(max_k + 1)]
            L_coeffs = (np.asarray(L_rows) / EV_TO_J_MOL) if L_rows else np.zeros((0, 6))
            surfaces[ph] = {
                "coeffs_A": cA / EV_TO_J_MOL,
                "coeffs_B": cB / EV_TO_J_MOL,
                "L_coeffs": L_coeffs,
                "rk_order": int(L_coeffs.shape[0]) if L_coeffs.size else rk_order,
            }

        compound_stoich_meta = metadata.get(
            "compound_stoichiometry",
            metadata.get("line_compound_stoichiometry", {}),
        )

        line_compound_fits = {}
        for ph, per_formula_coeffs in line_compound_params.items():
            stoich = compound_stoich_meta.get(ph)
            if stoich is not None:
                m_a, n_b = int(stoich[0]), int(stoich[1])
                N = m_a + n_b
                x_ref = n_b / N
            else:
                x_ref = float(np.median(composition_intervals.get(ph, (0.5, 0.5))))
                m_a, n_b = _stoichiometry_from_reference_composition(x_ref)
                N = m_a + n_b
            line_compound_fits[ph] = {
                "coeffs": per_formula_coeffs / (EV_TO_J_MOL * N),
                "x_reference": x_ref,
                "stoichiometry": (m_a, n_b),
            }

        compound_2sl_fits = {}
        for ph, params in compound_2sl_params.items():
            stoich = compound_stoich_meta.get(ph)
            if stoich is None:
                continue
            m_a, n_b = int(stoich[0]), int(stoich[1])
            compound_2sl_fits[ph] = {
                "site_ratios": (m_a, n_b),
                "stoichiometry": (m_a, n_b),
                "x_stoich": m_a / (m_a + n_b),
                "endmembers": {k: v for k, v in params.items()},
            }

        obj = object.__new__(PhaseDiagram)
        obj.reference_element = reference_element
        obj.phases = phases
        obj.composition_intervals = composition_intervals
        obj.df = pd.DataFrame(
            columns=["phase", "composition", "temperature", "free_energy"]
        )
        obj.tangents = None
        obj.temperatures = None
        obj.tangent_types = None
        obj._calc_kwargs = {}
        obj._calphad_surfaces = surfaces
        obj._line_compound_fits = line_compound_fits
        obj._compound_2sl_fits = compound_2sl_fits
        obj._tdb_metadata = metadata
        return obj

    @staticmethod
    def from_parquet(filename):
        """
        Load a :class:`PhaseDiagram` previously saved with :meth:`to_parquet`.

        The original folder structure is not required; only the processed
        DataFrame and metadata stored inside the Parquet file are used.

        Parameters
        ----------
        filename : str
            Path to the Parquet file.

        Returns
        -------
        PhaseDiagram
        """
        import json
        import pyarrow.parquet as pq

        table = pq.read_table(filename)
        meta = table.schema.metadata or {}

        reference_element = meta[b"calphy_reference_element"].decode()
        phases = json.loads(meta[b"calphy_phases"].decode())
        composition_intervals = {
            k: tuple(v)
            for k, v in json.loads(meta[b"calphy_composition_intervals"].decode()).items()
        }

        flat_df = table.to_pandas()

        # Reconstruct array columns by grouping on the integer row_id.
        scalar_cols = [
            c for c in flat_df.columns if c not in ("row_id", "temperature", "free_energy")
        ]
        rows = []
        for row_id, grp in flat_df.groupby("row_id", sort=True):
            grp = grp.sort_values("temperature")
            row = {c: grp[c].iloc[0] for c in scalar_cols}
            row["temperature"] = grp["temperature"].to_numpy()
            row["free_energy"] = grp["free_energy"].to_numpy()
            rows.append(row)

        reconstructed_df = pd.DataFrame(rows)

        obj = object.__new__(PhaseDiagram)
        obj.reference_element = reference_element
        obj.phases = phases
        obj.composition_intervals = composition_intervals
        obj.df = reconstructed_df
        obj.tangents = None
        obj.temperatures = None
        obj.tangent_types = None
        obj._calc_kwargs = {}
        obj._calphad_surfaces = {}
        return obj

    @classmethod
    def from_df(cls, df, reference_element, phases=None, composition_intervals=None):
        """
        Construct a :class:`PhaseDiagram` directly from a DataFrame.

        This bypasses folder reading and all pre-processing (``gather_results``,
        ``clean_df``, etc.).  Use it when you already have a tidy DataFrame
        — e.g. one previously obtained from :attr:`PhaseDiagram.df` or
        assembled manually.

        Parameters
        ----------
        df : pandas.DataFrame
            Must contain at least the columns:

            * ``phase`` — str, phase label (e.g. ``'fcc'``, ``'lqd'``).
            * ``composition`` — float, reference-element mole fraction.
            * ``temperature`` — array-like of floats (one per row).
            * ``free_energy`` — array-like of floats (one per row).

            Any additional columns are preserved unchanged.

        reference_element : str
            The element whose fraction defines the composition axis
            (e.g. ``'Ag'``).

        phases : list of str, optional
            Ordered list of phase names.  Defaults to the unique values of
            ``df['phase']`` in the order they first appear.

        composition_intervals : dict, optional
            ``{phase: (x_lo, x_hi)}`` bounds for each phase.  Phases not
            supplied are auto-detected from the data.

        Returns
        -------
        PhaseDiagram
        """
        import numpy as np

        df = df.copy()

        # Ensure array columns are numpy arrays
        for col in ("temperature", "free_energy"):
            if col in df.columns:
                df[col] = df[col].apply(np.asarray)

        if phases is None:
            # Preserve insertion order
            seen = {}
            for p in df["phase"]:
                seen[p] = None
            phases = list(seen)

        comp_intervals = dict(composition_intervals or {})
        for phase in phases:
            if phase not in comp_intervals:
                df_p = df.loc[df["phase"] == phase, "composition"]
                if len(df_p) > 0:
                    comp_intervals[phase] = (float(df_p.min()), float(df_p.max()))

        obj = object.__new__(cls)
        obj.reference_element = reference_element
        obj.phases = phases
        obj.composition_intervals = comp_intervals
        obj.df = df.reset_index(drop=True)
        obj.tangents = None
        obj.temperatures = None
        obj.tangent_types = None
        obj._calc_kwargs = {}
        obj._calphad_surfaces = {}
        return obj


def plot_pd(ax, pd_obj, phase_colors=None, two_phase_alpha=0.35,
            boundary_lw=1.2, boundary_color='k'):
    """
    Plot a phase diagram onto an existing matplotlib Axes.

    Only two-phase regions are colored; single-phase regions are uncolored.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    pd_obj : PhaseDiagram
        After calculate() has been called.
    phase_colors : dict, optional
        Mapping from tangent-type string ``'phaseA-phaseB'`` to a colour.
        Regions not listed fall back to ``'#cccccc'``.
    two_phase_alpha : float
        Alpha for the two-phase fill regions.
    boundary_lw : float
        Line width of the two-phase boundary lines.
    boundary_color : str
        Colour of the two-phase boundary lines.
    """
    if phase_colors is None:
        phase_colors = {}

    if pd_obj.temperatures is None:
        raise ValueError("pd_obj.calculate() must be called before plotting.")

    temps    = np.array(pd_obj.temperatures)
    tangents = pd_obj.tangents
    ttypes   = pd_obj.tangent_types

    two_phase = defaultdict(list)

    for T, tlist, typlist in zip(temps, tangents, ttypes):
        if len(tlist) == 0:
            continue
        order   = np.argsort(tlist[:, 0])
        tlist   = tlist[order]
        typlist = typlist[order]

        for (xl, xr), ttype in zip(tlist, typlist):
            two_phase[ttype].append((T, xl, xr))

    for ttype, entries in two_phase.items():
        arr = np.array(entries)
        T_vals, xl_vals, xr_vals = arr[:, 0], arr[:, 1], arr[:, 2]
        color = phase_colors.get(ttype, '#cccccc')
        ax.fill_betweenx(T_vals, xl_vals, xr_vals,
                         color=color, alpha=two_phase_alpha, label=ttype)
        ax.plot(xl_vals, T_vals, color=boundary_color, lw=boundary_lw)
        ax.plot(xr_vals, T_vals, color=boundary_color, lw=boundary_lw)
