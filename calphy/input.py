"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^2, Ralf Drautz^2
^1: Max Planck Institut für Eisenforschung, Dusseldorf, Germany
^2: Ruhr-University Bochum, Bochum, Germany

calphy is published and distributed under the Academic Software License v1.0 (ASL).
calphy is distributed in the hope that it will be useful for non-commercial academic research,
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
calphy API is published and distributed under the BSD 3-Clause "New" or "Revised" License
See the LICENSE FILE for more details.

More information about the program can be found in:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz.
“Automated Free Energy Calculation from Atomistic Simulations.” Physical Review Materials 5(10), 2021
DOI: 10.1103/PhysRevMaterials.5.103801

For more information contact:
sarath.menon@ruhr-uni-bochum.de/yury.lysogorskiy@icams.rub.de
"""

from typing_extensions import Annotated
from typing import Any, Callable, Dict, List, ClassVar, Literal, Optional, Union
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationError,
    model_validator,
    field_validator,
    conlist,
    PrivateAttr,
)
import difflib
from pydantic.functional_validators import AfterValidator, BeforeValidator
from annotated_types import Len
import mendeleev
from tqdm import tqdm

import yaml
import numpy as np
import copy
import datetime
import itertools
import os
import warnings
from pyscal3 import System
from pyscal3.core import structure_dict, element_dict, _make_crystal
from ase.io import read, write
import shutil

__version__ = "1.8.5"


def _check_equal(val):
    if not (val[0] == val[1] == val[2]):
        return False
    return True


def to_list(v: Any) -> List[Any]:
    return np.atleast_1d(v)


def _to_str(val):
    if np.isscalar(val):
        return str(val)
    else:
        return [str(x) for x in val]


def _to_int(val):
    if np.isscalar(val):
        return int(val)
    else:
        return [int(x) for x in val]


def _to_none(val):
    if val in [
        "none",
        "None",
    ]:
        return None
    return val


def _to_float(val):
    if np.isscalar(val):
        return float(val)
    else:
        return [float(x) for x in val]


def _extract_elements_from_pair_coeff(pair_coeff_string):
    """
    Extract element symbols from pair_coeff string.
    Returns None if pair_coeff doesn't contain element specifications.

    Parameters
    ----------
    pair_coeff_string : str
        The pair_coeff command string (e.g., "* * potential.eam.fs Cu Zr")

    Returns
    -------
    list or None
        List of element symbols in order, or None if no elements found
    """
    if pair_coeff_string is None:
        return None

    pcsplit = pair_coeff_string.strip().split()
    elements = []

    # Start collecting after we find element symbols
    # Elements are typically after the potential filename
    started = False

    for p in pcsplit:
        # Check if this looks like an element symbol
        # Element symbols are 1-2 characters, start with uppercase
        if len(p) <= 2 and p[0].isupper():
            try:
                # Verify it's a valid element using mendeleev
                _ = mendeleev.element(p)
                elements.append(p)
                started = True
            except Exception:
                # Not a valid element, might be done collecting
                if started:
                    # We already started collecting elements and hit a non-element
                    break

    return elements if len(elements) > 0 else None


# --------------------------------------------------------------------------- #
# Strict input models: unknown keys are hard errors, with hints.
# A typo in an input key must never silently fall back to a default.
# --------------------------------------------------------------------------- #

#: input keys that existed in calphy v1, with their migration message
_REMOVED_KEYS = {
    "savefile": "job-state pickling was removed in calphy v2; rerun from the input file",
    "save_job": "job-state pickling was removed in calphy v2; rerun from the input file",
    "load_job": "job-state pickling was removed in calphy v2; rerun from the input file",
    "seed": "use md.seed -- one master seed now controls every stochastic step "
            "(the old quantum_thermal_bath seed was never actually applied)",
}

#: lazily built: ({field name -> [locations]}, {model class -> location label})
_INPUT_LAYOUT = None


def _input_layout():
    """Map every input field to the block(s) it lives in.

    Built lazily on the first unknown-key error, once every model below is
    defined.  Nested blocks are discovered from Calculation's fields whose
    defaults are input models, so new blocks are picked up automatically.
    """
    global _INPUT_LAYOUT
    if _INPUT_LAYOUT is None:
        field_locs = {}
        model_locs = {Calculation: "the calculation block"}
        for name, field in Calculation.model_fields.items():
            field_locs.setdefault(name, []).append("the calculation block")
            if isinstance(field.default, _StrictInput):
                block_cls = type(field.default)
                model_locs[block_cls] = "the '%s:' block" % name
                for fname in block_cls.model_fields:
                    field_locs.setdefault(fname, []).append("the '%s:' block" % name)
        _INPUT_LAYOUT = (field_locs, model_locs)
    return _INPUT_LAYOUT


def _unknown_key_message(cls, key):
    """One error line for unknown ``key`` on model ``cls``, with the best hint."""
    if key in _REMOVED_KEYS:
        return "input key '%s': %s" % (key, _REMOVED_KEYS[key])
    field_locs, model_locs = _input_layout()
    here = model_locs.get(cls, "the '%s' block" % cls.__name__)
    msg = "unknown input key '%s' in %s" % (key, here)
    close = difflib.get_close_matches(key, cls.model_fields.keys(), n=1)
    if close:
        return "%s -- did you mean '%s'?" % (msg, close[0])
    if key in field_locs:  # right key, wrong block
        return "%s -- '%s' belongs in %s" % (msg, key, field_locs[key][0])
    close = difflib.get_close_matches(key, field_locs.keys(), n=1, cutoff=0.75)
    if close:
        return "%s -- did you mean '%s' (in %s)?" % (
            msg, close[0], field_locs[close[0]][0],
        )
    return msg


class _StrictInput(BaseModel):
    """Base class for every input block: unknown keys raise, with hints."""

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="before")
    @classmethod
    def _reject_unknown_keys(cls, data):
        if isinstance(data, dict):
            unknown = [
                k for k in data
                if isinstance(k, str) and k not in cls.model_fields
            ]
            if unknown:
                raise ValueError(
                    "; ".join(_unknown_key_message(cls, k) for k in unknown)
                )
        return data


class UFMP(_StrictInput, title="UFM potential input options"):
    p: Annotated[float, Field(default=50.0)]
    # sigma may be a scalar (single-component UFM reference, original behaviour) or a
    # mapping of per-element-pair length scales for the two-leg UFM reference path, e.g.
    #   sigma: {H_H: 0.9, O_O: 2.8, H_O: 0.9}
    # Keys are "<elementA>_<elementB>" (order-insensitive). Missing cross terms are
    # filled by LAMMPS geometric mixing.
    sigma: Annotated[Union[float, Dict[str, float]], Field(default=1.5)]
    # When set, activates the two-leg UFM path: the system is first switched from the
    # real potential to the (possibly multi-component) UFM at `sigma`, then from that
    # UFM to a single-component UFM at `single_sigma` whose absolute free energy is
    # known analytically (find_fe). Leaving this None reproduces the original
    # single-leg behaviour exactly.
    single_sigma: Annotated[Optional[float], Field(default=None)]
    # p for the single-component endpoint; defaults to `p` when not given.
    single_p: Annotated[Optional[float], Field(default=None)]


class MonteCarlo(
    _StrictInput, title="Options for Monte Carlo moves during particle swap moves"
):
    n_steps: Annotated[
        int, Field(default=1, gt=0, description="perform swap moves every n_step")
    ]
    n_swaps: Annotated[
        int, Field(default=0, ge=0, description="number of swap moves to perform")
    ]
    forward_swap_types: Annotated[
        List[int],
        Field(default=[], description="atom types to swap during forward integration"),
    ]
    reverse_swap_types: Annotated[
        List[int],
        Field(default=[], description="atom types to swap during reverse integration"),
    ]
    allow_all_swaps: Annotated[
        bool,
        Field(
            default=True,
            description="allow swapping between all atom types including fictitious ones",
        ),
    ]
    use_custom_lammps: Annotated[
        bool,
        Field(
            default=False,
            description="Whether to use the custom modified LAMMPS version",
        ),
    ]


class CompositionScaling(_StrictInput, title="Composition scaling input options"):
    _input_chemical_composition: PrivateAttr(default=None)
    output_chemical_composition: Annotated[dict, Field(default={}, required=False)]
    restrictions: Annotated[
        List[str], BeforeValidator(to_list), Field(default=[], required=False)
    ]


class MD(_StrictInput, title="MD specific input options"):
    timestep: Annotated[
        float,
        Field(
            default=0.001,
            description="timestep for md simulation",
            example="timestep: 0.001",
        ),
    ]
    n_small_steps: Annotated[int, Field(default=10000, gt=0)]
    n_every_steps: Annotated[int, Field(default=10, gt=0)]
    n_repeat_steps: Annotated[int, Field(default=10, gt=0)]
    n_cycles: Annotated[int, Field(default=100, gt=0)]
    # scalar in the input file (the integration-stage damping); Phase.__init__
    # rewrites these at runtime into [equilibration, integration] pairs, with
    # the equilibration value taken from the nose_hoover/berendsen/qtb block.
    # (The former Union[float, 2-list] input form never validated -- gt broke
    # the list branch -- and the drivers never handled a user-given list.)
    thermostat_damping: Annotated[float, Field(default=0.1, gt=0)]
    barostat_damping: Annotated[float, Field(default=0.1, gt=0)]
    cmdargs: Annotated[str, Field(default="")]
    init_commands: Annotated[List, Field(default=[])]
    # master seed for every stochastic choice a job makes: velocity creation,
    # langevin/atom-swap/qtb fix seeds, composition-scaling atom picks. None
    # (the default) draws a fresh seed per run; the seed actually used is
    # always backfilled into the simfolder copy of the input file, so any run
    # can be reproduced exactly by rerunning that file.
    seed: Annotated[Optional[int], Field(default=None, gt=0)]


class NoseHoover(_StrictInput, title="Specific input options for Nose-Hoover thermostat"):
    thermostat_damping: Annotated[float, Field(default=0.1, gt=0)]
    barostat_damping: Annotated[float, Field(default=0.1, gt=0)]


class Berendsen(_StrictInput, title="Specific input options for Berendsen thermostat"):
    thermostat_damping: Annotated[float, Field(default=100.0, gt=0)]
    barostat_damping: Annotated[float, Field(default=100.0, gt=0)]


class QuantumThermalBath(_StrictInput, title="Dammak quantum thermal bath (LAMMPS fix qtb)"):
    """
    Colored-noise Langevin thermostat that injects quantum statistics into
    classical MD (Dammak et al., Phys. Rev. Lett. 103, 190601, 2009).
    Active when ``mode: fe-qtb`` is set at the top level. The QTB
    thermostat is paired with ``fix nph`` (NPT) or ``fix nve`` (NVT)
    inside calphy; do NOT also pair with Nose-Hoover or Langevin.
    """
    thermostat_damping: Annotated[float, Field(default=0.1, gt=0)]
    barostat_damping: Annotated[float, Field(default=0.1, gt=0)]
    f_max: Annotated[float, Field(default=200.0, gt=0,
        description="Upper cutoff frequency for QTB power spectrum (THz). "
                    "Must exceed the highest phonon frequency of the system.")]
    n_f: Annotated[int, Field(default=100, gt=0,
        description="Number of frequency bins discretising the QTB spectrum.")]


class Queue(_StrictInput, title="Options for configuring queue"):
    scheduler: Annotated[str, Field(default="local")]
    cores: Annotated[int, Field(default=1, gt=0)]
    jobname: Annotated[str, Field(default="calphy")]
    walltime: Annotated[str, Field(default="23:59:00")]
    queuename: Annotated[str, Field(default="")]
    memory: Annotated[str, Field(default="3GB")]
    commands: Annotated[List, Field(default=[])]
    options: Annotated[Dict[str, str], Field(default={})]


class Tolerance(_StrictInput, title="Tolerance settings for convergence"):
    lattice_constant: Annotated[float, Field(default=0.0002, ge=0)]
    spring_constant: Annotated[float, Field(default=0.1, gt=0)]
    # Structural phase-stability checks during equilibration are OFF by
    # default: solid_fraction=0 means the melt check (solid_fraction < this)
    # never fires, and liquid_fraction=1.0 means the solidify check
    # (solid_fraction > this) never fires.  Set solid_fraction > 0 (e.g. 0.7)
    # to re-enable melt detection for a solid run, or liquid_fraction < 1
    # (e.g. 0.05) to re-enable solidification detection for a liquid run.
    solid_fraction: Annotated[float, Field(default=0.0, ge=0)]
    liquid_fraction: Annotated[float, Field(default=1.0, ge=0)]
    pressure: Annotated[float, Field(default=10.0, ge=0)]


class PhaseTransitionDetection(_StrictInput, title="Settings for the pre-flight temperature-range scan"):
    mode: Annotated[
        Literal["none", "adapt", "warn", "stop"],
        Field(
            default="none",
            description=(
                "Controls the pre-flight temperature-range scan that runs "
                "before a reversible-scaling (ts) sweep.  The scan performs a "
                "single fast real-thermostat temperature ramp (T0 -> Tf under "
                "NPT) and watches the fluctuation response functions for the "
                "onset of a phase transition.\n"
                "  'none'  — the scan is disabled; the ts sweep runs over the "
                "requested [T0, Tf] range as-is (default).\n"
                "  'adapt' — if the scan detects a transition, reduce the upper "
                "temperature to the detected clean onset and run the ts sweep "
                "over [T0, T_clean] (keeping the same number of switching "
                "steps).  If the scan is clean the range is unchanged.\n"
                "  'warn'  — run the scan and log the detected clean range, but "
                "do NOT modify the calculation; the ts sweep runs over the full "
                "requested range.  Use this to observe detection without "
                "changing the outcome.\n"
                "  'stop'  — if a transition is detected, raise "
                "PhaseTransitionError reporting the clean range so you can "
                "re-submit with a corrected temperature range."
            ),
        ),
    ]
    prescan_steps: Annotated[
        int,
        Field(
            default=20000,
            ge=1000,
            description=(
                "Number of MD steps for the pre-flight temperature ramp from "
                "T0 to Tf.  Typically a little shorter than the production "
                "switching length.  Note that a *faster* ramp superheats (or "
                "supercools) the metastable phase further and yields noisier "
                "response-function signals, pushing the detected onset closer "
                "to the collapse; if the adapted ts sweep melts/freezes during "
                "its equilibration, increase this (or lower onset_fraction).  "
                "Default 20000."
            ),
        ),
    ]
    onset_fraction: Annotated[
        float,
        Field(
            default=0.85,
            gt=0,
            le=1.0,
            description=(
                "Fractional safety margin applied to the detected onset: the "
                "adapted upper temperature is "
                "T_clean = T0 + onset_fraction * (T_onset - T0).  The detected "
                "onset is the foot of the deviation in a *fast* diagnostic ramp, "
                "but a reversible-scaling sweep then EQUILIBRATES at the boundary "
                "for many steps and the metastable phase has more time to "
                "nucleate the transition there — so its practical stability "
                "limit is somewhat below the ramp onset, and the exact onset is "
                "also noisy.  Backing the boundary off by a fraction of the "
                "super-heated/cooled span makes the adapted sweep robust to both "
                "effects.  1.0 disables the margin (cut exactly at the onset); "
                "lower values give more margin (and less usable range).  This is "
                "the main knob to turn if an adapted ts run still melts/freezes. "
                "Default 0.85."
            ),
        ),
    ]
    # Detector-internal calibration (peak_threshold, min_agreement, onset_sigma,
    # onset_level, the smoothing/fluctuation windows, slope-break parameters …)
    # are intentionally NOT exposed here.  They are stable defaults living in
    # calphy.range_scan.RangeScan; onset_sigma in particular had no effect on the
    # result (the slope-break onset always dominates the final boundary), and
    # onset_level overlapped with onset_fraction.  Users tune the scan through
    # mode, prescan_steps and onset_fraction only.

class MeltingTemperature(_StrictInput, title="Input options for melting temperature mode"):
    guess: Annotated[Union[float, None], Field(default=None, gt=0)]
    step: Annotated[int, Field(default=200, ge=20)]
    attempts: Annotated[int, Field(default=5, ge=1)]


class MaterialsProject(_StrictInput, title="Input options for materials project"):
    api_key: Annotated[str, Field(default="", exclude=True)]
    conventional: Annotated[bool, Field(default=True)]
    target_natoms: Annotated[
        int,
        Field(
            default=1500,
            description="The structure parsed from materials project would be repeated to approximately this value",
        ),
    ]

    @field_validator("api_key", mode="after")
    def resolve_api_key(cls, v: str) -> str:
        if not v:
            return v
        value = os.getenv(v)
        if not value:
            raise ValueError(
                f"Environment variable '{v}' not found or empty. "
                f"Set it before running, e.g.:\n  export {v}='your_api_key_here'"
            )
        return value


class Calculation(_StrictInput, title="Main input class"):
    monte_carlo: Optional[MonteCarlo] = MonteCarlo()
    composition_scaling: Optional[CompositionScaling] = CompositionScaling()
    md: Optional[MD] = MD()
    nose_hoover: Optional[NoseHoover] = NoseHoover()
    berendsen: Optional[Berendsen] = Berendsen()
    quantum_thermal_bath: Optional[QuantumThermalBath] = QuantumThermalBath()
    queue: Optional[Queue] = Queue()
    tolerance: Optional[Tolerance] = Tolerance()
    phase_transition_detection: Optional[PhaseTransitionDetection] = PhaseTransitionDetection()
    uhlenbeck_ford_model: Optional[UFMP] = UFMP()
    melting_temperature: Optional[MeltingTemperature] = MeltingTemperature()
    materials_project: Optional[MaterialsProject] = MaterialsProject()

    element: Annotated[List[str], BeforeValidator(to_list), Field(default=[])]
    n_elements: Annotated[int, Field(default=0)]
    mass: Annotated[List[float], BeforeValidator(to_list), Field(default=[])]
    _element_dict: dict = PrivateAttr(default={})
    kernel: Annotated[int, Field(default=0)]
    inputfile: Annotated[str, Field(default="")]

    mode: Annotated[Union[str, None], Field(default=None)]
    lattice: Annotated[str, Field(default="")]
    file_format: Annotated[str, Field(default="lammps-data")]

    # pressure properties
    pressure: Annotated[
        Union[
            None,
            float,
            conlist(float, min_length=1, max_length=3),  # Allow lists of 1-3 floats
            conlist(
                conlist(float, min_length=3, max_length=3), min_length=1, max_length=2
            ),
        ],
        Field(default=0),
    ]

    _pressure: float = PrivateAttr(default=None)
    _pressure_stop: float = PrivateAttr(default=None)
    _pressure_input: Any = PrivateAttr(default=None)
    _iso: bool = PrivateAttr(default=False)
    _pressure_coupling: str = PrivateAttr(default="iso")
    _fix_lattice: bool = PrivateAttr(default=False)

    pressure_coupling: Annotated[
        Union[str, None],
        Field(
            default=None,
            description="Pressure coupling style for LAMMPS barostat: 'iso', 'aniso', or 'tri'. "
            "If not set, 'iso' is used for scalar/1D pressure and 'aniso' for 3-component pressure.",
        ),
    ]

    temperature: Annotated[
        Union[float, conlist(float, min_length=1, max_length=2)], Field(default=0)
    ]
    temperature_high: Annotated[float, Field(default=0.0)]
    _temperature: float = PrivateAttr(default=None)
    _temperature_high: float = PrivateAttr(default=None)
    _temperature_stop: float = PrivateAttr(default=None)
    _temperature_input: float = PrivateAttr(default=None)

    melting_cycle: Annotated[bool, Field(default=True)]

    pair_style: Annotated[
        Union[List[str], None], BeforeValidator(to_list), Field(default=None)
    ]
    pair_coeff: Annotated[
        Union[List[str], None], BeforeValidator(to_list), Field(default=None)
    ]
    pair_mode: Annotated[Union[str, None], Field(default=None)]
    potential_file: Annotated[Union[str, None], Field(default=None)]
    fix_potential_path: Annotated[bool, Field(default=True)]
    _pair_style_with_options: List[str] = PrivateAttr(default=None)

    reference_phase: Annotated[str, Field(default="")]
    lattice_constant: Annotated[float, Field(default=0)]
    repeat: Annotated[
        conlist(int, min_length=3, max_length=3), Field(default=[1, 1, 1])
    ]

    script_mode: Annotated[bool, Field(default=False)]
    execution_mode: Annotated[str, Field(default="executable")]
    lammps_executable: Annotated[Union[str, None], Field(default=None)]
    mpi_executable: Annotated[Union[str, None], Field(default=None)]

    npt: Annotated[bool, Field(default=True)]
    n_equilibration_steps: Annotated[int, Field(default=25000)]
    n_switching_steps: Annotated[
        Union[int, conlist(int, min_length=2, max_length=2)],
        Field(default=[50000, 50000]),
    ]
    _n_switching_steps: int = PrivateAttr(default=50000)
    _n_sweep_steps: int = PrivateAttr(default=50000)
    n_print_steps: Annotated[int, Field(default=0)]
    n_iterations: Annotated[int, Field(default=1)]
    lambda_schedule: Annotated[str, Field(default="linear")]
    equilibration_control: Annotated[Union[str, None], Field(default=None)]
    folder_prefix: Annotated[Union[str, None], Field(default=None)]

    # add second level options; for example spring constants
    spring_constants: Annotated[Union[List[float], None], Field(default=None)]

    # some input keywords that will be used for the phase diagram mode
    phase_name: Annotated[str, Field(default="")]
    reference_composition: Annotated[float, Field(default=0.00)]

    # structure items
    _structure: Any = PrivateAttr(default=None)

    # just check for nlements in compscale
    _totalelements = PrivateAttr(default=0)

    # internal flag set when mode == "fe-qtb"; threaded through phase.py /
    # solid.py to switch the thermostat to QTB and the Einstein reference
    # to the quantum harmonic-oscillator form.
    _qtb: bool = PrivateAttr(default=False)

    @model_validator(mode="after")
    def _validate_all(self) -> "Input":
        if self.script_mode:
            raise ValueError(
                "script_mode was removed in calphy v2. calphy now always drives "
                "the LAMMPS executable directly; set lammps_executable (or use "
                "$CALPHY_LAMMPS_EXECUTABLE / lmp on PATH) and remove script_mode "
                "from the input file."
            )

        if self.execution_mode not in ("executable", "library"):
            raise ValueError(
                "execution_mode must be 'executable' (drive the lmp binary, "
                "the default) or 'library' (drive LAMMPS through pylammpsmpi; "
                "needs the optional dependency: pip install calphy[library]), "
                "got %r" % self.execution_mode
            )

        if not (len(self.element) == len(self.mass)):
            raise ValueError("mass and elements should have same length")

        # QTB-flavoured fe mode. We resolve it here so all downstream
        # dispatch sees mode=="fe" and the _qtb flag drives QTB behaviour.
        if self.mode == "fe-qtb":
            if self.reference_phase and self.reference_phase.lower() == "liquid":
                raise ValueError(
                    "mode=fe-qtb is solids-only. The Uhlenbeck-Ford liquid "
                    "reference is intrinsically classical and cannot be paired "
                    "with QTB sampling."
                )
            self._qtb = True
            self.mode = "fe"

        self.n_elements = len(self.element)

        if self.potential_file is not None:
            warnings.warn(
                "potential_file is deprecated and will be removed in a future version. "
                "Use pair_style/pair_coeff (with pair_mode: overlay for multi-potential setups) instead.",
                DeprecationWarning,
                stacklevel=2,
            )

        if self.pair_mode is not None:
            self.pair_mode = self.pair_mode.lower()
            if self.pair_mode not in ["overlay"]:
                raise ValueError("pair_mode should be one of: overlay")

        self.lambda_schedule = self.lambda_schedule.lower()
        if self.lambda_schedule not in ["linear", "uniform_temperature"]:
            raise ValueError(
                "lambda_schedule must be one of: 'linear', 'uniform_temperature'"
            )

        if self.pair_mode == "overlay":
            if self.pair_style is None or self.pair_coeff is None:
                raise ValueError("pair_mode overlay requires pair_style and pair_coeff")
            if len(self.pair_style) != len(self.pair_coeff):
                raise ValueError(
                    "pair_mode overlay requires pair_style and pair_coeff to have the same length"
                )
            for ps in self.pair_style:
                if ps.split()[0].startswith("hybrid"):
                    raise ValueError(
                        "pair_mode overlay expects component pair styles, not a hybrid pair_style"
                    )

        # Validate element/mass/pair_coeff ordering consistency
        # This is critical for multi-element systems where LAMMPS type numbers
        # are assigned based on element order: element[0]=Type1, element[1]=Type2, etc.
        if (
            len(self.element) > 1
            and self.pair_coeff is not None
            and len(self.pair_coeff) > 0
        ):
            extracted_elements = _extract_elements_from_pair_coeff(self.pair_coeff[0])

            if extracted_elements is not None:
                # pair_coeff specifies elements - check ordering
                if set(extracted_elements) != set(self.element):
                    raise ValueError(
                        f"Element mismatch between 'element' and 'pair_coeff'!\n"
                        f"  element:    {self.element}\n"
                        f"  pair_coeff: {extracted_elements}\n"
                        f"The elements specified must be the same."
                    )

                if list(extracted_elements) != list(self.element):
                    raise ValueError(
                        f"Element ordering mismatch detected!\n\n"
                        f"  element:    {self.element}\n"
                        f"  pair_coeff: {extracted_elements}\n"
                        f"  mass:       {self.mass}\n\n"
                        f"For multi-element systems, all three must be in the SAME order.\n\n"
                        f"Why this matters:\n"
                        f"  - Element order determines LAMMPS type numbers:\n"
                        f"      element[0] → Type 1, element[1] → Type 2, etc.\n"
                        f"  - The pair_coeff elements must match this type order\n"
                        f"  - The mass values must correspond to the same order\n"
                        f"  - Composition transformations depend on this ordering\n\n"
                        f"Please reorder your input so element, mass, and pair_coeff\n"
                        f"all use the same element ordering."
                    )

        if self.pressure_coupling is not None:
            self.pressure_coupling = self.pressure_coupling.lower()
            if self.pressure_coupling not in ("iso", "aniso", "tri"):
                raise ValueError(
                    "pressure_coupling must be one of 'iso', 'aniso', or 'tri'"
                )

        self._pressure_input = copy.copy(self.pressure)
        if self.pressure is None:
            self._iso = True
            self._fix_lattice = True
            self._pressure = None
            self._pressure_stop = None
        elif np.isscalar(self.pressure):
            self._pressure = self.pressure
            self._pressure_stop = self.pressure
            self._iso = True
            self._fix_lattice = False
        elif np.shape(self.pressure) == (1,):
            self._iso = True
            self._fix_lattice = False
            self._pressure = self.pressure[0]
            self._pressure_stop = self.pressure[0]
        elif np.shape(self.pressure) == (2,):
            self._iso = True
            self._fix_lattice = False
            self._pressure = self.pressure[0]
            self._pressure_stop = self.pressure[1]
        elif np.shape(self.pressure) == (1, 3):
            if not _check_equal(self.pressure[0]):
                raise ValueError("All pressure terms must be equal")
            self._iso = False
            self._fix_lattice = False
            self._pressure = self.pressure[0][0]
            self._pressure_stop = self.pressure[0][0]
        elif np.shape(self.pressure) == (2, 3):
            if not (_check_equal(self.pressure[0]) and _check_equal(self.pressure[1])):
                raise ValueError("All pressure terms must be equal")
            self._iso = False
            self._fix_lattice = False
            self._pressure = self.pressure[0][0]
            self._pressure_stop = self.pressure[1][0]
        else:
            raise ValueError("Unknown format for pressure")

        # Resolve the final LAMMPS barostat keyword.
        # An explicit pressure_coupling always takes priority; otherwise
        # fall back to the legacy _iso flag (True → "iso", False → "aniso").
        if self.pressure_coupling is not None:
            self._pressure_coupling = self.pressure_coupling
        else:
            self._pressure_coupling = "iso" if self._iso else "aniso"

        self._temperature_input = copy.copy(self.temperature)
        # guess a melting temp of the system, this will be mostly ignored
        # chem = mendeleev.element(self.element[0])
        # self._melting_temperature = chem.melting_point
        try:
            chem = mendeleev.element(self.element[0])
            self._melting_temperature = chem.melting_point
        except Exception:
            self._melting_temperature = None

        if self.temperature == 0:
            # the only situation in which it can be None is if mode is melting temp
            if len(self.element) > 1:
                raise ValueError(
                    "Cannot guess start temperature for more than one species, please specify"
                )
            # now try to guess
            if self._melting_temperature is None:
                raise ValueError(
                    "Could not guess start temperature for more than one species, please specify"
                )
            self._temperature = self._melting_temperature
            self._temperature_stop = self._melting_temperature
            if self.temperature_high == 0:
                self._temperature_high = 2 * self._melting_temperature
            else:
                self._temperature_high = self.temperature_high

        elif np.shape(np.atleast_1d(self.temperature)) == (1,):
            temp = np.atleast_1d(self.temperature)
            self._temperature = temp[0]
            self._temperature_stop = temp[0]
            if self.temperature_high == 0:
                self._temperature_high = 2 * temp[0]
            else:
                self._temperature_high = self.temperature_high

        elif np.shape(self.temperature) == (2,):
            temp = self.temperature
            self._temperature = temp[0]
            self._temperature_stop = temp[1]
            if self.temperature_high == 0:
                self._temperature_high = 2 * temp[1]
            else:
                self._temperature_high = self.temperature_high

        # fix pair styles
        # two main lists
        # _pair_style_with_options, read in as it is from file
        # _pair_style_names, just the names of the pair styles

        _pair_style_names = []

        for ps in self.pair_style:
            ps_split = ps.split()
            _pair_style_names.append(ps_split[0])

        # only set if its None
        self._pair_style_with_options = self.pair_style
        self._pair_style_names = _pair_style_names

        # now fix pair coeffs with path
        if self.fix_potential_path:
            self.pair_coeff = self.fix_paths(self.pair_coeff)

        if np.isscalar(self.n_switching_steps):
            self._n_sweep_steps = self.n_switching_steps
            self._n_switching_steps = self.n_switching_steps
        else:
            self._n_sweep_steps = self.n_switching_steps[1]
            self._n_switching_steps = self.n_switching_steps[0]

        # here we also prepare lattice dict
        for count, element in enumerate(self.element):
            self._element_dict[element] = {}
            self._element_dict[element]["mass"] = self.mass[count]
            self._element_dict[element]["count"] = 0
            self._element_dict[element]["composition"] = 0.0
            self._element_dict[element]["atomic_number"] = mendeleev.element(
                element
            ).atomic_number

        # generate temporary filename if needed
        write_structure_file = False
        rename_structure_file = False

        if self.lattice == "":
            # fetch from dict
            if len(self.element) > 1:
                raise ValueError(
                    "Cannot create lattice for more than one element, provide a lammps-data file explicitly"
                )
            if self.element[0] in element_dict.keys():
                self.lattice = element_dict[self.element[0]]["structure"]
                self.lattice_constant = element_dict[self.element[0]][
                    "lattice_constant"
                ]
            else:
                raise ValueError(
                    "Could not find structure, please provide lattice and lattice_constant explicitely"
                )

            if self.repeat == [1, 1, 1]:
                self.repeat = [5, 5, 5]

            structure = _make_crystal(
                self.lattice.lower(),
                lattice_constant=self.lattice_constant,
                repetitions=self.repeat,
                element=self.element,
            )
            structure = structure.write.ase()

            # extract composition
            types, typecounts = np.unique(
                structure.get_chemical_symbols(), return_counts=True
            )

            for c, t in enumerate(types):
                self._element_dict[t]["count"] = typecounts[c]
                self._element_dict[t]["composition"] = typecounts[c] / np.sum(
                    typecounts
                )

            self._natoms = len(structure)
            self._original_lattice = self.lattice.lower()
            write_structure_file = True

        elif self.lattice.lower() in structure_dict.keys():
            if len(self.element) > 1:
                raise ValueError(
                    "Cannot create lattice for more than one element, provide a lammps-data file explicitly"
                )

            # this is a valid structure
            if self.lattice_constant == 0:
                # we try try to get lattice_constant
                if self.element[0] in element_dict.keys():
                    self.lattice_constant = element_dict[self.element[0]][
                        "lattice_constant"
                    ]
                else:
                    raise ValueError("Please provide lattice_constant!")
            # now create lattice
            structure = _make_crystal(
                self.lattice.lower(),
                lattice_constant=self.lattice_constant,
                repetitions=self.repeat,
                element=self.element,
            )
            structure = structure.write.ase()

            # extract composition
            types, typecounts = np.unique(
                structure.get_chemical_symbols(), return_counts=True
            )

            for c, t in enumerate(types):
                self._element_dict[t]["count"] = typecounts[c]
                self._element_dict[t]["composition"] = typecounts[c] / np.sum(
                    typecounts
                )

            # concdict_counts = {str(t): typecounts[c] for c, t in enumerate(types)}
            # concdict_frac = {str(t): typecounts[c]/np.sum(typecounts) for c, t in enumerate(types)}
            # self._composition = concdict_frac
            # self._composition_counts = concdict_counts
            self._natoms = len(structure)
            self._original_lattice = self.lattice.lower()
            write_structure_file = True

        elif self.lattice.split("-")[0] == "mp":
            # confirm here that API key exists
            if not self.materials_project.api_key:
                raise ValueError("could not find API KEY, pls set it.")
            # now we need to fetch the structure
            try:
                from mp_api.client import MPRester
            except ImportError:
                raise ImportError(
                    "Could not import mp_api, make sure you install mp_api package!"
                )
            # now all good
            rest = {
                "use_document_model": False,
                "include_user_agent": True,
                "api_key": self.materials_project.api_key,
            }
            with MPRester(**rest) as mpr:
                docs = mpr.materials.summary.search(material_ids=[self.lattice])

            structures = []
            for doc in docs:
                struct = doc["structure"]
                if self.materials_project.conventional:
                    aseatoms = struct.to_conventional().to_ase_atoms()
                else:
                    aseatoms = struct.to_primitive().to_ase_atoms()
                structures.append(aseatoms)
            structure = structures[0]

            if np.prod(self.repeat) == 1:
                x = int(
                    np.ceil(
                        (self.materials_project.target_natoms / len(structure))
                        ** (1 / 3)
                    )
                )
                structure = structure.repeat(x)
            else:
                structure = structure.repeat(self.repeat)

            # extract composition
            types, typecounts = np.unique(
                structure.get_chemical_symbols(), return_counts=True
            )

            for c, t in enumerate(types):
                self._element_dict[t]["count"] = typecounts[c]
                self._element_dict[t]["composition"] = typecounts[c] / np.sum(
                    typecounts
                )

            self._natoms = len(structure)
            self._original_lattice = self.lattice.lower()
            write_structure_file = True

        else:
            # this is a file
            if not os.path.exists(self.lattice):
                raise ValueError(f"File {self.lattice} could not be found")
            if self.file_format == "lammps-data":
                # create atomic numbers for proper reading
                Z_of_type = dict(
                    [
                        (count + 1, self._element_dict[element]["atomic_number"])
                        for count, element in enumerate(self.element)
                    ]
                )
                structure = read(
                    self.lattice,
                    format="lammps-data",
                    style="atomic",
                    Z_of_type=Z_of_type,
                )
                # structure = System(aseobj, format='ase')
                rename_structure_file = True
            else:
                raise TypeError("Only lammps-data files are supported!")

            # extract composition
            # this is the types read in from the file
            types, typecounts = np.unique(
                structure.get_chemical_symbols(), return_counts=True
            )
            for c, t in enumerate(types):
                self._element_dict[t]["count"] = typecounts[c]
                self._element_dict[t]["composition"] = typecounts[c] / np.sum(
                    typecounts
                )

            self._natoms = len(structure)
            self._original_lattice = os.path.basename(self.lattice)
            self.lattice = os.path.abspath(self.lattice)

        # if needed, write the file out
        if write_structure_file:
            structure_filename = ".".join(
                [self.create_identifier(), str(self.kernel), "data"]
            )
            structure_filename = os.path.join(os.getcwd(), structure_filename)
            write(structure_filename, structure, format="lammps-data")
            self.lattice = structure_filename

        if rename_structure_file:
            structure_filename = ".".join(
                [self.create_identifier(), str(self.kernel), "data"]
            )
            structure_filename = os.path.join(os.getcwd(), structure_filename)
            shutil.copy(self.lattice, structure_filename)
            self.lattice = structure_filename

        if self.mode == "composition_scaling":
            # we also should check if actual contents are present
            input_chem_comp = {}
            for key, val in self._element_dict.items():
                input_chem_comp[key] = val["count"]
            self.composition_scaling._input_chemical_composition = input_chem_comp

            # now we should check output chem comp and see there are no keys extra
            for key in self.composition_scaling.output_chemical_composition.keys():
                if (
                    key
                    not in self.composition_scaling._input_chemical_composition.keys()
                ):
                    raise ValueError(
                        "An element in output composition is not possible with the given potential"
                    )

            natoms1 = np.sum(
                [
                    val
                    for key, val in self.composition_scaling._input_chemical_composition.items()
                ]
            )
            natoms2 = np.sum(
                [
                    val
                    for key, val in self.composition_scaling.output_chemical_composition.items()
                ]
            )
            if not (natoms1 == natoms2):
                raise ValueError(
                    f"Input and output number of atoms are not conserved! Input {self.composition_scaling._input_chemical_composition}, output {self.composition_scaling.output_chemical_composition}, total atoms in structure {len(structure)}"
                )
        return self

    def fix_paths(self, potlist):
        """
        Fix paths for potential files to complete ones.

        Also expands environment variables (e.g. ``$USER``, ``${USER}``) and
        the home-directory shortcut ``~`` in potential file paths before
        resolving them to absolute paths.  This allows input files to contain
        portable paths such as ``/home/$USER/potentials/Cu.eam``.
        """
        fixedpots = []
        pair_mode = getattr(self, "pair_mode", None)
        pair_style_names = getattr(self, "_pair_style_names", None) or []
        for pot in potlist:
            pcraw = pot.split()
            path_index = 2
            if (
                pair_mode == "overlay"
                and len(pcraw) >= 4
                and pcraw[2] in pair_style_names
            ):
                path_index = 3
            if len(pcraw) > path_index:
                raw_filename = pcraw[path_index]
                # Expand environment variables ($USER, ${USER}, $HOME, …) and ~
                expanded = os.path.expandvars(os.path.expanduser(raw_filename))
                abs_filename = os.path.abspath(expanded)
                if os.path.exists(abs_filename):
                    pcnew = " ".join(
                        [*pcraw[:path_index], abs_filename, *pcraw[path_index + 1 :]]
                    )
                    fixedpots.append(pcnew)
                elif expanded != raw_filename:
                    # A variable was expanded even though the file was not found
                    # at that location yet – still substitute so LAMMPS receives
                    # the resolved path rather than a literal "$USER" string.
                    pcnew = " ".join(
                        [*pcraw[:path_index], expanded, *pcraw[path_index + 1 :]]
                    )
                    fixedpots.append(pcnew)
                else:
                    fixedpots.append(pot)
            else:
                fixedpots.append(pot)
        return fixedpots

    def create_identifier(self):
        """
        Generate an identifier

        Parameters
        ----------
        calc: dict
            a calculation dict

        Returns
        -------
        identistring: string
            unique identification string
        """
        # lattice processed
        prefix = self.mode

        # Handle both validated and non-validated objects
        if self._temperature is not None:
            ts = int(self._temperature)
        else:
            # Fallback for non-validated objects
            temp = (
                self.temperature
                if np.isscalar(self.temperature)
                else self.temperature[0]
            )
            ts = int(temp)

        if hasattr(self, "_pressure") and self._pressure is not None:
            ps = "%d" % (int(self._pressure))
        elif self.pressure is None:
            ps = "None"
        elif np.isscalar(self.pressure):
            ps = "%d" % (int(self.pressure))
        else:
            if isinstance(self.pressure[0], list):
                ps = "%d" % (int(self.pressure[0][0]))
            else:
                ps = "%d" % (int(self.pressure[0]))

        if hasattr(self, "_original_lattice") and self._original_lattice:
            l = self._original_lattice
        else:
            l = self.lattice
        l = l.split("/")
        l = l[-1]

        if self.folder_prefix is None:
            identistring = "-".join(
                [prefix, l.lower(), self.reference_phase, str(ts), str(ps)]
            )
        else:
            identistring = "-".join(
                [
                    self.folder_prefix,
                    prefix,
                    l.lower(),
                    self.reference_phase,
                    str(ts),
                    str(ps),
                ]
            )
        return identistring

    def get_folder_name(self):
        identistring = self.create_identifier()
        simfolder = os.path.join(os.getcwd(), identistring)
        return simfolder

    def create_folders(self):
        """
        Create the necessary folder for calculation

        Parameters
        ----------
        calc : dict
            calculation block

        Returns
        -------
        folder : string
            create folder
        """
        simfolder = self.get_folder_name()

        # if folder exists, delete it -> then create
        if os.path.exists(simfolder):
            raise ValueError(
                f"Simulation folder {simfolder} exists. Please remove and run again!"
            )

        os.mkdir(simfolder)
        return simfolder


def read_inputfile(file, validate=True):
    """Read input file and parse calculations.

    Parameters
    ----------
    file : str
        Path to input YAML file
    validate : bool, optional
        If True, perform full Pydantic validation (structure creation, etc.).
        If False, skip expensive validation (faster, for job submission).
        Default is True for backward compatibility and to ensure tests pass.

    Returns
    -------
    list
        List of Calculation objects
    """
    if not os.path.exists(file):
        raise FileNotFoundError(f"Input file {file} not found.")

    return _read_inputfile(file, validate=validate)


def _read_inputfile(file, validate=False):
    with open(file, "r") as fin:
        data = yaml.safe_load(fin)
    if not isinstance(data, dict) or "calculations" not in data:
        if isinstance(data, dict) and "element" in data:
            raise ValueError(
                "%s uses the legacy calphy input format (top-level keys "
                "instead of a 'calculations:' list). Automatic conversion "
                "was removed in calphy v2; convert the file once with an "
                "older calphy (pip install 'calphy<2', then "
                "calphy_convert_input) or restructure it by hand." % file
            )
        raise ValueError("%s has no 'calculations:' list" % file)
    stray = [k for k in data if k != "calculations"]
    if stray:
        raise ValueError(
            "unknown top-level key(s) %s in %s -- everything except "
            "'calculations:' goes inside a calculation block" % (stray, file)
        )
    calculations = []
    for count, calc in enumerate(tqdm(data["calculations"])):
        calc["kernel"] = count
        calc["inputfile"] = file
        if "pressure" in calc.keys():
            calc["pressure"] = _to_none(calc["pressure"])
        if validate:
            # Full validation - used when actually running calculations
            calculations.append(Calculation(**calc))
        else:
            # Skip expensive validation - for fast parsing
            calculations.append(Calculation.model_construct(**calc))
    return calculations


def generate_metadata():
    metadata = {}
    metadata["software"] = {}
    metadata["software"]["name"] = "calphy"
    metadata["software"]["doi"] = "10.5281/zenodo.10527452"
    metadata["software"]["version"] = __version__
    metadata["software"]["repository"] = "https://github.com/ICAMS/calphy"
    metadata["software"]["webpage"] = "https://calphy.org/"

    metadata["files"] = {}
    metadata["files"]["input_file.yml"] = "input file"
    metadata["files"]["report.yaml"] = "results after thermodynamic integration"
    metadata["files"]["input_configuration.data"] = "input atomic configuration"

    return metadata
