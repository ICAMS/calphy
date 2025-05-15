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
from typing import Any, Callable, List, ClassVar, Optional, Union
from pydantic import (
    BaseModel,
    Field,
    ValidationError,
    model_validator,
    conlist,
    PrivateAttr,
)
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

__version__ = "1.4.2"


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


class UFMP(BaseModel, title="UFM potential input options"):
    p: Annotated[float, Field(default=50.0)]
    sigma: Annotated[float, Field(default=1.5)]


class MonteCarlo(
    BaseModel, title="Options for Monte Carlo moves during particle swap moves"
):
    n_steps: Annotated[
        int, Field(default=1, gt=0, description="perform swap moves every n_step")
    ]
    n_swaps: Annotated[
        int, Field(default=0, ge=0, description="number of swap moves to perform")
    ]
    reverse_swap: Annotated[bool, Field(default=True)]
    swap_types: Annotated[
        conlist(int, min_length=2, max_length=2),
        Field(default=[1, 2], description="which atoms to swap between"),
    ]


class CompositionScaling(BaseModel, title="Composition scaling input options"):
    _input_chemical_composition: PrivateAttr(default=None)
    output_chemical_composition: Annotated[dict, Field(default={}, required=False)]
    restrictions: Annotated[
        List[str], BeforeValidator(to_list), Field(default=[], required=False)
    ]


class MD(BaseModel, title="MD specific input options"):
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
    thermostat_damping: Annotated[
        Union[float, conlist(float, min_length=2, max_length=2)],
        Field(default=0.1, gt=0),
    ]
    barostat_damping: Annotated[
        Union[float, conlist(float, min_length=2, max_length=2)],
        Field(default=0.1, gt=0),
    ]
    cmdargs: Annotated[str, Field(default="")]
    init_commands: Annotated[List, Field(default=[])]


class NoseHoover(BaseModel, title="Specific input options for Nose-Hoover thermostat"):
    thermostat_damping: Annotated[float, Field(default=0.1, gt=0)]
    barostat_damping: Annotated[float, Field(default=0.1, gt=0)]


class Berendsen(BaseModel, title="Specific input options for Berendsen thermostat"):
    thermostat_damping: Annotated[float, Field(default=100.0, gt=0)]
    barostat_damping: Annotated[float, Field(default=100.0, gt=0)]


class Queue(BaseModel, title="Options for configuring queue"):
    scheduler: Annotated[str, Field(default="local")]
    cores: Annotated[int, Field(default=1, gt=0)]
    jobname: Annotated[str, Field(default="calphy")]
    walltime: Annotated[str, Field(default="23:59:00")]
    queuename: Annotated[str, Field(default="")]
    memory: Annotated[str, Field(default="3GB")]
    commands: Annotated[List, Field(default=[])]
    options: Annotated[List, Field(default=[])]
    modules: Annotated[List, Field(default=[])]


class Tolerance(BaseModel, title="Tolerance settings for convergence"):
    lattice_constant: Annotated[float, Field(default=0.0002, ge=0)]
    spring_constant: Annotated[float, Field(default=0.1, gt=0)]
    solid_fraction: Annotated[float, Field(default=0.7, ge=0)]
    liquid_fraction: Annotated[float, Field(default=0.05, ge=0)]
    pressure: Annotated[float, Field(default=0.5, ge=0)]


class MeltingTemperature(BaseModel, title="Input options for melting temperature mode"):
    guess: Annotated[Union[float, None], Field(default=None, gt=0)]
    step: Annotated[int, Field(default=200, ge=20)]
    attempts: Annotated[int, Field(default=5, ge=1)]


class Calculation(BaseModel, title="Main input class"):
    monte_carlo: Optional[MonteCarlo] = MonteCarlo()
    composition_scaling: Optional[CompositionScaling] = CompositionScaling()
    md: Optional[MD] = MD()
    nose_hoover: Optional[NoseHoover] = NoseHoover()
    berendsen: Optional[Berendsen] = Berendsen()
    queue: Optional[Queue] = Queue()
    tolerance: Optional[Tolerance] = Tolerance()
    uhlenbeck_ford_model: Optional[UFMP] = UFMP()
    melting_temperature: Optional[MeltingTemperature] = MeltingTemperature()
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
            conlist(float, min_length=1, max_length=2),
            conlist(
                conlist(float, min_length=3, max_length=3), min_length=1, max_length=2
            ),
        ],
        Field(default=0),
    ]

    _pressure_stop: float = PrivateAttr(default=None)
    _pressure: float = PrivateAttr(default=None)
    _pressure_stop: float = PrivateAttr(default=None)
    _pressure_input: Any = PrivateAttr(default=None)
    _iso: bool = PrivateAttr(default=False)
    _fix_lattice: bool = PrivateAttr(default=False)

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
    potential_file: Annotated[Union[str, None], Field(default=None)]
    fix_potential_path: Annotated[bool, Field(default=True)]
    _pair_style_with_options: List[str] = PrivateAttr(default=None)

    reference_phase: Annotated[str, Field(default="")]
    lattice_constant: Annotated[float, Field(default=0)]
    repeat: Annotated[
        conlist(int, min_length=3, max_length=3), Field(default=[1, 1, 1])
    ]

    script_mode: Annotated[bool, Field(default=False)]
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

    @model_validator(mode="after")
    def _validate_all(self) -> "Input":

        if not (len(self.element) == len(self.mass)):
            raise ValueError("mass and elements should have same length")

        self.n_elements = len(self.element)

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

        self._temperature_input = copy.copy(self.temperature)
        # guess a melting temp of the system, this will be mostly ignored
        # chem = mendeleev.element(self.element[0])
        # self._melting_temperature = chem.melting_point
        try:
            chem = mendeleev.element(self.element[0])
            self._melting_temperature = chem.melting_point
        except:
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
                    f"Input and output number of atoms are not conserved! Input {self.dict_to_string(self.input_chemical_composition)}, output {self.dict_to_string(self.output_chemical_composition)}, total atoms in structure {structure.natoms}"
                )
        return self

    def fix_paths(self, potlist):
        """
        Fix paths for potential files to complete ones
        """
        fixedpots = []
        for pot in potlist:
            pcraw = pot.split()
            if len(pcraw) >= 3:
                filename = pcraw[2]
                filename = os.path.abspath(filename)
                if os.path.exists(filename):
                    pcnew = " ".join([*pcraw[:2], filename, *pcraw[3:]])
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
        ts = int(self._temperature)
        if self._pressure is None:
            ps = "None"
        else:
            ps = "%d" % (int(self._pressure))
        l = self._original_lattice
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

    @property
    def savefile(self):
        simfolder = self.get_folder_name()
        return os.path.join(simfolder, "job.npy")


def save_job(job):
    filename = os.path.join(job.simfolder, "job.npy")
    np.save(filename, job)


def load_job(filename):
    job = np.load(filename, allow_pickle=True).flatten()[0]
    return job


def read_inputfile(file):
    if not os.path.exists(file):
        raise FileNotFoundError(f"Input file {file} not found.")

    with open(file, "r") as fin:
        data = yaml.safe_load(fin)

    if "element" in data.keys():
        # old format
        outfile = _convert_legacy_inputfile(file)
    else:
        outfile = file
    calculations = _read_inputfile(outfile)
    return calculations


def _read_inputfile(file):
    with open(file, "r") as fin:
        data = yaml.safe_load(fin)
    calculations = []
    for count, calc in enumerate(tqdm(data["calculations"])):
        calc["kernel"] = count
        calc["inputfile"] = file
        if "pressure" in calc.keys():
            calc["pressure"] = _to_none(calc["pressure"])
        calculations.append(Calculation(**calc))
    return calculations


def _convert_legacy_inputfile(file, return_calcs=False):
    with open(file, "r") as fin:
        data = yaml.safe_load(fin)
    if not "element" in data.keys():
        # new format
        raise ValueError("Not old format, exiting..")

    # prepare combos
    calculations = []
    for cc, ci in enumerate(data["calculations"]):
        mode = ci["mode"]
        if mode == "melting_temperature":
            calc = {}
            for key in [
                "md",
                "queue",
                "tolerance",
                "melting_temperature",
                "nose_hoover",
                "berendsen",
                "composition_scaling",
                "temperature_high",
            ]:
                if key in data.keys():
                    calc[key] = copy.copy(data[key])
            for key in [
                "element",
                "mass",
                "script_mode",
                "lammps_executable",
                "mpi_executable",
            ]:
                if key in data.keys():
                    calc[key] = data[key]
            for key in [
                "mode",
                "pair_style",
                "pair_coeff",
                "pair_style_options",
                "npt",
                "repeat",
                "n_equilibration_steps",
                "n_switching_steps",
                "n_print_steps",
                "n_iterations",
                "potential_file",
                "spring_constants",
                "melting_cycle",
                "equilibration_control",
                "folder_prefix",
                "temperature_high",
            ]:
                if key in ci.keys():
                    calc[key] = ci[key]

            # calc['pressure'] = float(np.atleast_1d(ci["pressure"]) if "pressure" in ci.keys() else np.atleast_1d(0))
            # calc['temperature'] = float(np.atleast_1d(ci["temperature"]) if "temperature" in ci.keys() else np.atleast_1d(0))
            # calc['lattice'] = str(ci["lattice"]) if "lattice" in ci.keys() else 'none'
            # calc['reference_phase'] = str(ci["reference_phase"]) if "reference_phase" in ci.keys() else 'none'
            # calc['lattice_constant'] = float(ci["lattice_constant"]) if "lattice_constant" in ci.keys() else 0
            calc["kernel"] = cc
            calc["inputfile"] = file
            calculations.append(calc)

        else:
            pressure = np.atleast_1d(ci["pressure"])
            temperature = np.atleast_1d(ci["temperature"])
            lattice = np.atleast_1d(ci["lattice"])
            reference_phase = np.atleast_1d(ci["reference_phase"])
            if "lattice_constant" in ci.keys():
                lattice_constant = np.atleast_1d(ci["lattice_constant"])
            else:
                lattice_constant = [0 for x in range(len(lattice))]

            lat_props = [
                {
                    "lattice": lattice[x],
                    "lattice_constant": lattice_constant[x],
                    "reference_phase": reference_phase[x],
                }
                for x in range(len(lattice))
            ]

            if (mode == "fe") or (mode == "alchemy") or (mode == "composition_scaling"):
                combos = itertools.product(lat_props, pressure, temperature)
            elif mode == "ts" or mode == "tscale" or mode == "mts":
                if not len(temperature) == 2:
                    raise ValueError("ts/tscale mode needs 2 temperature values")
                temperature = [temperature]
                combos = itertools.product(lat_props, pressure, temperature)
            elif mode == "pscale":
                if not len(pressure) == 2:
                    raise ValueError("pscale mode needs 2 pressure values")
                pressure = [pressure]
                combos = itertools.product(lat_props, pressure, temperature)

            cc = 0
            for combo in combos:
                calc = {}
                for key in [
                    "md",
                    "queue",
                    "tolerance",
                    "melting_temperature",
                    "nose_hoover",
                    "berendsen",
                    "composition_scaling",
                    "temperature_high",
                ]:
                    if key in data.keys():
                        calc[key] = copy.copy(data[key])
                for key in [
                    "element",
                    "mass",
                    "script_mode",
                    "lammps_executable",
                    "mpi_executable",
                ]:
                    if key in data.keys():
                        calc[key] = data[key]
                for key in [
                    "mode",
                    "pair_style",
                    "pair_coeff",
                    "pair_style_options",
                    "npt",
                    "repeat",
                    "n_equilibration_steps",
                    "n_switching_steps",
                    "n_print_steps",
                    "n_iterations",
                    "potential_file",
                    "spring_constants",
                    "melting_cycle",
                    "equilibration_control",
                    "folder_prefix",
                    "temperature_high",
                ]:
                    if key in ci.keys():
                        calc[key] = ci[key]
                # print(combo)
                calc["lattice"] = str(combo[0]["lattice"])
                calc["lattice_constant"] = float(combo[0]["lattice_constant"])
                calc["reference_phase"] = str(combo[0]["reference_phase"])
                calc["pressure"] = _to_float(combo[1])
                calc["temperature"] = _to_float(combo[2])
                calc["kernel"] = cc
                calc["inputfile"] = file
                cc += 1
                calculations.append(calc)

    if return_calcs:
        return calculations
    else:
        newdata = {}
        newdata["calculations"] = calculations
        # print(newdata)
        outfile = "new." + file
        warnings.warn(
            f"Old style input file calphy < v2 found. Converted input in {outfile}. Please check!"
        )
        with open(outfile, "w") as fout:
            yaml.safe_dump(newdata, fout)
        return outfile


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
