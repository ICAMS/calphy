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
from typing import Any, Callable, List
from pydantic import BaseModel, Field, ValidationError, model_validator, conlist
from pydantic.functional_validators import AfterValidator, BeforeValidator
from annotated_types import Len

import yaml
import numpy as np
import copy
import datetime
import itertools
import os
import warnings

def read_report(folder):
    """
    Read the finished calculation report
    """
    repfile = os.path.join(folder, "report.yaml")
    if not os.path.exists(repfile):
        raise FileNotFoundError(f"file {repfile} not found")

    with open(repfile, 'r') as fin:
        data = yaml.safe_load(fin)
    return data


def to_list(v: Any) -> List[Any]:
    return np.atleast_1d(v)

class CompositionScaling(BaseModel, title='Composition scaling input options'):
    input_chemical_composition: Annotated[dict, Field(default=None)]
    output_chemical_composition: Annotated[dict, Field(default=None)]
    restrictions: Annotated[List[str], BeforeValidator(to_list),
                                            Field(default=None)]

class MD(BaseModel, title='MD specific input options'):
    timestep: Annotated[float, Field(default=0.001, 
                                     description='timestep for md simulation', 
                                     example='timestep: 0.001'),]
    n_small_steps: Annotated[int, Field(default=10000)]
    n_every_steps: Annotated[int, Field(default=10)]
    n_repeat_steps: Annotated[int, Field(default=10)]
    n_cycles: Annotated[int, Field(default=100)]
    thermostat_damping: Annotated[str, Field(default=0.1)]
    barostat_damping: Annotated[str, Field(default=0.1)]
    cmdargs: Annotated[str, Field(default=None)]
    init_commands: Annotated[str, Field(default=None)]


class NoseHoover(BaseModel, title='MD specific input options'):
    thermostat_damping: Annotated[float, Field(default=0.1)]
    thermostat_damping: Annotated[float, Field(default=0.1)]


class Input(BaseModel, title='Main input class'):
    a_list: Annotated[List[int], BeforeValidator(to_list),
                                      Field(default=None,
                                      description='A test list',
                                      repr=True)]
    element: Annotated[List[str], BeforeValidator(to_list),
                                      Field(default=None)]
    n_elements: Annotated[int, Field(default=None)]
    mass: Annotated[List[int], BeforeValidator(to_list),
                                      Field(default=None)]
    mode: Annotated[str, Field(default=None)]
    lattice: Annotated[str, Field(default=None)]
    pressure: Annotated[ None | float | conlist(float, min_length=1, max_length=3) | conlist(conlist(float, min_length=3, max_length=3), min_length=1, max_length=2) , 
                                      Field(default=None)]
    #pressure_stop: = 0
    #pressure_input: = None
    #temperature: None
    #temperature_stop: None
    #temperature_high: None
    #temperature_input: None
    
    #self._iso = False
    #self._fix_lattice = False
    #self._melting_cycle = True
    #self._pair_style = None
    #self._pair_style_options = None
    #self._pair_coeff = None
    #self._potential_file = None
    #self._fix_potential_path = True
    #self._reference_phase = None
    #self._lattice_constant = 0
    #self._repeat = [1, 1, 1]
    #self._script_mode = False
    #self._lammps_executable = None
    #self._mpi_executable = None
    #self._npt = True
    #self._n_equilibration_steps = 25000
    #self._n_switching_steps = 50000
    #self._n_sweep_steps = 50000
    #self._n_print_steps = 0
    #self._n_iterations = 1
    #self._equilibration_control = None
    #self._folder_prefix = None

    #add second level options; for example spring constants
    #self._spring_constants = None
    #self._ghost_element_count = 0

    @model_validator(mode='after')
    def _validate_lengths(self) -> 'Input':
        if not (len(self.element) == len(self.mass)):
            raise ValueError('mass and elements should have same length')
        return self

    @model_validator(mode='after')
    def _validate_nelements(self) -> 'Input':
        self.n_elements = len(self.element)
        return self


    @model_validator(mode='after')
    def _validate_pressure(self) -> 'Input':
        self.input_pressure = copy.copy(self.pressure)
        if self.pressure is None:
            self.iso = True
            self.fix_lattice = True
            self.pressure = None
            self.pressure_stop = None
        elif np.isscalar(self.pressure):
            self.pressure_stop = self.pressure
            self.iso = True
            self.fix_lattice = False
        elif np.shape(self.pressure) == (1,):
            print(np.shape(self.pressure))
        else:
            print(np.shape(self.pressure))
        return self


