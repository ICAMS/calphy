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
from pydantic import BaseModel, Field, ValidationError, model_validator, conlist, ClassVar
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

def _check_equal(val):
    if not (val[0]==val[1]==val[2]):
        return False
    return True


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
    n_small_steps: Annotated[int, Field(default=10000, gt=0)]
    n_every_steps: Annotated[int, Field(default=10, gt=0)]
    n_repeat_steps: Annotated[int, Field(default=10, gt=0)]
    n_cycles: Annotated[int, Field(default=100, gt=0)]
    thermostat_damping: Annotated[str, Field(default=0.1, gt=0)]
    barostat_damping: Annotated[str, Field(default=0.1, gt=0)]
    cmdargs: Annotated[str, Field(default=None)]
    init_commands: Annotated[str, Field(default=None)]


class NoseHoover(BaseModel, title='Specific input options for Nose-Hoover thermostat'):
    thermostat_damping: Annotated[float, Field(default=0.1, gt=0)]
    thermostat_damping: Annotated[float, Field(default=0.1, gt=0)]

class Berendsen(BaseModel, title='Specific input options for Berendsen thermostat'):
    thermostat_damping: Annotated[float, Field(default=100.0, gt=0)]
    thermostat_damping: Annotated[float, Field(default=100.0, gt=0)]

class Queue(BaseModel, title='Options for configuring queue'):
    scheduler: Annotated[str, Field(default='local')]
    cores: Annotated[int, Field(default=1, gt=0)]
    jobname: Annotated[str, Field(default='calphy')]
    walltime: Annotated[str, Field(default=None)]
    queuename: Annotated[str, Field(default=None)]
    memory: Annotated[str, Field(default="3GB")]
    commands: Annotated[List[str], Field(default=None)]
    options: Annotated[List[str], Field(default=None)]
    modules: Annotated[List[str], Field(default=None)]

class Tolerance(BaseModel, title='Tolerance settings for convergence'):
    lattice_constant: Annotated[float, Field(default=0.0002, ge=0)]
    spring_constant: Annotated[float, Field(default=0.1, gt=0)]
    solid_fraction: Annotated[float, Field(default=0.7, ge=0)]
    liquid_fraction: Annotated[float, Field(default=0.05, ge=0)]
    pressure: Annotated[float, Field(default=0.5, ge=0)]

class MeltingTemperature(BaseModel, title='Input options for melting temperature mode'):
    guess: Annotated[float, Field(default=None, gt=0)]
    step: Annotated[int, Field(default=200, ge=20)]
    attempts: Annotated[int, Field(default=5, ge=1)]

class Input(BaseModel, title='Main input class'):
    element: Annotated[List[str], BeforeValidator(to_list),
                                      Field(default=None)]
    n_elements: Annotated[int, Field(default=None)]
    mass: Annotated[List[int], BeforeValidator(to_list),
                                      Field(default=None)]
    mode: Annotated[str, Field(default=None)]
    lattice: Annotated[str, Field(default=None)]
    
    #pressure properties
    pressure: Annotated[ None | float | conlist(float, min_length=1, max_length=2) | conlist(conlist(float, min_length=3, max_length=3), min_length=1, max_length=2) , 
                                      Field(default=None)]
    pressure_stop: ClassVar[float] = None
    pressure_input: ClassVar[Any] = None
    iso: ClassVar[bool] = False
    fix_lattice: ClassVar[bool] = False

    temperature: Annotated[ float | conlist(float, min_length=1, max_length=2),
                            Field(default=None)]
    temperature_high: Annotated[float, Field(default=None)]
    temperature_stop: ClassVar[float] = None
    temperature_input: ClassVar[Any] = None
    
    melting_cycle: Annotated[ bool, Field(default=True)]
    
    pair_style: Annotated[ List(str), BeforeValidator(to_list),
                            Field(default=None)]
    pair_coeff: Annotated[ List(str), BeforeValidator(to_list),
                            Field(default=None)]
    potential_file: Annotated[ str, Field(default=None)]
    pair_style_options: ClassVar[List(str)] = None
    fix_potential_path: ClassVar[bool] = True
    
    reference_phase: Annotated[ str, Field(default = None)]
    lattice_constant: Annotated[int, Field(default = 0)]
    repeat: Annotated[conlist(int, min_length=3, max_length=3), 
                            Field(default=[1,1,1])]
    
    script_mode: Annotated[ bool, Field(default = False)]
    lammps_executable: Annotated[ str, Field(default = None)]
    mpi_executable: Annotated[ str, Field(default = None)]
    
    npt: Annotated[ bool, Field(default = True)]
    n_equilibration_steps: Annotated[ int, Field(default = 25000)]
    n_switching_steps: Annotated[ int | conlist(int, min_length=2,
                max_length=2), Field(default = [50000, 50000])]
    n_sweep_steps: ClassVar[int] = None
    n_print_steps: Annotated[int, Field(default = 0)]
    n_iterations: Annotated[int, Field(default = 1)]
    equilibration_control: Annotated[str, Field(default = None)]
    folder_prefix: Annotated[str, Field(default = None)]

    #add second level options; for example spring constants
    spring_constants: Annotated[List(float), Field(default = None)]
    ghost_element_count: ClassVar[int] = 0

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
        self.pressure_input = copy.copy(self.pressure)
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
            self.iso = True
            self.fix_lattice = False
            self.pressure = self.pressure_input[0]
            self.pressure_stop = self.pressure_input[0]
        elif np.shape(self.pressure) == (2,):
            self.iso = True
            self.fix_lattice = False
            self.pressure = self.pressure_input[0]
            self.pressure_stop = self.pressure_input[1]
        elif np.shape(self.pressure) == (1, 3):
            if not _check_equal(self.pressure[0]):
                raise ValueError('All pressure terms must be equal')
            self.iso = False
            self.fix_lattice = False
            self.pressure = self.pressure_input[0][0]
            self.pressure_stop = self.pressure_input[0][0]                
        elif np.shape(self.pressure) == (2, 3):
            if not (_check_equal(self.pressure[0]) and _check_equal(self.pressure[1])):
                raise ValueError('All pressure terms must be equal')
            self.iso = False
            self.fix_lattice = False
            self.pressure = self.pressure_input[0][0]
            self.pressure_stop = self.pressure_input[1][0]                                
        else:
            raise ValueError('Unknown format for pressure')
        return self


    @model_validator(mode='after')
    def _validate_temperature(self) -> 'Input':
        self.temperature_input = copy.copy(self.temperature)
        if self.temperature is None:
            pass
        elif shape(np.atleast_1d(self.temperature)) == (1,):
            temp = self.atleast_1d(self.temperature_input)
            self.temperature = temp[0]
            self.temperature_stop = temp[0]
            self.temperature_high is None:
                self.temperature_high = 2*temp[0]
        elif shape(self.temperature) == (2,):
            temp = self.temperature_input
            self.temperature = temp[0]
            self.temperature_stop = temp[1]
            self.temperature_high is None:
                self.temperature_high = 2*temp[1]
        return self

    @model_validator(mode='after')
    def _validate_pair_style(self) -> 'Input':
        ps_lst = []
        ps_options_lst = []
        for ps in self.pair_style:
            ps_split = ps.split()
            ps_lst.append(ps_split[0])
            if len(ps) > 1:
                ps_options_lst.append(" ".join(ps_split[1:]))
            else:
                ps_options_lst.append("")

        #val = self.fix_paths(val)
        self.pair_style = ps_lst

        #only set if its None
        if self.pair_style_options is None:
            self.pair_style_options = ps_options_lst
        return self

    @model_validator(mode='after')
    def _validate_time(self) -> 'Input':
        if np.isscalar(self.n_switching_steps):
            self.n_sweep_steps = self.n_switching_steps
        else:
            self.n_sweep_steps = self.n_switching_steps[1]
            self.n_switching_steps = self.n_switching_steps[0]
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
                pcnew = " ".join([*pcraw[:2], filename, *pcraw[3:]])
                fixedpots.append(pcnew)
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
        #lattice processed
        prefix = self.mode
        if prefix == 'melting_temperature':
            ts = int(0)
            ps = int(0)
            l = 'tm'
        else:
            ts = int(self._temperature)
            if self._pressure is None:
                ps = "None"
            else:
                ps = "%d"%(int(self._pressure))
            l = self.lattice
            l = l.split('/')
            l = l[-1]
        
        if self.folder_prefix is None:
            identistring = "-".join([prefix, l, str(ts), str(ps)])
        else:
            identistring = "-".join([self.folder_prefix, prefix, l, str(ts), str(ps)])
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

        #if folder exists, delete it -> then create
        try:
            if os.path.exists(simfolder):
                shutil.rmtree(simfolder)
        except OSError:
            newstr = '-'.join(str(datetime.datetime.now()).split())
            newstr = '-'.join([simfolder, newstr])
            shutil.move(simfolder, newstr)
        
        os.mkdir(simfolder)
        return simfolder
        
    @property
    def savefile(self):
        simfolder = self.get_folder_name()
        return os.path.join(simfolder, 'job.npy')

def save_job(job):
    filename = os.path.join(job.simfolder, 'job.npy')
    np.save(filename, job)

def load_job(filename):
    job = np.load(filename, allow_pickle=True).flatten()[0]
    return job

def check_dict(indict, key, retval=None):
    if key in indict.items():
        return indict[key]
    else:
        return retval
