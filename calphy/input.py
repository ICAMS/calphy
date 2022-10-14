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

import os
import yaml
import warnings
import copy
import yaml
import itertools
import shutil
import numpy as np

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

class InputTemplate:
    def __init__(self):
        pass
    
    def add_from_dict(self, indict, keys=None):
        if keys is None:
            for key, val in indict.items():
                setattr(self, key, val) 
        elif isinstance(keys, list):
            for key, val in indict.items():
                if key in keys:
                    setattr(self, key, val)

    def from_dict(self, indict):
        for key, val in indict.items():
            if isinstance(val, dict):
                setattr(self, key, InputTemplate())
                getattr(self, key).from_dict(val)
            else:
                setattr(self, key, val)

    def to_dict(self):
        tdict = copy.deepcopy(self.__dict__)
        for key, val in tdict.items():
            if isinstance(val, InputTemplate):
                tdict[key] = val.to_dict()
        return tdict

    def check_and_convert_to_list(self, data, check_none=False):
        """
        Check if the given item is a list, if not convert to a single item list
        """
        if not isinstance(data, list):
            data = [data]

        if check_none:
            data = [None if x=="None" else x for x in data]

        return data
    
    @staticmethod
    def convert_to_list(data, check_none=False):
        """
        Check if the given item is a list, if not convert to a single item list
        """
        if not isinstance(data, list):
            data = [data]

        if check_none:
            data = [None if x=="None" else x for x in data]

        return data

    def merge_dicts(self, dicts):
        """
        Merge dicts from a given list
        """
        merged_dict = {}
        for d in dicts:
            for key, val in d.items():
                merged_dict[key] = val
        return merged_dict

class CompositionScaling(InputTemplate):
    def __init__(self):
        self._input_chemical_composition = None
        self._output_chemical_composition = None
        self._restrictions = None

    @property
    def input_chemical_composition(self):
        return self._input_chemical_composition        

    @input_chemical_composition.setter
    def input_chemical_composition(self, val):
        if isinstance(val, list):
            val = self.merge_dicts(val) 
        self._input_chemical_composition = val         

    @property
    def output_chemical_composition(self):
        return self._output_chemical_composition

    @output_chemical_composition.setter
    def output_chemical_composition(self, val):
        if isinstance(val, list):
            val = self.merge_dicts(val) 
        self._output_chemical_composition = val

    @property
    def restrictions(self):
        return self._restrictions

    @restrictions.setter
    def restrictions(self, val):
        self._restrictions = self.check_and_convert_to_list(val)       


class Calculation(InputTemplate):
    def __init__(self):
        super(InputTemplate, self).__init__()
        self._element = None
        self._n_elements = 0
        self._mass = 1.0
        self._mode = None
        self._lattice = None
        self._pressure = 0
        self._pressure_stop = 0
        self._pressure_input = None
        self._temperature = None
        self._temperature_stop = None
        self._temperature_high = None
        self._temperature_input = None
        self._iso = False
        self._fix_lattice = False
        self._melting_cycle = True
        self._pair_style = None
        self._pair_coeff = None
        self._potential_file = None
        self._fix_potential_path = True
        self._reference_phase = None
        self._lattice_constant = 0
        self._repeat = [1, 1, 1]
        self._npt = True
        self._n_equilibration_steps = 25000
        self._n_switching_steps = 50000
        self._n_sweep_steps = 50000
        self._n_print_steps = 0
        self._n_iterations = 1
        self._equilibration_control = None
        self._folder_prefix = None

        #add second level options; for example spring constants
        self._spring_constants = None
        self._ghost_element_count = 0
        
        self.md = InputTemplate()
        self.md.timestep = 0.001
        self.md.n_small_steps = 10000
        self.md.n_every_steps = 10
        self.md.n_repeat_steps = 10
        self.md.n_cycles = 100
        self.md.thermostat_damping = 0.1
        self.md.barostat_damping = 0.1
        self.md.cmdargs = None

        self.nose_hoover = InputTemplate()
        self.nose_hoover.thermostat_damping = 0.1
        self.nose_hoover.barostat_damping = 0.1

        self.berendsen = InputTemplate()
        self.berendsen.thermostat_damping = 100.0
        self.berendsen.barostat_damping = 100.0

        self.queue = InputTemplate()
        self.queue.scheduler = "local"
        self.queue.cores = 1
        self.queue.jobname = "calphy"
        self.queue.walltime = None
        self.queue.queuename = None
        self.queue.memory = "3GB"
        self.queue.commands = None
        self.queue.options = None
        self.queue.modules = None
        
        self.tolerance = InputTemplate()
        self.tolerance.lattice_constant = 0.0002
        self.tolerance.spring_constant = 0.01
        self.tolerance.solid_fraction = 0.7
        self.tolerance.liquid_fraction = 0.05
        self.tolerance.pressure = 0.5
        
        #specific input options
        self.melting_temperature = InputTemplate()
        self.melting_temperature.guess = None
        self.melting_temperature.step = 200
        self.melting_temperature.attempts = 5

        #new mode for composition trf
        self.composition_scaling = CompositionScaling()
    
    def __repr__(self):
        """
        String of the class
        """
        data = "%s system with T=%s, P=%s in %s lattice for mode %s"%(self.to_string(self.reference_phase),
            self.to_string(self._temperature), self.to_string(self._pressure), self.to_string(self.lattice), self.to_string(self.mode)) 
        return data
    
    def to_string(self, val):
        if val is None:
            return "None"
        else:
            return str(val)

    def to_bool(self, val):        
        if val == "True":
            val = True
        elif val == "False":
            val = False
        elif val == 1:
            val = True
        elif val == 0:
            val = False
        elif val == "1":
            val = True
        elif val == "0":
            val = False
        return val

    @property
    def element(self):
        return self._element

    @element.setter
    def element(self, val):
        val = self.check_and_convert_to_list(val)
        self._n_elements = len(val)
        self._element = val
    
    @property
    def n_elements(self):
        return self._n_elements
    
    @property
    def mass(self):
        return self._mass
    
    @mass.setter
    def mass(self, val):
        val = self.check_and_convert_to_list(val)
        if self.element is not None:
            if not len(self.element) == len(val):
                raise ValueError("Elements and mass must have same length")
        self._mass = val

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, val):
        #add check
        self._mode = val
        if val == "melting_temperature":
            self.temperature = 0
            self.pressure = 0

    @property
    def lattice(self):
        return self._lattice

    @lattice.setter
    def lattice(self, val):
        self._lattice = val

    @property
    def pressure(self):
        return self._pressure_input

    @pressure.setter
    def pressure(self, val):
        """
        Pressure input can be of many diff types:
        Input                      iso           fix_lattice  Mode  add. constraints
        1. None                    True          True         any   No
        2. scalar                  True          False        any   No
        3. list of scalar len 1    True          False        any   No
        4. list of scalar len 2    True          False        ps    No
        5. list of scalar len 3    False         False        any   All terms equal
        6. list of list   len 1x3  False         False        any   Each inner term equal
        7. list of list   len 2x3  False         False        ps    Each inner term equal
        """
        def _check_equal(val):
            if not (val[0]==val[1]==val[2]):
                return False
            return True

        self._pressure_input = val

        error_message = "Available pressure types are of shape: None, 1, 3, (1x1), (2x1), (3x1), (2x3)"
        # Case: 1
        if val is None:
            self._iso = True
            self._fix_lattice = True
            self._pressure = None
            self._pressure_stop = None

        # Cases: 3,4,5,6,7
        elif isinstance(val, list):
            shape = np.shape(val)
            indx = shape[0]
            indy = shape[1] if len(shape)==2 else None

            if indx == 1:
                if indy is None:
                    # Case: 3
                    self._iso = True
                    self._fix_lattice = False
                    self._pressure = val[0]
                    self._pressure_stop = val[0]
                elif indy == 3:
                    # Case: 6
                    if _check_equal(val[0]):
                        self._iso = False
                        self._fix_lattice = False
                        self._pressure = val[0][0]
                        self._pressure_stop = val[0][0]
                    else:
                        raise NotImplementedError("Pressure should have px=py=pz")
                else:
                    raise NotImplementedError(error_message)
            elif indx == 2:
                if indy is None:
                    # Case: 4
                    self._iso = True
                    self._fix_lattice = False
                    self._pressure = val[0]
                    self._pressure_stop = val[1]
                elif indy == 3:
                    # Case: 7
                    self._iso = False
                    self._fix_lattice = False                    
                    if _check_equal(val[0]) and _check_equal(val[1]):
                        self._pressure = val[0][0]
                        self._pressure_stop = val[1][0]
                    else:
                        raise NotImplementedError("Pressure should have px=py=pz")
                else:
                    raise NotImplementedError(error_message)

            elif indx == 3:
                if indy is None:
                    # Case: 5
                    if _check_equal(val):
                        self._iso = False
                        self._fix_lattice = False                    
                        self._pressure = val[0]
                        self._pressure_stop = val[0]                                
                    else:
                        raise NotImplementedError("Pressure should have px=py=pz")
                else:
                    raise NotImplementedError(error_message)
            else:
                raise NotImplementedError()
        
        # Case: 2
        else:
            self._pressure = val
            self._pressure_stop = val
            self._iso = True
            self._fix_lattice = False

    @property
    def temperature(self):
        return self._temperature_input

    @temperature.setter
    def temperature(self, val):
        self._temperature_input = val
        if val is None:
                self._temperature = val
                self._temperature_stop = val
                self._temperature_high = val
        elif isinstance(val, list):
            if len(val) == 2:
                self._temperature = val[0]
                self._temperature_stop = val[1]
                if self._temperature_high is None:
                    self._temperature_high = 2*val[1]
            else:
                raise ValueError("Temperature cannot have len>2")
        else:
            self._temperature = val
            self._temperature_stop = val
            if self._temperature_high is None:
                self._temperature_high = 2*val

    @property
    def temperature_high(self):
        return self._temperature_high

    @temperature_high.setter
    def temperature_high(self, val):
        self._temperature_high = val
    
    @property
    def pair_style(self):
        return self._pair_style
    
    @pair_style.setter
    def pair_style(self, val):
        val = self.check_and_convert_to_list(val)
        #val = self.fix_paths(val)
        self._pair_style = val

    @property
    def pair_coeff(self):
        return self._pair_coeff
    
    @pair_coeff.setter
    def pair_coeff(self, val):
        val = self.check_and_convert_to_list(val)
        if self._fix_potential_path:
            val = self.fix_paths(val)
        self._pair_coeff = val

    @property
    def potential_file(self):
        return self._potential_file

    @potential_file.setter
    def potential_file(self, val):
        if os.path.exists(val):
            self._potential_file = val
        else:
            raise FileNotFoundError("File %s not found"%val)

    @property
    def reference_phase(self):
        return self._reference_phase
    
    @reference_phase.setter
    def reference_phase(self, val):
        self._reference_phase = val

    @property
    def lattice_constant(self):
        return self._lattice_constant
    
    @lattice_constant.setter
    def lattice_constant(self, val):
        self._lattice_constant = val

    @property
    def repeat(self):
        return self._repeat
    
    @repeat.setter
    def repeat(self, val):
        if isinstance(val, list):
            if not len(val) == 3:
                raise ValueError("repeat should be three")
        else:
            val = [val, val, val]
        self._repeat = val
    
    @property
    def npt(self):
        return self._npt
    
    @npt.setter
    def npt(self, val):
        self._npt = val
    
    @property
    def n_equilibration_steps(self):
        return self._n_equilibration_steps
    
    @n_equilibration_steps.setter
    def n_equilibration_steps(self, val):
        self._n_equilibration_steps = val

    @property
    def n_switching_steps(self):
        return self._n_switching_steps
    
    @n_switching_steps.setter
    def n_switching_steps(self, val):
        if isinstance(val, list):
            if len(val) == 2:
                self._n_switching_steps = val[0]
                self._n_sweep_steps = val[1]
            else:
                raise TypeError("n_switching_steps should be len 1 or 2")
        else:
            self._n_switching_steps = val
            self._n_sweep_steps = val

    @property
    def n_print_steps(self):
        return self._n_print_steps
    
    @n_print_steps.setter
    def n_print_steps(self, val):
        self._n_print_steps = val

    @property
    def n_iterations(self):
        return self._n_iterations
    
    @n_iterations.setter
    def n_iterations(self, val):
        self._n_iterations = val

    @property
    def equilibration_control(self):
        return self._equilibration_control

    @equilibration_control.setter
    def equilibration_control(self, val):
        if val not in ["berendsen", "nose-hoover", None]:
            raise ValueError("Equilibration control should be either berendsen or nose-hoover or None")
        self._equilibration_control = val

    @property
    def melting_cycle(self):
        return self._melting_cycle
    
    @melting_cycle.setter
    def melting_cycle(self, val):
        val = self.to_bool(val)
        if isinstance(val, bool):
            self._melting_cycle = val
        else:
            raise TypeError("Melting cycle should be either True/False")

    @property
    def spring_constants(self):
        return self._spring_constants

    @spring_constants.setter
    def spring_constants(self, val):
        val = self.check_and_convert_to_list(val, check_none=True)
        self._spring_constants = val

    @property
    def folder_prefix(self):
        return self._folder_prefix

    @folder_prefix.setter
    def folder_prefix(self, val):
        self._folder_prefix = val

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

    def create_folders(self, prefix=None):
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
        identistring = self.create_identifier()
        if prefix is None:
            simfolder = os.path.join(os.getcwd(), identistring)
        else:
            simfolder = os.path.join(prefix, identistring)

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
    
    @classmethod
    def generate(cls, indata):
        if not isinstance(indata, dict):
            if os.path.exists(indata):
                with open(indata) as file:
                    indata = yaml.load(file, Loader=yaml.FullLoader)            
        if isinstance(indata, dict):
            calc = cls()
            calc.element = indata["element"]
            calc.mass = indata["mass"]
            if "md" in indata.keys():
                calc.md.add_from_dict(indata["md"])
            if "queue" in indata.keys():
                calc.queue.add_from_dict(indata["queue"])
            if "tolerance" in indata.keys():
                calc.tolerance.add_from_dict(indata["tolerance"])
            if "melting_temperature" in indata.keys():
                calc.melting_temperature.add_from_dict(indata["melting_temperature"])
            if "nose_hoover" in indata.keys():
                calc.nose_hoover.add_from_dict(indata["nose_hoover"])
            if "berendsen" in indata.keys():
                calc.berendsen.add_from_dict(indata["berendsen"])
            if "composition_scaling" in indata.keys():
                calc.composition_scaling.add_from_dict(indata["composition_scaling"])
            #if temperature_high is present, set it
            if "temperature_high" in indata.keys():
                calc.temperature_high = indata["temperature_high"]
            return calc
        else:
            raise FileNotFoundError('%s input file not found'% indata)
    
    @staticmethod
    def read_file(file):
        if os.path.exists(file):
            with open(file) as file:
                indata = yaml.load(file, Loader=yaml.FullLoader)
        else:
            raise FileNotFoundError('%s input file not found'% file)
        return indata
    
def read_inputfile(file):
    """
    Read calphy inputfile
    
    Parameters
    ----------
    file : string
        input file
    
    Returns
    -------
    options : dict
        dictionary containing input options
    """
    indata = Calculation.read_file(file)
    calculations = []
    for ci in indata["calculations"]:
        #get main variables
        mode = ci["mode"]
        if mode == "melting_temperature":
            calc = Calculation.generate(indata)
            calc.add_from_dict(ci, keys=["mode", "pair_style", "pair_coeff", "repeat", "n_equilibration_steps",
                                "n_switching_steps", "n_print_steps", "n_iterations", "spring_constants", "equilibration_control",
                                "folder_prefix"])
            calc.pressure = Calculation.convert_to_list(ci["pressure"], check_none=True) if "pressure" in ci.keys() else 0
            calc.temperature = Calculation.convert_to_list(ci["temperature"]) if "temperature" in ci.keys() else None
            calc.lattice = Calculation.convert_to_list(ci["lattice"]) if "lattice" in ci.keys() else None
            calc.reference_phase = Calculation.convert_to_list(ci["reference_phase"]) if "reference_phase" in ci.keys() else None
            calc.lattice_constant = Calculation.convert_to_list(ci["lattice_constant"]) if "lattice_constant" in ci.keys() else 0 
            calculations.append(calc)
        else:
            pressure = Calculation.convert_to_list(ci["pressure"], check_none=True) if "pressure" in ci.keys() else []
            temperature = Calculation.convert_to_list(ci["temperature"]) if "temperature" in ci.keys() else []
            lattice = Calculation.convert_to_list(ci["lattice"])
            reference_phase = Calculation.convert_to_list(ci["reference_phase"])
            if "lattice_constant" in ci.keys():
                lattice_constant = Calculation.convert_to_list(ci["lattice_constant"])
            else:
                lattice_constant = [0 for x in range(len(lattice))]
            if not len(lattice_constant)==len(reference_phase)==len(lattice):
                raise ValueError("lattice constant, reference phase and lattice should have same length")
            lat_props = [{"lattice": lattice[x], "lattice_constant":lattice_constant[x], "reference_phase":reference_phase[x]} for x in range(len(lattice))]
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
            #create calculations
            for combo in combos:
                calc = Calculation.generate(indata)
                calc.add_from_dict(ci, keys=["mode", "pair_style", "pair_coeff", "repeat", "n_equilibration_steps",
                                "n_switching_steps", "n_print_steps", "n_iterations", "potential_file", "spring_constants",
                                "melting_cycle", "equilibration_control", "folder_prefix", "temperature_high"])
                calc.lattice = combo[0]["lattice"]
                calc.lattice_constant = combo[0]["lattice_constant"]
                calc.reference_phase = combo[0]["reference_phase"]
                calc.pressure = combo[1]
                calc.temperature = combo[2]
                calculations.append(calc)
    return calculations
