"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^2, Ralf Drautz^2
^1: Max Planck Institut für Eisenforschung, Dusseldorf, Germany 
^2: Ruhr-University Bochum, Bochum, Germany

calphy is published and distributed under the Academic Software License v1.0 (ASL). 
calphy is distributed in the hope that it will be useful for non-commercial academic research, 
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the ASL for more details. 

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

class InputTemplate:
    def __init__(self):
        pass
    
    def from_dict(self, indict, keys=None):
        if keys is None:
            for key, val in indict.items():
                setattr(self, key, val) 
        elif isinstance(keys, list):
            for key, val in indict.items():
                if key in keys:
                    setattr(self, key, val)
                    
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
        self._pair_style = None
        self._pair_coeff = None
        self._reference_phase = None
        self._lattice_constant = 0
        self._repeat = [1, 1, 1]
        self._npt = True
        self._n_equilibration_steps = 25000
        self._n_switching_steps = 50000
        self._n_sweep_steps = 50000
        self._n_print_steps = 0
        self._n_iterations = 1
        
        self.md = InputTemplate()
        self.md.timestep = 0.001
        self.md.n_small_steps = 10000
        self.md.n_every_steps = 10
        self.md.n_repeat_steps = 10
        self.md.n_cycles = 100
        self.md.thermostat_damping = 0.1
        self.md.barostat_damping = 0.1

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
    
    def __repr__(self):
        """
        String of the class
        """
        data = "%s system with T=%d, P=%d in %s lattice for mode %s"%(self.reference_phase,
            int(self._temperature), int(self._pressure), self.lattice, self.mode) 
        return data
    
    def apply_check(self, key1, key2):
        key2 = self.check_and_convert_to_list(key2)
        for key in key2:
            val_ref = getattr(self, key1)
            val_tar = getattr(self, key)
            if val_ref is not None:
                if val_tar is not None:
                    if not len(val_ref) == len(val_tar):
                        raise ValueError("Length of %s and %s should be same"%(val_ref, val_tar))
            
    @property
    def element(self):
        return self._element

    @element.setter
    def element(self, val):
        val = self.check_and_convert_to_list(val)
        self._nelements = len(val)
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
        self._pressure_input = val
        if val is None:
            self._fix_lattice = True
        elif isinstance(val, list):
            if len(val) == 3:
                if (val[0]==val[1]==val[2]):
                    self._pressure = val[0]
                    self._iso = False
                else:
                    raise NotImplementedError()
            elif len(value) == 2:
                self._pressure = val[0]
                self._pressure_stop = val[1]
            else:
                raise NotImplementedError()
        else:
            self._pressure = val
            self._iso = True
            self._fix_lattice = False

    @property
    def temperature(self):
        return self._temperature_input

    @temperature.setter
    def temperature(self, val):
        self._temperature_input = val
        if isinstance(val, list):
            if len(val) == 2:
                self._temperature = val[0]
                self._temperature_stop = val[1]
                self._temperature_high = 2*val[1]
            else:
                raise ValueError("Temperature cannot have len>2")
        else:
            self._temperature = val
            self._temperature_stop = val
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
        self._pair_style = val

    @property
    def pair_coeff(self):
        return self._pair_coeff
    
    @pair_coeff.setter
    def pair_coeff(self, val):
        val = self.check_and_convert_to_list(val)
        self._pair_coeff = val

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
    
    def check_and_convert_to_list(self, data):
        """
        Check if the given item is a list, if not convert to a single item list
        """
        if not isinstance(data, list):
            return [data]
        else:
            return data
    
    @staticmethod
    def convert_to_list(data):
        """
        Check if the given item is a list, if not convert to a single item list
        """
        if not isinstance(data, list):
            return [data]
        else:
            return data

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
                fixedpots.append(pcraw)
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
            ts = int(self._calc.temperature)
            ps = int(self._calc.pressure)
            l = self.lattice
            l = l.split('/')
            l = l[-1]
        identistring = "-".join([prefix, l, str(ts), str(ps)])
        return identistring

    def create_folders(prefix=None):
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
                calc.md.from_dict(indata["md"])
            if "queue" in indata.keys():
                calc.queue.from_dict(indata["queue"])
            if "tolerance" in indata.keys():
                calc.tolerance.from_dict(indata["tolerance"])
            if "melting_temperature" in indata.keys():
                calc.melting_temperature.from_dict(indata["melting_temperature"])
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
            calc.from_dict(ci, keys=["mode", "pair_style", "pair_coeff", "repeat", "n_equilibration_steps",
                                "n_switching_steps", "n_print_steps", "n_iterations"])
            calc.pressure = Calculation.convert_to_list(ci["pressure"]) if "pressure" in ci.keys() else []
            calc.temperature = Calculation.convert_to_list(ci["temperature"]) if "temperature" in ci.keys() else []
            calc.lattice = Calculation.convert_to_list(ci["lattice"]) if "lattice" in ci.keys() else None
            calc.reference_phase = Calculation.convert_to_list(ci["reference_phase"]) if "reference_phase" in ci.keys() else None
            calc.lattice_constant = Calculation.convert_to_list(ci["lattice_constant"]) if "lattice_constant" in ci.keys() else 0 
            calculations.append(calc)
        else:
            pressure = Calculation.convert_to_list(ci["pressure"]) if "pressure" in ci.keys() else []
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
            if (mode == "fe") or (mode == "alchemy"):
                combos = itertools.product(lat_props, pressure, temperature)
            elif mode == "ts" or mode == "tscale":
                if not len(temperature) == 2:
                    raise ValueError("ts/tscale mode needs 2 temperature values")
                temperature = [temperature for x in range(len(lattice))]
                combos = itertools.product(lat_props, pressure, temperature)
            elif mode == "pscale":
                if not len(pressure) == 2:
                    raise ValueError("pscale mode needs 2 pressure values")
                pressure = [pressure for x in range(len(lattice))]
                combos = itertools.product(lat_props, pressure, temperature)
            #create calculations
            for combo in combos:
                calc = Calculation.generate(indata)
                calc.from_dict(ci, keys=["mode", "pair_style", "pair_coeff", "repeat", "n_equilibration_steps",
                                "n_switching_steps", "n_print_steps", "n_iterations"])
                calc.lattice = combo[0]["lattice"]
                calc.lattice_constant = combo[0]["lattice_constant"]
                calc.reference_phase = combo[0]["reference_phase"]
                calc.pressure = combo[1]
                calc.temperature = combo[2]
                calculations.append(calc)
    return calculations