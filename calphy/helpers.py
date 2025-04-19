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
import shutil
import warnings
import logging
import numpy as np

from pylammpsmpi import LammpsLibrary
from lammps import lammps
from ase.io import read, write

import pyscal3.core as pc
from pyscal3.trajectory import Trajectory

class LammpsScript:
    def __init__(self):
        self.script = []

    def command(self, command_str):
        self.script.append(command_str)

    def write(self, infile):
        with open(infile, 'w') as fout:
            for line in self.script:
                fout.write(f'{line}\n')

def create_object(cores, directory, timestep, cmdargs="", 
    init_commands=(), script_mode=False):
    """
    Create LAMMPS object

    Parameters
    ----------
    cores : int
        number of cores

    directory: string
        location of the work directory

    timestep: float
        timestep for the simulation

    Returns
    -------
    lmp : LammpsLibrary object
    """
    if script_mode:
        lmp = LammpsScript()
    else:
        if cmdargs == "":
            cmdargs = None
        lmp = LammpsLibrary(
            cores=cores, working_directory=directory, cmdargs=cmdargs
        )

    commands = [["units", "metal"], 
    ["boundary", "p p p"],
    ["atom_style", "atomic"],
    ["timestep", str(timestep)],
    ["box", "tilt large"]]

    if len(init_commands) > 0:
        #we need to replace some initial commands
        for rc in init_commands:
            #split the command
            raw = rc.split()
            for x in range(len(commands)):
                if raw[0] == commands[x][0]:
                    #we found a matching command
                    commands[x] = [rc]
                    break
            else:
                #its a new command, add it to the list
                commands.append([rc])

    for command in commands:
        lmp.command(" ".join(command))

    return lmp


def create_structure(lmp, calc):
    """
    Create structure using LAMMPS

    Parameters
    ----------
    lmp: LammpsLibrary object

    calc: dict
        calculation dict with the necessary input

    Returns
    -------
    lmp : LammpsLibrary object
    """
    lmp.command("read_data      %s" % calc.lattice)
    return lmp


def set_mass(lmp, options):
    if options.mode == 'composition_scaling':
        lmp.command(f'mass * {options.mass[-1]}')

    else:
        for i in range(options.n_elements):
            lmp.command(f'mass {i+1} {options.mass[i]}')
    return lmp


def set_potential(lmp, options):
    """
    Set the interatomic potential

    Parameters
    ----------
    lmp : LammpsLibrary object

    options : dict

    Returns
    -------
    lmp : LammpsLibrary object
    """
    #lmp.pair_style(options.pair_style_with_options[0])
    #lmp.pair_coeff(options.pair_coeff[0])
    lmp.command(f'pair_style {options._pair_style_with_options[0]}')
    lmp.command(f'pair_coeff {options.pair_coeff[0]}')

    lmp = set_mass(lmp, options)

    return lmp

def read_data(lmp, file):
    lmp.command(f"read_data {file}")
    return lmp

def get_structures(file, species, index=None):
    traj = Trajectory(file)
    if index is None:
        aseobjs = traj[:].to_ase(species=species)
    else:
        aseobjs = traj[index].to_ase(species=species)
    return aseobjs

def remap_box(lmp, x, y, z):
    lmp.command("run 0")
    lmp.command(
        "change_box     all x final 0.0 %f y final 0.0 %f z final 0.0 %f remap units box"
        % (x, y, z)
    )
    return lmp


def compute_msd(lmp, options):
    elements = options.element
    str1 = "fix  4 all ave/time %d %d %d " % (
        int(options.md.n_every_steps),
        int(options.md.n_repeat_steps),
        int(options.md.n_every_steps * options.md.n_repeat_steps),
    )

    # set groups
    for i in range(len(elements)):
        lmp.command("group  g%d type %d" % (i + 1, i + 1))

    str2 = []
    for i in range(len(elements)):
        lmp.command("compute          c%d g%d msd com yes" % (i + 1, i + 1))
        lmp.command("variable         msd%d equal c_c%d[4]" % (i + 1, i + 1))
        str2.append("v_msd%d" % (i + 1))
    str2.append("file")
    str2.append("msd.dat")
    str2 = " ".join(str2)
    command = str1 + str2
    lmp.command(command)
    return lmp


"""
PYSCAL helper routines
---------------------------------------------------------------------
"""


def find_solid_fraction(file):
    sys = pc.System(file)
    try:
        sys.find.neighbors(method="cutoff", cutoff=0)
    except RuntimeError:
        sys.find.neighbors(
            method="cutoff", cutoff=5.0
        )  # Maybe add value as convergence param?
    sys.find.solids(cluster=False)
    solids = np.sum(sys.atoms.solid)
    return solids


def write_data(lmp, file):
    lmp.command(f"write_data {file}")
    return lmp

def prepare_log(file, screen=False):
    logger = logging.getLogger(__name__)
    handler = logging.FileHandler(file)
    formatter = logging.Formatter("%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    if screen:
        scr = logging.StreamHandler()
        scr.setLevel(logging.INFO)
        scr.setFormatter(formatter)
        logger.addHandler(scr)
    return logger

def check_if_any_is_none(data):
    """
    Check if any elements of a list is None, if so return True
    """
    if not isinstance(data, list):
        data = [data]

    for d in data:
        if d is None:
            return True

    return False


def check_if_any_is_not_none(data):
    """
    Check if any element is not None
    """
    if not isinstance(data, list):
        data = [data]

    for d in data:
        if d is not None:
            return True

    return False


def replace_nones(data, replace_data, logger=None):
    """
    Replace Nones in the given array
    """
    if not len(data) == len(replace_data):
        raise ValueError("both arrays must have same length")

    for count, d in enumerate(data):
        if d is None:
            data[count] = replace_data[count]
            if logger is not None:
                logger.info(
                    "Replacing input spring constant None with %f" % replace_data[count]
                )

    return data


def validate_spring_constants(data, klo=0.0001, khi=1000.0, logger=None):
    """
    Validate spring constants and replace them if needed
    """
    # first find a sane value
    sane_k = 0.1
    found = False

    for d in data:
        if klo <= d <= khi:
            sane_k = d
            found = True
            break

    if not found:
        raise ValueError("No spring constant values are between %f and %f" % (klo, khi))

    for count, d in enumerate(data):
        if not (klo <= d <= khi):
            data[count] = sane_k
            if logger is not None:
                logger.info(
                    "Replace insane k %s for element %d with %f"
                    % (str(d), count, sane_k)
                )

    return data
