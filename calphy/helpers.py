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
from pylammpsmpi import LammpsLibrary
import logging
import numpy as np
from lammps import lammps
import calphy.lattice as pl
import pyscal.core as pc


def create_object(cores, directory, timestep):
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
    lmp = LammpsLibrary(mode="local", cores=cores, 
        working_directory=directory)
    
    lmp.units("metal")
    lmp.boundary("p p p")
    lmp.atom_style("atomic")
    lmp.timestep(timestep)
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
    l, alat, apc, conc = pl.prepare_lattice(calc)

    if l == "file":
        lmp.command("read_data      %s"%calc.lattice)
    else:
        lmp.lattice(l, alat)
        lmp.region("box block", 0, calc.repeat[0], 
            0, calc.repeat[1], 
            0, calc.repeat[2])
        lmp.create_box("1 box")
        lmp.create_atoms("1 box")
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
    lmp.pair_style(options.pair_style[0])
    lmp.pair_coeff(options.pair_coeff[0])

    for i in range(options.n_elements):
        lmp.mass(i+1, options.mass[i])
    return lmp


def read_dump(lmp, file, species=1):
    # Read atoms positions, velocities and box parameters.
    lmp.command("lattice          fcc 4.0")
    lmp.command("region           box block 0 2 0 2 0 2")
    lmp.command("create_box       %d box"%species)
    lmp.command("read_dump        %s 0 x y z vx vy vz scaled no box yes add keep"%file)
    return lmp


def set_hybrid_potential(lmp, options, eps):
    pc =  options.pair_coeff[0]
    pcraw = pc.split()
    pcnew = " ".join([*pcraw[:2], *[options.pair_style[0],], *pcraw[2:]])
    
    lmp.command("pair_style       hybrid/overlay %s ufm 7.5"%options.pair_style[0])
    lmp.command("pair_coeff       %s"%pcnew)
    lmp.command("pair_coeff       * * ufm %f 1.5"%eps) 

    for i in range(options.n_elements):
        lmp.mass(i+1, options.mass[i])
    return lmp

def set_double_hybrid_potential(lmp, options, pair_style, pair_coeff):
    
    pc1 =  pair_coeff[0]
    pcraw1 = pc1.split()
    
    pc2 =  pair_coeff[1]
    pcraw2 = pc2.split()

    if pair_style[0] == pair_style[1]:
        pcnew1 = " ".join([*pcraw1[:2], *[pair_style[0],], "1", *pcraw1[2:]])
        pcnew2 = " ".join([*pcraw2[:2], *[pair_style[1],], "2", *pcraw2[2:]])
    else:
        pcnew1 = " ".join([*pcraw1[:2], *[pair_style[0],], *pcraw1[2:]])
        pcnew2 = " ".join([*pcraw2[:2], *[pair_style[1],], *pcraw2[2:]])

    lmp.command("pair_style       hybrid/overlay %s %s"%(pair_style[0], pair_style[1]))
    
    lmp.command("pair_coeff       %s"%pcnew1)
    lmp.command("pair_coeff       %s"%pcnew2) 

    for i in range(options.n_elements):
        lmp.mass(i+1, options.mass[i])
    return lmp

def remap_box(lmp, x, y, z):
    lmp.command("run 0")
    lmp.command("change_box     all x final 0.0 %f y final 0.0 %f z final 0.0 %f remap units box"%(x, y, z))
    return lmp


def compute_msd(lmp, options):
    elements = options.element
    str1 = "fix  4 all ave/time %d %d %d "%(int(options.md.n_every_steps),
                                           int(options.md.n_repeat_steps), 
              int(options.md.n_every_steps*options.md.n_repeat_steps))

    #set groups
    for i in range(len(elements)):
        lmp.command("group  g%d type %d"%(i+1, i+1))

    str2 = []
    for i in range(len(elements)):
        lmp.command("compute          c%d g%d msd com yes"%(i+1, i+1))
        lmp.command("variable         msd%d equal c_c%d[4]"%(i+1, i+1))
        str2.append("v_msd%d"%(i+1))
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
    sys = pc.System()
    sys.read_inputfile(file)
    sys.find_neighbors(method="cutoff", cutoff=0)
    solids = sys.find_solids()
    return solids

def reset_timestep(file, conf):
    sys = pc.System()
    sys.read_inputfile(file, customkeys=["vx", "vy", "vz", "mass"])
    sys.to_file(conf, customkeys=["vx", "vy", "vz", "mass"])


"""
NOrmal helper routines
---------------------------------------------------------------------
"""
def prepare_log(file):
    logger = logging.getLogger(__name__)
    handler = logging.FileHandler(file)
    formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    return logger