"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^1, Ralf Drautz^1
^1: Ruhr-University Bochum, Bochum, Germany

More information about the program can be found in:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz. 
“Automated Free Energy Calculation from Atomistic Simulations.” 
ArXiv:2107.08980 [Cond-Mat], July 19, 2021. 
http://arxiv.org/abs/2107.08980.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

See the LICENSE file.

For more information contact:
sarath.menon@ruhr-uni-bochum.de
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
        lmp.command("read_data      %s"%calc["lattice"])
    else:
        lmp.lattice(l, alat)
        lmp.region("box block", 0, calc["repeat"][0], 
            0, calc["repeat"][1], 
            0, calc["repeat"][2])
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
    lmp.pair_style(options["md"]["pair_style"])
    lmp.pair_coeff(options["md"]["pair_coeff"])

    for i in range(options["nelements"]):
        lmp.mass(i+1, options["mass"][i])
    return lmp


def read_dump(lmp, file, species=1):
    # Read atoms positions, velocities and box parameters.
    lmp.command("lattice          fcc 4.0")
    lmp.command("region           box block 0 2 0 2 0 2")
    lmp.command("create_box       %d box"%species)
    lmp.command("read_dump        %s 0 x y z vx vy vz scaled no box yes add keep"%file)
    return lmp


def set_hybrid_potential(lmp, options, eps):
    pc =  options["md"]["pair_coeff"]
    pcraw = pc.split()
    pcnew = " ".join([*pcraw[:2], *[options["md"]["pair_style"],], *pcraw[2:]])
    
    lmp.command("pair_style       hybrid/overlay %s ufm 7.5"%options["md"]["pair_style"])
    lmp.command("pair_coeff       %s"%pcnew)
    lmp.command("pair_coeff       * * ufm %f 1.5"%eps) 

    for i in range(options["nelements"]):
        lmp.mass(i+1, options["mass"][i])
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

    for i in range(options["nelements"]):
        lmp.mass(i+1, options["mass"][i])
    return lmp

def remap_box(lmp, x, y, z):
    lmp.command("run 0")
    lmp.command("change_box     all x final 0.0 %f y final 0.0 %f z final 0.0 %f remap units box"%(x, y, z))
    return lmp


def compute_msd(lmp, options):
    elements = options["element"]
    str1 = "fix  4 all ave/time %d %d %d "%(int(options["md"]["nevery"]),
                                           int(options["md"]["nrepeat"]), 
              int(options["md"]["nevery"]*options["md"]["nrepeat"]))

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