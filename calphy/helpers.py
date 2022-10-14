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
from ase.io import read, write
from pyscal.trajectory import Trajectory


def create_object(cores, directory, timestep, cmdargs=None):
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
        working_directory=directory, cmdargs=cmdargs)
    
    lmp.units("metal")
    lmp.boundary("p p p")
    lmp.atom_style("atomic")
    lmp.timestep(timestep)
    return lmp

def create_structure(lmp, calc, species=None):
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
    l, alat, apc, conc, dumpfile = pl.prepare_lattice(calc)

    if species is None:
        species = len(conc)

    if l == "file":
        if dumpfile:
            reset_timestep(calc.lattice, calc.lattice, keys=None)
            lmp.command("lattice          fcc 4.0")
            lmp.command("region           box block 0 2 0 2 0 2")
            lmp.command("create_box       %d box"%species)
            lmp.command("read_dump        %s 0 x y z scaled no box yes add keep"%calc.lattice)
        else:
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

def convert_to_data_file(inputfile, outputfile, ghost_elements=0):
    atoms = read(inputfile, format="lammps-dump-text")
    write(outputfile, atoms, format="lammps-data")
    
    if ghost_elements > 0:
        lines = []
        with open(outputfile, "r") as fin:
            for line in fin:
                raw = line.strip().split()
                if (len(raw) == 3) and (raw[2] == "types"):
                    raw[0] = "%d"%ghost_elements
                    raw.append("\n")
                    rline = " ".join(raw)
                    lines.append(rline)
                else:
                    lines.append(line)

        with open(outputfile, "w") as fout:
            for line in lines:
                fout.write(line)

def get_structures(file, species, index=None):
    traj = Trajectory(file)
    if index is None:
        aseobjs = traj[:].to_ase(species=species)
    else:
        aseobjs = traj[index].to_ase(species=species)
    return aseobjs

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

def reset_timestep(file, conf, keys=["vx", "vy", "vz", "mass"]):
    sys = pc.System()
    sys.read_inputfile(file, customkeys=keys)
    sys.to_file(conf, customkeys=keys)


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
                logger.info("Replacing input spring constant None with %f"%replace_data[count])

    return data

def validate_spring_constants(data, klo=0.0001, khi=1000.0, logger=None):
    """
    Validate spring constants and replace them if needed
    """
    #first find a sane value
    sane_k = 0.1
    found = False

    for d in data:
        if (klo <= d <= khi):
            sane_k = d
            found = True
            break

    if not found:
        raise ValueError("No spring constant values are between %f and %f"%(klo, khi))

    for count, d in enumerate(data):
        if not (klo <= d <= khi):
            data[count] = sane_k
            if logger is not None:
                logger.info("Replace insane k %s for element %d with %f"%(str(d), count, sane_k))

    return data

