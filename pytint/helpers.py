"""
Helper methods for pytint
"""
from pylammpsmpi import LammpsLibrary
import logging
import numpy as np
from lammps import lammps
import pytint.lattice as pl
import pyscal.core as pc
"""
LAMMPS helper routines
--------------------------------------------------------------------
"""
def check_data_file(file):
    try:
        lmp = lammps()
        lmp.command("read_data %s"%file)
        natoms = lmp.natoms
        lmp.close()
    except:
        raise TypeError("LAMMPS could not read in the data file. Please check!")
    return natoms

def create_object(cores, directory, timestep):
    """
    """
    lmp = LammpsLibrary(mode="local", cores=cores, 
        working_directory=directory)
    
    lmp.units("metal")
    lmp.boundary("p p p")
    lmp.atom_style("atomic")
    lmp.timestep(timestep)
    return lmp

def create_structure(lmp, calc):
    l, alat, apc = pl.prepare_lattice(calc)

    if l == "file":
        lmp.command("read_data      %s"%file)
    else:
        lmp.lattice(l, alat)
        lmp.region("box block", 0, calc["repeat"][0], 
            0, calc["repeat"][1], 
            0, calc["repeat"][2])
        lmp.create_box("1 box")
        lmp.create_atoms("1 box")
    return lmp


def set_potential(lmp, options):
    lmp.pair_style(options["md"]["pair_style"])
    lmp.pair_coeff(options["md"]["pair_coeff"])

    for i in range(options["nelements"]):
        lmp.mass(i+1, options["mass"][i])
    return lmp


def read_dump(lmp, file):
    # Read atoms positions, velocities and box parameters.
    lmp.command("lattice          fcc 4.0")
    lmp.command("region           box block 0 2 0 2 0 2")
    lmp.command("create_box       1 box")
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

def calculate_concentrations(file):
    sys = pc.System()
    sys.read_inputfile(file)
    atoms = sys.atoms
    types = [atom.type for atom in atoms]
    unq_types, unq_counts = np.unique(types, return_counts=True)
    return unq_counts/np.sum(unq_counts)

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