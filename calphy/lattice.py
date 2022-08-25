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

from mendeleev import element
import os
from pylammpsmpi import LammpsLibrary
import numpy as np
import pyscal.core as pc

"""
Conversion factors for creating initial lattices
"""
latticedict = {
    "BCC" :{"LQD": 1.00000, "BCC":1.00000, "FCC":0.79370, "HCP":1.12246, "DIA":0.62996, "SC":1.25992, "N":2},
    "FCC" :{"LQD": 1.00000, "BCC":1.25992, "FCC":1.00000, "HCP":1.78179, "DIA":0.79370, "SC":1.58740, "N":4},
    "HCP" :{"LQD": 1.00000, "BCC":0.89090, "FCC":0.79370, "HCP":1.00000, "DIA":0.62996, "SC":0.89089, "N":4},
    "DIA" :{"LQD": 1.00000, "BCC":1.58740, "FCC":0.79370, "HCP":1.25992, "DIA":1.00000, "SC":2.00000, "N":8},
    "SC"  :{"LQD": 1.00000, "BCC":0.79370, "FCC":0.62996, "HCP":1.12247, "DIA":0.50000, "SC":1.00000, "N":1},
}

def get_lattice(symbol, lat):
    """
    Find lattice constants of an element

    Parameters
    ----------
    symbol : string
        symbol of chemical element

    lattice_list : list of strings
        list of lattices

    Returns
    -------
    lattice_constants : list of floats
        list of lattice constant values

    atoms_per_cell : list of ints
        number of atoms per cell

    lammps_lattice : list of strings
        the main lattice to be used in lammps
    """

    chem = element(symbol)
    
    mainlat = chem.lattice_structure
    
    if mainlat == "HEX":
        mainlat = "HCP"

    mainalat = chem.lattice_constant

    #print(mainlat, lat)
    newlat = latticedict[mainlat][lat]*mainalat
    lattice_constant = newlat

    if lat == "LQD":
        atoms_per_cell = latticedict[mainlat]["N"]
        lammps_lattice = mainlat.lower()    
    else:
        atoms_per_cell = latticedict[lat]["N"]
        lammps_lattice = lat.lower()

    return lattice_constant, atoms_per_cell, lammps_lattice

def check_dump_file(infile):
    try:
        #now use pyscal to read it in,
        sys = pc.System()
        sys.read_inputfile(infile)
        atoms = sys.atoms
        natoms = len(atoms)
        types = [atom.type for atom in atoms]
        xx, xxcounts = np.unique(types, return_counts=True)
        conc = xxcounts/np.sum(xxcounts)
        return (natoms, conc)
    except:
        return None

def check_data_file(infile):
    try:
        lmp = LammpsLibrary(mode="local", cores=1, 
            working_directory=os.getcwd())
        lmp.units("metal")
        lmp.boundary("p p p")
        lmp.atom_style("atomic")
        lmp.timestep(0.001)            
        lmp.read_data(infile)
        natoms = lmp.natoms
        #now we convert to a dump file and read the concentration
        trajfile = ".".join([infile, "dump"])
        lmp.command("mass * 1.0")
        lmp.dump("2 all custom", 1, trajfile,"id type x y z")
        lmp.run(0)
        lmp.undump(2)
        #now use pyscal to read it in,
        sys = pc.System()
        sys.read_inputfile(trajfile)
        atoms = sys.atoms
        types = [atom.type for atom in atoms]
        xx, xxcounts = np.unique(types, return_counts=True)
        conc = xxcounts/np.sum(xxcounts)
        lmp.close()
        return (natoms, conc)
    except:
        return None

def prepare_lattice(calc):
    #process lattice
    lattice = calc.lattice.upper()
    dumpfile = False

    if lattice in ["BCC", "FCC", "HCP", "DIA", "SC", "LQD"]:
        #process lattice
        #throw error for multicomponent
        if calc.n_elements > 1:
            raise ValueError("Only files supported for multicomponent")

        alat, apc, l = get_lattice(calc.element[0], calc.lattice)

        #replace lattice constant
        if calc.lattice_constant != 0:
            alat = calc.lattice_constant
        
        conc = [1,]

    elif os.path.exists(calc.lattice):
        calc.lattice = os.path.abspath(calc.lattice)

        res = check_dump_file(calc.lattice)
        dumpfile = True

        if res is None:
            res = check_data_file(calc.lattice)
            dumpfile = False

        if res is not None:
            natoms = res[0]
            conc = res[1]
            l = "file"
            alat = 1.00
            apc = natoms
        else:
            raise ValueError("An input file was provided but it was neither data or dump file")

    else:
        raise ValueError("Unknown lattice found. Allowed options are BCC, FCC, HCP, DIA, SC or LQD; or an input file.")
    
    if l == "dia":
        l = "diamond"

    return l, alat, apc, conc, dumpfile 

