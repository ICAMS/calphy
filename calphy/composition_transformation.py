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

import re
import numpy as np
import os
import pyscal.core as pc
from mendeleev import element
from ase.io import read, write
from ase.atoms import Atoms

class CompositionTransformation:
    """
    Class for performing composition transformations and 
    generating necessary pair styles for such transformations.

    Parameters
    ----------
    input_structure: ASE object, LAMMPS Data file or LAMMPS dump file
        input structure which is used for composition transformation

    input_chemical_formula: dict
        dictionary of input chemical 

    output_chemical_formula: string
        the required chemical composition string
    """
    def __init__(self, input_structure, input_chemical_composition, output_chemical_composition):
        self.structure = self.prepare_structure(input_structure)
        self.input_chemical_composition = input_chemical_composition
        self.output_chemical_composition = output_chemical_composition
        self.actual_species = None
        self.new_species = None
        self.maxtype = None        
        self.atom_mark = None
        self.atom_species = None
        self.mappings = None
        self.unique_mappings = None
        self.prepare_mappings()
    
    def dict_to_string(self, inputdict):
        strlst = []
        for key, val in inputdict.items():
            strlst.append(str(key))
            strlst.append(str(val))
        return "".join(strlst)
    
    def is_data_file(self, filename):
        try:
            atoms = read(filename, format="lammps-data", style="atomic")
            return atoms
        except:
            return None
    
    def is_dump_file(self, filename):
        try:
            atoms = read(filename, format="lammps-dump-text")
            return atoms
        except:
            return None
        
    def prepare_structure(self, input_structure):
        """
        Check the format of a given input file and validate it.
        """
        if isinstance(input_structure, Atoms):
            return input_structure
        elif os.path.exists(input_structure):
            atoms = self.is_data_file(input_structure)
            if atoms is None:
                atoms = self.is_dump_file(input_structure)
            return atoms
    
    def convert_to_pyscal(self):
        """
        Convert a given system to pyscal and give a dict of type mappings
        """
        pstruct = pc.System()
        pstruct.read_inputfile(self.structure, format="ase")
        
        #here we have to validate the input composition dict; and map it
        composition = pstruct.get_concentration()
        atomsymbols = []
        atomtypes = []
        for key, val in self.input_chemical_composition.items():
            if not val==0:
                found = False
                for key1, val1 in composition.items():
                    if val1==val:
                        found = True
                        atomsymbols.append(key)
                        atomtypes.append(int(key1))
                        del composition[key1]
                        break
                if not found:
                    raise ValueError("Input structure and composition do not match!")
        
        self.pyscal_structure = pstruct
        self.typedict = dict(zip(atomsymbols, atomtypes))
        self.reversetypedict = dict(zip(atomtypes, atomsymbols))
        self.natoms = self.pyscal_structure.natoms
        
        self.actual_species = len(self.typedict)
        self.new_species = len(self.output_chemical_composition) - len(self.typedict)
        self.maxtype = self.actual_species + 1 #+ self.new_species
        print(self.typedict)

    def check_if_atoms_conserved(self):
        """
        Check if a given transformation is possible by checking if number of atoms are conserved
        """
        natoms1 = np.sum([val for key, val in self.input_chemical_composition.items()])
        natoms2 = np.sum([val for key, val in self.output_chemical_composition.items()])
        if not (natoms1==natoms2==self.natoms):
            raise ValueError(f"Input and output number of atoms are not conserved! Input {self.dict_to_string(self.input_chemical_composition)}, output {self.dict_to_string(self.output_chemical_composition)}, total atoms in structure {self.natoms}")

    def get_composition_transformation(self):
        """
        From the two given composition transformation, find the transformation dict
        """
        fdiff = {}
        for key, val in self.output_chemical_composition.items():
            if key in self.input_chemical_composition.keys():
                fdiff[key] = val - self.input_chemical_composition[key]
            else:
                fdiff[key] = val - 0
        to_remove = {}
        to_add = {}

        for key, val in fdiff.items():
            if val < 0:
                to_remove[key] = np.abs(val)
            else:
                to_add[key] = val
        
        self.to_remove = to_remove
        self.to_add = to_add

    def get_random_index_of_species(self, species):
        """
        Get a random index of a given species
        """
        ids = [count for count, x in enumerate(self.atom_type) if x==species]
        return ids[np.random.randint(0, len(ids))]
    
    def mark_atoms(self):
        for i in range(self.natoms):
            self.atom_mark.append(False)
            self.atom_type = [atom.type for atom in self.pyscal_structure.iter_atoms()]
            self.mappings = [f"{x}-{x}" for x in self.atom_type]
            
    def update_mark_atoms(self):
        self.marked_atoms = []
        for key, val in self.to_remove.items():
            #print(f"Element {key}, count {val}")
            for i in range(100000):
                rint = self.get_random_index_of_species(self.typedict[key])
                self.atom_mark[rint] = True
                self.marked_atoms.append(rint)
                val -= 1
                if (val <= 0):
                    break 
    
    def get_mappings(self):
        #in a cycle add things to the typedict
        for key, val in self.to_add.items():
            #print(f"Element {key}, count {val}")
            if key in self.typedict.keys():
                newtype = self.typedict[key]
            else:
                newtype = self.maxtype
                self.typedict[key] = self.maxtype
                self.reversetypedict[self.maxtype] = key
                self.maxtype += 1
                #print(f"Element {key}, newtype {newtype}")
        
        #now we go through
        for key, val in self.to_add.items():
            #now get all
            for x in range(val):
                index = np.random.randint(0, len(self.marked_atoms))
                int_to_change = self.marked_atoms[index]
                del self.marked_atoms[index]
                self.mappings[int_to_change] = f"{self.atom_type[int_to_change]}-{self.typedict[key]}"
        
        self.unique_mappings, self.unique_mapping_counts = np.unique(self.mappings, return_counts=True)
    
    def prepare_pair_lists(self):
        self.pair_list_old = []
        self.pair_list_new = []
        for mapping in self.unique_mappings:
            map_split = mapping.split("-")
            #conserved atom
            if (map_split[0]==map_split[1]):
                self.pair_list_old.append(self.reversetypedict[int(map_split[0])])
                self.pair_list_new.append(self.reversetypedict[int(map_split[0])])
            else:
                self.pair_list_old.append(self.reversetypedict[int(map_split[0])])
                self.pair_list_new.append(self.reversetypedict[int(map_split[1])]) 
        self.new_atomtype = np.array(range(len(self.unique_mappings)))+1
        self.mappingdict = dict(zip(self.unique_mappings, self.new_atomtype))
    
    def update_types(self):
        for x in range(len(self.atom_type)):
            self.atom_type[x] = self.mappingdict[self.mappings[x]]
        atoms = self.pyscal_structure.atoms
        for count, atom in enumerate(atoms):
            atom.type = self.atom_type[count]
        self.pyscal_structure.atoms = atoms
            
    def iselement(self, symbol):
        try:
            _ = element(symbol)
            return True
        except:
            return False
    
    def update_pair_coeff(self, pair_coeff):
        pcsplit = pair_coeff.strip().split()
        pc_before = []
        pc_after = []

        started = False
        stopped = False

        for p in pcsplit:
            if ((not self.iselement(p)) and (not started)):
                pc_before.append(p)
            elif (self.iselement(p) and (not started)):
                started = True
            elif ((not self.iselement(p)) and started):
                stopped = True
            elif ((not self.iselement(p)) and stopped):
                pc_after.append(p)
                
        pc_old = " ".join([*pc_before, *self.pair_list_old, *pc_after])
        pc_new = " ".join([*pc_before, *self.pair_list_new, *pc_after])
        return pc_old, pc_new
    
    def write_structure(self, outfilename):
        self.pyscal_structure.to_file(outfilename)
        
    def prepare_mappings(self):
        self.atom_mark = []
        self.atom_species = []
        self.mappings = []
        self.unique_mappings = []

        self.convert_to_pyscal()
        self.check_if_atoms_conserved()
        self.get_composition_transformation()

        self.mark_atoms()
        self.update_mark_atoms()
        self.get_mappings()
        self.prepare_pair_lists()
        self.update_types()