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
import pyscal.core as pc
from mendeleev import element

class CompositionTransformation:
    def __init__(self, input_structure, output_chemical_formula):
        self.structure = input_structure
        self.input_chemical_formula = input_structure.get_chemical_formula()
        self.output_chemical_formula = output_chemical_formula
        self.input_formula_dict = self.chemical_formula_to_dict(self.input_chemical_formula)
        self.output_formula_dict = self.chemical_formula_to_dict(self.output_chemical_formula)
        self.check_if_atoms_conserved()
        self.get_composition_transformation()
        self.actual_species = len(self.input_formula_dict)
        self.new_species = len(self.output_formula_dict) - len(self.input_formula_dict)
        self.maxtype = self.actual_species + self.new_species        
        self.atom_mark = None
        self.atom_species = None
        self.mappings = None
        self.unique_mappings = None
        self.prepare_mappings()
    
    def chemical_formula_to_dict(self, formula):
        """
        Convert a chemical formula to dict
        """
        formula_split = re.findall('(\d+|[A-Za-z]+)', formula)
        formula_dict = {}
        for x in range(0, len(formula_split), 2):
            formula_dict[formula_split[x]] = int(formula_split[x+1])
        return formula_dict

    def check_if_atoms_conserved(self):
        """
        Check if a given transformation is possible by checking if number of atoms are conserved
        """
        natoms1 = np.sum([val for key, val in self.input_formula_dict.items()])
        natoms2 = np.sum([val for key, val in self.output_formula_dict.items()])
        if not (natoms1==natoms2):
            raise ValueError(f"Input and output number of atoms are not conserved! Input {self.input_chemical_formula}, output {self.output_chemical_formula}")

    def get_composition_transformation(self):
        """
        From the two given composition transformation, find the transformation dict
        """
        fdiff = {}
        for key, val in self.output_formula_dict.items():
            if key in self.input_formula_dict.keys():
                fdiff[key] = val - self.input_formula_dict[key]
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

    def convert_to_pyscal(self):
        """
        Convert a given system to pyscal and give a dict of type mappings
        """
        pstruct = pc.System()
        pstruct.read_inputfile(self.structure, format="ase")
        atomsymbols = np.unique(self.structure.get_chemical_symbols())
        atomtypes = np.array(range(1, len(atomsymbols)+1))
        self.pyscal_structure = pstruct
        self.typedict = dict(zip(atomsymbols, atomtypes))
        self.reversetypedict = dict(zip(atomtypes, atomsymbols))
        self.natoms = self.pyscal_structure.natoms

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
            
    def update_mark_atoms(self):
        for key, val in self.to_remove.items():
            print(f"Element {key}, count {val}")
            for i in range(100000):
                rint = self.get_random_index_of_species(self.typedict[key])
                self.atom_mark[rint] = True
                val -= 1
                if (val <= 0):
                    break 
    
    def get_mappings(self):
        for key, val in self.to_add.items():
            print(f"Element {key}, count {val}")
            if key in self.typedict.keys():
                newtype = self.typedict[key]
            else:
                newtype = self.maxtype
                self.typedict[key] = self.maxtype
                self.reversetypedict[self.maxtype] = key
                self.maxtype += 1
        for i in range(100000):
            for x in range(self.natoms):
                if self.atom_mark[x]:
                    self.mappings.append(f"{self.atom_type[x]}-{newtype}")
                    val -= 1
                else:
                    self.mappings.append(f"{self.atom_type[x]}-{self.atom_type[x]}")
            if (val <= 0):
                break 
        self.unique_mappings = np.unique(self.mappings)
    
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
            if ((not iselement(p)) and (not started)):
                pc_before.append(p)
            elif (iselement(p) and (not started)):
                started = True
            elif ((not iselement(p)) and started):
                stopped = True
            elif ((not iselement(p)) and stopped):
                pc_after.append(p)
                
        pc_old = " ".join([*pc_before, *self.pair_list_old, *pc_after])
        pc_new = " ".join([*pc_before, *self.pair_list_new, *pc_after])
        return pc_old, pc_new
        
    def prepare_mappings(self):
        self.atom_mark = []
        self.atom_species = []
        self.mappings = []
        self.unique_mappings = []

        self.convert_to_pyscal()
        self.mark_atoms()
        self.update_mark_atoms()
        self.get_mappings()
        self.prepare_pair_lists()
        self.update_types()