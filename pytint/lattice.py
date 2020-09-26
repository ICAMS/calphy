
from mendeleev import element

latticedict = {
	"BCC" :{"LQD": 1.00000, "BCC":1.00000, "FCC":0.79370, "HCP":1.12246, "DIA":0.62996, "SC":1.25992, "N":2},
	"FCC" :{"LQD": 1.00000, "BCC":1.25992, "FCC":1.00000, "HCP":1.78179, "DIA":0.79370, "SC":1.58740, "N":4},
	"HCP" :{"LQD": 1.00000, "BCC":0.89090, "FCC":0.79370, "HCP":1.00000, "DIA":0.62996, "SC":0.89089, "N":4},
	"DIA" :{"LQD": 1.00000, "BCC":1.58740, "FCC":0.79370, "HCP":1.25992, "DIA":1.00000, "SC":2.00000, "N":8},
	"SC"  :{"LQD": 1.00000, "BCC":0.79370, "FCC":0.62996, "HCP":1.12247, "DIA":0.50000, "SC":1.00000, "N":1},
}

def get_lattice(symbol, lattice_list):
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

	lattice_constants = []
	atoms_per_cell = []
	lammps_lattice = []
	#print(lattice_list)
	for lat in lattice_list:
		#print(mainlat, lat)
		newlat = latticedict[mainlat][lat]*mainalat
		lattice_constants.append(newlat)
		if lat == "LQD":
			atoms_per_cell.append(latticedict[mainlat]["N"])
			lammps_lattice.append(mainlat.lower())	
		else:
			atoms_per_cell.append(latticedict[lat]["N"])
			lammps_lattice.append(lat.lower())

	return lattice_constants, atoms_per_cell, lammps_lattice
