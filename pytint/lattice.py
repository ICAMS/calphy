
from mendeleev import element

latticedict = {
	"BCC" :{"BCC":1.00000, "FCC":0.79370, "HCP":1.12246, "DIA":0.62996, "SC":1.25992, "N":2},
	"FCC" :{"BCC":1.25992, "FCC":1.00000, "HCP":1.78179, "DIA":0.79370, "SC":1.58740, "N":4},
	"HCP" :{"BCC":0.89090, "FCC":0.79370, "HCP":1.00000, "DIA":0.62996, "SC":0.89089, "N":4},
	"DIA" :{"BCC":1.58740, "FCC":0.79370, "HCP":1.25992, "DIA":1.00000, "SC":2.00000, "N":8},
	"SC"  :{"BCC":0.79370, "FCC":0.62996, "HCP":1.12247, "DIA":0.50000, "SC":1.00000, "N":1},
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
	"""

	chem = element(symbol)
	
	mainlat = si.lattice_structure
	
	if mainlat == "HEX":
		mainlat = "HCP"

	mainalat = si.lattice_constant

	lattice_constants = []
	atoms_per_cell = []
	
	for lat in lattice_list:
		newlat = latticedict[mainlat][lat]*mainalat
		lattice_list.append(newlat)
		atoms_per_cell.append(latticedict[lat]["N"])

	return lattice_constants, atoms_per_cell
