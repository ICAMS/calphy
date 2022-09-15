# calphy

`calphy`(pronounced _cal-phee_) is a Python library and command line tool for calculation of free energies. It uses [LAMMPS](https://www.lammps.org/) as the molecular dynamics driver to enable calculation of free energies using thermodynamic integration in a completely automated and efficient manner. The various methods that calphy can perform are:

- <img src="https://render.githubusercontent.com/render/math?math=F(V_i,T_i)"> and <img src="https://render.githubusercontent.com/render/math?math=G(P_i,T_i)"> for both solid and liquid phases at the given thermodynamic conditions using [non-equilibrium Hamiltonian interpolation](https://linkinghub.elsevier.com/retrieve/pii/S0927025615007089).
- <img src="https://render.githubusercontent.com/render/math?math=F(T_i \to T_f)"> and <img src="https://render.githubusercontent.com/render/math?math=G(T_i \to T_f)">, temperature dependence of Gibbs and Helmholtz free energies using the [reversible scaling](https://link.aps.org/doi/10.1103/PhysRevLett.83.3973) approach.
- Calculation of solid-solid or solid-liquid phase transition temperatures.
- Calculation of coexistence lines using [dynamic Clausius-Clapeyron integration](http://aip.scitation.org/doi/10.1063/1.1420486).
- Calculation of specific heat <img src="https://render.githubusercontent.com/render/math?math=c_P(T)"> as a function of temperature.
- Calculation of <img src="https://render.githubusercontent.com/render/math?math=F(x, T)"> and <img src="https://render.githubusercontent.com/render/math?math=G(x, T)"> for multicomponent systems using [alchemical transformations](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.103801).
- [Upsampling calculations](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.103801).

Calphy works with all [interatomic potentials implemented in LAMMPS](https://docs.lammps.org/pairs.html) and can obtain reliable results with low error bars in as low as 50 ps of simulation time with less than 3000 atoms. More information about the methods in calphy can be found in the [associated publication](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.103801).

## Installation, examples and documentation

- Read the `calphy` installation instructions and documentation including examples [here](https://calphy.readthedocs.io/en/latest/index.html).

- For a description of the methods and algorithms implemented in `calphy`, please see the [associated publication](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.103801).

- For a set of examples presented in the `calphy` publication, please see [here](https://github.com/srmnitc/calphy-publication-examples).
 
## Citing calphy

If you find calphy useful, please consider citing:  
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz.   
[“Automated Free Energy Calculation from Atomistic Simulations.”](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.103801) 
Physical Review Materials 5(10), 2021  
DOI: 10.1103/PhysRevMaterials.5.103801  

Download bibtex [here](https://journals.aps.org/prmaterials/export/10.1103/PhysRevMaterials.5.103801?type=bibtex&download=true)


[![Anaconda-Server Badge](https://anaconda.org/conda-forge/calphy/badges/installer/conda.svg)](https://conda.anaconda.org/conda-forge)

[![PyPI version](https://badge.fury.io/py/calphy.svg)](https://badge.fury.io/py/calphy)