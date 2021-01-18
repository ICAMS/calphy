
# pytint
Python library for free energy calculations using thermodynamic integration

### Installing LAMMPS

`pytint` works with the standard version of LAMMPS. Currently supported pair styles are `pace`(works with [lammps-ace](https://git.noc.ruhr-uni-bochum.de/atomicclusterexpansion/lammps-ace)) and `eam` (with standard LAMMPS). If you want use other pair styles such as `snap` or `sw`, please [contact](mailto:sarath.menon@rub.de).

`pytint` needs LAMMPS compiled as a library with Python support. It can be done by the following instructions-

```
cd lammps
mkdir build_lib
cd build_lib
cmake -D BUILD_LIB=ON -D BUILD_SHARED_LIBS=ON -D BUILD_MPI=ON -D PKG_MANYBODY=ON -D PKG_USER-MISC=ON -D PKG_USER-PACE=ON ../cmake
make # -j${NUM_CPUS}
cp liblammps${SHLIB_EXT}* ../src
cd ../src
make install-python 
```

The `include` files and compiled files should be available in the paths. A full set of instructions can be found [here](https://lammps.sandia.gov/doc/Python_install.html).

### Installing `pytint`

Requirements for `pytint` is given in `requirements.txt`.

After installing the requirements, `pytint` can be installed by,

```
python setup.py install
```

### Running a calculation

A commented input file can be found in `example/input.yaml`.

The cli has a help command `tint --help`.

To run a calculation:

```
tint -i input.yaml
```

### Publications

Freitas, Rodrigo, Mark Asta, and Maurice de Koning. “Nonequilibrium Free-Energy Calculation of Solids Using LAMMPS.” Computational Materials Science 112 (February 2016): 333–41. https://doi.org/10.1016/j.commatsci.2015.10.050.  

Paula Leite, Rodolfo, and Maurice de Koning. “Nonequilibrium Free-Energy Calculations of Fluids Using LAMMPS.” Computational Materials Science 159 (March 2019): 316–26. https://doi.org/10.1016/j.commatsci.2018.12.029.  

Koning, Maurice de, A. Antonelli, and Sidney Yip. “Optimized Free-Energy Evaluation Using a Single Reversible-Scaling Simulation.” Physical Review Letters 83, no. 20 (November 15, 1999): 3973–77. https://doi.org/10.1103/PhysRevLett.83.3973.

