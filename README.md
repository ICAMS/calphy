
# pytint
Python library for free energy calculations using thermodynamic integration

### Installing LAMMPS

pytint needs a custom version of lammps compiled as a library. A custom version with fix adapt implemented for eam, snap, pace and sw potentials is available [here](https://git.noc.ruhr-uni-bochum.de/sarath/lammps-ace).

It can be compiled by-
```
mkdir build_lib
cd build_lib
cmake -D BUILD_LIB=ON -D BUILD_SHARED_LIBS=ON -D BUILD_MPI=ON -D PKG_MANYBODY=ON -D PKG_USER-MISC=ON -D PKG_USER-PACE=ON ../cmake
make # -j${NUM_CPUS}
cp liblammps${SHLIB_EXT}* ../src  # For compatibility with the original make system.
cd ../src
make install-python 
```

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

