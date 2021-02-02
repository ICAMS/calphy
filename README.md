# pytint
Python library for free energy calculations using thermodynamic integration

> `pytint` is still in alpha phase. The documentation and features might change abruptly! Significant changes are coming up in the new release. See `update/input` branch.

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

#### Install dependencies

The following packages need to be installed. 

- numpy (`conda install -c conda-forge numpy`)
- scipy (`conda install -c conda-forge scipy`)
- pyyaml (`conda install -c conda-forge pyyaml`)
- mendeleev (`conda install -c conda-forge mendeleev`)
- pylammpsmpi (`conda install -c conda-forge pylammpsmpi`)
- pyscal (`conda install -c conda-forge pyscal`)

### Install `pytint`

After installing the requirements, `pytint` can be installed by,

```
git clone https://git.noc.ruhr-uni-bochum.de/atomicclusterexpansion/pytint.git
cd pytint
python setup.py install
```

## Building the documentation

```
cd pytint/docs
pip install -r requirements.txt
make html
```

The files will be in `pytint/docs/build/html`.


## Running a calculation

`pytint` can be run as both a Python library and as a command line tool. The recommended way to use `pytint` is through the command line. After installation, `pytint` can be accessed from the terminal using,

```
tint --help
```

The main option one needs to specify is the `--input` or `-i`. This keyword species the location of the input file. The format of the inputfile is discussed in detail [here](inputfile.md).

```
tint -i input.yaml
```

Such a command will read the input file and start NEHI calculations **for each temperature** mentioned in the input file. Alternatively, one can use the `--mode` option to launch a reversible scaling calculation.

```
tint -i input.yaml -m rs
```

In this case, **one** NEHI calculation is done for the first temperature mentioned in the input file, and then a reversible scaling calculation is done to extend the free energy up to the last temperature specified in the input file. 

### Publications

Freitas, Rodrigo, Mark Asta, and Maurice de Koning. “Nonequilibrium Free-Energy Calculation of Solids Using LAMMPS.” Computational Materials Science 112 (February 2016): 333–41. https://doi.org/10.1016/j.commatsci.2015.10.050.  

Paula Leite, Rodolfo, and Maurice de Koning. “Nonequilibrium Free-Energy Calculations of Fluids Using LAMMPS.” Computational Materials Science 159 (March 2019): 316–26. https://doi.org/10.1016/j.commatsci.2018.12.029.  

Koning, Maurice de, A. Antonelli, and Sidney Yip. “Optimized Free-Energy Evaluation Using a Single Reversible-Scaling Simulation.” Physical Review Letters 83, no. 20 (November 15, 1999): 3973–77. https://doi.org/10.1103/PhysRevLett.83.3973.

