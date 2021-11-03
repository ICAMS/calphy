# calphy

`calphy`(pronounced _cal-phee_) is a Python library and command line tool for calculation of free energies. It uses [LAMMPS](https://www.lammps.org/) as the molecular dynamics driver to enable calculation of free energies using thermodynamic integration in a completely automated and efficient manner. The various methods that calphy can perform are:

- $F(V_i,T_i)$ and $G(P_i,T_i)$ for both solid and liquid phases at the given thermodynamic conditions using [non-equilibrium Hamiltonian interpolation](https://linkinghub.elsevier.com/retrieve/pii/S0927025615007089).
- $F(T_i \to T_f)$ and $G(T_i \to T_f)$, temperature dependence of Gibbs and Helmholtz free energies using the [reversible scaling](https://link.aps.org/doi/10.1103/PhysRevLett.83.3973) approach.
- Calculation of solid-solid or solid-liquid phase transition temperatures.
- Calculation of coexistence lines using [dynamic Clausius-Clapeyron integration](http://aip.scitation.org/doi/10.1063/1.1420486).
- Calculation of specific heat $c_P(T)$ as a function of temperature.
- Calculation of $F(x, T)$ and $G(x, T)$ for multicomponent systems using [alchemical transformations](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.103801).
- [Upsampling calculations](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.103801).

Calphy works with all [interatomic potentials implemented in LAMMPS](https://docs.lammps.org/pairs.html) and can obtain reliable results with low error bars in as low as 50 ps of simulation time with less than 3000 atoms. More information about the methods in calphy can be found in the [associated publication](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.103801).


## Installation

### Supported operating systems

`calphy` can be installed on Linux and Mac OS based systems. On Windows systems, it is recommended to use  [Windows subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install).

### Installation using [Conda](https://anaconda.org/)

calphy can be installed directly using [Conda](https://docs.conda.io/en/latest/) from the [conda-forge channel](https://conda-forge.org/) by the following statement:

```
conda install -c conda-forge calphy
```

### Installation from the repository

calphy can be built from the repository by-

```
git clone https://github.com/ICAMS/calphy.git
cd calphy
python setup.py install --user
```

Please note that the list of dependencies mentioned below has to installed manually by the user.

### Using a conda environment

It is **strongly** recommended to install and use `calphy` within a conda environment. To see how you can install conda see [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Once a conda distribution is available, the following steps will help set up an environment to use `calphy`. First step is to clone the repository.

```
git clone https://github.com/ICAMS/calphy.git
```

After cloning, an environment can be created from the included file-

```
cd calphy
conda env create -f environment.yml
```

Note that the conda-forge distribution of LAMMPS will be automatically installed. Alternatively, you can use an existing version of LAMMPS (compiled as library). If a LAMMPS distribution need not be installed, use the  `calphy/environment-nolammps.yml` file instead to create the environment.
This will install the necessary packages and create an environment called calphy. It can be activated by,

```
conda activate calphy
```

then, install `calphy` using,

```
python setup.py install
```
The environment is now set up to run calphy.

### Dependencies

- lammps                            `conda install -c conda-forge lammps`
- mendeleev           >=0.7.0       `pip install mendeleev`
- pylammpsmpi         >=0.0.8       `pip install pylammpsmpi`
- pyscal              >=2.10.14     `pip install git+https://github.com/srmnitc/pyscal`
- pyyaml              >=5.4.1       `pip install pyyaml`
- scipy               >=1.7.0       `pip install scipy`
- tqdm                >=4.61.2      `pip install tqdm`

#### Optional

- matplotlib          >=3.4.2       `pip install matplotlib`
- pytest              >=6.2.4       `pip install pytest`

## About [LAMMPS](https://www.lammps.org/) for `calphy`

calphy uses LAMMPS as the driver for molecular dynamics simulations. For calphy to work, LAMMPS needs to be compiled as a library along with the Python interface. The easiest way to do this is to install LAMMPS through the conda-forge channel using:

```
conda install -c conda-forge lammps
```

Alternatively, when interatomic potentials with special compilation needs are to be used, LAMMPS (Stable release 29 Sept 2021 and above) can be compiled manually using the following set of instructions:

First step is to obtain the stable version from [here](https://github.com/lammps/lammps/archive/refs/tags/stable_29Sep2021.tar.gz) and extract the archive. From the extracted archive, the following steps, used in the [conda-forge recipe](https://github.com/conda-forge/lammps-feedstock/blob/master/recipe/build.sh) can be run:

```
mkdir build_lib
cd build_lib
cmake -D BUILD_LIB=ON -D BUILD_SHARED_LIBS=ON -D BUILD_MPI=ON -D PKG_MPIIO=ON -D LAMMPS_EXCEPTIONS=yes -D PKG_MANYBODY=ON -D PKG_MISC=ON -D PKG_MISC=ON -D PKG_EXTRA-COMPUTE=ON -D PKG_EXTRA-DUMP=ON -D PKG_EXTRA-FIX=ON -D PKG_EXTRA-PAIR=ON ../cmake
make
cp liblammps${SHLIB_EXT}* ../src 
cd ../src
make install-python 
mkdir -p $PREFIX/include/lammps
cp library.h $PREFIX/include/lammps
cp liblammps${SHLIB_EXT}* "${PREFIX}"/lib/
cd ..
```

The above commands only builds the [MANYBODY](https://docs.lammps.org/Packages_details.html#pkg-manybody) package. To use some of the other potentials, the following commands could be added to the `cmake` call.

- `-D PKG_ML-PACE` for performant [Atomic Cluster Expansion](https://docs.lammps.org/Packages_details.html#pkg-ml-pace) potential.
- `-D PKG_ML-SNAP=ON`for [SNAP potential](https://docs.lammps.org/Packages_details.html#pkg-ml-snap).
- `-D PKG_MEAM=ON` for [MEAM potential](https://docs.lammps.org/Packages_details.html#meam-package).
- `-D PKG_KIM=ON` for [KIM support](https://docs.lammps.org/Packages_details.html#pkg-kim).

Now the python wrapper can be installed:

```
cd ../src
make install-python
```
In the case of a conda environment, the following commands can be used to copy the compiled libraries to an accessible path:

```
mkdir -p $CONDA_PREFIX/include/lammps
cp library.h $CONDA_PREFIX/include/lammps
cp liblammps${SHLIB_EXT}* $CONDA_PREFIX/lib/
```

The combined libraries **should be available** on the system or the environment paths. Please see here for more details.

Once LAMMPS is compiled and the libraries are available in an accessible location, the following commands can be used within python to test the installation:

```
from lammps import lammps
lmp = lammps()
```
Now, we install pylammpsmpi using,

```
pip install pylammpsmpi
```
and finally calphy:

```
cd calphy
python setup.py install --user
```









