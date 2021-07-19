# calphy

Python library and command line tool for calculation of free energies.


### Dependencies

- lammps              2021.05.27  
- mendeleev           0.7.0       `pip install mendeleev`
- pylammpsmpi         0.0.8       `pip install pylammpsmpi`
- pyscal              2.10.14     `pip install git+https://github.com/srmnitc/pyscal`
- pyyaml              5.4.1       `pip install pyyaml`
- scipy               1.7.0       `pip install scipy`
- tqdm                4.61.2      `pip install tqdm`

#### Optional

- matplotlib          3.4.2       `pip install matplotlib`
- pytest              6.2.4       `pip install pytest`

### Installation

> **NOTE**: If you are planning to use a custom version of LAMMPS, read the 'Notes on the LAMMPS installation' section first.
    

It is **strongly** recommended to install and use `calphy` within a conda environment. To see how you can install conda see [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Once a conda distribution is available, the following steps will help set up an environment to use `calphy`. First step is to clone the repository.

```
git clone https://git.noc.ruhr-uni-bochum.de/atomicclusterexpansion/calphy.git
```

After cloning, an environment can be created from the included file-

```
cd calphy
conda env create -f environment.yml
```

This will install the necessary packages and create an environment called calphy. It can be activated by,

```
conda activate calphy
```

then, install `calphy` using,

```
python setup.py install --user
```
The environment is now set up to run calphy.

### Notes on the LAMMPS installation

> **NOTE**: 08 July 2021: The upcoming release of [LAMMPS](https://github.com/lammps/lammps/releases) will significantly change the compilation procedure and render the below information outdated.

- The above commands will install LAMMPS from the conda-forge channel. With this version of LAMMPS, only `eam` pair style is supported!

- For other pair styles such as Stillinger-Weber, MEAM, ADP, and SNAP, a custom version of LAMMPS is available. Please contact [here](mailto:sarath.menon@rub.de).

- For using with pair style PACE, the LAMMPS distribution has to be compiled manually:

```
wget https://github.com/lammps/lammps/archive/refs/tags/stable_29Oct2020.tar.gz
tar -xvf stable_29Oct2020.tar.gz
cd lammps-stable_29Oct2020
mkdir build_lib
cd build_lib
cmake -D BUILD_LIB=ON -D BUILD_SHARED_LIBS=ON -D BUILD_MPI=ON -D PKG_MANYBODY=ON -D PKG_USER-MISC=ON -D PKG_USER-PACE=ON ../cmake
make # -j${NUM_CPUS}
cp liblammps${SHLIB_EXT}* ../src
```

If the commands are run within a conda environment and errors during compilation are observed, try installing the following packages:

```
conda install -c conda-forge cmake gcc_linux-64 gxx_linux-64 gfortran_linux-64
```
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

The combined libraries **should be available** on the system or the environment paths. Please see [here](https://lammps.sandia.gov/doc/Python_install.html) for more details. 

Once LAMMPS is compiled and the libraries are available in an accessible location, the following commands can be used within python to test the installation:

```
from lammps import lammps
lmp = lammps()
```