# calphy

Python library and command line tool for calculation of free energies.

## Installation

In this section, the installation method for `calphy` is explained. Depending on which potential is used, a different version of LAMMPS might be required, leading to a different installation method. The table below lists the pair styles, and the corresponding installation methods.



| pair style                                              	| installation 	|
|---------------------------------------------------------	|--------------	|
| eam, eam/fs, eam/alloy                                  	| [method 1](#id1)     	|
| pace                                                    	| [method 2](#id2)     	|
| adp, snap, sw, tersoff, tersoff/mod, tersoff/modc, meam 	| [method 3](#id3)     	|

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

### Method 1 

> This method only supports the `eam`, `eam/fs` and `eam/alloy` pair styles. However, direct free energy calculations for solid state can be carried out for any pair style using this installation.

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

### Method 2

> This method supports `pace` pair style [1]. Additionally it can also be used if a custom version of LAMMPS is required. Drect free energy calculations for solid state can be carried out for any pair style using this installation.

> **NOTE**: 08 July 2021: The upcoming release of [LAMMPS](https://github.com/lammps/lammps/releases) will significantly change the compilation procedure and render the below information outdated.


For using with pair style PACE, the LAMMPS distribution has to be compiled manually:

```
git clone https://github.com/ICAMS/calphy.git
conda env create -f calphy/environment-nolammps.yml
```

Now proceed with LAMMPS installation

```
wget https://github.com/lammps/lammps/archive/refs/tags/patch_27May2021.tar.gz
tar -xvf patch_27May2021.tar.gz
cd lammps-patch_27May2021
mkdir build_lib
cd build_lib
cmake -D BUILD_LIB=ON -D BUILD_SHARED_LIBS=ON -D BUILD_MPI=ON -D PKG_MANYBODY=ON -D PKG_USER-MISC=ON -D PKG_USER-PACE=ON ../cmake
make 
cp liblammps${SHLIB_EXT}* ../src
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

Now, we install `pylammpsmpi` using,

```
pip install pylammpsmpi
```

and finally `calphy`:

```
cd calphy
python setup.py install --user
```

### Method 3

> This method supports adp, meam, snap, sw, tersoff, tersoff/mod and tersoff/modc pair styles. Additionally it can also be used if a custom version of LAMMPS is required. Drect free energy calculations for solid state can be carried out for any pair style using this installation.

> **NOTE**: 08 July 2021: The upcoming release of [LAMMPS](https://github.com/lammps/lammps/releases) will significantly change the compilation procedure and render the below information outdated.

```
git clone https://github.com/ICAMS/calphy.git
conda env create -f calphy/environment-nolammps.yml
```

Now download LAMMPS,

```
wget https://github.com/lammps/lammps/archive/refs/tags/patch_27May2021.tar.gz
tar -xvf patch_27May2021.tar.gz
```

Get the modified pair styles

```
git clone https://github.com/srmnitc/lammps-calphy.git
```

The cloned repository contains the necessary files for pair styles. For example, to use sw pair style, replace the original LAMMPS files with the modified ones:

```
cp lammps-calphy/src/MANYBODY/pair_sw.* lammps-patch_27May2021/src/MANYBODY/
```

After copying the necessary files, LAMMPS can be installed:

```
cd lammps-patch_27May2021
mkdir build_lib
cd build_lib
cmake -D BUILD_LIB=ON -D BUILD_SHARED_LIBS=ON -D BUILD_MPI=ON -D PKG_MANYBODY=ON -D PKG_USER-MISC=ON ../cmake
make
cp liblammps${SHLIB_EXT}* ../src
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

Now, we install `pylammpsmpi` using,

```
pip install pylammpsmpi
```

and finally `calphy`:

```
cd calphy
python setup.py install --user
```

[1]  Lysogorskiy, Yury, Cas van der Oord, Anton Bochkarev, Sarath Menon, Matteo Rinaldi, Thomas Hammerschmidt, Matous Mrovec, et al. “Performant Implementation of the Atomic Cluster Expansion (PACE) and Application to Copper and Silicon.” Npj Computational Materials 7, no. 1 (December 2021): 97. https://doi.org/10.1038/s41524-021-00559-9.
