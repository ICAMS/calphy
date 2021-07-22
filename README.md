# calphy

Python library and command line tool for calculation of free energies. Read the documentation [here](https://calphy.readthedocs.io/en/latest/index.html).

## Dependencies

- lammps              2021.05.27  
- mendeleev           0.7.0       `pip install mendeleev`
- pylammpsmpi         0.0.8       `pip install pylammpsmpi`
- pyscal              2.10.14     `pip install git+https://github.com/srmnitc/pyscal`
- pyyaml              5.4.1       `pip install pyyaml`
- scipy               1.7.0       `pip install scipy`
- tqdm                4.61.2      `pip install tqdm`

### Optional

- matplotlib          3.4.2       `pip install matplotlib`
- pytest              6.2.4       `pip install pytest`


## Installing

### Standard installation

Standard installation uses the conda forge LAMMPS package. In this version, only `pair_style eam` is compatible with calphy.
    
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

This will install the necessary packages and create an environment called calphy. It can be activated by,

```
conda activate calphy
```

then, install `calphy` using,

```
python setup.py install --user
```
The environment is now set up to run calphy.

### PACE support

> **NOTE**: 08 July 2021: The upcoming release of [LAMMPS](https://github.com/lammps/lammps/releases) will significantly change the compilation procedure and render the below information outdated.

If you want support for `pair_style pace` [1], or want to use a custom version of LAMMPS, it is recommended to follow these steps:

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
make # -j${NUM_CPUS}
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


### Support for other pair styles

`calphy` also supports Stillinger-Weber, Tersoff, MEAM, SNAP, and ADP pair styles using a custom version. If you are interested in using these pair styles, please contact [us](mailto:sarath.menon@rub.de).

## Documentation

- [Inputfile](examples/inputfile.md)
- [Command line interface](examples/cli.md)
- Examples  
    > Note that the examples use small system sizes and switching times so that the functionality can be illustrated. For better results, it is recommended to increase both of these quantities.
    - [Free energy calculation](examples/example_01)
    - [BCC to FCC transition in Fe](examples/example_02): Use free energy calculation and temperature sweep to calculate phase transition temperature in Fe.
    - [Melting temperature of Cu](examples/example_03): Use free energy calculation and temperature sweep to calculate melting temperature.
    - [Pressure-temperature phase diagram of Cu](examples/example_04): Calculate Gibbs free energy and reversible scaling to calculate the pressure-temperature phase diagram of Cu.

## Citing calphy

If you find calphy useful, please consider citing:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz. 
“Automated Free Energy Calculation from Atomistic Simulations.” 
ArXiv:2107.08980 [Cond-Mat], July 19, 2021. 
http://arxiv.org/abs/2107.08980.
