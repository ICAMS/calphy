# pytint

Python library for thermodynamic integration

## Installing

It is **strongly** recommended to install and use `pytint` within a conda environment. To see how you can install conda see [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Once a conda distribution is available, the following steps will help set up an environment to use `pytint`. First step is to clone the repository.

```
git clone https://git.noc.ruhr-uni-bochum.de/atomicclusterexpansion/pytint.git
```

After cloning, an environment can be created from the included file-

```
cd pytint
conda env create -f environment.yml
```

This will install the necessary packages and create an environment called pytint. It can be activated by,

```
conda activate pytint
```

then, install `pytint` using,

```
python setup.py install
```
The environment is now set up to run pytint.

### Notes on the LAMMPS installation

- The above commands will install LAMMPS from the conda-forge channel. With this version of LAMMPS, only `eam` pair style is supported!

- For other pair styles such as Stillinger-Weber, MEAM, and SNAP, a custom version of LAMMPS will soon be available.

- For using with pair style PACE, the lammps-ace distribution has to compiled manually (for the moment). The following steps can help:

```
git clone https://git.noc.ruhr-uni-bochum.de/atomicclusterexpansion/lammps-ace
cd lammps-ace
mkdir build_lib
cd build_lib
cmake -D BUILD_LIB=ON -D BUILD_SHARED_LIBS=ON -D BUILD_MPI=ON -D PKG_MANYBODY=ON -D PKG_USER-MISC=ON -D PKG_USER-PACE=ON ../cmake
make # -j${NUM_CPUS}
cp liblammps${SHLIB_EXT}* ../src
cd ../src
make install-python 
```
The combined libraries **should be available** on the system or the environment paths. Please see [here](https://lammps.sandia.gov/doc/Python_install.html) for more details.

## Documentation

- [Inputfile](examples/inputfile.md)
- [Command line interface](examples/cli.md)
- Examples
    - [Cu](examples/Cu_EAM/cu_example.ipynb): Example demonstrating the calculation of melting temperature for Copper.
    - [Fe](examples/Fe_EAM/fe_example.ipynb): Calculation of free energy of BCC Fe and calculating the phase transformation temperature for bcc->fcc.
    - [Ti](examples/Ti_EAM/ti_example.ipynb): Example demonstrating the calculation of temperature for the hcp->bcc transformation in Ti.
    - [CuZr](examples/ZrCu_EAM/zrcu_example.ipynb): Example demonstarting the calculation of free energy for pure Cu and CuZr B2 structure)
