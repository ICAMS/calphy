
# Installation

### Supported operating systems

`calphy` can be installed on Linux and Mac OS based systems. On Windows systems, it is recommended to use  [Windows subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install).

### Normal installation

```{tab} Conda
`conda install -c conda-forge calphy`
```

```{tab} pip
`pip install calphy`
```

```{tab} from repository
`git clone https://github.com/ICAMS/calphy.git`  
`cd calphy`  
`pip install .`
```

```{tab} Singularity
A singularity containter can be used for running calphy locally or on HPC machines.
The containerised environment contains all of the packages required to run calphy.

### Downloading the container

The containerised environemnt can be pulled from the repository with:       
`singularity pull --arch amd64 library://sebastianhavens/calphy/calphy:latest`


[Calphy image repository](https://cloud.sylabs.io/library/sebastianhavens/calphy/calphy)

### Ruinning jobs using the container

On HPC machines you can usually load the singularity module with:
`module load singularity` if it is not already available.

You can initiate calculations using this container with the following line:

`singularity exec --bind $PWD --pwd $PWD {location_of_.sif_image} calphy_kernel -i {input_file} -k 0`

where `{location_of_.sif_image}` is the file location of the containerlised environment you just pulled and `{input_file}` is the name of your input file.
This line can be placed within a slurm script.

In the calphy input file, the scheduler should be set to local.<br>
For parallel calculations to run effectively, the OpenMPI module on the host system must be at least 4.1.2.

---
```


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
pip install .
```
The environment is now set up to run calphy.

### Dependencies

calphy requires Python ≥ 3.10 and the following packages (all installed automatically when using `pip install` or the supplied conda environment file):

- [LAMMPS](https://www.lammps.org/) compiled as a Python library (`conda install -c conda-forge lammps`)
- `numpy >= 2`
- `scipy`
- `matplotlib`
- `pyyaml`
- `tqdm`
- `pydantic`
- `mendeleev`
- [`pyscal3`](https://pyscal.org/)

#### Optional

- `pytest >= 7` for running the test-suite
- `mp_api` for fetching structures from the [Materials Project](https://materialsproject.org/)

### About [LAMMPS](https://www.lammps.org/) for `calphy`

calphy uses LAMMPS as the driver for molecular dynamics simulations. For calphy to work, LAMMPS needs to be compiled as a library along with the Python interface. The easiest way to do this is to install LAMMPS through the conda-forge channel using:

```
conda install -c conda-forge lammps
```

Alternatively, when interatomic potentials with special compilation needs are to be used, a recent stable release of LAMMPS can be compiled manually using the following set of instructions.

In order to help with installing all the prerequisites, an environment file which does not include the LAMMPS distribution is also provided. **This is only required if you want use a conda environment.** This environment can be installed using:

```
cd calphy
conda env create -f environment-nolammps.yml
```

Activate the environment using:

```
conda activate calphy2
```

Obtain the stable version from [here](https://github.com/lammps/lammps/releases) and extract the archive. From the extracted archive, the following steps, used in the [conda-forge recipe](https://github.com/conda-forge/lammps-feedstock/blob/master/recipe/build.sh) can be run:

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

(**Optional**) The above commands only builds the [MANYBODY](https://docs.lammps.org/Packages_details.html#pkg-manybody) package. To use some of the other potentials, the following commands could be added to the `cmake` call.

- `-D PKG_ML-PACE=ON` for performant [Atomic Cluster Expansion](https://docs.lammps.org/Packages_details.html#pkg-ml-pace) potential (from [June 2022 version](https://github.com/lammps/lammps/releases/tag/patch_27Oct2021)).
- `-D PKG_ML-SNAP=ON`for [SNAP potential](https://docs.lammps.org/Packages_details.html#pkg-ml-snap).
- `-D PKG_MEAM=ON` for [MEAM potential](https://docs.lammps.org/Packages_details.html#meam-package).
- `-D PKG_KIM=ON` for [KIM support](https://docs.lammps.org/Packages_details.html#pkg-kim).

Install the python wrapper:

```
cd ../src
make install-python
```

**In the case of a conda environment**, the following commands can be used to copy the compiled libraries to an accessible path (sometimes `PREFIX` needs to be used instead of `CONDA_PREFIX`):

```
mkdir -p $CONDA_PREFIX/include/lammps
cp library.h $CONDA_PREFIX/include/lammps
cp liblammps${SHLIB_EXT}* $CONDA_PREFIX/lib/
```

Once LAMMPS is compiled and the libraries are available in an accessible location, the following commands can be used within python to test the installation:

```
from lammps import lammps
lmp = lammps()
```

Finally install calphy:

```
cd calphy
pip install .
```









