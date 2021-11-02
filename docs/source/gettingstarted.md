# calphy

Python library and command line tool for calculation of free energies.

## Installation

### Supported operating systems

`calphy` can be installed on Linux and Mac OS based systems. On Windows systems, it is recommended to use  [Windows subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install).

### Installation using [Conda](https://anaconda.org/)

calphy can be installed directly using [Conda](https://docs.conda.io/en/latest/) from the [conda-forge channel](https://conda-forge.org/) by the following statement:

```
conda install -c conda-forge calphy
```


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

### 