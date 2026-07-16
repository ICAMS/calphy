
# Installation

### Supported operating systems

`calphy` can be installed on Linux and Mac OS based systems. On Windows systems, it is recommended to use  [Windows subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install).

### How calphy uses LAMMPS

calphy drives molecular dynamics with [LAMMPS](https://www.lammps.org/). As of calphy v2 it does this by default by **running an external `lmp` executable** as a subprocess — it no longer needs LAMMPS compiled as a Python library. (A live-session library backend through `pylammpsmpi` is still available as an opt-in: install `pip install calphy[library]` and set `execution_mode: library` in the input file; see [`execution_mode`](execution_mode).) In practice the default means:

- You **bring your own `lmp` binary** (from conda-forge, an HPC module, a container, or a manual build).
- calphy locates it at run time in this order:
  1. the `lammps_executable` key in the input file,
  2. the environment variable `$CALPHY_LAMMPS_EXECUTABLE`,
  3. `lmp` on your `PATH`.
- For parallel runs (`queue.cores > 1`) the MPI launcher is resolved the same way: `mpi_executable` → `$CALPHY_MPI_EXECUTABLE` → `mpirun` on `PATH`.

Before the first simulation, calphy runs a quick **preflight** check (`lmp -h`) and reports, with a clear message, any LAMMPS package a calculation needs but the binary does not provide. Set `CALPHY_SKIP_PREFLIGHT=1` to bypass it.

### Normal installation

```{tab} Conda
`conda install -c conda-forge calphy`

The conda-forge package pulls in a `lammps` build that provides the `lmp` binary, so a full setup works out of the box.
```

```{tab} pip
`pip install calphy`

`pip` installs the Python package only — **you still need an `lmp` binary** on your `PATH` (see *Installing LAMMPS* below).
```

```{tab} from repository
`git clone https://github.com/ICAMS/calphy.git`
`cd calphy`
`pip install .`

As with `pip`, provide an `lmp` binary separately.
```

```{tab} Singularity
A singularity container can be used for running calphy locally or on HPC machines.
The containerised environment contains all of the packages required to run calphy, including the `lmp` binary.

### Downloading the container

The containerised environment can be pulled from the repository with:
`singularity pull --arch amd64 library://sebastianhavens/calphy/calphy:latest`

[Calphy image repository](https://cloud.sylabs.io/library/sebastianhavens/calphy/calphy)

### Running jobs using the container

On HPC machines you can usually load the singularity module with:
`module load singularity` if it is not already available.

You can initiate calculations using this container with the following line:

`singularity exec --bind $PWD --pwd $PWD {location_of_.sif_image} calphy_kernel -i {input_file} -k 0`

where `{location_of_.sif_image}` is the file location of the containerised environment you just pulled and `{input_file}` is the name of your input file.
This line can be placed within a slurm script.

In the calphy input file, the scheduler should be set to local.<br>
For parallel calculations to run effectively, the OpenMPI module on the host system must be at least 4.1.2.

---
```

### Using a conda environment

It is **strongly** recommended to install and use `calphy` within a conda environment. To see how you can install conda see [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Once a conda distribution is available, clone the repository and create the environment from the included file:

```
git clone https://github.com/ICAMS/calphy.git
cd calphy
conda env create -f environment.yml
conda activate calphy
pip install .
```

`environment.yml` installs the conda-forge `lammps` package (which ships the `lmp` binary calphy drives) together with `openmpi` (for `mpirun`). If you would rather supply your own `lmp` binary — e.g. an HPC module or a custom build — use `calphy/environment-nolammps.yml` instead, which creates the same environment without LAMMPS.

### Dependencies

calphy requires Python ≥ 3.10 and the following packages (all installed automatically when using `pip install` or the supplied conda environment file):

- an external [LAMMPS](https://www.lammps.org/) **`lmp` binary** — *not* a pip/conda dependency; provide it yourself (see below)
- `numpy >= 2`
- `scipy`
- `matplotlib`
- `pyyaml`
- `tqdm`
- `pydantic >= 2`
- `mendeleev`
- [`pyscal3`](https://pyscal.org/)

#### Optional

- `pytest >= 7` for running the test-suite
- `mp_api` for fetching structures from the [Materials Project](https://materialsproject.org/)

## Installing LAMMPS

calphy needs an `lmp` executable that includes the LAMMPS packages its methods rely on. You do **not** need to build LAMMPS as a Python library.

### Which LAMMPS packages calphy needs

| calphy feature | LAMMPS package | provides |
|---|---|---|
| EAM / MEAM / most `pair_style`s | `MANYBODY` (and the relevant potential package) | your interatomic potential |
| solid free energy, `ts`, `tscale` | `EXTRA-FIX` | `fix ti/spring` (Frenkel–Ladd spring) |
| liquid free energy, `melting_temperature` | `EXTRA-PAIR` | `pair_style ufm`, `pair_style hybrid/scaled` |
| Monte-Carlo swaps (`monte_carlo.n_swaps > 0`) | `MC` | `fix atom/swap` |
| `mode: fe-qtb` (quantum thermal bath) | `QTB` | `fix qtb` |

The preflight check verifies these for each calculation and tells you exactly which package to add if one is missing.

### Option 1 — conda-forge (recommended)

```
conda install -c conda-forge lammps
```

This provides an `lmp` binary with the common packages (MANYBODY, EXTRA-FIX, EXTRA-PAIR, MC). If you need `mode: fe-qtb`, check that the build includes the QTB package (`lmp -h | grep qtb`); if not, use a build/module that has it or compile one (below).

### Option 2 — an HPC module or existing binary

If your cluster already provides LAMMPS:

```
module load lammps          # or your site's module name
```

and point calphy at it (any one of these):

```
# in the input file
lammps_executable: /path/to/lmp

# or in the environment
export CALPHY_LAMMPS_EXECUTABLE=/path/to/lmp
export CALPHY_MPI_EXECUTABLE=/path/to/mpirun   # only needed for cores > 1
```

### Option 3 — compile it yourself

Obtain a recent stable release from [the LAMMPS releases page](https://github.com/lammps/lammps/releases), extract it, and build the **executable** with the packages calphy needs:

```
cd lammps-*/            # extracted source
mkdir build && cd build
cmake -D BUILD_MPI=ON \
      -D PKG_MANYBODY=ON \
      -D PKG_EXTRA-FIX=ON \
      -D PKG_EXTRA-PAIR=ON \
      -D PKG_MC=ON \
      -D PKG_QTB=ON \
      ../cmake
make -j
```

This produces an `lmp` binary in the `build` directory. Put it on your `PATH` (or point `lammps_executable`/`$CALPHY_LAMMPS_EXECUTABLE` at it). Add further package flags for special potentials, for example:

- `-D PKG_ML-PACE=ON` for the [Atomic Cluster Expansion](https://docs.lammps.org/Packages_details.html#pkg-ml-pace) potential.
- `-D PKG_ML-SNAP=ON` for the [SNAP potential](https://docs.lammps.org/Packages_details.html#pkg-ml-snap).
- `-D PKG_MEAM=ON` for the [MEAM potential](https://docs.lammps.org/Packages_details.html#meam-package).
- `-D PKG_KIM=ON` for [KIM support](https://docs.lammps.org/Packages_details.html#pkg-kim).

### Checking the LAMMPS setup

Confirm calphy can find and use the binary:

```
which lmp          # or: echo $CALPHY_LAMMPS_EXECUTABLE
lmp -h | head      # should print the LAMMPS help / style listing
```

If `lmp -h` lists the packages in the table above, you are ready to run calphy.
