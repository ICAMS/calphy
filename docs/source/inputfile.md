
# `calphy` input keywords


The inputfile is `yaml` formatted. In this section the possible keys in the inputfile is discussed. The table below gives a quick overview of all available keywords in calphy.

A sample file looks like this:

```
calculations:
- key1: value1
  key2: value2
  md:
    key1: value1
    key2: value2
```


## `calculations`


````{grid} 1 2 3 4
:outline:
```{grid-item} [](element)
```
```{grid-item} [](mass)
```
```{grid-item} [](mode)
```
```{grid-item} [](lattice)
```
```{grid-item} [](reference_phase)
```
```{grid-item} [](temperature)
```
```{grid-item} [](pressure)
```
```{grid-item} [](pressure_coupling)
```
```{grid-item} [](temperature_high)
```
```{grid-item} [](lattice_constant)
```
```{grid-item} [](repeat)
```
```{grid-item} [](n_iterations)
```
```{grid-item} [](n_switching_steps)
```
```{grid-item} [](n_equilibration_steps)
```
```{grid-item} [](pair_style)
```
```{grid-item} [](pair_coeff)
```
```{grid-item} [](pair_mode)
```
```{grid-item} [](n_print_steps)
```
```{grid-item} [](potential_file)
```
```{grid-item} [](fix_potential_path)
```
```{grid-item} [](file_format)
```
```{grid-item} [](spring_constants)
```
```{grid-item} [](equilibration_control)
```
```{grid-item} [](melting_cycle)
```
```{grid-item} [](folder_prefix)
```
```{grid-item} [](script_mode)
```
```{grid-item} [](lammps_executable)
```
```{grid-item} [](mpi_executable)
```
```{grid-item} [](npt)
```
```{grid-item} [](phase_name)
```
```{grid-item} [](reference_composition)
```
````

### `md`

````{grid} 1 2 3 4
:outline:
```{grid-item} [](timestep)
```
```{grid-item} [](n_small_steps)
```
```{grid-item} [](n_every_steps)
```
```{grid-item} [](n_repeat_steps)
```
```{grid-item} [](n_cycles)
```
```{grid-item} [](thermostat_damping)
```
```{grid-item} [](barostat_damping)
```
```{grid-item} [](cmdargs)
```
```{grid-item} [](init_commands)
```
```` 

### `queue` 

````{grid} 1 2 3 4
:outline:
```{grid-item} [](scheduler)
```
```{grid-item} [](cores)
```
```{grid-item} [](jobname)
```
```{grid-item} [](walltime)
```
```{grid-item} [](queuename)
```
```{grid-item} [](memory)
```
```{grid-item} [](commands)
```
```{grid-item} [](options)
```
````

### `tolerance` 

````{grid} 1 2 3 4
:outline:
```{grid-item} [](tol_lattice_constant)
```
```{grid-item} [](tol_spring_constant)
```
```{grid-item} [](tol_solid_fraction)
```
```{grid-item} [](tol_liquid_fraction)
```
```{grid-item} [](tol_pressure)
```
````

### `melting_temperature` 

````{grid} 1 2 3 4
:outline:
```{grid-item} [](mt_guess)
```
```{grid-item} [](step)
```
```{grid-item} [](attempts)
```
````

### `composition_scaling` 

````{grid} 1 2 3 4
:outline:
```{grid-item} [](output_chemical_composition)
```
```{grid-item} [](compscaling_restrictions)
```
````

### `monte_carlo` 

````{grid} 1 2 3 4
:outline:
```{grid-item} [](mc_n_steps)
```
```{grid-item} [](mc_n_swaps)
```
```{grid-item} [](mc_forward_swap_types)
```
```{grid-item} [](mc_reverse_swap_types)
```
```{grid-item} [](mc_allow_all_swaps)
```
```{grid-item} [](mc_use_custom_lammps)
```
````

### `uhlenbeck_ford_model` 

````{grid} 1 2 3 4
:outline:
```{grid-item} [](ufm_p)
```
```{grid-item} [](ufm_sigma)
```
````

### `materials_project` 

````{grid} 1 2 3 4
:outline:
```{grid-item} [](mp_api_key)
```
```{grid-item} [](mp_conventional)
```
```{grid-item} [](mp_target_natoms)
```
````

### `nose_hoover` 

````{grid} 1 2 3 4
:outline:
```{grid-item} [](nose_hoover_thermostat_damping)
```
```{grid-item} [](nose_hoover_barostat_damping)
```
````

### `berendsen` 

````{grid} 1 2 3 4
:outline:
```{grid-item} [](berendsen_thermostat_damping)
```
```{grid-item} [](berendsen_barostat_damping)
```
````

---
---

## `calculations` block

All other options are specified in blocks. The main block is the `calculations` block. This block can list all the calculations that `calphy` should perform and it can have more than one entry. A sample calculation block is shown below.

```
calculations:
- mode: ts 
  temperature: [1300, 1400]
  pressure: [0]
  lattice: [fcc]
  repeat: [2, 2, 2]
  reference_phase: solid
  n_iterations: 1
```

The various keys are-

---

(element)=
#### `element`

_type_: string/list of strings \
_default_: None \
_example_:
```
element: 'Cu'
element: ['Cu', 'Zr']
```
Chemical symbol(s) of the element(s) in the simulation.   

---

(mass)=
#### `mass`

_type_: float/list of floats \
_default_: 1.0 \
_example_:
```
mass: 63.546
mass: [63.546, 91.224]
```

Mass of the element(s) in the simulation. It should follow the same order as that of `element`, and should be of the same length.


(mode)=
#### `mode`  

_type_: string, `fe` or `ts` or `alchemy` or `melting_temperature` or `tscale` or `pscale` or `composition_scaling` \
_default_: None \
_example_:
```
mode: fe
mode: ts
```

Calculation mode. A small description of the different modes are given below.   
- `fe` performs a direct free energy calculation
- `ts` performs a direct free energy calculation followed by reversible scaling to find temperature dependence.
- `alchemy` is used for switching between two different interatomic potentials, or for integration over concentration. 
- `melting_temperature` can be used for automated calculation of melting temperature. 
- `tscale` is used for similar purpose as `ts`, but scales the temperature directly.
- `pscale` calculates the free energy as a function of the pressure.
- `composition_scaling` performs integration over composition.

---

(temperature)=
#### `temperature`  

_type_: float or list of floats \
_default_: None \
_example_:
```
temperature: 1200
temperature: [1200, 1300]
```

Temperatures for the simulation in Kelvin. The way temperature is used in `calphy` depends on the selected mode of calculation. 
- mode `ts` or `tscale`, a temperature sweep is carried out. In that case, only two values of temperature should be specified.
- mode `melting_temperature`, if provided it is used as an initial guess temperature. If not, the experimental melting temperature is used as a guess value.
- all other modes, a calculation is launched for each temperature on the list.

---

(pressure)=
#### `pressure`  

_type_: None or float or list of floats \
_default_: 0 \
_example_:
```
pressure: None
pressure: 0
pressure: [0,0,0]
pressure: [0, 10000]
pressure: [[1,1,1], [1000, 1000, 1000]]
```

Pressure for the simulation in bars. Depending on the pressure input, other options are automatically set. These options are `iso` and `fix_lattice`. `iso` decides if the barostat relaxes the system isotropically or anisotropically. `fix_lattice` does not relax the lattice at all. A table consistning of possible pressure input options and their effect is below:

| `pressure` input options  | | | | | |
| :-: | :-: | :-: | :-: | :-: | :-: |
| Type | Example | `iso` | `fix_lattice` | `mode` | Restrictions |
| `None` | `pressure: None` | True | True | `fe`, `alchemy` | None |
| scalar | `pressure: 100` | True | False | all except `pscale` | None |
| `(1,)` | `pressure: [100]` | True | False | all except `pscale` | None |
| `(2,)` | `pressure: [100, 200]` | True | False | `pscale` | None |
| `(3,)` | `pressure: [100, 100, 100]` | False | False | all except `pscale` | px=py=pz |
| `(1,3)` | `pressure: [[100, 100, 100]]` | False | False | all except `pscale` | px=py=pz |
| `(2,3)` | `pressure: [[100, 100, 100], [200, 200, 200]]` | False | False | `pscale` | px=py=pz |

---

(pressure_coupling)=
#### `pressure_coupling`

_type_: string \
_default_: `iso` for scalar/1-D pressure, `aniso` for 3-component pressure \
_example_:
```
pressure_coupling: iso
pressure_coupling: aniso
pressure_coupling: tri
```

Explicitly sets the LAMMPS barostat coupling style. When provided, this overrides the value that calphy would otherwise infer from the shape of the `pressure` input. The accepted values map directly to the LAMMPS `fix npt` keywords:

- `iso` — all three box dimensions are scaled equally (isotropic coupling).
- `aniso` — each box dimension is scaled independently, but off-diagonal (tilt) components are held fixed.
- `tri` — all six cell components (three lengths and three tilt factors) are allowed to change independently. **Required for triclinic crystal systems.** The simulation box must already be triclinic (i.e. the input structure must have non-zero tilt factors).

```{note}
For triclinic systems, provide a LAMMPS data file as the `lattice` input and set `pressure_coupling: tri`.
```

.. versionadded:: TBD

---

(lattice)=
#### `lattice`

_type_: string or list of strings \
_default_: None \
_example_:
```
lattice: FCC
lattice: [FCC, conf.data]
```
   
Lattice to be used for the calculations. The `lattice` option supports three forms:

- A built-in unit cell name: one of `bcc`, `fcc`, `hcp`, `diamond`, or `sc`. calphy will build the cell with [pyscal3](https://pyscal.org/) using `lattice_constant` and `repeat`. **Single-species only.**
- A path to a LAMMPS data file (`file_format: lammps-data`). This is the preferred form for multi-element systems and for triclinic cells.
- A [Materials Project](https://materialsproject.org/) ID, e.g. `mp-30`. The structure is fetched via `mp_api` and replicated to roughly `materials_project.target_natoms`. Requires the [`materials_project`](mp_api_key) block to be configured.

---

(reference_phase)=
#### `reference_phase`         
   
_type_: string or list of strings \
_default_: None \
_example_:
```
state: solid
state: [solid, liquid]
```

The protocol to be used for the calculation. The `reference_phase` input is closely related to the `lattice` command. It should be of the same length as the `lattice` input. For each of the lattice specified, `reference_phase` command specifies which reference state to use.

---

(lattice_constant)=
#### `lattice_constant`

_type_: float or list of floats \
_default_: experimental lattice constant \
_example_:
```
lattice_constant: 3.68
lattice_constant: [3.68, 5.43]
```

Lattice constant for input structures. Lattice constant values to be used for initial structure creation. Only required if the structure lattice is specified. If not specified, experimental lattice constant is used.

---

(repeat)=
#### `repeat`

_type_: list of ints of length 3 \
_default_: [1, 1, 1] \
_example_:
```
repeat: [3,3,3]
```

`repeat` command specifies the number of unit cells required in each x, y and z directions. This is only used if structure lattice is specified.

---

(n_iterations)=
#### `n_iterations`         

_type_: int \
_default_: 1 \
_example_:
```
n_iterations: 3
```

The number of backward and forward integrations to be carried out in all modes. If more than one integration cycle is used, the errors can also be evaluated.

---

(temperature_high)=
#### `temperature_high`

_type_: float \
_default_: 2*temperature \
_example_:
```
temperature_high: 1600
```

The temperature used to melt a solid structure to create a liquid. If `reference_phase` is chosen as `liquid`, calphy performs a melting cycle to create an equilibrated liquid structure. calphy starts from the given input structure, and heats it using 2 times the highest temperature provided in the `temperature` option. If the structure is not melted, the temperature is increased progressively. `temperature_high` keyword can be used to set the temperature to overheat the structure and melt it.   

---

(npt)=
#### `npt`

_type_: bool \
_default_: True \  
_example_:
```
npt: True
```

`npt` determines if calculations are carried out in the NPT ensemble. This option is used with modes `alchemy`, and `ts`. The effect is described below:

- for mode `ts`: the reversible scaling approach is carried out in NPT if `npt` is True, otherwise the `NVT` ensemble is used.

- for mode `alchemy`: If `npt` is False, the NVT ensemble is used, meaning that the calculated work during alchemical transformation is calculated on the equilibrated volume of the first pair style. 

---

(pair_style)=
#### `pair_style`

_type_: string or list of strings \
_default_: None \
_example_:
```
pair_style: eam/alloy
pair_style: [eam/alloy, eam/alloy]
pair_style: eam/fs
pair_style: pace
pair_style: ["coul/long 11.0", "eam/alloy"]
```

The [pair style](https://lammps.sandia.gov/doc/pair_style.html) command for LAMMPS. All styles supported in LAMMPS can be used with calphy. Except for mode `alchemy` and `pair_mode: overlay`, only the first value in the list will be used. For mode `alchemy`, there should be two pair styles specified, and the alchemical transformation is carried out between the two. For `pair_mode: overlay`, list the component pair styles and include any pair-style options with the corresponding component.

---

(pair_coeff)=
#### `pair_coeff`

_type_: string or list of strings \
_default_: None \
_example_:
```
pair_coeff: "* * Cu_EAM/Cu01.eam.alloy Cu"
pair_coeff: ["* * Cu_EAM/Cu01.eam.alloy Cu", "* * Cu_EAM/Cu02.eam.alloy Cu"]
pair_coeff: "* * CuZr_EAM/CuZr.eam.fs Cu Zr"
```

The [pair coeff](https://lammps.sandia.gov/doc/pair_coeff.html) command for LAMMPS. It should be the same length as `pair_style`. Except for mode `alchemy`, only the first value in the list will be used. For mode `alchemy`, there should be two pair styles specified, and the alchemical transformation is carried out between the two.

---

(pair_mode)=
#### `pair_mode`

_type_: string \
_default_: None \
_example_:
```
pair_mode: overlay
pair_style: ["coul/long 11.0", "eam/alloy"]
pair_coeff:
  - "* * coul/long"
  - "* * eam/alloy /path/to/myfile elt1 elt2"
```

Set `pair_mode: overlay` to combine multiple component pair styles as a single physical potential using LAMMPS `hybrid/overlay`. The `pair_style` values should be the component styles, not the full `hybrid/overlay` command. The `pair_coeff` values may either include the component style name, as in native LAMMPS hybrid syntax, or omit it when the list order matches `pair_style`. In liquid free-energy integration and reversible scaling, calphy flattens this overlay into `hybrid/scaled` so all real-potential components are scaled together.

---

(potential_file)=
#### `potential_file`

_type_: string \
_default_: None \
_example_:
```
potential_file: "/home/calc/potential.inp"
```

```{deprecated}
`potential_file` is deprecated and will be removed in a future release. Use `pair_style` and `pair_coeff` (with [`pair_mode: overlay`](pair_mode) for multi-potential setups) instead.
```

If specified, the `pair_style` and `pair_coeff` commands are not used, but rather the potential is read in from the provided input file using the `include` command in LAMMPS. Due to the `hybrid/scaled` styles employed in calphy, **this option only works with mode `fe` and `reference_phase: solid`.**

---

(fix_potential_path)=
#### `fix_potential_path`

_type_: bool \
_default_: True \
_example_:
```
fix_potential_path: True
```

If True (the default), calphy expands `~`, environment variables (e.g. `$USER`, `${HOME}`), and converts the potential filename in each `pair_coeff` entry to an absolute path before it is passed to LAMMPS. Set this to `False` if you want to pass `pair_coeff` strings to LAMMPS verbatim.

---

(file_format)=
#### `file_format`

_type_: string \
_default_: `lammps-data` \
_example_:
```
file_format: lammps-data
```

File format used to read the input structure when `lattice` is a path to a file on disk. Currently only `lammps-data` is supported.

---

(n_equilibration_steps)=
#### <a name="n_equilibration_steps"></a>`n_equilibration_steps`

_type_: int \
_default_: 25000 \
_example_:
```
n_equilibration_steps: 10000
```

The number of time steps for equilibrating the system.

---

(n_switching_steps)=
#### `n_switching_steps`

_type_: int or list of ints \
_default_: 50000 \
_example_:
```
n_switching_steps: 10000
n_switching_steps: [10000, 20000]
```

The number of switching steps. If a list of two integers is provided, the first value is used for mode `fe` while the second value will be used for all other modes.

---

(n_print_steps)=
#### `n_print_steps`        

_type_: int \
_default_: 0 \
_example_:
```
n_print_steps: 100
```

Record MD trajectory during temperature sweep runs in the given interval of time steps. Default 0, trajectories are never recorded.


---

(spring_constants)=
#### `spring_constants`        

_type_: list of floats \
_default_: None \
_example_:
```
spring_constants: 1.2
spring_constants: [1.2, 1.3]
```

Spring constants for Einstein crystal. If specified, the automatic calculation is not performed. Should be equal to the number of species in the system.


---

(equilibration_control)=
#### `equilibration_control`        

_type_: str \
_default_: None \
_example_:
```
equilibration_control: berendsen
equilibration_control: nose-hoover
```

The barostat and thermostat combination to be used for the equilibration stage. By default, Berendsen will be used for solid reference and Nose-Hoover will be used for liquid. The damping parameters can be tuned using the [nose_hoover](nose_hoover_block) block or the [berendsen](berendsen_block) block. Used only for equilibration stage.

---

(melting_cycle)=
#### `melting_cycle`        

_type_: bool \
_default_: False \
_example_:
```
melting_cycle: True
melting_cycle: False
```

If True, a melting cycle is carried out to melt the given input structure. Only used if the `reference_phase` is `"liquid"`.


---

(folder_prefix)=
#### `folder_prefix`        

_type_: string \
_default_: None \
_example_:
```
folder_prefix: set1
```  

Prefix string to be added to folder names for calculation. Folders for calculations in calphy are named as `mode-lattice-temperature-pressure`. Therefore, if more than one calculation is run with the same parameters, they will be overwritten. To prevent this, `folder_prefix` can be used. If `folder_prefix` is provided, the folders will be named as `folder_prefix-mode-lattice-temperature-pressure`.


---

(script_mode)=
#### `script_mode`        

_type_: bool \
_default_: False \
_example_:
```
script_mode: False
```  

If True, a LAMMPS executable script is written and executed instead of the library interface of LAMMPS.
Works only with `reference_phase: solid`, and `mode: fe`.
Needs specification of [`lammps_executable`](lammps_executable) and [`mpi_executable`](mpi_executable).


---

(lammps_executable)=
#### `lammps_executable`        

_type_: string \
_default_: None \
_example_:
```
lammps_executable: lmp_mpi
```  

LAMMPS executable to run the calculations with.
Works only with `reference_phase: solid`, and `mode: fe`.
Works only if [`script_mode`](script_mode) is `True`.



---

(mpi_executable)=
#### `mpi_executable`        

_type_: string \
_default_: None \
_example_:
```
mpi_executable: mpiexec
```  

MPI executable to run the LAMMPS with.
Works only with `reference_phase: solid`, and `mode: fe`.
Works only if [`script_mode`](script_mode) is `True`.


---

(phase_name)=
#### `phase_name`

_type_: string \
_default_: "" \
_example_:
```
phase_name: bcc_AB
```

Label attached to a calculation when it is part of a phase-diagram workflow (see the [`phase_diagram`](phase_diagram_block) module). Calculations belonging to the same physical phase share a `phase_name`, and post-processing utilities (`gather_results`, `find_transition_temperature`, …) group rows by this column.

---

(reference_composition)=
#### `reference_composition`

_type_: float \
_default_: 0.0 \
_example_:
```
reference_composition: 0.25
```

For composition sweeps in the phase-diagram workflow this records the *target* composition that the structure has been transformed to (used by composition-scaling integration). Normally set automatically by the phase-diagram driver and rarely supplied by hand.

---
---

## `md` block

MD block consists of the various options required for the MD runs. An example block is given below and the keys are discussed.

```
md:
  timestep: 0.001
  thermostat_damping: 0.1
  barostat_damping: 0.1
```

---

(timestep)=
#### `timestep`        

_type_: float \
_default_: 0.001 \
_example_:
```
timestep: 0.001
```

The timestep for MD in picoseconds. 

---

(thermostat_damping)=
#### `thermostat_damping`

_type_: float \
_default_: 0.1 \
_example_:
```
thermostat_damping: 0.1
```

The thermostat damping constant for MD. For damping during equilibration stage see [`equilibration_control`](equilibration_control).

---

(barostat_damping)=
#### <a name="barostat_damping"></a>`barostat_damping`

_type_: float \
_default_: 0.1 \
_example_:
```
barostat_damping: 0.1
```

Pressure damping for MD. For damping during equilibration stage see [`equilibration_control`](equilibration_control).

---

(n_small_steps)=
#### `n_small_steps`

_type_: int \
_default_: 10000 \
_example_:
```
n_small_steps: 10000
```

The number of time steps for equilibration cycles to calculate spring constant and average volume.

---

(n_every_steps)=
#### `n_every_steps`         

_type_: int \
_default_: 10 \
_example_:
```
n_every_steps: 100
```

Keywords to tune how often average values are recorded in LAMMPS. Please see [here](https://docs.lammps.org/fix_ave_time.html) for more details.

---

(n_repeat_steps)=
#### `n_repeat_steps`         

_type_: int \
_default_: 10 \
_example_:
```
n_repeat_steps: 10
```

Keywords to tune how often average values are recorded in LAMMPS. Please see [here](https://docs.lammps.org/fix_ave_time.html) for more details.

---

(n_cycles)=
#### `n_cycles`        

_type_: int \
_default_: 100 \
_example_:
```
n_cycles: 100
```

Number of cycles to try converging the pressure of the system. If the pressure is not converged after `n_cycles`, an error will be raised. In each `n_cycle`, `n_small_steps` MD steps will be run.


---

(cmdargs)=
#### `cmdargs`

_type_: string \
_default_: "" \
_example_:
```
md:
  cmdargs: "-k on g 1 -sf kk -pk kokkos newton on neigh half"
```

Extra command-line arguments passed verbatim to the LAMMPS Python library when the LAMMPS object is created. Most commonly used to enable accelerator packages such as `KOKKOS`, `GPU` or `INTEL`.

---

(init_commands)=
#### `init_commands`

_type_: list of strings \
_default_: None \
_example_:
```
init_commands:
  - timestep 0.002
  - atom_style charge
  - neighbor 0.6 bin
```

Provides the possibility to replace or add initial commands when the LAMMPS object is initialised. If the command is already used in calphy, for example `timestep` or `atom_style` they will be replaced. If it is a new command, it will be added. This commands receive higher priority than the ones that already exist. For examples if you provide `timestep: 0.002` in the `md` block, and `timestep 0.004` in `init_commands`, the timestep used would be 0.004.

---
---

## `queue` block

This block consists of the options for specifying the scheduler for carrying out the calculations. An example block is given below-

```
queue:
  scheduler: local
  cores: 2
  jobname: ti
  walltime: "23:50:00"
  queuename: shorttime
  memory: 3GB
  commands:
    - conda activate env
```

---

(scheduler)=
#### `scheduler`

_type_: string \
_default_: local \
_example_:
```
scheduler: slurm
```

The scheduler to be used for the job. Can be `local`, `slurm` or `sge`. The code has been tested only on local and slurm.

---

(cores)=
#### `cores`

_type_: int \
_default_: 1 \
_example_:
```
cores: 4
```

The number of cores to be used for the job.

---

(jobname)=
#### `jobname`         

_type_: string \
_default_: calphy \
_example_:
```
jobname: cu
```

Name of the job. Not used for `local`.

---

(walltime)=
#### `walltime`         

_type_: string \
_example_:
```
walltime: "23:50:00"
```

Walltime for the job. Not used for `local`.

---

(queuename)=
#### `queuename`         

_type_: string \
_example_:
```
queuename: "shorttime"
```

Name of the queue. Not used for `local`.

---

(memory)=
#### `memory`         

_type_: string \
_example_:
```
memory: 3GB
```

Memory to be used per core. Not used for `local`.

---

(commands)=
#### `commands`         

_type_: list of strings \
_example_:
```
command:
 - source .bashrc
 - conda activate ace
 - module load lammps
```

Command that will be run **before** the actual calculations are carried out. This section can be used to specify commands that need to be run before the actual calculation. If the calculations are run within a conda environment, the activate command for conda should be specified here. If additional modules need to be loaded, that can also be specified here.

---

List of commands to be run before the main calculation command. These commands will be executed in the submission script.

---

(options)=
#### `options`         

_type_: string \
_example_:
```
options:
  memory: 3GB
```

Extra options to be added to the submission script.

---
---

## `tolerance` block

This block helps to tune the internal convergence parameters that `calphy` uses. Generally, tuning these parameters are not required.

```
tolerance:
   lattice_constant: 0.0002
   spring_constant: 0.1
   solid_fraction: 0.7
   liquid_fraction: 0.05
   pressure: 10.0
```

---

(tol_lattice_constant)=
#### `lattice_constant`

_type_: float \
_default_: 0.0002 \
_example_:
```
lattice_constant: 0.0002
```

Convergence tolerance (in Å) on the average lattice constant during the equilibration / volume-averaging stage.

---

(tol_spring_constant)=
#### `spring_constant`

_type_: float \
_default_: 0.1 \
_example_:
```
spring_constant: 0.1
```

tolerance for the convergence of spring constant calculation.

---

(tol_solid_fraction)=
#### `solid_fraction`

_type_: float \
_default_: 0.7 \
_example_:
```
solid_fraction: 0.7
```

The minimum amount of solid particles that should be there in solid.

---

(tol_liquid_fraction)=
#### `liquid_fraction`

_type_: float \
_default_: 0.05 \
_example_:
```
liquid_fraction: 0.05
```

Maximum fraction of solid atoms allowed in liquid after melting.

---

(tol_pressure)=
#### `pressure`

_type_: float \
_default_: 10.0 \
_example_:
```
pressure: 10.0
```

Tolerance (in bars) for the convergence of the average pressure during the equilibration / volume-averaging stage.

---
---

## `melting_temperature` block

This block contains keywords that are used only for the mode `melting_temperature`.

```
melting_temperature:
   guess: 1300
   step: 200
   attempts: 5
```

---

(mt_guess)=
#### `guess`

_type_: float \
_default_: None (uses experimental melting point of `element[0]` from `mendeleev`) \
_example_:
```
guess: 1300
```

Initial guess (in Kelvin) for the melting temperature search. If not provided, calphy uses the experimental melting point of the first element from the [mendeleev](https://mendeleev.readthedocs.io/) database. Only used if mode is `melting_temperature`.

---

(step)=
#### `step`

_type_: int \
_example_:
```
step: 200
```

Temperature interval for search of melting temperature. Only used if mode is `melting_temperature`.

---

(attempts)=
#### `attempts`

_type_: int \
_example_:
```
attempts: 5
```

The number of maximum attempts to try find the melting temperature in a automated manner. Only used if mode is `melting_temperature`.

---
---

## `composition_scaling` block

This block contains keywords that are used only for the mode `composition_scaling`.

```
composition_scaling:
  output_chemical_composition:
     Cu: 513
     Zr: 511
```

---

(output_chemical_composition)=
#### `output_chemical_composition`

_type_: list \
_example_:
```
output_chemical_composition:
   Cu: 513
   Zr: 511
```

The output chemical composition in number of atoms. The total number of atoms should be equal to the input provided.

---

(compscaling_restrictions)=
#### `restrictions`

_type_: list of strings \
_default_: `[]` \
_example_:
```
composition_scaling:
  output_chemical_composition:
    Al: 494
    Li: 2
    O: 3
    C: 1
  restrictions:
    - "Al-O"
```

List of `"A-B"` pair strings that **must not be transformed into one another** during composition scaling. In the example above, the algorithm is forbidden from converting any Al atom directly into an O atom (or vice-versa). If no valid transformation chain can be found that respects all restrictions, an error is raised.

---
---

(nose_hoover_block)=
## <a name="nose_hoover"></a>`nose_hoover` block

This block contains keywords that are used only for the equilibration stage if [`equilibration_control`](equilibration_control) is `nose_hoover`.

```
nose_hoover:
   thermostat_damping: 0.1
   barostat_damping: 0.1
```

---

(nose_hoover_thermostat_damping)=
#### `thermostat_damping`

_type_: float \
_default_: 0.1 \
_example_:
```
thermostat_damping: 0.1
```

The thermostat damping constant for equilibration MD. 

---

(nose_hoover_barostat_damping)=
#### `barostat_damping`

_type_: float \
_default_: 0.1 \
_example_:
```
barostat_damping: 0.1
```

Pressure damping for equilibration MD. 

---
---

(berendsen_block)=
## <a name="berendsen"></a>`berendsen` block

This block contains keywords that are used only for the equilibration stage if [`equilibration_control`](equilibration_control) is `berendsen`.

```
berendsen:
   thermostat_damping: 100.0
   barostat_damping: 100.0
```

---

(berendsen_thermostat_damping)=
#### `thermostat_damping`

_type_: float \
_default_: 100.0 \
_example_:
```
thermostat_damping: 100.0 
```

The thermostat damping constant for equilibration MD. 

---

(berendsen_barostat_damping)=
#### `barostat_damping`

_type_: float \
_default_: 100.0 \
_example_:
```
barostat_damping: 100.0  
```

Pressure damping for equilibration MD. 

---
---

## `monte_carlo` block

This block configures Monte Carlo particle-swap moves that can be interleaved with the molecular-dynamics integration during alchemical / composition-scaling calculations. By default no swaps are performed (`n_swaps: 0`).

```
monte_carlo:
  n_steps: 100
  n_swaps: 10
  forward_swap_types: [1, 2]
  reverse_swap_types: [1, 2]
  allow_all_swaps: True
  use_custom_lammps: False
```

---

(mc_n_steps)=
#### `n_steps`

_type_: int \
_default_: 1 \
_example_:
```
n_steps: 100
```

Attempt a batch of swap moves every `n_steps` MD steps.

---

(mc_n_swaps)=
#### `n_swaps`

_type_: int \
_default_: 0 \
_example_:
```
n_swaps: 10
```

Number of swap moves to attempt in each batch. The default of `0` disables Monte Carlo swaps entirely.

---

(mc_forward_swap_types)=
#### `forward_swap_types`

_type_: list of ints \
_default_: `[]` \
_example_:
```
forward_swap_types: [1, 2]
```

LAMMPS atom-type IDs that are eligible for swapping during the forward integration leg.

---

(mc_reverse_swap_types)=
#### `reverse_swap_types`

_type_: list of ints \
_default_: `[]` \
_example_:
```
reverse_swap_types: [1, 2]
```

LAMMPS atom-type IDs eligible for swapping during the backward integration leg.

---

(mc_allow_all_swaps)=
#### `allow_all_swaps`

_type_: bool \
_default_: True \
_example_:
```
allow_all_swaps: True
```

If True, swaps between any pair of allowed types are accepted, including the auxiliary "fictitious" types introduced for alchemical interpolation. Set to False to restrict swaps to only the real species.

---

(mc_use_custom_lammps)=
#### `use_custom_lammps`

_type_: bool \
_default_: False \
_example_:
```
use_custom_lammps: True
```

Use the custom LAMMPS build that ships with the modified `fix atom/swap` patch required for some swap topologies. Leave at the default unless you have built that LAMMPS variant.

---
---

## `uhlenbeck_ford_model` block

Parameters of the Uhlenbeck–Ford reference fluid used as the reference state for **liquid** free-energy calculations. The defaults reproduce the values used in the [original calphy paper](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.103801) and rarely need to be touched.

```
uhlenbeck_ford_model:
  p: 50.0
  sigma: 1.5
```

---

(ufm_p)=
#### `p`

_type_: float \
_default_: 50.0 \
_example_:
```
p: 50.0
```

Strength parameter of the Uhlenbeck–Ford pair potential. Larger values make the reference fluid more repulsive.

---

(ufm_sigma)=
#### `sigma`

_type_: float \
_default_: 1.5 \
_example_:
```
sigma: 1.5
```

Length parameter (in Å) of the Uhlenbeck–Ford pair potential.

---
---

## `materials_project` block

Configuration for fetching input structures directly from the [Materials Project](https://materialsproject.org/) using `mp_api`. Only used when `lattice` is given as a Materials Project ID such as `mp-30`.

```
materials_project:
  api_key: MP_API_KEY
  conventional: True
  target_natoms: 1500
```

---

(mp_api_key)=
#### `api_key`

_type_: string \
_default_: "" \
_example_:
```
api_key: MP_API_KEY
```

Name of an environment variable holding your Materials Project API key (for example `MP_API_KEY`). calphy reads `os.environ[api_key]` at runtime so the key itself is never written into input files or log output. Set the variable in your shell before running:

```
export MP_API_KEY="your_api_key_here"
```

---

(mp_conventional)=
#### `conventional`

_type_: bool \
_default_: True \
_example_:
```
conventional: True
```

If True, the conventional cell is used; if False, the primitive cell is used.

---

(mp_target_natoms)=
#### `target_natoms`

_type_: int \
_default_: 1500 \
_example_:
```
target_natoms: 1500
```

The structure parsed from Materials Project is replicated isotropically until it contains approximately this many atoms (ignored if `repeat` is set to anything other than `[1, 1, 1]`).

---
