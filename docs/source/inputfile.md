
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
```{grid-item} [](n_print_steps)
```
```{grid-item} [](potential_file)
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
```{grid-item} [](modules)
```
```{grid-item} [](options)
```
````

### `tolerance` 

````{grid} 1 2 3 4
:outline:
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

(lattice)=
#### `lattice`

_type_: string or list of strings \
_default_: None \
_example_:
```
lattice: FCC
lattice: [FCC, LQD]
lattice: [FCC, conf.data]
```
   
Lattice to be used for the calculations. The `lattice` option can use either LAMMPS for creation of input structure or use an input file in the LAMMPS data format. To use LAMMPS to create the structure, the keyword specified should be from the following: `bcc`, `fcc`, `hcp`, `diamond`, and `sc`. Lattice creation can **only be used for single species**. The preferred method is to provide a LAMMPS data file can be specified which contains the configuration.

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
```

The [pair style](https://lammps.sandia.gov/doc/pair_style.html) command for LAMMPS. All styles supported in LAMMPS can be used with calphy. Except for mode `alchemy`, only the first value in the list will be used. For mode `alchemy`, there should be two pair styles specified, and the alchemical transformation is carried out between the two.

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

(potential_file)=
#### `potential_file`

_type_: string \
_default_: None \
_example_:
```
pair_coeff: "/home/calc/potential.inp"
```

If specified, the `pair_style` and `pair_coeff` commands are not used, but rather the potential is read in from the provided input file using `include` command in LAMMPS. This allows the use of more complex or multiple potential files. Due to the `hybrid/scaled` styles employed in calphy, **this option only works with mode `fe` and `reference_phase` solid.**

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
_default_: True \
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
  modules:
    - anaconda/4
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

(modules)=
#### `modules`         

_type_: list of strings \
_example_:
```
modules:
  - anaconda
  - lammps
```

List of modules to be loaded before running the calculations. The given module names will be prefixed by `module load`.

---

(options)=
#### `options`         

_type_: string \
_example_:
```
options:
  - memory: 3GB
```

Extra options to be added to the submission script.

---
---

## `tolerance` block

This block helps to tune the internal convergence parameters that `calphy` uses. Generally, tuning these parameters are not required.

```
tolerance:
   spring_constant: 0.01
   solid_fraction: 0.7
   liquid_fraction: 0.05
   pressure: 0.5
```

---

(tol_spring_constant)=
#### `spring_constant`

_type_: float \
_default_: 0.01 \
_example_:
```
spring_constant: 0.01
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
_default_: 0.5 \
_example_:
```
pressure: 0.5
```

tolerance for the convergence of pressure.

---
---

## `melting_temperature` block

This block contains keywords that are used only for the mode `melting_temperature`.

```
melting_temperature:
   step: 200
   attempts: 5
```

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
     - Cu: 513
     - Zr: 511
```

---

(output_chemical_composition)=
#### `output_chemical_composition`

_type_: list \
_example_:
```
output_chemical_composition:
   - Cu: 513
   - Zr: 511
```

The output chemical composition in number of atoms. The total number of atoms should be equal to the input provided.

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