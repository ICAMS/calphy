
# `calphy` input file

The inputfile is `yaml` formatted. In this section the possible keys in the inputfile is discussed. The input file consists of two main keys, and four separate blocks. For a sample of the inputfile see the end of this document. The table below gives a quick overview of all available keywords in calphy.

| Main keys  | |
| :-: | :-: |
| [element](#element) | [mass](#mass) | 

| `calculations` block | | | | |
| :-: | :-: | :-: | :-: | :-: |
| [mode](#mode) | [lattice](#lattice) | [reference_phase](#reference_phase) | [temperature](#temperature) | [pressure](#pressure) |
| [temperature_high](#temperature_high) | [lattice_constant](#lattice_constant) | [repeat](#repeat) | [n_iterations](#n_iterations) | [n_switching_steps](#n_switching_steps) | 
| [n_equilibration_steps](#n_equilibration_steps) | [pair_style](#pair_style) | [pair_coeff](#pair_coeff) | [n_print_steps](#n_print_steps) | [potential_file](#potential_file) | 

| `md` block  | | | | |
| :-: | :-: | :-: | :-: | :-: |
| [timestep](#timestep) | [n_small_steps](#n_small_steps) | [n_every_steps](#n_every_steps) | [n_repeat_steps](#n_repeat_steps) | [n_cycles](#n_cycles) |
| [thermostat_damping](#thermostat_damping) | [barostat_damping](#barostat_damping) |

| `queue` block  | | | | |
| :-: | :-: | :-: | :-: | :-: |
| [scheduler](#scheduler) | [cores](#cores) | [jobname](#jobname) | [walltime](#walltime) | [queuename](#queuename) |
| [memory](#memory) | [commands](#commands) | [modules](#modules) | [options](#options) |

| `tolerance` block | | | | |
| :-: | :-: | :-: | :-: | :-: |
| [lattice_constant](#tol_lattice_constant) | [spring_constant](#tol_spring_constant) | [solid_fraction](#tol_solid_frac) | [liquid_fraction](#tol_liquid_frac) | [pressure](#tol_pressure) |

| `melting_temperature` block | | | | |
| :-: | :-: | :-: | :-: | :-: |
| [step](#step) | [attempts](#attempts) |

---
---

## main keys

---

#### <a name="element"></a>`element`

_type_: string/list of strings  
_default_: None  
_example_:
```
element: 'Cu'
element: ['Cu', 'Zr']
```
Chemical symbol(s) of the element(s) in the simulation.   

---

#### <a name="mass"></a>`mass`

_type_: float/list of floats
_default_: 1.0
_example_:
```
mass: 63.546
mass: [63.546, 91.224]
```

Mass of the element(s) in the simulation. It should follow the same order as that of `element`, and should be of the same length.

---
---

## `calculations` block

Other than these two main keys, all other options are specified in blocks. The first block is the `calculations` block. This block can list all the calculations that `calphy` should perform and it can have more than one entry. A sample calculation block is shown below.

```
calculations:
- mode: ts 
  temperature: [1300, 1400]
  pressure: [0]
  lattice: [FCC]
  repeat: [2, 2, 2]
  reference_phase: solid
  n_iterations: 1
```

The various keys are-

---

#### <a name="mode"></a>`mode`  

_type_: string, `fe` or `ts` or `mts` or `alchemy` or `melting_temperature` or `tscale` or `pscale`  
_default_: None  
_example_:
```
mode: fe
mode: ts
```

Calculation mode. A small description of the different modes are given below.   
- `fe` performs a direct free energy calculation
- `ts` performs a direct free energy calculation followed by reversible scaling to find temperature dependence.
- `mts` performs only reversible scaling part and can be used for dynamic Clausius-Clapeyron integration.
- `alchemy` is used for switching between two different interatomic potentials, or for integration over concentration. 
- `melting_temperature` can be used for automated calculation of melting temperature. 
- `tscale` is used for similar purpose as `ts`, but scales the temperature directly.
- `pscale` calculates the free energy as a function of the pressure.

---

#### <a name="temperature"></a>`temperature`  

_type_: float or list of floats  
_default_: None  
_example_:
```
temperature: 1200
temperature: [1200, 1300]
```

Temperatures for the simulation in Kelvin. The way temperature is used in `calphy` depends on the selected mode of calculation. 
- mode `ts` or `tscale` or `mts`, a temperature sweep is carried out. In that case, only two values of temperature should be specified.
- mode `melting_temperature`, if provided it is used as an initial guess temperature. If not, the experimental melting temperature is used as a guess value.
- all other modes, a calculation is launched for each temperature on the list.

---

#### <a name="pressure"></a>`pressure`  

_type_: None or float or list of floats  
_default_: None  
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

#### <a name="lattice"></a>`lattice`

_type_: string or list of strings  
_default_: None  
_example_:
```
lattice: FCC
lattice: [FCC, LQD]
lattice: [FCC, conf.data]
```
   
Lattice to be used for the calculations. The `lattice` option can use either LAMMPS for creation of input structure or use an input file in the LAMMPS data format. To use LAMMPS to create the structure, the keyword specified should be from the following: `BCC`, `FCC`, `HCP`, `DIA`, `SC` and `LQD`. LAMMPS lattice creation can **only be used for single species**. If `LQD` is specified, a solid structure will be first created and melted within the MD run. Alternatively, a LAMMPS data file can be specified which contains the configuration.

---

#### <a name="reference_phase"></a>`reference_phase`         
   
_type_: string or list of strings  
_default_: None  
_example_:
```
state: solid
state: [solid, liquid]
```

The protocol to be used for the calculation. The `reference_phase` input is closely related to the `lattice` command. It should be of the same length as the `lattice` input. For each of the lattice specified, `reference_phase` command specifies which reference state to use.

---

#### <a name="lattice_constant"></a>`lattice_constant`

_type_: float or list of floats  
_default_: Experimental lattice constant  
_example_:
```
lattice_constant: 3.68
lattice_constant: [3.68, 5.43]
```

Lattice constant for input structures. Lattice constant values to be used for initial structure creation. Only required if the structure is created through LAMMPS. If not specified, the experimental lattice constant will be used.

---

#### <a name="repeat"></a>`repeat`

_type_: list of ints of length 3   
_default_: [1, 1, 1]  
_example_:
```
repeat: [3,3,3]
```

`repeat` command specifies the number of unit cells required in each x, y and z directions. This is only used if `lattice` command uses one of the LAMMPS structure keywords.

---

#### <a name="n_iterations"></a>`n_iterations`         

_type_: int     
_example_:
```
n_iterations: 3
```

The number of backward and forward integrations to be carried out in all modes. If more than one integration cycle is used, the errors can also be evaluated.

---

#### <a name="temperature_high"></a>`temperature_high`

_type_: float  
_default_: 2*temperature  
_example_:
```
temperature_high: 1600
```

The temperature used to melt a solid structure to create a liquid. If `reference_phase` is chosen as `liquid`, calphy performs a melting cycle to create an equilibrated liquid structure. calphy starts from the given input structure, and heats it using 2 times the highest temperature provided in the `temperature` option. If the structure is not melted, the temperature is increased progressively. `temperature_high` keyword can be used to set the temperature to overheat the structure and melt it.   

---

#### <a name="npt"></a>`npt`

_type_: bool  
_default_: True  
_example_:
```
npt: True
```

`npt` determines if calculations are carried out in the NPT ensemble. This option is used with modes `alchemy`, `ts` and `mts`. The effect is described below:

- for mode `ts` and `mts`: the reversible scaling approach is carried out in NPT if `npt` is True, otherwise the `NVT` ensemble is used.

- for mode `alchemy`: If `npt` is False, the NVT ensemble is used, meaning that the calculated work during alchemical transformation is calculated on the equilibrated volume of the first pair style. 

---

#### <a name="pair_style"></a>`pair_style`

_type_: string or list of strings       
_example_:
```
pair_style: eam/alloy
pair_style: [eam/alloy, eam/alloy]
pair_style: eam/fs
pair_style: pace
```

The [pair style](https://lammps.sandia.gov/doc/pair_style.html) command for LAMMPS. All styles supported in LAMMPS can be used with calphy. Except for mode `alchemy`, only the first value in the list will be used. For mode `alchemy`, there should be two pair styles specified, and the alchemical transformation is carried out between the two.

---

#### <a name="pair_coeff"></a>`pair_coeff`

_type_: string or list of strings    
_default_: None  
_example_:
```
pair_coeff: "* * Cu_EAM/Cu01.eam.alloy Cu"
pair_coeff: ["* * Cu_EAM/Cu01.eam.alloy Cu", "* * Cu_EAM/Cu02.eam.alloy Cu"]
pair_coeff: "* * CuZr_EAM/CuZr.eam.fs Cu Zr"
```

The [pair coeff](https://lammps.sandia.gov/doc/pair_coeff.html) command for LAMMPS. It should be the same length as `pair_style`. Except for mode `alchemy`, only the first value in the list will be used. For mode `alchemy`, there should be two pair styles specified, and the alchemical transformation is carried out between the two.

---

#### <a name="potential_file"></a>`potential_file`

_type_: string      
_default_: None  
_example_:
```
pair_coeff: "/home/calc/potential.inp"
```

If specified, the `pair_style` and `pair_coeff` commands are not used, but rather the potential is read in from the provided input file using `include` command in LAMMPS. This allows the use of more complex or multiple potential files. Due to the `hybrid/scaled` styles employed in calphy, **this option only works with mode `fe` and `reference_phase` solid.**

---

#### <a name="n_equilibration_steps"></a>`n_equilibration_steps`

_type_: int   
_default_: 25000  
_example_:
```
n_equilibration_steps: 10000
```

---

The number of time steps for equilibrating the system.

#### <a name="n_switching_steps"></a>`n_switching_steps`

_type_: int or list of ints  
_default_: 50000  
_example_:
```
n_switching_steps: 10000
n_switching_steps: [10000, 20000]
```

The number of switching steps. If a list of two integers is provided, the first value is used for mode `fe` while the second value will be used for all other modes.

---

#### <a name="n_print_steps"></a>`n_print_steps`        

_type_: int  
_default_: 0  
_example_:
```
n_print_steps: 100
```

Record MD trajectory during temperature sweep runs in the given interval of time steps. Default 0, trajectories are never recorded.

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

#### <a name="timestep"></a>`timestep`        

_type_: float 
_default_: 0.001  
_example_:
```
timestep: 0.001
```

The timestep for MD in picoseconds. 

---

#### <a name="thermostat_damping"></a>`thermostat_damping`

_type_: float   
_default_: 0.1   
_example_:
```
thermostat_damping: 0.1
```

The thermostat damping constant for MD. 

---

#### <a name="barostat_damping"></a>`barostat_damping`

_type_: float  
_default_: 0.1  
_example_:
```
barostat_damping: 0.1
```

Pressure damping for MD. 

---

#### <a name="n_small_steps"></a>`n_small_steps`

_type_: int  
_default_: 10000  
_example_:
```
n_small_steps: 10000
```

The number of time steps for equilibration cycles to calculate spring constant and average volume.

---

#### <a name="n_every_steps"></a>`n_every_steps`         

_type_: int  
_default_: 10  
_example_:
```
n_every_steps: 100
```

Keywords to tune how often average values are recorded in LAMMPS. Please see [here](https://docs.lammps.org/fix_ave_time.html) for more details.

---

#### <a name="n_repeat_steps"></a>`n_repeat_steps`         

_type_: int  
_default_: 10    
_example_:
```
n_repeat_steps: 10
```

Keywords to tune how often average values are recorded in LAMMPS. Please see [here](https://docs.lammps.org/fix_ave_time.html) for more details.

---

#### <a name="n_cycles"></a>`n_cycles`        

_type_: int  
_default_: 100    
_example_:
```
n_cycles: 100
```

Number of cycles to try converging the pressure of the system. If the pressure is not converged after `n_cycles`, an error will be raised. In each `n_cycle`, `n_small_steps` MD steps will be run.

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

#### <a name="scheduler"></a>`scheduler`

_type_: string         
_example_:
```
scheduler: slurm
```

The scheduler to be used for the job. Can be `local`, `slurm` or `sge`. The code has been tested only on local and slurm.

---

#### <a name="cores"></a>`cores`

_type_: int           
_example_:
```
cores: 4
```

The number of cores to be used for the job.

---

#### <a name="jobname"></a>`jobname`         

_type_: string         
_example_:
```
jobname: cu
```

Name of the job. Not used for `local`.

---

#### <a name="walltime"></a>`walltime`         

_type_: string         
_example_:
```
walltime: "23:50:00"
```

Walltime for the job. Not used for `local`.

---

#### <a name="queuename"></a>`queuename`         

_type_: string         
_example_:
```
queuename: "shorttime"
```

Name of the queue. Not used for `local`.

---

#### <a name="memory"></a>`memory`         

_type_: string         
_example_:
```
memory: 3GB
```

Memory to be used per core. Not used for `local`.

---

#### <a name="commands"></a>`commands`         

_type_: list of strings         
_example_:
```
command:
 - source .bashrc
 - conda activate ace
 - module load lammps
```

Command that will be run **before** the actual calculations are carried out. This section can be used to specify commands that need to be run before the actual calculation. If the calculations are run within a conda environment, the activate command for conda should be specified here. If additional modules need to be loaded, that can also be specified here.

---

#### <a name="modules"></a>`modules`         

_type_: list of strings         
_example_:
```
modules:
  - anaconda
  - lammps
```

List of modules to be loaded before running the calculations. The given module names will be prefixed by `module load`.

---

#### <a name="options"></a>`options`         

_type_: string         
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

#### <a name="tol_spring_constant"></a>`spring_constant`

_type_: float  
_default_: 0.01  
_example_:
```
spring_constant: 0.01
```

tolerance for the convergence of spring constant calculation.

---

#### <a name="tol_solid_fraction"></a>`solid_fraction`

_type_: float         
_default_: 0.7  
_example_:
```
solid_fraction: 0.7
```

The minimum amount of solid particles that should be there in solid.

---

#### <a name="tol_liquid_fraction"></a>`liquid_fraction`

_type_: float         
_default_: 0.05  
_example_:
```
liquid_fraction: 0.05
```

Maximum fraction of solid atoms allowed in liquid after melting.

---

#### <a name="tol_pressure"></a>`pressure`

_type_: float         
_default_: 0.5  
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

#### <a name="step"></a>`step`

_type_: int         
_example_:
```
step: 200
```

Temperature interval for search of melting temperature. Only used if mode is `melting_temperature`.

---

#### <a name="attempts"></a>`attempts`

_type_: int         
_example_:
```
attempts: 5
```

The number of maximum attempts to try find the melting temperature in a automated manner. Only used if mode is `melting_temperature`.

---


