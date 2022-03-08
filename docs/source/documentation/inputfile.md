
# `calphy` input file

The inputfile is `yaml` formatted. In this section the possible keys in the inputfile is discussed. The input file consists of two main keys, and four separate blocks. For a sample of the inputfile see the end of this document. The table below gives a quick overview of all available keywords in calphy.

| Main keys  | |
| :-: | :-: |
| [element](#element) | [mass](#mass) | 

| `calculations` block  | | | | |
| :-: | :-: | :-: | :-: | :-: |
| [mode](#mode) | [lattice](#lattice) | [state](#state) | [temperature](#temperature) | [pressure](#pressure) |
| [lattice_constant](#lattice_constant) | [iso](#iso) | [fix_lattice](#fix_lattice) | [repeat](#repeat) | [nsims](#nsims) | 

| `md` block  | | | | |
| :-: | :-: | :-: | :-: | :-: |
| [pair_style](#pair_style) | [pair_coeff](#pair_coeff) | [timestep](#timestep) | [nsmall](#nsmall) | [nevery](#nevery) |
| [nrepeat](#nrepeat) | [ncycles](#ncycles) | [tdamp](#tdamp) | [pdamp](#pdamp) | [te](#te) |
| [ts](#ts) | [tguess](#tguess) | [dtemp](#dtemp) | [maxattempts](#maxattempts) | [traj_interval](#traj_interval) |  

| `queue` block  | | | | |
| :-: | :-: | :-: | :-: | :-: |
| [scheduler](#scheduler) | [cores](#cores) | [jobname](#jobname) | [walltime](#walltime) | [queuename](#queuename) |
| [memory](#memory) | [commands](#commands) | [modules](#modules) | [options](#options) |

| `conv` block  | | | | |
| :-: | :-: | :-: | :-: | :-: |
| [alat_tol](#alat_tol) | [k_tol](#k_tol) | [solid_frac](#solid_frac) | [liquid_frac](#liquid_frac) | [p_tol](#p_tol) |

## main keys

#### <a name="element"></a>`element`

_type_: string/list of strings  
_example_:
```
element: 'Cu'
element: ['Cu', 'Zr']
```
Chemical symbol(s) of the element(s) in the simulation.   
   
#### <a name="mass"></a>`mass`

_type_: float/list of floats    
_example_:
```
mass: 63.546
mass: [63.546, 91.224]
```

Mass of the element(s) in the simulation. It should follow the same order as that of `element`, and should be of the same length.

## `calculations` block

Other than these two main keys, all other options are specified in blocks. The first block is the `calculations` block. This block can list all the calculations that `calphy` should perform and it can have more than one entry. A sample calculation block is shown below.

```
calculations:
- mode: ts 
  temperature: [1300, 1400]
  pressure: [0]
  lattice: [FCC]
  repeat: [2, 2, 2]
  state: solid
  nsims: 1
```

The various keys are-

#### <a name="mode"></a>`mode`  

_type_: string, `fe` or `ts` or `mts` or `alchemy` or `melting_temperature` or `tscale` or `pscale`
_example_:
```
mode: fe
mode: ts
```

Calculation mode. The modes can be chosen from `fe`, `ts`, `mts`, `alchemy`, `melting_temperature`, `tscale` or `pscale`. `fe` performs a direct free energy calculation, while `ts` performs a direct free energy calculation followed by reversible scaling to find temperature dependence. `mts` performs only reversible scaling part and can be used for dynamic Clausius-Clapeyron integration. Mode `alchemy` is used for switching between two different interatomic potentials, or for integration over concentration.  `melting_temperature` can be used for automated calculation of melting temperature. `tscale` is used for similar purpose as `ts`, but scales the temperature directly. `pscale` calculates the free energy as a function of the pressure.
   
#### <a name="temperature"></a>`temperature`  

_type_: float or list of floats   
_example_:
```
temperature: [1200, 1300]
```

Temperatures for the simulation in Kelvin. The way temperature is used in `calphy` depends on the selected mode of calculation. If the `mode` is `fe` or `alchemy`, a calculation is launched for each temperature on the list. If the mode is `ts` or `mts`, a temperature sweep is carried out. In that case, only two values of temperature should be specified.


#### <a name="pressure"></a>`pressure`  

_type_: list of floats   
_example_:
```
pressure: [0, 10000]
```

Pressure for the simulation in bars. Calculations are performed at each pressure specified for all modes except `pscale`. For each pressure specified in the list, one calculation will be started. For example, in the following combination-

```
mode: fe
temperature: [1000, 1200]
pressure: [0, 10000]
```  

A total of 2 temperatures * 2 pressures, 4 calculations will be started. If the mode is `ts` for the same configuration above, 2 calculations will be started. However if mode is `pscale`, 2 temperatures*1 pressure range, total of two calculations will be started.

#### <a name="lattice"></a>`lattice`

_type_: string of list of strings   
_example_:
```
lattice: FCC
lattice: [FCC, LQD]
lattice: [FCC, conf.data]
```
   
Lattice to be used for the calculations. The `lattice` option can use either LAMMPS for creation of input structure or use an input file in the LAMMPS data format. To use LAMMPS to create the structure, the keyword specified should be from the following: `BCC`, `FCC`, `HCP`, `DIA`, `SC` and `LQD`. LAMMPS lattice creation can **only be used for single species**. If `LQD` is specified, a solid structure will be first created and melted within the MD run. Alternatively, a LAMMPS data file can be specified which contains the configuration.
    
#### <a name="state"></a>`state`         
   
_type_: string of list of strings   
_example_:
```
state: solid
state: [solid, liquid]
```

The protocol to be used for the calculation. The `state` input is closely related to the `lattice` command. It should be of the same length as the `lattice` input. For each of the lattice specified, `state` command specifies which reference state to use.

#### <a name="lattice_constant"></a>`lattice_constant`

_type_: list of floats   
_example_:
```
lattice_constant: 3.68
lattice_constant: [3.68, 5.43]
```

Lattice constant for input structures. Lattice constant values to be used for initial structure creation. Only required if the structure is created through LAMMPS. If not specified, the experimental lattice constant will be used.

#### <a name="iso"></a>`iso`            

_type_: list of bools   
_example_:
```
iso: True
iso: [True, False]
```

Specify whether the barostat will control the pressure isotropically or anisotropically.
   
#### <a name="fix_lattice"></a>`fix_lattice`

_type_: list of bools   
_example_:
```
fix_lattice: True
fix_lattice: [True, False]
```

This option can be used to skip the equilibriation cycle. For example, if you need to calculate the free energy of a strained lattice, this option can be set to True. Then, either a data file with the strained structure can be provided, or `lattice_constant` option can be used to specify a strained lattice constant.

#### <a name="repeat"></a>`repeat`

_type_: list of ints of length 3     
_example_:
```
repeat: [3,3,3]
```

`repeat` command specifies the number of unit cells required in each x, y and z directions. This is only used if `lattice` command uses one of the LAMMPS structure keywords.

#### <a name="nsims"></a>`nsims`         

_type_: int     
_example_:
```
nsims: 3
```

The number of backward and forward integrations to be carried out in the non-equilibrium method. If more than one integration cycle is used, the errors can also be evaluated.
   
#### <a name="thigh"></a>`thigh`

_type_: float     
_example_:
```
thigh: 1600
```

The temperature used to melt a solid structure to create a liquid. If `state` is chosen as `liquid`, calphy performs a melting cycle to create an equilibrated liquid structure. calphy starts from the given input structure, and heats it using 2 times the highest temperature provided in the `temperature` option. If the structure is not melted, the temperature is increased progressively. `thigh` keyword can be used to set the temperature to overheat the structure and melt it.   

#### <a name="npt"></a>`npt`

_type_: bool    
_example_:
```
npt: True
```

`npt` determines if alchemical transformations are carried out in the NPT ensemble. If False, the calculations are carried out in the NVT ensemble. This means that the calculated work during alchemical transformation is calculated on the equilibrated volume of the first pair style. Only used for mode `alchemy`. 


## `md` block

MD block consists of the various options required for the MD runs. An example block is given below and the keys are discussed.

```
md:
  pair_style: eam/alloy
  pair_coeff: "* * Cu_EAM/Cu01.eam.alloy Cu"
  timestep: 0.001
  tdamp: 0.1
  pdamp: 0.1
  te: 10000
  ts: 15000
```

#### <a name="pair_style"></a>`pair_style`

_type_: string or list of strings       
_example_:
```
pair_style: eam/alloy
pair_style: eam/fs
pair_style: pace
```

The [pair style](https://lammps.sandia.gov/doc/pair_style.html) command for LAMMPS. All styles supported in LAMMPS can be used with calphy. Except for mode `alchemy`, only the first value in the list will be used. For mode `alchemy`, there should be two pair styles specified, and the alchemical transformation is carried out between the two.

#### <a name="pair_coeff"></a>`pair_coeff`

_type_: string or list of strings      
_example_:
```
pair_coeff: "* * Cu_EAM/Cu01.eam.alloy Cu"
pair_coeff: "* * CuZr_EAM/CuZr.eam.fs Cu Zr"
```

The [pair coeff](https://lammps.sandia.gov/doc/pair_coeff.html) command for LAMMPS. It should be the same length as `pair_style`. Except for mode `alchemy`, only the first value in the list will be used. For mode `alchemy`, there should be two pair styles specified, and the alchemical transformation is carried out between the two.

#### <a name="timestep"></a>`timestep`        

_type_: float         
_example_:
```
timestep: 0.001
```

The timestep for MD in picoseconds. 

#### <a name="tdamp"></a>`tdamp`

_type_: float         
_example_:
```
tdamp: 0.1
```

The thermostat damping constant for MD. A [Langevin thermostat](https://docs.lammps.org/fix_langevin.html) is used for calculations in calphy.

#### <a name="pdamp"></a>`pdamp`

_type_: float         
_example_:
```
pdamp: 0.1
```

Pressure damping for MD. A [Parrinello-Rahman](https://docs.lammps.org/fix_nh.html) barostat is used for calculations in calphy.

#### <a name="te"></a>`te`

_type_: int         
_example_:
```
te: 10000
```

The number of time steps for equilibrating the system.

#### <a name="ts"></a>`ts`

_type_: int         
_example_:
```
ts: 10000
```

The number of time steps for switching between the system of interest and reference system.

#### <a name="tguess"></a>`tguess`

_type_: int         
_example_:
```
tguess: 1300
```

An initial starting guess for melting temperature. Only used if mode is `melting_temperature`.

#### <a name="dtemp"></a>`dtemp`

_type_: int         
_example_:
```
dtemp: 200
```

Temperature interval for search of melting temperature. Only used if mode is `melting_temperature`.

#### <a name="maxattempts"></a>`maxattempts`

_type_: int         
_example_:
```
maxattempts: 5
```

The number of maximum attempts to try find the melting temperature in a automated manner. Only used if mode is `melting_temperature`.

#### <a name="nsmall"></a>`nsmall`

_type_: int         
_example_:
```
nsmall: 10000
```

The number of time steps for equilibration cycles to calculate spring constant and average volume.
   
#### <a name="nevery"></a>`nevery`         

_type_: int         
_example_:
```
nevery: 100
```

Keywords to tune how often average values are recorded in LAMMPS. Please see [here](https://docs.lammps.org/fix_ave_time.html) for more details.
   
#### <a name="nrepeat"></a>`nrepeat`         

_type_: int         
_example_:
```
nrepeat: 10
```

Keywords to tune how often average values are recorded in LAMMPS. Please see [here](https://docs.lammps.org/fix_ave_time.html) for more details.

#### <a name="ncycles"></a>`ncycles`        

_type_: int         
_example_:
```
ncycles: 100
```

Number of cycles to try converging the pressure of the system. If the pressure is not converged after `ncycles`, an error will be raised. In each `ncycle`, `nsmall` MD steps will be run.

#### <a name="traj_interval"></a>`traj_interval`        

_type_: int         
_example_:
```
traj_interval: 100
```

Record MD trajectory during temperature sweep runs in the given interval of time steps. Default 0, trajectories are never recorded.


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

#### <a name="scheduler"></a>`scheduler`

_type_: string         
_example_:
```
scheduler: slurm
```

The scheduler to be used for the job. Can be `local`, `slurm` or `sge`. The code has been tested only on local and slurm.

#### <a name="cores"></a>`cores`

_type_: int           
_example_:
```
cores: 4
```

The number of cores to be used for the job.
   
#### <a name="jobname"></a>`jobname`         

_type_: string         
_example_:
```
jobname: cu
```

Name of the job. Not used for `local`.

#### <a name="walltime"></a>`walltime`         

_type_: string         
_example_:
```
walltime: "23:50:00"
```

Walltime for the job. Not used for `local`.
  
#### <a name="queuename"></a>`queuename`         

_type_: string         
_example_:
```
queuename: "shorttime"
```

Name of the queue. Not used for `local`.

#### <a name="memory"></a>`memory`         

_type_: string         
_example_:
```
memory: 3GB
```

Memory to be used per core. Not used for `local`.

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
   
#### <a name="modules"></a>`modules`         

_type_: list of strings         
_example_:
```
modules:
  - anaconda
  - lammps
```

List of modules to be loaded before running the calculations. The given module names will be prefixed by `module load`.

#### <a name="options"></a>`options`         

_type_: string         
_example_:
```
options:
  - memory: 3GB
```

Extra options to be added to the submission script.

## `conv` block

This block helps to tune the internal convergence parameters that `calphy` uses. Generally, tuning these parameters are not required.

```
conv:
   k_tol: 0.01
   solid_frac: 0.7
   liquid_frac: 0.05
   p_tol: 0.5
```

#### <a name="k_tol"></a>`k_tol`

_type_: float         
_example_:
```
ktol: 0.01
```

tolerance for the convergence of spring constant calculation.

#### <a name="solid_frac"></a>`solid_frac`

_type_: float         
_example_:
```
solid_frac: 0.7
```

The minimum amount of solid particles that should be there in solid.

#### <a name="liquid_frac"></a>`liquid_frac`

_type_: float         
_example_:
```
liquid_frac: 0.05
```

Maximum fraction of solid atoms allowed in liquid after melting.

#### <a name="p_tol"></a>`p_tol`

_type_: float         
_example_:
```
p_tol: 0.5
```

tolerance for the convergence of pressure.





