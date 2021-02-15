
# `pytint` input file

The inputfile is `yaml` formatted. In this section the possible keys in the inputfile is discussed. The input file consists of two main keys, and four separate blocks. For a sample of the inputfile see the end of this document. The two main keys are as shown below-

## main keys

- `element` : Chemical symbol(s) of the element(s) in the simulation.   
   _type_: string/list of strings  
   _example_:
   ```
   element: 'Cu'
   element: ['Cu', 'Zr']
   ```

- `mass` : Mass of the element(s) in the simulation. It should follow the same order as that of `element`.     
   _type_: float/list of floats    
   _example_:
   ```
   mass: 63.546
   mass: [63.546, 91.224]
   ```

## `calculation` block

Other than these two main keys, all other options are specified in blocks. The first block is the `calculations` block. This block can list all the calculations that `pytint` should perform and it can have more than one entry. A sample calculation block is shown below.

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

- `mode` : Calculation mode. Either `fe` for free energy calculations or `ts` for temperature sweep.     
   _type_: string, `fe` or `ts`  
   _example_:
   ```
   mode: fe
   mode: ts
   ```
- `temperature` : Temperatures for the simulation in Kelvin.         
   _type_: list of floats   
   _example_:
   ```
   temperature: [1200, 1300]
   ```
   The way temperature is used `pytint` depends on the selected mode of calculation. If the `mode` is `fe`, a free energy calculation is launched for each temperature on the list. If the mode is `ts`, a temperature sweep is carried out. In that case, only two values of temperature should be specified.
   
- `pressure` : Pressure for the simulation in bars.         
   _type_: list of floats   
   _example_:
   ```
   pressure: [0, 10000]
   ```
   For each pressure specified in the list, one calculation will be started. For example, in the following combination-
   ```
   mode: fe
   temperature: [1000, 1200]
   pressure: [0, 10000]
   ```  
   A total of 2 temperatures * 2 pressures, 4 calculations will be started. If the mode is `ts` for the same configuration above, 2 calculations will be started. 

- `lattice` : Lattice to be used for the calculations.         
   _type_: string of list of strings   
   _example_:
   ```
   lattice: FCC
   lattice: [FCC, LQD]
   lattice: [FCC, conf.data]
   ```
   The `lattice` option can use either LAMMPS for creation of input structure or use an input file in the LAMMPS data format. To use LAMMPS to create the structure, the keyword specified should be from the following: `BCC`, `FCC`, `HCP`, `DIA`, `SC` and `LQD`. LAMMPS lattice creation can **only be used for single species**. If `LQD` is specified, a solid structure will be first created and melted within the MD run. Alternatively, a LAMMPS data file can be specified which contains the configuration.

- `state` : The protocol to be used for the calculation.         
   _type_: string of list of strings   
   _example_:
   ```
   state: solid
   state: [solid, liquid]
   ```
   The `state` input is closely related to the `lattice` command. It should be of the same length as the `lattice` input. For each of the lattice specified, `state` command specifies which reference state to use.

- `lattice_constant` : lattice constant for input structures          
   _type_: list of floats   
   _example_:
   ```
   lattice_constant: 3.68
   lattice_constant: [3.68, 5.43]
   ```
   Lattice constant values to be used for initial structure creation. Only required if the structure is created through LAMMPS. If not specified, the experimental lattice constant will be used.

- `iso` : Specify if the barostat is isotropic or anisotropic.            
   _type_: list of bools   
   _example_:
   ```
   iso: True
   iso: [True, False]
   ```
   Specify whether the barostat will control the pressure isotropically or anisotropically.

   
- `repeat` : The number of unit cells to be replicated in each direction.         
   _type_: list of ints of length 3     
   _example_:
   ```
   repeat: [3,3,3]
   ```
   `repeat` command specifies the number of unit cells required in each x, y and z directions. This is only used if `lattice` command uses one of the LAMMPS structure keywords.
   
- `nsims` : The number of backward and forward interactions to be carried out.         
   _type_: int     
   _example_:
   ```
   nsims: 3
   ```

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

- `pair_style` : The [pair style](https://lammps.sandia.gov/doc/pair_style.html) command for LAMMPS.         
   _type_: string       
   _example_:
   ```
   pair_style: eam/alloy
   pair_style: eam/fs
   pair_style: pace
   ```

- `pair_coeff` : The [pair coeff](https://lammps.sandia.gov/doc/pair_coeff.html) command for LAMMPS.         
   _type_: string       
   _example_:
   ```
   pair_coeff: "* * Cu_EAM/Cu01.eam.alloy Cu"
   pair_coeff: "* * CuZr_EAM/CuZr.eam.fs Cu Zr"
   ```

- `timestep` : The timestep for MD in picoseconds.         
   _type_: float         
   _example_:
   ```
   timestep: 0.001
   ```

- `tdamp` : The thermostat damping for MD in time units.         
   _type_: float         
   _example_:
   ```
   tdamp: 0.1
   ```
- `pdamp` : The pressure damping for MD in time units.         
   _type_: float         
   _example_:
   ```
   pdamp: 0.1
   ```
- `te` : The number of time steps for equilibrating the system.         
   _type_: int         
   _example_:
   ```
   te: 10000
   ```
- `ts` : The number of time steps for switching between the system of interest and reference system.         
   _type_: int         
   _example_:
   ```
   ts: 10000
   ```

## `queue` block

This block consints of the options for specifying the scheduler for carrying out the calculations. An example block is given below-

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

- `scheduler` : The scheduler to be used for the job. Can be `local`, `slurm` or `sge`.         
   _type_: string         
   _example_:
   ```
   scheduler: slurm
   ```
   The code has been tested only on local and slurm.

- `cores` : The number of cores to be used for the job.  
   _type_: int           
   _example_:
   ```
   cores: 4
   ```
- `jobname` : Name of the job.         
   _type_: string         
   _example_:
   ```
   jobname: cu
   ```
   Not used for `local`.

- `walltime` : Walltime for the job.         
   _type_: string         
   _example_:
   ```
   walltime: "23:50:00"
   ```
   Not used for `local`.

- `queuename` : Name of the queue.         
   _type_: string         
   _example_:
   ```
   queuename: "shorttime"
   ```
   Not used for `local`.

- `memory` : memory to be used per core.         
   _type_: string         
   _example_:
   ```
   memory: 3GB
   ```
   Not used for `local`.

- `commands` : Command that will be run **before** the actual calculations are carried out.         
   _type_: list of strings         
   _example_:
   ```
   command:
     - source .bashrc
     - conda activate ace
     - module load lammps
   ```
   This section can be used to specify commands that need to be run before the actual calculation. If the calculations are run within a conda environment, the activate command for conda should be specified here. If additional modules need to be loaded, that can also be specified here.
   
## `conv` block

This block helps to tune the internal convergence parameters that `pytint` uses. Generally, tuning these parameters are not required.

```
conv:
   k_tol: 0.01
   solid_frac: 0.7
   liquid_frac: 0.05
   p_tol: 0.5
```

- `k_tol` : tolerance for the convergence of spring constant calculation.         
   _type_: float         
   _example_:
   ```
   ktol: 0.01
   ```

- `solid_frac` : The minimum amount of solid particles that should be there in solid.         
   _type_: float         
   _example_:
   ```
   solid_frac: 0.7
   ```

- `liquid_frac` : Maximum fraction of solid atoms allowed in liquid after melting.         
   _type_: float         
   _example_:
   ```
   liquid_frac: 0.05
   ```
- `p_tol` : tolerance for the convergence of pressure.         
   _type_: float         
   _example_:
   ```
   p_tol: 0.5
   ```






