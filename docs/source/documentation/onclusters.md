# Running `calphy` on computing clusters

The easiest way to run `calphy` on computing clusters is by using `queue` block in the input file. The `queue` block looks like this:

```
queue:
  scheduler: slurm
  cores: 20
  jobname: fec
  walltime: "23:50:00"
  queuename: shorttime
  memory: 3GB
  commands:
    - source ~/.bashrc
    - conda activate calphy
```

The default `scheduler` is local, which means that the calculations are run on the local machine. Instead `slurm` or `sge` can be specified to run on computing clusters. The number of cores, walltime, queue name etc can be set during the respective keywords. Note that if you are using a conda environment, it needs to be activated. This can be done using the `commands` argument. Anything listed within the `commands` argument is copied directly to the submission script. 

## Adding a new scheduler

```
calculations:
- mode: ts 
  temperature: [1200, 1400]
  pressure: [0]
  lattice: [FCC, LQD]
  repeat: [5, 5, 5]
  state: [solid, liquid]
  nsims: 1
```

In the current version only `slurm` and `sge` schedulers are supported. There are two ways in which a new scheduler can be added, they are described below:

### Editing `scheduler.py`

The scheduler are implemented in `calphy` within `calphy/scheduler.py` file. An existing class can be copied and modified to support the new submission script.

### Adding calphy to a submission script

An easier method would be to add `calphy` manually within a submission script. `calphy` offers a `calphy_kernel` command line interface which can be used for direct submission. The above `queue` block can be translated to an equivalent submission script as follows:

```
#!/bin/bash
#SBATCH --job-name=fec
#SBATCH --time=23:50:00
#SBATCH --partition=shorttime
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=3GB

source ~/.bashrc
conda activate calphy

calphy_kernel -i input.yaml -k 0
```

Note that `calphy_kernel` is called instead of `calphy`. Additionally an extra argument `-k` or `--kernel` is required, which takes an integer argument. This argument refers to the which calculation needs to be run from the input file. For example if the `calculations`block of the input file is as follows:

```
calculations:
- mode: ts 
  temperature: [1200, 1400]
  pressure: [0]
  lattice: [FCC, LQD]
  repeat: [5, 5, 5]
  state: [solid, liquid]
  nsims: 1
```

Specify `-k 0` will run a calculation block equivalent to:

```
- mode: ts 
  temperature: [1200, 1400]
  pressure: [0]
  lattice: [FCC]
  repeat: [5, 5, 5]
  state: [solid]
  nsims: 1
```
and `-k 1` will run:

```
- mode: ts 
  temperature: [1200, 1400]
  pressure: [0]
  lattice: [LQD]
  repeat: [5, 5, 5]
  state: [liquid]
  nsims: 1
```