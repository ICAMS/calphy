## Running `calphy`

After the calculations are specified in the input file, `calphy` can be run by:

```
calphy -i inputfilename
```

The `calphy` command will wrap each calculation defined in the `calculations` block into a script and submit it based on the scheduler specified. 

Alternatively, if you want to run a calculation directly from the terminal, the `calphy_kernel` command can be used.

```
calphy_kernel -i input.yaml -k 0
```

An extra argument `-k` or `--kernel` is required, which takes an integer argument. This argument refers to the which calculation needs to be run from the input file. For example if the `calculations`block of the input file is as follows:

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