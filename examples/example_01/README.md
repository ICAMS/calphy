# Example 01: Calculation of free energy

In this example, a simple calculation of the Helmholtz free energy $F(NVT)$ is illustrated. We use a Fe BCC structure at 100K.

The EAM potential that will be used: [Meyer, R, and P Entel. “Martensite-austenite transition and phonon dispersion curves of Fe1−xNix studied by molecular-dynamics simulations.” Phys. Rev. B 57, 5140.](https://doi.org/10.1103/PhysRevB.57.5140)

The reference data is from: [Freitas, Rodrigo, Mark Asta, and Maurice de Koning. “Nonequilibrium Free-Energy Calculation of Solids Using LAMMPS.” Computational Materials Science 112 (February 2016): 333–41.](https://doi.org/10.1016/j.commatsci.2015.10.050)

## Solid free energy

The [input file](input.yaml) for calculation of solid BCC Fe at 100 K is provided in the folder. A detailed description of the input file is available [here](../inputfile.md). The calculation can be started from the terminal using:

```
calphy -i input.yaml
```

This should give a message:

```
Total number of 1 calculations found
```

`calphy` is running in the background executing the calculation. A new folder called `fe-BCC-100-0`, and files called `fe-BCC-100-0.sub.err` and `fe-BCC-100-0.sub.out` will also be created. 

After the calculation is over, the results are available in the `report.yaml` file. The file is shown below:

```
average:
  spring_constant: '3.32'
  vol/atom: 11.994539429692749
input:
  concentration: '1'
  element: Fe
  lattice: bcc
  pressure: 0.0
  temperature: 100
results:
  error: 0.0
  free_energy: -4.263568357143783
  pv: 0.0
  reference_system: 0.015029873513789175
  work: -4.278598230657573
```

The calculated free energy is $-4.2636$ eV/atom, the value reported in the publication is $-4.2631147(1)$ eV/atom. The calculated value can be further improved by increasing the system size, and by increasing the switching time. 





