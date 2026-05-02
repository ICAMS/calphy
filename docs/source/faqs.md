# FAQs

**My LAMMPS run dies immediately with a `pair_style` / `pair_coeff` error.**

The most common cause is that the `pair_coeff` element list does not match the `element` list (or is in a different order). For multi-element systems calphy enforces strict ordering: `element[i]` is mapped to LAMMPS atom type `i+1`, and the elements at the end of `pair_coeff` must appear in the same order. Reorder both lists so they agree.

**`calphy_kernel` complains that the simulation folder already exists.**

calphy refuses to overwrite an existing run folder so that data is never silently lost. Either remove the offending `mode-lattice-temp-pressure` directory or set [`folder_prefix`](folder_prefix) to a unique tag.

**My liquid free-energy calculation reports MSD = 0 / spring constant warnings.**

This happens when the structure has not actually melted during the equilibration cycle. Increase [`temperature_high`](temperature_high) (calphy will overheat to this value before quenching) and/or set [`melting_cycle: True`](melting_cycle).

**How do I run calphy with the GPU / KOKKOS package?**

Use the [`md.cmdargs`](cmdargs) keyword to pass the standard LAMMPS accelerator command-line arguments:

```
md:
  cmdargs: "-k on g 1 -sf kk -pk kokkos newton on neigh half"
```

