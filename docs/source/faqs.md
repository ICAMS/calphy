# FAQs

**My LAMMPS run dies immediately with a `pair_style` / `pair_coeff` error.**

The most common cause is that the `pair_coeff` element list does not match the `element` list (or is in a different order). For multi-element systems calphy enforces strict ordering: `element[i]` is mapped to LAMMPS atom type `i+1`, and the elements at the end of `pair_coeff` must appear in the same order. Reorder both lists so they agree.

**`calphy_kernel` complains that the simulation folder already exists.**

calphy refuses to overwrite an existing run folder so that data is never silently lost. Either remove the offending `mode-lattice-temp-pressure` directory or set [`folder_prefix`](folder_prefix) to a unique tag.

**My liquid free-energy calculation reports MSD = 0 / spring constant warnings.**

This happens when the structure has not actually melted during the equilibration cycle. Increase [`temperature_high`](temperature_high) (calphy will overheat to this value before quenching) and/or set [`melting_cycle: True`](melting_cycle).

**My liquid free-energy calculation loses atoms during the integration (`ERROR: Lost atoms`).**

This is most common for **molecular liquids** (e.g. water) where the single-component Uhlenbeck–Ford reference cannot represent two well-separated length scales. There are two distinct failure modes:

- *Detonation* — the UFM energy spikes just before the crash. The single `sigma` is larger than a bonded distance, so bonded pairs sit inside the UFM repulsive core and are blown apart. Use the two-leg path with a per-element-pair [`sigma`](ufm_sigma) dict and a small cross-term (e.g. `H_O: 0.9`) so bonded pairs stay outside the core.
- *Drift* — atoms wander off slowly with no energy spike. The UFM is purely repulsive, so once the real potential is switched off the reference must confine the atoms on its own; a too-small [`p`](ufm_p) (the confinement is `p` in units of `kB·T`) lets light atoms escape. Increase `p` (the default `50` is usually sufficient).

The two effects are independent: `sigma` controls detonation, `p` controls drift. See the [`uhlenbeck_ford_model`](ufm_sigma) block for the two-leg setup.

**How do I run calphy with the GPU / KOKKOS package?**

Use the [`md.cmdargs`](cmdargs) keyword to pass the standard LAMMPS accelerator command-line arguments:

```
md:
  cmdargs: "-k on g 1 -sf kk -pk kokkos newton on neigh half"
```

