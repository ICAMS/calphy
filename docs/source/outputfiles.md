(outputfiles)=
# Output files

Every calphy calculation writes its results into a dedicated **simulation
folder** created in the working directory. The folder is named

```
<mode>-<lattice>-<phase>-<temperature>-<pressure>
```

for example `fe-fcc-solid-500-0` (an [](mode) `fe` run on FCC, solid
[](reference_phase), 500 K, 0 bar). A custom prefix can be prepended with the
`folder_prefix` input key.

All plain-text data files (`*.dat`) carry a `#`-commented header line naming
their columns and units, so they can be read directly with
`numpy.loadtxt(..., comments="#")`, `pandas`, or any plotting tool. Unless
noted otherwise the units follow LAMMPS `metal` units:

| Quantity | Unit |
|---|---|
| energy | eV/atom |
| length | Å |
| volume | Å³ |
| pressure | bar |
| temperature | K |
| mean-squared displacement | Å² |

---

## `report.yaml` — the primary result

The machine-readable summary of the calculation. It has three blocks:

**`input`** — the thermodynamic state the run targeted:

| Key | Meaning |
|---|---|
| `temperature` | Target temperature [K] ([](temperature)) |
| `pressure` | Target pressure [bar] ([](pressure)) |
| `lattice` | Input lattice / structure |
| `element` | Element(s) ([](element)) |
| `concentration` | Composition of each element |

**`average`** — quantities measured during equilibration:

| Key | Meaning |
|---|---|
| `vol_atom` | Equilibrium volume per atom [Å³] |
| `spring_constant` | Fitted Einstein spring constant(s) [eV/Å²] (solid reference only) |
| `density` | Number density [1/Å³] (liquid reference only) |

**`results`** — the free energy and its components (all in `unit`, eV/atom):

| Key | Meaning |
|---|---|
| `free_energy` | **The result** — Helmholtz/Gibbs free energy per atom |
| `error` | Statistical uncertainty on `free_energy`: the standard error of the mean over the [](n_iterations) independent switching runs |
| `reference_system` | Free energy of the analytic reference (Einstein crystal or Uhlenbeck–Ford model) |
| `einstein_crystal` | Einstein-crystal contribution (solid) |
| `com_correction` | Fixed-centre-of-mass finite-size correction (solid) |
| `work` | Reversible switching work between reference and system of interest |
| `dissipation` | Mean switching dissipation, `0.5·(W_forward + W_backward)`. A measure of the irreversibility of the switching — ideally close to `0`; large values indicate the switching was too fast or a structural change occurred |
| `ts_dissipation` | *(ts / tscale only)* maximum energy dissipation along the reversible-scaling temperature sweep. Clean sweeps give ~1e-4 eV/atom; much larger values flag a hidden phase transition and a contaminated `temperature_sweep.dat` |
| `pv` | Pressure–volume contribution |
| `unit` | Units of the above (`eV/atom`) |

---

## Free-energy curves

Produced by the temperature/pressure-scaling modes. These are the files you
plot to get a free-energy curve rather than a single point.

### `temperature_sweep.dat` — [](mode) `ts` / `tscale`

Free energy as a function of temperature from a reversible-scaling sweep.

| Column | Meaning |
|---|---|
| 1 | `temperature` [K] |
| 2 | `free_energy` [eV/atom] |
| 3 | `error` [eV/atom] — standard error of the mean across iterations |

The spacing of the temperature points is controlled by [](lambda_schedule)
(`linear` clusters samples at the low-temperature end; `uniform_temperature`
spaces them evenly in temperature).

### `pressure_sweep.dat` — [](mode) `pscale`

Free energy as a function of pressure.

| Column | Meaning |
|---|---|
| 1 | `pressure` [bar] |
| 2 | `free_energy` [eV/atom] |
| 3 | `error` [eV/atom] |

---

## Raw switching data

The per-step data recorded during each non-equilibrium switching run. There is
one forward/backward pair **per iteration** ([](n_iterations)); calphy
integrates these to obtain `work` and `dissipation`. They are mainly of
interest for diagnostics (checking hysteresis, plotting the switching path).

### `forward_<i>.dat` / `backward_<i>.dat` — [](mode) `fe`

The Frenkel–Ladd / non-equilibrium switching between the system and its
analytic reference for iteration `i`.

*Solid reference*

| Column | Meaning |
|---|---|
| 1 | `dU_sys` — system potential energy [eV/atom] |
| 2 … *N*+1 | `dU_ref` — Einstein reference energy of each of the *N* elements [eV/atom] |
| last | `lambda` — switching parameter (0 → 1) |

*Liquid reference*

| Column | Meaning |
|---|---|
| 1 | `dU_sys` [eV/atom] |
| 2 | `dU_ref` — Uhlenbeck–Ford reference energy [eV/atom] |
| 3 | `lambda` |

For the two-leg liquid path the files are split as `forward_leg1_<i>.dat`,
`forward_leg2_<i>.dat` (and the `backward_` equivalents).

### `forward_<i>.dat` / `backward_<i>.dat` — [](mode) `alchemy`

| Column | Meaning |
|---|---|
| 1 | `dU_1` — potential energy of system 1 [eV/atom] |
| 2 | `dU_2` — potential energy of system 2 [eV/atom] |
| 3 | `lambda` |

### `ts.forward_<i>.dat` / `ts.backward_<i>.dat` — [](mode) `ts` / `tscale`

The reversible-scaling sweep for iteration `i`.

| Column | Meaning |
|---|---|
| 1 | `dU` — potential energy [eV/atom] |
| 2 | `press` [bar] |
| 3 | `vol` [Å³] |
| 4 | `lambda` (mapped to temperature `T₀/lambda`) |

### `ps.forward_<i>.dat` / `ps.backward_<i>.dat` — [](mode) `pscale`

Same four columns as `ts.*`, but `lambda` maps to pressure.

---

## Averaging data

Time series written during the NPT/NVT equilibration and volume-convergence
stage.

### `avg.dat`

| Column | Meaning |
|---|---|
| 1 | `TimeStep` |
| 2–4 | `lx`, `ly`, `lz` — box lengths [Å] |
| 5 | `press` [bar] |
| 6 | `pe` — potential energy [eV/atom] |
| 7 | `etotal` — total energy [eV/atom] |
| 8 | `temp` [K] |

### `msd.dat`

| Column | Meaning |
|---|---|
| 1 | `TimeStep` |
| 2 … | `msd<k>` — mean-squared displacement of element group *k* [Å²] |

The MSD is used to fit the Einstein spring constant for the solid reference.

---

## Trajectories and configurations

These are written in native LAMMPS formats and are only produced when the
relevant input options are set (`n_print_steps` for trajectories):

| File | Format | Content |
|---|---|---|
| `traj.*.dat` | LAMMPS dump | Atomic trajectory snapshots (equilibration, switching, melting checks) |
| `conf.*.data` | LAMMPS data | Equilibrated configurations (restartable structures) |

---

## Logs

| File | Content |
|---|---|
| `calphy.log` | The human-readable run log — stage timing, convergence messages, warnings, and the citations to use |
| `calphy.seg<k>.lmp` | The generated LAMMPS input script for execution segment `k` |
| `*.seglog` | The LAMMPS log for each segment |

The `*.seg*` files are produced by the executable runner (calphy drives the
`lmp` binary segment by segment) and are primarily useful for debugging a
failed run — the LAMMPS error message will be at the end of the relevant
`*.seglog`.
