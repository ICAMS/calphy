# Validation of solid `fe` with Monte Carlo swaps (configurational entropy)

Runs for the cluster. They validate the Frenkel-Ladd path with species swaps
(`mode: fe`, `monte_carlo.n_swaps > 0`, see the `monte_carlo` block of the
input documentation) beyond what the unit tests and the short fork-gated test
in `tests/test_fe_swaps_lammps.py` can do on a laptop.

## Prerequisites

- calphy from branch `fe-swaps` (or a release that contains it).
- The thermoatoms LAMMPS fork (`github.com/thermoatoms/lammps`, branch
  `develop`), built with at least `MANYBODY`, `MC`, `EXTRA-FIX` and, for ACE
  potentials, `ML-PACE` (enables the `localE` fast path of `fix atom/swap`).
  Only the fork accepts an equal-style variable as the swap temperature of
  `fix atom/swap`; stock LAMMPS stops at the first switching run.
- `lammps_executable` in every input must point at that build.

Generate the structure and potential files first:

```bash
python make_degenerate_cu.py --repeat 6          # cu2.data, cu2.eam.alloy (T1)
python make_alloy_inputs.py --elements Cu Ni --masses 63.546 58.693 \
    --potential /path/to/CuNi.eam.alloy --lattice-file /path/to/cuni_fcc_random.data \
    --temperature 800 --cores 16                  # T3 and T4 inputs
```

Then run each input with `calphy -i inputs/<name>.yaml` from this directory,
and evaluate with `analyse.py` (see the end of this file).

## Test matrix

| id | system | what it checks | pass criterion |
|---|---|---|---|
| T1 | fcc Cu, 864 atoms, half of the atoms relabelled type 2, two-element copy of the Cu EAM with identical functions, equal masses, 500 K | bookkeeping and dissipation of the swap machinery in a system where every arrangement is exactly degenerate at every lambda | `F_swaps - F_noswaps = k_B T ln(1/2) = -0.029865 eV/atom` within 2 sigma of the combined error and within 0.5 meV/atom; acceptance 1.0 everywhere; `dissipation` of the same size in both runs |
| T2 | B2 CuZr, 1024 atoms, Mendelev ZrCu EAM, 800 K | an ordered compound: the swap acceptance at lambda = 0 is small and the arrangement stays ordered; the ordering work cancels the ideal mixing term | first shell Warren-Cowley parameter alpha_CuZr stays close to -1 in the dumps at lambda = 0; `F_swaps - F_noswaps` is small compared with `k_B T ln(1/2)` (expect a few meV/atom at most); if the B2 order is lost instead, report alpha against lambda and treat the run as a measurement, not a failure |
| T3 | random fcc solid solution (e.g. Cu-Ni), 800 K, switching lengths 25k, 50k, 100k, 200k steps | convergence against the switching rate: the arrangement must follow the switching near lambda = 0, where the acceptance is lowest | `W_fwd + W_bwd` decreases with the switching length; `free_energy` of the two longest runs agree within 1 meV/atom and 2 sigma |
| T4 | the T3 alloy by two routes | cross-check of the configurational free energy: Frenkel-Ladd with swaps at composition x against composition scaling with swaps from the pure metal | `F_FL+swaps(x) = F_fe(pure) + dF_cs(pure -> x)` within 2 sigma (both routes include the ideal mixing term and the short range order) |

T3 and T4 need an alloy potential that is not in the repository; the
generator writes the inputs with the paths given on the command line.  A
random solid solution data file can be made with `make_degenerate_cu.py
--relabel-only` applied to any fcc data file, or with ASE.

## What to look at

- `swap.equilibration.dat`: the configurational equilibration must converge
  (`calphy.log` says so); if the acceptance is below about 1 percent, raise
  `monte_carlo.n_swaps` or `md.n_small_steps` rather than the tolerance.
- `swap.fe.forward_<i>.dat` / `swap.fe.backward_<i>.dat`: acceptance against
  lambda (`analyse.py acceptance`).  It is 1 at lambda = 1 and falls to the
  alloy's value at lambda = 0.  The lambda range where it is below a few
  percent is where the arrangement lags, and the forward/backward mismatch
  comes from there.
- `report.yaml`: `free_energy`, `error`, `work`, `dissipation`, `ideal_mixing`.
  Compare free energies between runs, never the work alone: the measured
  spring constants differ between runs and shift the reference and the work
  against each other (set `spring_constants` explicitly to avoid that).

## Analysis

```bash
python analyse.py summary  <run folder> [<run folder> ...]     # report table
python analyse.py degenerate <swaps folder> <noswaps folder>   # T1 criterion
python analyse.py acceptance <run folder> [--bins 20]          # acceptance vs lambda
python analyse.py sro <dump file> --cutoff 3.0 [--types 1 2]   # Warren-Cowley alpha per frame
python analyse.py rate <run folder> [<run folder> ...]         # dissipation vs n_switching_steps
python analyse.py compare <FL+swaps folder> <pure fe folder> <composition scaling folder>  # T4
```

`sro` reads LAMMPS dumps with ASE (`traj.fe.forward_<i>.dat` when
`n_print_steps > 0`).  The cutoff is the midpoint between the first and second
neighbour shell of the lattice (fcc Cu: 3.0 A at a = 3.615 A; B2 CuZr: 3.3 A
for the first shell at a = 3.26 A).
