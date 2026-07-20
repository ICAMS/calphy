# Baseline capture info (Part 0)

Reference data frozen from **current `main` (pylammpsmpi path)** before the ExecutableRunner
refactor. These files are an **immutable contract**: after the refactor the code must reproduce
them within the Part 6 tolerances. Never edit baselines to make code pass — fix the code.

## Environment

| item | value |
|---|---|
| date | 2026-07-08 |
| calphy version | 1.8.4 |
| calphy commit | 412cefd7f79fd97ff97882a089c52f8842ca0302 |
| calphy branch | only_exec |
| LAMMPS version | 20260211 |
| lammps package | 2026.2.11 (pypi wheel; ships `lmp` binary) |
| lmp binary | $CONDA_PREFIX/bin/lmp |
| pylammpsmpi | 0.3.12 |
| mpi4py | 4.0.1 |
| openmpi | 5.0.8 |
| mpirun | $CONDA_PREFIX/bin/mpirun |
| python | 3.12.11 |
| platform | Darwin 25.5.0 arm64 (macOS 26.5, build 25F71) |

## Test suite status (baseline)

Command: `pytest tests -q` (from repo root, conda env `calphy`).

```
105 passed, 2 skipped, 34 warnings in 23.72s
```

No failures. The 2 skips are pre-existing (`ss` in the first test file). Warnings are all
third-party Pydantic v2 deprecation / ASE FutureWarning noise, not calphy errors.

## Baseline scenario settings (Part 0.2)

Common settings for every scenario (unless noted in the scenario table):

- structure: **generated FCC Cu, 5×5×5 = 500 atoms** at the EAM equilibrium lattice constant
  (`lattice: FCC`, `repeat: [5,5,5]`), potential `Cu01.eam.alloy` — chosen over the stretched
  `tests/conf1.data` box so melt/solidify/prescan paths start from an equilibrated crystal.
- `queue.cores: 4`, `queue.scheduler: local`
- `n_equilibration_steps: 2500`, `n_switching_steps: 5000`, `md.n_small_steps: 1000`
- `n_iterations: 2`, `spring_constants: null` (forces the spring-constant convergence path)

See `inputs/B*.yaml` for the exact per-scenario files and `<id>/report.yaml` for the frozen results.

## Capture results (Part 0.2)

10 of 12 scenarios captured and committed to `tests/baselines/<id>/` (`report.yaml` + `calphy.log`
+ input). Two scenarios are **deferred** (see below).

| id | mode / phase | free_energy (eV/atom) | notes |
|----|--------------|-----------------------|-------|
| B1 | fe / solid, 500 K, 0 bar | −3.6439 | zero-pressure convergence |
| B2 | fe / solid, 500 K, 10000 bar | −3.5694 | finite-pressure path |
| B3 | fe / solid, 500 K, fix_lattice | −3.4021 | constrained (pressure: null); thermal P≈32.6 kbar |
| B4 | fe / liquid, 1400 K, melt cycle | −4.2259 | heat-melt-quench |
| B5 | fe / liquid, 1400 K, rattle | −4.2378 | melting_cycle: False |
| B6 | ts / solid, 400→700 K | −3.6005 @400 K | reversible scaling |
| B7 | ts / liquid, 1400→1600 K | −4.2266 @1400 K | RS, liquid ref |
| B8 | tscale / solid, 400→700 K | −3.6005 @400 K | matches B6 |
| B9 | pscale / solid, 500 K, 0→50000 bar | −3.6445 @0 bar | matches B1 |
| B10 | alchemy, Cu eam/alloy → eam/fs, 500 K | see report | potential-scaling path |

Internal consistency checks (frozen as contract): **B1≈B9** (0 bar/500 K), **B6≈B8** (400 K solid),
**B4≈B5≈B7** (1400 K liquid).

## Deferred scenarios (blocked; not captured)

- **B11 — fe-qtb / solid, 300 K.** BLOCKED: the LAMMPS build in this env has no QTB package
  (`ERROR: Unrecognized fix style 'qtb'`). Current calphy drives LAMMPS through the `lammps` python
  module (pylammpsmpi), same wheel, so QTB would require swapping the whole build — which would
  change the LAMMPS under the 10 captured baselines. Input kept at `inputs/B11.yaml`. Capture later
  with a QTB-enabled binary; the QTB command stream is still covered by the Part 2 golden and Part 4
  preflight (no binary needed).

- **B12 — ts + prescan / solid, 400→700 K, `phase_transition_detection.mode: warn`.** RESOLVED.
  `calphy/range_scan.py` (the pure-numpy detector behind `scan_temperature_range`) was ported from
  the `main` branch; its `RangeScan` / `plot_scan` / `ScanResult` API matches what `phase.py` already
  expected, and the LAMMPS driving stays in `phase.py` (executable-runner native), so no backend
  changes were needed. The `prescan.txt` golden is generated and end-to-end validated against real
  LAMMPS (ramp 400→700 K, `prescan.forward.dat` analysed, `prescan_signals.png` written, correctly
  reports no transition for solid Cu in-range). Input kept at `inputs/B12.yaml`.

## Environment / install note

calphy is installed **non-editable** in site-packages (a byte-identical copy of the repo `only_exec`
checkout at capture time — verified with `diff -rq`). Baselines therefore reflect the repo code.
Before Part 1 the package should be reinstalled editable (`pip install -e .`) so `calphy_kernel`
invoked from any directory uses the repo source rather than the stale copy.
