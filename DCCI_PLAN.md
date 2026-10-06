# Implementation plan: dynamic Clausius-Clapeyron integration (mode `dcci`)

**Status:** proposal, awaiting owner review
**Date:** 2026-09-10
**Reference:** M. de Koning, A. Antonelli, S. Yip, *Single-simulation determination of phase
boundaries: A dynamic Clausius-Clapeyron integration method*, J. Chem. Phys. 115, 11025 (2001),
DOI 10.1063/1.1420486 (`reference_paper/11025_1_online.pdf`).

**Goal:** given one known coexistence point (T_m, P_i) for a solid and a liquid, trace the
whole coexistence line P_coex(T) between P_i and P_f in one nonequilibrium run per direction,
with two LAMMPS cells (solid, liquid) driven side by side, using the reversible-scaling
formulation of the paper. Output a `coexistence_line.dat` with a forward/backward hysteresis
estimate, in the same spirit as `ts` writes `temperature_sweep.dat`.

**Non-goals (v1):** multi-component systems with different compositions in the two phases,
solid-solid boundaries (works in principle, not validated), Monte Carlo swaps, restart/resume
of an interrupted sweep, the paper's Woon/argon reproduction.

---

## 1. The method, in calphy's variables

Scaled system: `H_RS = K + lambda * U`, thermostatted at the *fixed* kinetic temperature `T0 = T_m`,
barostatted at the *scaled* pressure `P_RS`. Mapping to the real system (paper Eqs. 10, 11):

```
T      = T0 / lambda
P_real = P_RS / lambda
```

Coexistence is preserved under a perturbation `(d lambda, d P_RS)` when the reversible work is
equal in both cells (Eqs. 18-20):

```
(u_s - u_l) d lambda + (v_s - v_l) d P_RS = 0
```

with `u = <U_real>/N = <pe>/(N lambda)` in eV/atom (LAMMPS reports `lambda * U_real` under
`hybrid/scaled`, exactly as `integrate_rs` already divides by lambda) and `v = <V>/N` in
A^3/atom, `P_RS` converted bar -> eV/A^3 with the existing `EV_A3_TO_BAR`. Per-atom quantities
make the two cells free to have different N. The 3/2 k_B T ln(lambda) kinetic term of Eq. 12 is
identical per atom in both phases and cancels, so it never appears.

**Independent variable: `P_RS`, ramped linearly in time. Dependent variable: `lambda`.**

```
d lambda / d P_RS = -(v_s - v_l) / (u_s - u_l)              (Eq. 20 inverted)
```

This is the paper's own choice for the LJ test ("the scaled pressure P_RS was chosen as the
independent variable"), and it is forced on us by LAMMPS, see Section 2. The singular case is
`u_s = u_l` (never at a melting line); the alternative formulation (lambda independent, Eq. 23) is
singular at `v_s = v_l`, i.e. at a melting-curve maximum, and needs a time-dependent barostat
target that `fix npt` cannot provide. Section 6 lists it as an extension.

Diagnostic written per block, for the user and for the tests: the real-space slope
`dP/dT = dH / (T dV)` with `dH = (u_l - u_s) + P_real (v_l - v_s)`, which must agree with the
finite-difference slope of the produced curve.

### Discrete integration (one "block" = `n_block_steps` MD steps)

Block k spans `P_RS,k -> P_RS,k+1` (known in advance, linear ramp). Let
`f_k = -(v_s - v_l)/(u_s - u_l)` from the block-k averages, `h = P_RS,k+1 - P_RS,k`.

* Inside block k, lambda is ramped linearly from `lambda_k` to a *predicted* end value
  `lambda_k + h f_(k-1)` (held constant at `lambda_0 = 1` in the first block). This keeps the
  scaled Hamiltonian continuous in time instead of jumping at block boundaries.
* After the block, the corrected `lambda_(k+1) = lambda_k + h f_k` (Euler with the block average;
  optional trapezoid `h (f_k + f_(k-1))/2` behind `dcci.integrator`). The O(h^2) mismatch between
  predicted and corrected end value is the block-boundary perturbation; it vanishes with block size.
* Real-space point k: `T_k = T0/lambda_k`, `P_k = P_RS,k/lambda_k`.

Per-step updating (the paper) is the `n_block_steps: 1` limit and is only sensible in
`execution_mode: library`.

---

## 2. LAMMPS mechanics and the design constraints they impose (verified in the checkout)

| Fact (source) | Consequence for the design |
|---|---|
| `fix npt` accepts only numeric `Pstart Pstop` (`src/fix_nh.cpp` lines 141-240, `utils::numeric`); no `v_` variable support. Its ramp is `(step - beginstep)/(endstep - beginstep)` (lines 2367-2381). | A time-dependent barostat target must be a linear ramp fixed at fix creation. Hence `P_RS` is the independent, linear variable and the fix is created once per sweep with `iso P_RS,i P_RS,f`. Blocks are run with `run n_block start 0 stop N` so the ramp follows global progress, not the block. Re-creating the fix per block would reset the Nose-Hoover chain (eta, omega) every block and is ruled out. |
| `pair_style hybrid/scaled` looks its scale variable up *by name every compute* (`src/pair_hybrid_scaled.cpp` lines 73-95, 437-450). | lambda can be an equal-style variable redefined per block by the driver (`variable lam equal <a>+<b>*(step-<s>)/<n>`, only literal numbers, no `$(...)` immediate evaluation so it survives `ExecutableRunner` segment boundaries). |
| `fix ave/time` flushes its file at every output (`fix_ave_time.cpp` lines 683, 904); `ExecutableRunner` suffixes replayed `file` targets with `.seg<k>` and `read_timeseries` concatenates them. | Block averages of `pe/atoms`, `vol/atoms`, `v_lam` are read back through the existing runner contract (`fix fav all ave/time 1 n_block n_block ... file dcci.<cell>.<dir>_<i>.dat`, one row per block). No new data accessor is needed; `sync()` per block then read the last row. Works identically for `LibraryRunner`. |
| `fix nh` writes/reads its state to restart files (`fix_nh.cpp` 1247-1340) and `run ... start/stop` exists (`run.cpp` 68-73). | Executable mode: one block = one segment, barostat state carries over (this is already what `test_t2_nose_hoover_restoration` checks for a constant target). Needs one new check for a *ramped* target, Part 0. |
| `reset_timestep` is not in the runner vocabulary (`runner.py` `KNOWN_TOKENS`). | Add it to `ONE_SHOT_TOKENS` (restart files carry the timestep, so replay is unaffected), or have the driver track the step counter and pass `start S0 stop S0+N`. Prefer `reset_timestep 0` before each sweep for readable logs. |

Cost in executable mode: one `lmp` process launch per block per cell. For a 1e5-step sweep with
`n_block_steps: 1000` that is 100 launches per cell per direction, well under the MD time. Library
mode has no per-block cost and can use `n_block_steps` of 10-100.

---

## 3. Code changes

### 3.1 `calphy/input.py`

* New block `dcci: DynamicClausiusClapeyron(_StrictInput)`:

  | key | default | meaning |
  |---|---|---|
  | `n_block_steps` | 1000 | MD steps per integration block (the `h` of Section 1). |
  | `integrator` | `"euler"` | `"euler"` or `"trapezoid"` corrector. |
  | `stop_at_target_pressure` | `true` | End the sweep at the first block whose *real* pressure reaches `pressure[1]`. With `dT/dP > 0` (lambda < 1) the real pressure runs ahead of `P_RS`, so the ramp target `P_RS,f = pressure[1]` is conservative; with a negative slope the sweep ends at `P_RS = pressure[1]` and the reached real range is logged. |
  | `n_check_blocks` | 0 | Run the solid-fraction melt/freeze checks every this many blocks (0 = only at the end of each sweep). |
  | `hysteresis_tolerance` | 5.0 | K. Warn (and flag in `report.yaml`) when the backward sweep misses `T0` at `P_i` by more than this. |
  | `parallel_cells` | `false` | Run the two cells concurrently (two threads, each cell's LAMMPS with `queue.cores/2`) instead of one after the other with all cores. Same results, half the sweep wall time; falls back to sequential with one core. |

* `mode: dcci` validation: `pressure` must be a 2-list (`_pressure` = P_i, `_pressure_stop` = P_f;
  already parsed), `temperature` a scalar > 0 (`_temperature` is `T0 = T_m`),
  `n_switching_steps` gives `_n_sweep_steps` = N (existing semantics), `reference_phase` is
  ignored and forced to `"solid"` so `create_identifier` yields
  `dcci-<lattice>-solid-<T>-<Pi>`. `monte_carlo.n_swaps > 0`, `pair_mode: overlay` with
  composition changes, and `fe-qtb` are rejected with clear messages.

### 3.2 `calphy/dcci.py` (new): `class DynamicCCI`

Mirrors `MeltingTemp` in shape (a coordinator owning a `Solid` and a `Liquid`), but runs them in
lockstep. Simfolder layout: `dcci-.../` with `solid/` and `liquid/` sub-folders (each sub-job gets
its own `calphy.log`, `input_file.yaml`, `conf.equilibration.data`, `averaging.log.lammps`), the
coupled-sweep files and `report.yaml` at the top level.

1. `prepare_cells()`: build two sub-`Calculation`s from the parent dict (mode `dcci`,
   `reference_phase` solid/liquid, `temperature: T0`, `pressure: P_i`), instantiate `Solid` /
   `Liquid`, attach the parent log handlers as `MeltingTemp.run_jobs` does, enable the phase
   checks with `MeltingTemp._enable_phase_detection` semantics (factor that helper out of
   `MeltingTemp` into a module-level function).
2. `equilibrate_cells()`: `soljob.run_averaging()`, `lqdjob.run_averaging()`. This reuses the
   pressure-convergence machinery unchanged and yields `conf.equilibration.data`, `lx/ly/lz`,
   `natoms` for each cell at (T0, P_i). Melted/frozen cells raise the existing errors.
3. `run_sweep(direction, iteration)`: opens two runners, one per cell:
   * `create_object`, `hybrid/scaled` with `v_lam` (initialised to the sweep's starting lambda),
     `read_data` (equilibration conf for forward; `conf.dcci.forward_<i>.<cell>.data` for
     backward), `pair_coeff`, `mass`, `remap_box`;
   * warm start + equilibration under `fix npt ... iso P_RS,start P_RS,start` (reuse
     `_warm_start_steps`), `unfix`, `reset_timestep 0`;
   * `fix f1 all npt temp T0 T0 tdamp iso P_RS,start P_RS,end pdamp`, thermo, the `ave/time`
     fix, optional `dump` on `n_print_steps`;
   * block loop: redefine `v_lam` in both cells, `run n_block start 0 stop N` in both, `sync()`
     both, read the last `ave/time` row of each, `cce_step(...)`, append a row to
     `dcci.<dir>_<i>.dat`, log `T_k, P_k, dP/dT`, apply the stop rule and the periodic checks;
   * end: snapshot + melt/freeze check, `write_data`, `close`, `rotate_logs`.
   By default the two runners are stepped one after the other in the same Python process, each
   cell's `lmp` getting the full core count in turn; `dcci.parallel_cells` runs the two per-block
   sequences (ramp, run, sync) in two threads with half the cores each. There is no MPI coupling:
   the cells exchange information only through the driver, once per block.
4. `integrate(iteration)`: pure numpy (`integrators.integrate_dcci`), see 3.4.
5. `submit_report()` / `clean_up()`: `report.yaml` at the top level with
   `results: {t0, p_start, p_stop_requested, p_stop_reached, t_stop, n_blocks_used,
   hysteresis_t_at_p_start, hysteresis_high, coexistence_line: coexistence_line.dat}`,
   `metadata.yaml` with the de Koning DOI appended to `publications`.

Pure functions kept free of LAMMPS for unit testing: `cce_step(u_s, u_l, v_s, v_l, dp_rs_bar)`,
`predict_lambda(...)`, `plan_blocks(n_sweep, n_block)`, `lambda_ramp_command(name, l0, l1, step0, n)`.

### 3.3 Dispatch and plumbing

* `queuekernel.py` (`setup_calculation`, `main`) and `routines.py`: `mode == "dcci"` ->
  `DynamicCCI(...).calculate_coexistence_line()`; extend the "Mode should be either ..." messages.
* `runner.py`: `required_styles` adds `hybrid/scaled`, `npt`, `ave/time`, `temp/com` for `dcci`
  (liquid equilibration already pulls `ufm`); `reset_timestep` into `ONE_SHOT_TOKENS`
  (plus a row in `tests/test_command_vocabulary.py`).
* `.gitignore`: `dcci-*`.
* `postprocessing.py`: `read_coexistence_line(folder)` returning a DataFrame, and a small
  `plot_coexistence_line(folder)`; `gather_results` picks up `report.yaml` of `dcci` folders.

### 3.4 Output files and the integration step

* `dcci.<cell>.<forward|backward>_<i>.dat` (raw `ave/time`, one row per block):
  `step  pe_atom  vol_atom  lambda`.
* `dcci.<forward|backward>_<i>.dat` (driver, one row per block):
  `block step lambda T[K] P_RS[bar] P[bar] u_s u_l v_s v_l dPdT[bar/K]`.
* `coexistence_line.dat`: forward and backward curves interpolated onto a common real-pressure
  grid; columns `P[bar] T_forward[K] T_backward[K] T_mean[K] T_err[K]` where `T_err` is
  half the local hysteresis combined with the standard error over `n_iterations` replicas
  (same convention as `integrate_rs`).
* `conf.dcci.forward_<i>.<cell>.data`, `conf.dcci.backward_<i>.<cell>.data`.

### 3.5 Docs and examples

* `docs/source/inputfile.md`: `mode: dcci` under `mode`; new `dcci` block section; note on
  `pressure` as a 2-list and `temperature` as the known melting temperature at `pressure[0]`.
* `docs/source/outputfiles.md`: the files above.
* `examples/example_13`: Cu (Mishin `Cu01.eam.alloy`, already in `tests/`) melting line
  0 -> 100 kbar starting from a `melting_temperature` run at 0 bar; notebook plotting
  `coexistence_line.dat` against two independent `melting_temperature` points.
* Recommended two-step workflow documented: `melting_temperature` at `P_i`, then `dcci` with
  `temperature: <Tm>`. Chaining the two inside one calculation is Section 6.

---

## 4. Tests

Unit (no LAMMPS):

* `tests/test_dcci_integrator.py`: `cce_step` sign and units (bar -> eV/A^3); constant
  `dU, dV` synthetic data reproduces the analytic `lambda(P_RS)`; trapezoid beats Euler at
  second order; `coexistence_line` regridding and hysteresis; `n_iterations` error combination.
* `tests/test_dcci_input.py`: scalar `pressure` rejected, `temperature: 0` rejected, identifier
  and folder name, block defaults, rejection of swaps / qtb / composition-changing overlays.
* Golden command stream (`tests/conftest.py` `RecordingRunner` extended with a
  `read_timeseries` stub returning synthetic block rows): the exact per-cell command sequence for
  a 3-block forward + backward sweep, frozen under `tests/golden/`.

Integration (`-m lammps`, run with the `lammps-bin` env recipe in the project memory):

* `test_dcci_ramp_continuity`: T2-style check that `fix npt iso P0 P1` under
  `run n start 0 stop N` across three segments reproduces a single-run trajectory's pressure
  target (compare volume and `press` time series, not just energy).
* `test_dcci_smoke_cu`: Cu01, `repeat [3,3,3]`, tiny equilibration, `n_switching_steps: 2000`,
  `n_block_steps: 500`: all files present, no NaN, `T` monotonic with `P`, forward and backward
  `dP/dT` positive and within 30 % of `dH/(T dV)`, `report.yaml` complete. Marked `slow`.

Physics validation (manual, results go into the example notebook, not CI):

1. Cu01: `dcci` 0 -> 100 kbar with N = 2e4, 5e4, 1e5; convergence of `T(P_f)` with N
   (paper Fig. 4) and hysteresis at `P_i` versus N (paper Fig. 5). Compare against independent
   `melting_temperature` runs at 0 and 100 kbar (must agree within their error bars).
2. Stretch: LJ melting line `P* = 1 -> 170` against Agrawal & Kofke (paper Fig. 6) via
   `md.init_commands: ["units lj"]`. Blocked on the metal-unit assumptions in equilibration
   (`tolerance.pressure` in bar, `EV_A3_TO_BAR`), so only if a units-aware constant is cheap.

---

## 5. Parts, in order (one PR per part, no Claude attribution in commits)

| Part | Content | Acceptance |
|---|---|---|
| 0 | LAMMPS behaviour checks, kept as `tests/test_dcci_lammps_semantics.py` (`-m lammps`): ramped `fix npt` under `run n start 0 stop N` blocks reproduces the single-run pressure trajectory exactly (a control without start/stop collapses after the first boundary); `hybrid/scaled` picks up a redefined scale variable within a segment (`pre no post no`) and across a boundary; `ave/time 1 n n` writes exactly one row per block after a restart. | Done 2026-09-10, all three green: Section 2 assumptions hold. |
| 1 | `input.py` block + validation, dispatch, `required_styles`, `reset_timestep` token, `.gitignore`. | `test_dcci_input.py`, vocabulary test, full suite green. |
| 2 | `dcci.py`: cells, equilibration, forward sweep, `cce_step`, raw and driver `.dat` files. | Integrator unit tests; golden stream; Cu smoke test forward-only. |
| 3 | Backward sweep, `integrate_dcci`, `coexistence_line.dat`, `report.yaml`, hysteresis flag, `n_iterations`, postprocessing reader/plot. | Full smoke test; unit tests for regridding. |
| 4 | Docs, example_13 with validation notebook. | Done 2026-09-10: Cu01 lines to 5, 10 and 158 GPa (53/155/29 blocks) agree with each other to ~10 K, close to <0.01-2.4 K, track the published Simon fits within a few percent, and match direct fe crossings of both phases at 5 and 10 GPa to 2 K once the 13 K offset of the starting `melting_temperature` (1340 K vs 1353 K direct) is accounted for. Version bump left to the owner. |

---

## 6. Extensions (not in v1, listed so the v1 design does not block them)

* **lambda as independent variable** (Eq. 23, needed near a melting-curve maximum where
  `dV -> 0`): requires a per-block barostat retarget without state loss. Route: make each block a
  restart cycle in both runners and re-specify `fix f1` with the same ID after `read_restart`
  (executable mode does this already; `LibraryRunner` would need a `write_restart/clear/
  read_restart` + `SessionState` replay). Verify `Modify::add_fix` re-applies restart state when
  the fix is re-specified before the first `run`.
* **Chained initial condition**: `dcci.initial_condition: melting_temperature` runs `MeltingTemp`
  at `P_i` first and feeds `T_m` in.
* **Per-step updating in library mode** via `run 1 pre no post no` (the paper's scheme).
* **Solid-solid boundaries**: two `Solid` cells with different lattices; only the equilibration
  differs.

---

## 7. Risks, assumptions, and one pre-existing finding

* **Accuracy is set by the starting point.** The integration itself reproduces the pressure
  dependence to a few kelvin; an error `dT0` in the starting coexistence temperature propagates
  as roughly `dT0 * T/T0`. At 0 bar `melting_temperature` (30000-step sweeps, one iteration)
  gave 1340 K where direct fe crossings of both phases give 1353 K, and the d-CCI line carries
  that 1 % offset unchanged. Spend the effort on the starting point (n_iterations, longer
  sweeps, or a direct fe crossing) rather than on the sweep.
* **Hysteresis is not a discretisation check.** Forward and backward Euler errors mirror each
  other, so a round trip can close while each direction is biased; averaging the two directions
  cancels the leading error. Convergence in `n_block_steps` / `n_switching_steps` has to be
  checked by rerunning (Cu01: 1000-step blocks to 158 GPa are 11 K below 500-step blocks at
  10 GPa; 53 and 155 blocks of 500 steps agree at 5 GPa).
* **A runner fragility the mode exposes.** Two of the validation runs (a `melting_temperature`
  at 200 kbar and the 155-block d-CCI) died with `LammpsExecutionError: exited with return code
  1` on a segment whose LAMMPS log ends cleanly (`Total wall time`, no `ERROR`), both under
  several concurrent `mpirun -np 4` jobs; the d-CCI lost only its report, the sweep files were
  complete and `integrate_dcci` recovered the line. A mode that launches hundreds of segments
  needs the ExecutableRunner to treat a clean log with a nonzero `mpirun` exit as a warning, or
  to retry the segment from its restart; worth a follow-up in `runner.py`.

* **Noise in the slope.** Small cells give `dV/atom` comparable to its fluctuation; the error
  integrates as a random walk that shrinks with N and is exposed by the hysteresis. Default
  `n_block_steps: 1000` plus 500+ atoms per cell matches the paper's regime.
* **Cells changing phase mid-sweep.** If the solid melts (or liquid freezes) the slope collapses
  silently; the `n_check_blocks` checks and the end-of-sweep checks raise `MeltedError` /
  `SolidifiedError`, with the block index in the log. No automatic recovery in v1.
* **Executable-mode overhead** scales with blocks, not steps; documented in the `n_block_steps`
  entry.
* **Assumption:** both cells contain the same composition (congruent melting). Multi-component
  input is accepted only when every element appears in both cells with the same fraction.
* **Finding in the existing `ts` / `tscale` sweeps (fixed on branch `sweep-barostat`).**
  `_reversible_scaling_forward`, `_reversible_scaling_backward` computed `pf = lf * pi` and
  logged "P pi -> pf", but the barostat that stayed active through the sweep was created with
  `iso pi pi`. The reversible-scaling relation (paper Eq. 11, and the `U + P V` integrand that
  `integrate_rs` uses) requires `P_RS = lambda * P` during the sweep. With a constant `P_RS = P`
  the sweep samples the real system on the path `(T0/lambda, P/lambda)` while the integration
  labels it as the isobar `P`. Because the `P V` term in the integrand uses the fixed target
  pressure, the leading error in the integrand is `-(alpha T) V (P/lambda - P)` plus a
  second-order compressibility term, i.e. small for a stiff solid at tens of kbar (Cu, 50 kbar,
  400 -> 700 K: about 1 meV/atom, invisible in a quick fe cross-check) and growing with
  `alpha T`, `V`, `P` and the temperature span; liquids at 100+ kbar are where it matters.
  At `pressure: 0` it is exactly zero. `tscale` and the pre-flight scan had the opposite
  mismatch (a scaled-pressure ramp applied to a real-temperature ramp, since v1.0.0), and the
  `tscale` backward sweep had been quenching to `T0` instead of ramping `Tf -> T0` since
  "switch to stepped ts" (b286426, released in 1.7.4), which corrupts every `tscale` result
  since then at any pressure. The `dcci` sweep must use the ramped `P_RS`, as in Section 2.
