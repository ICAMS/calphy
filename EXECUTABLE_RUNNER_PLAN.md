# Implementation plan: replace pylammpsmpi with an executable-driven runner

**Status:** approved design, ready for implementation
**Date:** 2026-07-08
**Owner decisions (locked, do not revisit without asking):**
1. The staged CLI (`calphy_run_averaging`, `calphy_process_averaging`, `calphy_run_integration`,
   `calphy_process_integration`) and the `script_mode` input key are **removed**. `script_mode: True`
   in an input file must raise a clear validation error. The new runner keeps a `dry_run` mode that
   emits segment scripts for debugging.
2. The `lmp=` injection parameter on `Phase`/`Solid`/`Liquid`/`Alchemy`/`create_object` is **removed**.
3. LAMMPS binary resolution order: input key `lammps_executable` → env var `CALPHY_LAMMPS_EXECUTABLE`
   → `lmp` on `PATH`. Fail at startup with a clear message. `mpi_executable` resolves the same way
   (input key → `CALPHY_MPI_EXECUTABLE` → `mpirun` on PATH), and is only used when `queue.cores > 1`.
4. Numerical baselines are produced by **running current `main` (with pylammpsmpi) once** before the
   refactor (Part 0), and are committed to the repo.

**Goal:** calphy drives an external LAMMPS binary for every calculation. All physics capability is
retained exactly as-is: every mode (`fe`, `ts`, `tscale`, `pscale`, `alchemy`, `composition_scaling`,
`melting_temperature`), both reference phases, melting cycle, iterative pressure / spring-constant
convergence, melt/solidify checks, phase-transition pre-scan, QTB, overlay potentials, MC swaps,
`n_iterations > 1`.

**Non-goals:** no physics changes, no new features, no performance work, no Windows support.

---

## Rules for the implementing agent

- Work strictly part by part, in order. Each part ends with its acceptance criteria green.
  Do not start part N+1 with part N red.
- One commit (or PR) per part, message `exec-runner part N: <title>`.
- Run the full test suite after every part: `pytest calphy/tests -x -q` (plus markers noted per part).
- Never edit files under `tests/golden/` or `tests/baselines/` after they are first committed
  (Parts 0 and 2) unless the plan explicitly says so. If a golden mismatch appears later, the code
  is wrong, not the golden — investigate the code.
- If LAMMPS semantics are in doubt (what a restart file preserves, command ordering constraints),
  **stop and ask the user** — do not guess. The doc to consult first: LAMMPS manual pages for
  `read_restart`, `write_restart`, and the specific fix.
- Any command string whose first token is not in the classification table (Appendix A) must make the
  runner raise — never silently pass through an unknown command.
- Keep the existing code style: `%`-formatting for LAMMPS command strings, existing logging patterns.

---

## Architecture in one page (read before coding)

Today, calphy talks to LAMMPS through `pylammpsmpi.LammpsLibrary`, created in
`calphy/helpers.py::create_object`. The audit finding this plan rests on: **calphy only ever sends
command strings and reads results back from files** (`avg.dat`, `msd.dat`, `forward_*.dat`,
`ts.*.dat`, dump files, `conf.*.data`). The live library instance is used for exactly one thing:
keeping LAMMPS state alive between the points where Python reads a file and decides what to do next.

The replacement, `ExecutableRunner`, preserves that one property with segments:

- `command(str)` appends to a buffer (and updates a `SessionState` tracker).
- `sync()` = flush: write the buffered commands to `calphy.seg<k>.lmp`, prepend restart/replay
  header (segment > 0), append `write_restart`, run
  `[mpirun -np N] lmp -in calphy.seg<k>.lmp -log calphy.seg<k>.log -screen none <cmdargs>`
  as a subprocess in the sim folder, check for errors, continue.
- Continuity between segments: LAMMPS restart files preserve box, atoms, velocities, groups,
  timestep, and the internal state of fixes **that are redefined with identical IDs** (Nose-Hoover
  chains). Everything else — pair_style/pair_coeff/mass, variables, computes, fixes, thermo — is
  replayed from `SessionState` at the top of each segment.
- `sync()` is called at every site where Python reads a file produced by the current session
  (Appendix B); `close()` does a final flush.
- Runs that must be continuous (ti/spring switching legs, RS/tscale/pscale sweeps) contain no reads,
  therefore are never split. This is an invariant, not a hope — the runner enforces the assertions
  in Part 3.

The driver process (`calphy_kernel`) stays alive for the whole calculation, exactly like today's
interactive mode — so all Python-side loops, checks, and orchestration code keep working unchanged.

---

## Part 0 — Baseline capture on current `main` (needs working pylammpsmpi env)

Purpose: freeze today's behavior as reference data before touching anything.

### 0.1 Record environment and test status

- Create `tests/baselines/BASELINE_INFO.md` recording: calphy commit hash, LAMMPS version
  (`python -c "from pylammpsmpi import LammpsLibrary; ..."` or `lmp -h | head`), pylammpsmpi
  version, platform, date.
- Run `pytest calphy/tests -q`; record pass/fail counts in the same file. Pre-existing failures are
  noted, not fixed.

### 0.2 Reference calculations (Cu EAM, tiny but real)

Use `tests/Cu01.eam.alloy` and `tests/conf1.data`. Create input files under
`tests/baselines/inputs/`, one per scenario, all with:
`queue.cores: 2` (or 1 if no MPI available), `n_equilibration_steps: 2500`,
`n_switching_steps: 5000`, `md.n_small_steps: 1000`, `n_iterations: 2`, fixed
`spring_constants: null` (so the convergence path runs).

| id | mode | reference_phase | key settings |
|----|------|-----------------|--------------|
| B1 | fe | solid | fcc Cu, 500 K, 0 bar |
| B2 | fe | solid | 500 K, 10000 bar (finite-pressure path) |
| B3 | fe | solid | 500 K, `fix_lattice: True` (constrained path) |
| B4 | fe | liquid | 1400 K, `melting_cycle: True` |
| B5 | fe | liquid | 1400 K, `melting_cycle: False` (rattle path) |
| B6 | ts | solid | 400 → 700 K |
| B7 | ts | liquid | 1400 → 1600 K |
| B8 | tscale | solid | 400 → 700 K |
| B9 | pscale | solid | 500 K, 0 → 50000 bar |
| B10 | alchemy | — | Cu EAM → Cu EAM with slightly perturbed second potential (any second `pair_style eam/alloy` file; scaling one potential into another is enough to exercise the path) |
| B11 | fe + QTB | solid | 300 K, `mode: fe-qtb` equivalent settings |
| B12 | ts + prescan | solid | B6 plus `phase_transition_detection.mode: warn` |

Skip `melting_temperature` (too slow for a baseline; it is a composition of B6+B7 mechanics) and
`composition_scaling` (no alloy fixture in repo; its MD command stream is covered by alchemy — its
Python pre/post-processing already has dedicated unit tests).

For each scenario, run on current main (`calphy_kernel -i input.yaml -k 0`) and commit to
`tests/baselines/<id>/`:
- `report.yaml` (the primary reference: fe, w, spring constants, vol, error bars),
- `calphy.log`,
- the input yaml.

### 0.3 Fixture capture for mock tests

From the B1, B4, B6, B12 run folders, copy into `tests/fixtures/`:
- `avg.dat` (from B1), `msd.dat` (from B1),
- `traj.equilibration_stage1.dat` and one melted + one solid dump snapshot
  (`traj.melt` from B4, `traj.equilibration_stage2.dat` from B1),
- `forward_1.dat`, `backward_1.dat` (from B1),
- `ts.forward_1.dat`, `ts.backward_1.dat` (from B6),
- `prescan.forward.dat` (from B12),
- `conf.equilibration.data` (from B1).

### Acceptance criteria (Part 0)

- `tests/baselines/` contains all 12 scenario folders with `report.yaml` + inputs + `BASELINE_INFO.md`.
- `tests/fixtures/` contains the files listed above.
- Nothing under `calphy/` (source) has changed.

---

## Part 1 — Decoupling micro-changes (safe, independent)

### 1.1 Kill the only library reads

- `calphy/phase.py::check_if_melted` (~line 311): replace `lmp.natoms` with `self.natoms`.
- `calphy/phase.py::check_if_solidfied` (~line 334): same.
- `self.natoms` is already set in `Phase.__init__` from `calc._natoms`; verify it is an int > 0 for
  all code paths that call these checks (it is — it comes from the input structure).

### 1.2 Kill the magic-attribute interface

The runner interface must be exactly `command()` + lifecycle methods. Remove every non-`command`
LAMMPS call:
- `calphy/liquid.py::rattle_structure`: `lmp.velocity("all create", X, Y)` →
  `lmp.command("velocity all create %f %d" % (X, Y))`; `lmp.run(n)` → `lmp.command("run %d" % n)`.
- `calphy/liquid.py::melt_structure`: same two conversions (~lines 117, 123).
- `calphy/liquid.py::run_averaging`: `lmp.velocity(...)` (~line 198) → `lmp.command(...)`.
- Verify with: `grep -rnE "lmp\.(velocity|run)\(" calphy/` → no hits.

### 1.3 Dead import

- `calphy/helpers.py` line 32: delete `from lammps import lammps` (unused).

### Acceptance criteria (Part 1)

- `grep -rn "lmp\.natoms\|lmp\.velocity(\|lmp\.run(" calphy/` → no hits.
- Full test suite status identical to Part 0 record.

---

## Part 2 — RecordingRunner + golden command streams (still on old code)

Purpose: freeze the exact sequence of LAMMPS commands each driver path produces. These goldens are
the contract that "capability is retained as-is": after the refactor, the same drivers must emit the
same logical command stream.

### 2.1 Test infrastructure (`tests/conftest.py` additions)

- `class RecordingRunner`: `self.commands = []`; `command(s)` appends `" ".join(s.split())`
  (normalize whitespace so cosmetic padding differences don't break diffs); `close()`, `clear()` are
  no-ops; `write(path)` writes commands to path (needed because old script-mode branches call it —
  they are not exercised here, but keep it harmless).
- Fixture `recorded_job(tmp_path, monkeypatch)`:
  - `np.random.seed(42)` — all seeds in commands become deterministic.
  - monkeypatch `calphy.helpers.create_object` to return a shared `RecordingRunner` (signature-
    compatible: accept and ignore all current kwargs).
  - copy needed fixture files from `tests/fixtures/` into the sim folder before the driver call, so
    every `np.loadtxt` / `Trajectory` read succeeds.
  - monkeypatch `calphy.helpers.find_solid_fraction` with a test-supplied sequence
    (e.g. `[N, N, 0]` to force the melting loop to iterate 3 times); default: returns "all solid"
    for solid checks and 0 for liquid checks so no error paths trigger.
- Helper `assert_golden(commands, "name")`: compares against `tests/golden/name.txt`;
  on first run (file missing) writes it (then commit); afterwards read-only comparison with a
  readable diff on failure.

### 2.2 Golden scenarios (each is one test in `tests/test_golden_streams.py`)

Build a minimal `Calculation` object per scenario (reuse baseline input yamls from Part 0 via
`read_inputfile`). Call the driver methods directly, then snapshot `runner.commands`:

| golden file | driver calls |
|---|---|
| `solid_fe_averaging.txt` | `Solid.run_averaging()` (zero-pressure path) |
| `solid_fe_averaging_finite_p.txt` | same, 10000 bar input |
| `solid_fe_averaging_fixlattice.txt` | same, `fix_lattice: True` |
| `solid_fe_integration.txt` | `Solid.run_integration(iteration=1)` after setting `job.lx/ly/lz/k` from fixture values |
| `liquid_fe_averaging_meltcycle.txt` | `Liquid.run_averaging()`, solid-fraction sequence forces 3 melt attempts |
| `liquid_fe_averaging_rattle.txt` | `Liquid.run_averaging()`, `melting_cycle: False` |
| `liquid_fe_integration.txt` | `Liquid.run_integration(iteration=1)` |
| `ts_forward_solid.txt` | `Phase._reversible_scaling_forward(iteration=1)` |
| `ts_backward_solid.txt` | `Phase._reversible_scaling_backward(iteration=1)` |
| `ts_uniform_temperature.txt` | forward with `lambda_schedule: uniform_temperature` |
| `ts_mc_swaps.txt` | forward with `monte_carlo.n_swaps > 0` and 2 swap types |
| `tscale.txt` | `Phase.temperature_scaling(iteration=1)` |
| `pscale.txt` | `Phase.pressure_scaling(iteration=1)` |
| `alchemy_averaging.txt` / `alchemy_integration.txt` | `Alchemy.run_averaging()` / `run_integration(1)` |
| `prescan.txt` | `Phase.scan_temperature_range()` (fixture prescan.forward.dat, mode `warn`) |
| `qtb_averaging.txt` / `qtb_integration.txt` | Solid QTB variants |
| `overlay_potential.txt` | Solid averaging with `pair_mode: overlay` (two pair styles) |

Notes:
- All scenarios run with `script_mode` **False** (library path) — that is the behavior contract.
- Loops that converge on data (pressure, spring constant): the static fixtures make them converge on
  cycle 1; that is fine — the multi-cycle repetition pattern is exercised with real LAMMPS in Part 6.
- Also add `tests/test_command_vocabulary.py`: iterate over every golden file, split each line's
  first token, assert it is in the Appendix A vocabulary list (import the table from
  `calphy.runner` once Part 3 lands; until then keep a literal copy in the test and reconcile in
  Part 3).

### Acceptance criteria (Part 2)

- All golden tests pass twice in a row (determinism check — run pytest twice).
- `tests/golden/*.txt` committed; every file non-empty; no file contains the substring
  `object at 0x` (sanity against repr leakage).

---

## Part 3 — `calphy/runner.py`: ExecutableRunner + SessionState

New module, no changes to existing modules yet (wire-in is Part 5). Everything here is unit-testable
without LAMMPS.

### 3.1 Errors (`calphy/errors.py`)

Add:
```python
class LammpsExecutionError(RuntimeError):
    """LAMMPS binary exited abnormally. Carries segment script path, log path, log excerpt."""
class RunnerStateError(RuntimeError):
    """Invalid command sequence for segmented execution (unknown command, dump across sync, ...)."""
```

### 3.2 Binary resolution (`calphy/runner.py`)

```python
def resolve_lammps_executable(explicit: str | None) -> str   # order: explicit -> $CALPHY_LAMMPS_EXECUTABLE -> shutil.which("lmp")
def resolve_mpi_executable(explicit: str | None) -> str      # order: explicit -> $CALPHY_MPI_EXECUTABLE -> shutil.which("mpirun")
```
- Resolved path must exist and be executable (`shutil.which` on the value; also accept absolute
  paths). On failure raise `ValueError` with a message that names all three lookup steps and shows
  how to set each.
- MPI resolution only runs when `cores > 1`.

### 3.3 SessionState

A pure-Python tracker fed every command string. It parses **only the first token** (plus id/name
tokens where noted) and maintains:

- `init_commands: list[str]` — the fixed header block (`units`, `boundary`, `atom_style`,
  `timestep`, `box`) captured verbatim from the first segment, in order.
- `pair_block: list[str]` — most recent contiguous potential definition: all `pair_style`,
  `pair_coeff`, `mass` commands. Rule: a new `pair_style` command **clears** the previous pair_block
  (calphy redefines the full block each time it switches potentials); `pair_coeff`/`mass` append.
- `groups: dict[name, command]` (insertion-ordered; redefinition overwrites).
- `computes: dict[id, command]`; `uncompute id` removes.
- `variables: dict[name, command]` (insertion-ordered; redefinition overwrites in place —
  **keep original insertion position**, calphy redefines `flambda`/`blambda`/`lambda` between runs).
- `fixes: dict[id, {"command": str, "fix_modify": list[str]}]`; `unfix id` removes;
  `fix_modify id ...` appends to the entry (error if id unknown).
- `dumps: dict[id, command]`; `undump id` removes.
- `thermo: str | None`, `thermo_style: str | None`, `echo: str | None` — last value wins.
- `timestep: str | None` — last value wins (also part of init; replay uses the latest).
- `has_box: bool` — set True after `read_data`/`read_restart` seen.

Classification (authoritative table in Appendix A): every incoming command is either
*init*, *sticky* (tracked above), *one-shot* (buffered for this segment, never replayed), or
*unknown* → raise `RunnerStateError` naming the command.

Immediate-evaluation guard: a sticky `variable` command whose text contains `$(` captures a value at
definition time and cannot be faithfully replayed. Track a flag; if such a variable is still live
when a replay header is generated → raise `RunnerStateError` ("immediate-evaluation variable X
cannot cross a segment boundary"). calphy's only such variable (`step0`) never crosses one.

### 3.4 ExecutableRunner

```python
class ExecutableRunner:
    def __init__(self, *, binary, mpi_command, cores, cmdargs, directory,
                 restart_name="calphy.restart", dry_run=False, timeout=None): ...
    def command(self, s: str) -> None      # normalize whitespace, classify, buffer
    def sync(self) -> None                 # flush buffer as one segment and execute it
    def close(self) -> None                # final sync (if buffer non-empty), keep logs
    def rotate_logs(self, stage_name: str) -> None   # see 3.6
    @property
    def logical_commands(self) -> list[str]  # every command ever passed to command(), verbatim order
```

Segment emission on `sync()` (skip entirely if buffer is empty):

- **Segment 0** (`has_box` was set by buffered commands): emit buffered commands verbatim, append
  `write_restart <restart_name>`.
- **Segment k>0**: emit, in this order:
  1. `init_commands` **minus** the `box` command and minus `boundary` (both stored in the restart;
     re-issuing `box` after a box exists is an error). Keep `units`, `atom_style`, `timestep`.
  2. `read_restart <restart_name>`
  3. `pair_block` (style, coeffs, mass — required after read_restart for eam/alloy etc.)
  4. `groups` values (idempotent re-definition of type-based groups)
  5. `computes` values
  6. `variables` values (original insertion order)
  7. `fixes` values in original insertion order, each followed by its `fix_modify` lines
  8. `thermo_style`, then `thermo`, then `echo` (if set)
  9. the buffered one-shot/new commands for this segment
  10. `write_restart <restart_name>`
- File-output redirection (see 3.5) is applied to replayed fix commands at step 7.
- Assertions before emitting a replay header: no live `dumps` (raise `RunnerStateError` — calphy
  always undumps before reading), no live immediate-eval variables.

Execution:

- Write segment to `<directory>/calphy.seg<k>.lmp`.
- argv: `[mpi, "-np", str(cores)] if cores > 1 else []` + `[binary, "-in", seg, "-log", seglog,
  "-screen", "none"]` + `cmdargs` (normalize `cmdargs` from str/list exactly as today's
  `create_object` does).
- `subprocess.run(argv, cwd=directory, capture_output=True, text=True, timeout=timeout)`.
- Failure = nonzero returncode, or timeout, or the string `ERROR` in the last 50 lines of the
  segment log. On failure raise `LammpsExecutionError` whose message includes: segment script path,
  log path, and the last 15 lines of the log (or stderr if the log is missing).
- `dry_run=True`: write the segment file, do not execute, do not check anything.

### 3.5 File-append semantics across segments

`fix ave/time ... file X` and `fix print ... file X` truncate `X` when the fix is (re)defined, but
calphy's convergence loops read the **cumulative** history of `avg.dat`/`msd.dat` across cycles.

- When a sticky fix command containing ` file <name>` is **replayed** in segment k, rewrite that
  token pair to ` file <name>.seg<k>`.
- Provide module-level helper:
  ```python
  def read_timeseries(directory, name, usecols, skiprows-comment-aware) -> np.ndarray
  ```
  which loads `<name>` plus every `<name>.seg*` in ascending segment order, strips each file's
  header comments, and vstacks. It must reproduce exactly what `np.loadtxt(name, usecols=...)`
  returns today when only one file exists.
- `fix print` never crosses a boundary in calphy (defined and unfixed around each sweep); if a
  replayed `fix print ... file` is encountered, raising `RunnerStateError` is acceptable and
  preferred over silent redirection. (`fix print` also supports `append` — do not use it; error
  instead, it means an invariant broke.)

### 3.6 Log rotation

Phase code today renames `log.lammps` to stage-named logs (`averaging.log.lammps`,
`integration.log.lammps`, `melted_error.log.lammps`, ...). With per-segment logs:

- `rotate_logs(stage_name)`: concatenate all segment logs produced since the last rotation (with
  `# --- segment k: calphy.seg<k>.lmp ---` separators) into `<stage_name>.log.lammps` in the sim
  folder, then reset the internal "since last rotation" list. Segment `.lmp` scripts are kept
  (they are small and are the debugging story); the rolling restart file is overwritten each
  segment and deleted on `close()`.

### 3.7 Unit tests (`tests/test_runner.py`, no LAMMPS binary needed)

State machine:
- every vocabulary token classifies without error; an unknown token (`min_style foo`) raises.
- fix defined → header of next segment contains it; unfix → absent; redefined fix keeps position.
- `fix_modify` attaches to its fix and replays immediately after it; `fix_modify` for unknown id raises.
- variable redefinition keeps original position; `$( )` variable crossing a boundary raises.
- new `pair_style` clears old pair block; overlay block (2 styles) replays completely.
- live dump at sync raises; dump+undump within a segment is fine.
- `file avg.dat` replay rewritten to `file avg.dat.seg1`, `seg2`, ...; non-replayed occurrence
  untouched in segment 0.
- segment 0 has no `read_restart` and ends with `write_restart`; segment 1 header order matches 3.4
  exactly (assert full expected text for one synthetic sequence).

Execution layer (monkeypatch `subprocess.run`):
- argv construction with/without MPI, cmdargs as str and list, cwd, log naming.
- nonzero rc → `LammpsExecutionError` includes log excerpt (write a fake seg log in tmp dir).
- `ERROR:` in log with rc 0 → still raises.
- dry_run writes files, calls nothing (assert the monkeypatched run was not invoked).
- resolution: explicit > env > PATH (use monkeypatched env and a tmp dir with a fake `lmp`);
  helpful error message when nothing resolves.

`read_timeseries`:
- single-file equivalence with `np.loadtxt`; multi-file concat ordering; comment-header stripping
  (use the real `avg.dat` fixture from Part 0 plus a synthetic `.seg1`).

### Acceptance criteria (Part 3)

- `tests/test_runner.py` green; module imports without pylammpsmpi installed
  (`python -c "import calphy.runner"` in an env without pylammpsmpi — simulate with a test that
  asserts `"pylammpsmpi" not in sys.modules` after import).
- Coverage of `calphy/runner.py` ≥ 90% (`pytest --cov=calphy.runner`).

---

## Part 4 — Preflight capability check

### 4.1 Implementation (`calphy/runner.py`)

```python
def get_binary_styles(binary) -> dict[str, set[str]]   # {"pair": {...}, "fix": {...}, "compute": {...}}
def preflight(calc, binary) -> None                    # raises ValueError listing missing styles
```
- `get_binary_styles`: run `[binary, "-h"]` (timeout 60 s), parse the style listing sections of the
  help output (`* Pair styles:`, `* Fix styles:`, `* Compute styles:` — verify exact headings
  against a real binary's output during implementation and keep the parser tolerant of column
  wrapping). Cache per (path, mtime) in a module dict.
- Required styles by calculation (derive from `calc`):
  - always: user pair styles (`calc._pair_style_names`), fixes `nve nvt npt langevin ave/time print momentum`, computes `msd temp/com pair`;
  - `reference_phase == "solid"` and mode in (fe, ts, tscale): fix `ti/spring` (package EXTRA-FIX);
  - `reference_phase == "liquid"` or mode `melting_temperature`: pair `ufm` (EXTRA-PAIR) and pair `hybrid/scaled`;
  - mode in (ts, tscale) or overlay potentials: pair `hybrid/scaled`;
  - `monte_carlo.n_swaps > 0`: fix `atom/swap` (MC);
  - QTB modes: fix `qtb` (QTB).
- Error message format: one line per missing style with the LAMMPS package that provides it and a
  closing hint ("rebuild LAMMPS with `-D PKG_EXTRA-FIX=yes` or use a full binary").
- Escape hatch: env var `CALPHY_SKIP_PREFLIGHT=1` skips (log a warning).
- Call site: wired in Part 5 (once per calculation, before the first segment executes).

### 4.2 Tests (`tests/test_preflight.py`)

- Canned `-h` outputs (full build / minimal build) as fixture strings; assert parsing, assert the
  missing-style error message content for a liquid calc against a minimal build.
- Cache behavior: second call does not re-invoke subprocess (monkeypatch counter).
- `CALPHY_SKIP_PREFLIGHT` honored.

### Acceptance criteria (Part 4)

- `tests/test_preflight.py` green. Still zero changes to phase/solid/liquid/alchemy.

---

## Part 5 — Wire-in and deletion (the big one)

Follow this order exactly. Run the golden tests (Part 2) after step 5.8 and keep them green through
the rest.

### 5.1 `helpers.py::create_object` — new construction

Replace the body: signature becomes
`create_object(calc, directory)` and it returns a configured `ExecutableRunner`:
- `binary = resolve_lammps_executable(calc.lammps_executable)`
- `mpi = resolve_mpi_executable(calc.mpi_executable) if calc.queue.cores > 1 else None`
- `cores = calc.queue.cores`, `cmdargs = calc.md.cmdargs`, `dry_run = False`
- run `preflight(calc, binary)` (Part 4) — once per calculation; cache makes repeats free.
- then emit the same five init commands as today (`units metal`, `boundary p p p`,
  `atom_style atomic`, `timestep X`, `box tilt large`) with the same `init_commands` override merge
  logic (keep that block verbatim).
- Delete the `script_mode` and `lmp` parameters and the `LammpsLibrary` branch; delete the
  `-screen none` cmdargs fiddling (runner owns `-screen`).
- Update every caller (`solid.py` ×3, `liquid.py` ×2, `alchemy.py` ×2, `phase.py` ×5) to the new
  signature: `lmp = ph.create_object(self.calc, self.simfolder)`.
- Delete `class LammpsScript` from helpers.py and the `from pylammpsmpi import LammpsLibrary` import.

### 5.2 Insert sync points (Appendix B is the authoritative list)

- `Phase.dump_current_snapshot`: append `lmp.sync()` after the `undump` — this single insertion
  covers every melt/solidify check and `melt_structure`.
- `Phase.run_iterative_pressure_convergence`: `lmp.sync()` immediately after the
  `run n_small_steps` inside the cycle loop, before reading `avg.dat`.
- `Phase.run_iterative_constrained_pressure_convergence`: same position, before `process_pressure()`.
- `Solid.run_iterative_spring_constant_convergence`: `lmp.sync()` after the in-loop `run`, before
  `analyse_spring_constants()`.
- `Phase.scan_temperature_range`: `lmp.sync()` after `unfix f2`, before reading `prescan.forward.dat`
  (then `lmp.close()` where `lammps_close` is today).
- No other sites. Do **not** add syncs inside integration/sweep methods — their reads happen after
  `close()`.

### 5.3 Cumulative-file reads → `read_timeseries`

Replace `np.loadtxt` on session-produced averaging files (and only those):
- `Phase.run_iterative_pressure_convergence` (`avg.dat`),
- `Phase.process_pressure`, `Phase.finalise_pressure` (`avg.dat`),
- `Solid.analyse_spring_constants` (`msd.dat`).
Keep plain `np.loadtxt` everywhere else (`forward_*.dat`, `ts.*.dat`, `prescan.forward.dat` are
written by a single `fix print` within one segment).

### 5.4 Delete the `minimal` twins and dispatchers

- `phase.py`: delete `run_minimal_pressure_convergence`,
  `run_minimal_constrained_pressure_convergence`; `run_pressure_convergence` and
  `run_constrained_pressure_convergence` become direct calls (or rename the `iterative` methods to
  the dispatcher names and delete the dispatchers — pick the rename, it keeps call sites stable).
- `solid.py`: delete `run_minimal_spring_constant_convergence`, `run_minimal_averaging`,
  `process_averaging_results`; `run_spring_constant_convergence` and `run_averaging` lose their
  dispatch (keep the interactive bodies).
- `solid.py::analyse_spring_constants` and `phase.py::process_pressure`/`finalise_pressure`: delete
  the `if self.calc.script_mode:` ncount branch — keep the interactive branch. **Note the existing
  branches are swapped in `analyse_spring_constants` (the script branch uses `n_small_steps`, the
  interactive one `n_equilibration_steps`) — keep exactly the behavior of `script_mode=False`,
  i.e. `ncount = n_equilibration_steps // (n_every_steps * n_repeat_steps)`.**

### 5.5 Delete every remaining `script_mode` branch

`grep -rn "script_mode" calphy/` and resolve each:
- `liquid.py`: melting-cycle guard (delete), averaging dump/check guard (keep the check code,
  remove the condition), script-write blocks at the end of `run_averaging`/`run_integration`
  (delete), `create_object` kwargs (gone in 5.1).
- `solid.py`: `run_averaging` dispatch (5.4), end-of-`run_integration` script-write else-branch
  (delete; keep close + `rotate_logs("integration")`).
- `alchemy.py`: three branches (averaging check guard, script-write blocks) — same treatment.
- `phase.py`: `_reversible_scaling_forward`/`_reversible_scaling_backward` — remove the
  `script_mode=...` kwargs and the `if self.calc.script_mode:` write-blocks (the backward one at
  ~line 1634 and the stability-check guard at ~1544 — keep the check, drop the condition);
  `temperature_scaling` write-block (~2107, delete); `scan_temperature_range` kwarg (gone).
- `pressure_scaling`: no script branch existed, but it also never closed the lmp — add
  `self.lammps_close(lmp)` + `rotate_logs("pressure_scaling")` at the end (bug fix, matches
  tscale's structure).
- `kernel.py::run_jobs`: delete the whole `if script_mode:` arm (lines ~80–139); keep the scheduler
  arm as the only path.

### 5.6 Remove `lmp=` injection and simplify close

- Remove the `lmp` parameter from `Phase.__init__`, `Solid.__init__`, `Liquid.__init__`,
  `Alchemy.__init__`; delete `self._lmp` and every `lmp=self._lmp` kwarg.
- `Phase.lammps_close(lmp)` → body is `lmp.close()`; keep the method (call sites are many) or
  inline it — keep it, one line.
- Log renames: replace every `os.rename(log.lammps → <stage>.log.lammps)` block (averaging,
  integration, RS forward/backward, tscale, pscale, prescan, melting, and the error paths in
  `check_if_melted`/`check_if_solidfied`/`melt_structure`/pressure-convergence failures) with
  `lmp.rotate_logs("<same stage name>")`. Error paths call it with the same error-log names used
  today (`melted_error`, `solidified_error`, `pressure_convergence_error`, ...).

### 5.7 Input schema (`input.py`)

- `script_mode` field: keep the key but add a validator: if truthy → raise
  `ValueError("script_mode was removed in calphy vX. calphy now always drives the LAMMPS "
  "executable directly; set lammps_executable (or PATH/CALPHY_LAMMPS_EXECUTABLE) and remove "
  "script_mode from the input file.")`. If falsy → accept silently (old files keep working).
- `lammps_executable`, `mpi_executable`: update field docs (now used by every run; optional with
  resolution order).
- Delete `savefile` property, `save_job`, `load_job`.

### 5.8 CLI and packaging

- `clitools.py`: delete `run_averaging`, `process_averaging`, `run_integration`,
  `process_integration` (keep `convert_legacy_inputfile`, `phase_diagram`).
- `pyproject.toml`: remove the four entry points; remove `"pylammpsmpi"` from dependencies.
- `environment.yml`: keep the `lammps` conda package (it ships the `lmp` binary — add a comment
  saying that is why it is there), drop `pylammpsmpi`.
- Run the golden tests now (see 5.9), then continue.

### 5.9 Golden re-verification

- Update the Part 2 conftest: `recorded_job` now monkeypatches `create_object` to return an
  `ExecutableRunner(dry_run=True, ...)` wrapper whose `logical_commands` is compared instead of
  `RecordingRunner.commands`. The `sync()` calls added in 5.2 are interface no-ops for the
  comparison (logical_commands excludes runner-injected replay/restart lines by construction).
- Every golden must match **byte-for-byte**. A mismatch means the refactor changed the physics
  command stream — fix the code, never the golden. (Exception: if a golden legitimately must change
  — e.g. it recorded a script-mode-only artifact — stop and ask the user.)

### 5.10 Final sweeps

- `grep -rn "pylammpsmpi" calphy/ pyproject.toml` → zero hits (docs handled in Part 7).
- `grep -rn "script_mode" calphy/` → only the input-validator error message.
- `grep -rn "LammpsScript\|save_job\|load_job\|_lmp" calphy/` → zero hits.
- `pip install -e .` in a fresh venv **without** pylammpsmpi/lammps-python: `import calphy`,
  `calphy --version`, `calphy_kernel --help` all work.

### Acceptance criteria (Part 5)

- All Part 2 goldens byte-identical. All unit tests green. Sweeps in 5.10 clean.

---

## Part 6 — Integration tests with real LAMMPS

New file `tests/test_integration_lammps.py`, everything marked `@pytest.mark.lammps` and skipped
with a clear reason when no binary resolves. Use the Cu EAM fixtures. Configure marker in
`pyproject.toml` (`[tool.pytest.ini_options] markers`).

### 6.1 Continuity tests (the physics-safety core — write these first)

- **T1 exact-restart energy continuity (the strongest single test):** build a session
  (Cu fcc, eam/alloy, nvt fix), `run 100`, force a segment boundary (`sync()`), then `run 0` and a
  `print`/thermo read of PE via `fix print` or log parsing. Assert PE in the last thermo line of
  segment k equals PE in the first thermo line of segment k+1 to ≤ 1e-6 relative (restart +
  pair replay is exact; any drift = replay bug).
- **T2 Nose-Hoover state restoration:** NPT at 500 K/0 bar, 2000 steps in one segment vs
  1000 + sync + 1000; parse both segment logs; assert no volume/pressure discontinuity at the
  boundary beyond normal per-step fluctuation (compare boundary jump against the run's step-to-step
  std; factor < 3).
- **T3 ave/time redirection:** drive `run_pressure_convergence` on a deliberately bad initial
  lattice constant so it takes ≥ 2 cycles; assert `avg.dat` and `avg.dat.seg*` exist and
  `read_timeseries` returns strictly increasing step column across the concatenation.
- **T4 group/compute/variable replay:** define the MSD computes (via `compute_msd` helper), sync,
  run, assert `msd.dat(.seg*)` columns are finite and non-zero.

### 6.2 End-to-end matrix

Re-run the Part 0 scenarios B1–B12 through the new code (same inputs, `calphy_kernel`-equivalent
call in-process). For each:
- run completes, `report.yaml` exists;
- free energy is finite, not NaN;
- `abs(fe_new - fe_baseline) < tol` where `tol = max(3 * (err_new + err_baseline), 2e-3 eV/atom)`
  — document the achieved deltas in the test output; if any scenario exceeds tol, stop and ask the
  user (do not widen tol yourself);
- spring constants (solid fe) within 15% of baseline;
- `thermodynamic_integration()` + `integrate_reversible_scaling()` etc. run without error.
Mark B6–B12 `@pytest.mark.slow` in addition; CI runs B1/B4/B5 on every PR, the full matrix nightly
or on demand.

### 6.3 Error paths

- fe solid at 5000 K → `MeltedError` raised AND `melted_error.log.lammps` exists.
- fe liquid at 300 K with melting_cycle → `SolidifiedError` (or the melt loop error) raised.
- broken `pair_coeff` path → `LammpsExecutionError`, message contains the LAMMPS `ERROR` line.
- nonexistent binary → resolution `ValueError` naming the three lookup steps.
- unconverged pressure (1 cycle allowed, absurd tolerance) → the existing `ValueError` with
  `pressure_convergence_error.log.lammps` present.

### 6.4 MPI variant

- Re-run T1 and B1 with `queue.cores: 2` (skip when `mpirun` unavailable). Assert identical
  pass criteria.

### Acceptance criteria (Part 6)

- All `lammps`-marked tests green locally with a conda-forge `lmp` binary, serial and `-np 2`.
- End-to-end deltas vs baselines recorded in the PR description.

---

## Part 7 — CI, docs, cleanup, release

### 7.1 CI

- Inspect `.github/workflows/` and extend/add:
  - **unit job** (no LAMMPS): install package, run pytest without the `lammps` marker; also assert
    `pip install .` succeeds without pylammpsmpi.
  - **integration job**: micromamba env with conda-forge `lammps` (binary) + `openmpi`; run
    `-m lammps` excluding `slow`.
  - **nightly/dispatch job**: full matrix including `slow`.

### 7.2 Documentation

- `docs/source/inputfile.md`: rewrite the `script_mode` entry as a removal notice; rewrite
  `lammps_executable`/`mpi_executable` (resolution order, env vars, now applies to all modes);
  remove the staged-CLI references; document `md.cmdargs` now reaching the binary argv (GPU/kokkos
  flags work).
- Installation docs + README: pip install + "bring your own `lmp` binary" (conda-forge lammps,
  HPC module, container); remove pylammpsmpi instructions.
- Grep `docs/`, `examples/`, `notebooks/` for `script_mode`/`pylammpsmpi` and update every hit.
- CHANGELOG: breaking-changes section — script_mode removal, staged CLI removal, `lmp=` hook
  removal (note for pyiron/external embedders: pin previous version until migrated), new binary
  resolution, new env vars, `CALPHY_SKIP_PREFLIGHT`.

### 7.3 Final acceptance (whole project)

- Fresh venv, `pip install -e ".[test]"`, **no** pylammpsmpi/mpi4py/lammps-python present:
  unit tests green.
- With a conda-forge `lmp` on PATH: integration tests green, serial + MPI.
- `grep -rn "pylammpsmpi" .` (excluding CHANGELOG/baselines/git history) → zero hits.
- One end-to-end manual smoke: `calphy -i tests/baselines/inputs/B1.yaml` from a clean folder
  completes and the report matches the baseline within tolerance.
- Version bump (major or clearly-flagged minor — breaking release).

---

## Appendix A — Command classification table (authoritative)

Complete vocabulary emitted by calphy (verified by grep over `calphy/*.py`). Unknown first token →
`RunnerStateError`.

| first token | class | replay rule |
|---|---|---|
| `units`, `atom_style` | init | replayed before `read_restart` in every segment |
| `boundary`, `box` | init | segment 0 only (stored in restart; `box` illegal after box exists) |
| `timestep` | init/sticky | replay latest value after `read_restart` |
| `pair_style` | sticky | starts a new pair_block (clears previous); replayed after `read_restart` |
| `pair_coeff`, `mass` | sticky | appended to pair_block; replayed with it |
| `group` | sticky | replayed (idempotent type-based groups) |
| `compute` | sticky | tracked by id; removed by `uncompute` |
| `variable` | sticky | tracked by name, redefinition keeps position; `$(`-immediate crossing a boundary → error |
| `fix` | sticky | tracked by id, original order; removed by `unfix`; ` file X` rewritten to ` file X.seg<k>` on replay |
| `fix_modify` | sticky | attached to its fix id, replayed after it |
| `thermo`, `thermo_style`, `echo` | sticky | last value wins, replayed near end of header |
| `run` | one-shot | never replayed |
| `velocity`, `displace_atoms`, `change_box` | one-shot | never replayed |
| `read_data`, `read_restart` | one-shot | sets `has_box` |
| `write_data`, `write_restart`, `print` | one-shot | — |
| `dump` | sticky-tracked | live dump at a segment boundary → error (calphy always undumps first) |
| `undump`, `unfix`, `uncompute` | removal | delete from state |

## Appendix B — Sync-point map (only these five)

| location | reason |
|---|---|
| `Phase.dump_current_snapshot` (after `undump`) | dump file read by pyscal (melt/solidify checks, melt_structure, get_structures) |
| `Phase.run_iterative_pressure_convergence` loop | reads `avg.dat` each cycle |
| `Phase.run_iterative_constrained_pressure_convergence` loop | reads `avg.dat` via `process_pressure` |
| `Solid.run_iterative_spring_constant_convergence` loop | reads `msd.dat` each cycle |
| `Phase.scan_temperature_range` (after `unfix f2`) | reads `prescan.forward.dat` |

Everything else reads files only after `close()`.

## Appendix C — Deletion inventory

- `helpers.py`: `LammpsScript`, `pylammpsmpi` import, `lammps` import (Part 1), old `create_object` body.
- `phase.py`: `run_minimal_pressure_convergence`, `run_minimal_constrained_pressure_convergence`,
  all `script_mode` branches, `lmp=`/`self._lmp`, `os.rename(log.lammps ...)` blocks (→ `rotate_logs`).
- `solid.py`: `run_minimal_spring_constant_convergence`, `run_minimal_averaging`,
  `process_averaging_results`, dispatchers' branches.
- `liquid.py`, `alchemy.py`: all `script_mode` branches and script-write blocks.
- `kernel.py`: entire script-mode arm of `run_jobs`.
- `clitools.py`: `run_averaging`, `process_averaging`, `run_integration`, `process_integration`.
- `input.py`: `savefile`, `save_job`, `load_job`; `script_mode` becomes a tombstone validator.
- `pyproject.toml`: 4 entry points, `pylammpsmpi` dependency.

## Appendix D — Open questions (ask the user if hit)

- Exact `-h` style-listing format of the target LAMMPS versions (Part 4 parser).
- Any golden that appears to require modification after Part 2.
- Any end-to-end delta exceeding the Part 6 tolerances.
- Minimum supported LAMMPS version (affects preflight messaging; suggest documenting the version
  used in Part 0 as the reference).

## Addendum (2026-07-16) — Decision 5: optional pylammpsmpi library backend

Owner decision, superseding the "pylammpsmpi is fully removed" consequence of the
original scope (the four locked decisions above stand unchanged):

- calphy gets **two backends** behind one contract, `BaseRunner` in `runner.py`:
  `command`, `sync`, `close`, `rotate_logs`, `read_timeseries`, `logical_commands`.
  A runner owns *all* traffic across the LAMMPS boundary — commands in, data out.
  Driver code (phase/solid/liquid/alchemy) never touches files or a LAMMPS handle
  directly, so future backend-specific accessors (direct thermo/structure reads via
  the library) can be added to the contract with a file-based default and a
  library-mode override, without an `if backend:` ever appearing in driver code.
- `ExecutableRunner` stays the **default**. `LibraryRunner` (`library_runner.py`)
  drives a live `pylammpsmpi.LammpsLibrary` session: commands forward immediately,
  `sync()` is a no-op (no segments, no restart replay, no `.seg*` files),
  `rotate_logs` rotates a scratch `-log` target, tolerant of the drivers'
  close-before-rotate ordering.
- Selection: input key `execution_mode: executable | library` (default `executable`).
  pylammpsmpi is an **optional dependency** (`pip install calphy[library]`), imported
  lazily inside `LibraryRunner.__init__` only — `import calphy` never needs it.
  Preflight and binary/MPI resolution are executable-mode only.
- Unified-backend guarantee: goldens are asserted against **both** backends
  (`tests/test_library_runner.py` runs full driver paths over a fake pylammpsmpi and
  compares to the same immutable `tests/golden/*.txt`). Goldens/baselines remain
  immutable contracts.

Implemented as parts 8 (BaseRunner extraction + runner-mediated reads),
9 (LibraryRunner + unit tests), 10 (wire-in: input key, create_object branch,
pyproject extra, golden-equivalence tests, docs).
