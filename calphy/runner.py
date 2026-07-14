"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^2, Ralf Drautz^2
^1: Max Planck Institut fuer Eisenforschung, Dusseldorf, Germany
^2: Ruhr-University Bochum, Bochum, Germany

calphy is published and distributed under the Academic Software License v1.0 (ASL).
See the LICENSE FILE for more details.
"""
# ruff: noqa: E501
"""Executable-driven LAMMPS runner.

calphy talks to LAMMPS only through ``command()`` strings and reads everything
back from files.  ``ExecutableRunner`` preserves that contract without a live
library: it buffers commands, and at each ``sync()`` flushes the buffer as one
segment script that it runs with the ``lmp`` binary.  Continuity between segments
is provided by ``write_restart``/``read_restart`` for box/atoms/velocities/
timestep plus a replay of the sticky session state (potential, groups, computes,
variables, fixes, thermo) tracked in :class:`SessionState`.

This module is pure-Python and imports no LAMMPS library -- it drives the binary
as a subprocess.
"""
import os
import glob
import shutil
import logging
import subprocess

import numpy as np

from calphy.errors import LammpsExecutionError, RunnerStateError

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Command classification (Appendix A of EXECUTABLE_RUNNER_PLAN.md)
# --------------------------------------------------------------------------- #
INIT_TOKENS = frozenset({"units", "atom_style", "boundary", "box", "timestep"})
STICKY_TOKENS = frozenset({
    "pair_style", "pair_coeff", "mass", "group", "compute", "variable",
    "fix", "fix_modify", "thermo", "thermo_style", "echo", "dump",
})
ONE_SHOT_TOKENS = frozenset({
    "run", "velocity", "displace_atoms", "change_box",
    "read_data", "read_restart", "write_data", "write_restart", "print",
})
REMOVAL_TOKENS = frozenset({"undump", "unfix", "uncompute"})

#: Every first token calphy is allowed to emit.  Anything else -> RunnerStateError.
KNOWN_TOKENS = INIT_TOKENS | STICKY_TOKENS | ONE_SHOT_TOKENS | REMOVAL_TOKENS


# --------------------------------------------------------------------------- #
# Binary resolution
# --------------------------------------------------------------------------- #
def _resolve_executable(explicit, env_var, default_name, label):
    """Resolve an executable: explicit arg -> env var -> ``default_name`` on PATH."""
    candidates = [
        (explicit, "the input key"),
        (os.environ.get(env_var), "$" + env_var),
        (default_name, "PATH"),
    ]
    for value, _ in candidates:
        if not value:
            continue
        resolved = shutil.which(value)
        if resolved is None and os.path.isabs(value) and os.access(value, os.X_OK):
            resolved = value
        if resolved is not None:
            return resolved
    raise ValueError(
        "Could not resolve the {label} executable. Tried, in order:\n"
        "  1. the input key ({explicit!r})\n"
        "  2. the environment variable {env} ({envval!r})\n"
        "  3. {default!r} on PATH\n"
        "Set one of these, e.g. `lammps_executable: /path/to/{default}` in the "
        "input file, `export {env}=/path/to/{default}`, or put {default} on PATH."
        .format(
            label=label, explicit=explicit, env=env_var,
            envval=os.environ.get(env_var), default=default_name,
        )
    )


def resolve_lammps_executable(explicit=None):
    """input key -> $CALPHY_LAMMPS_EXECUTABLE -> ``lmp`` on PATH."""
    return _resolve_executable(explicit, "CALPHY_LAMMPS_EXECUTABLE", "lmp", "LAMMPS")


def resolve_mpi_executable(explicit=None):
    """input key -> $CALPHY_MPI_EXECUTABLE -> ``mpirun`` on PATH (cores > 1 only)."""
    return _resolve_executable(explicit, "CALPHY_MPI_EXECUTABLE", "mpirun", "MPI")


# --------------------------------------------------------------------------- #
# Cumulative time-series reads across segments
# --------------------------------------------------------------------------- #
def _ordered_timeseries_files(directory, name):
    """``name`` (if present) followed by every ``name.seg<k>`` in ascending k."""
    base = os.path.join(directory, name)
    parts = [base] if os.path.exists(base) else []
    segs = glob.glob(base + ".seg*")

    def _segidx(path):
        try:
            return int(path.rsplit(".seg", 1)[1])
        except (IndexError, ValueError):
            return -1

    parts.extend(sorted(segs, key=_segidx))
    return parts


def read_timeseries(directory, name, usecols=None):
    """Load a session-produced time series that may be split across segments.

    Loads ``<directory>/<name>`` plus every ``<name>.seg<k>`` (ascending k) and
    concatenates them along the sample axis.  For a single file this reproduces
    exactly ``np.loadtxt(name, usecols=usecols)`` (comment headers stripped by
    numpy).  ``usecols`` may be a tuple of column indices as in ``np.loadtxt``.
    """
    parts = _ordered_timeseries_files(directory, name)
    if not parts:
        raise FileNotFoundError(
            "no time-series files for %r in %s" % (name, directory)
        )
    if len(parts) == 1:
        return np.loadtxt(parts[0], usecols=usecols)

    # Load each file as 2D (ndmin=2 normalizes single-row / single-column files
    # that np.loadtxt would otherwise squeeze), then vstack.
    blocks = [np.loadtxt(p, usecols=usecols, ndmin=2) for p in parts]
    blocks = [b for b in blocks if b.size]
    if not blocks:
        return np.loadtxt(parts[0], usecols=usecols)
    result = np.vstack(blocks)
    # Match np.loadtxt's shape: a single selected column comes back 1D.
    if result.shape[1] == 1:
        return result[:, 0]
    return result


# --------------------------------------------------------------------------- #
# Session state tracker
# --------------------------------------------------------------------------- #
def _rewrite_replayed_fix(cmd, segidx):
    """Rewrite a replayed fix command's ``file`` target to ``file.seg<segidx>``.

    A replayed ``fix print`` (or any ``append`` fix) would corrupt cumulative
    output, so those raise instead -- calphy defines/unfixes ``fix print`` within
    a single segment and never lets it cross a boundary.
    """
    tokens = cmd.split()
    style = tokens[3] if len(tokens) >= 4 else ""
    if style == "print":
        raise RunnerStateError(
            "fix print (id %r) cannot be replayed across a segment boundary"
            % (tokens[1] if len(tokens) > 1 else "?")
        )
    if "append" in tokens:
        raise RunnerStateError(
            "fix ... append (id %r) cannot cross a segment boundary"
            % (tokens[1] if len(tokens) > 1 else "?")
        )
    if "file" in tokens:
        i = tokens.index("file")
        if i + 1 < len(tokens):
            tokens[i + 1] = tokens[i + 1] + ".seg%d" % segidx
    return " ".join(tokens)


class SessionState:
    """Tracks the sticky LAMMPS state that must be replayed after each restart.

    Fed every command via :meth:`apply` (in emission order).  Only the first
    token (plus id/name tokens) is parsed.  Classification follows Appendix A.
    """

    def __init__(self):
        self._init = {}                 # token -> command (units/atom_style/boundary/box)
        self.timestep = None            # last value wins
        self.pair_block = []            # contiguous pair_style/pair_coeff/mass
        self.groups = {}                # name -> command (insertion-ordered)
        self.computes = {}              # id -> command
        self.variables = {}             # name -> command (position kept on redefine)
        self.fixes = {}                 # id -> {"command": str, "fix_modify": [str]}
        self.dumps = {}                 # id -> command
        self.thermo = None
        self.thermo_style = None
        self.echo = None
        self.has_box = False

    # -- ingestion ---------------------------------------------------------- #
    def apply(self, cmd):
        """Mutate the tracked state for one (already whitespace-normalized) command."""
        tokens = cmd.split()
        if not tokens:
            return
        t = tokens[0]
        if t in ("units", "atom_style", "boundary", "box"):
            self._init[t] = cmd
        elif t == "timestep":
            self.timestep = cmd
        elif t == "pair_style":
            self.pair_block = [cmd]         # a new style clears the old block
        elif t in ("pair_coeff", "mass"):
            self.pair_block.append(cmd)
        elif t == "group":
            self.groups[tokens[1]] = cmd
        elif t == "compute":
            self.computes[tokens[1]] = cmd
        elif t == "uncompute":
            self.computes.pop(tokens[1], None)
        elif t == "variable":
            self.variables[tokens[1]] = cmd   # dict keeps original position on overwrite
        elif t == "fix":
            self.fixes[tokens[1]] = {"command": cmd, "fix_modify": []}
        elif t == "fix_modify":
            fid = tokens[1]
            if fid not in self.fixes:
                raise RunnerStateError("fix_modify references unknown fix id %r" % fid)
            self.fixes[fid]["fix_modify"].append(cmd)
        elif t == "unfix":
            self.fixes.pop(tokens[1], None)
        elif t == "dump":
            self.dumps[tokens[1]] = cmd
        elif t == "undump":
            self.dumps.pop(tokens[1], None)
        elif t == "thermo":
            self.thermo = cmd
        elif t == "thermo_style":
            self.thermo_style = cmd
        elif t == "echo":
            self.echo = cmd
        elif t in ("read_data", "read_restart"):
            self.has_box = True
        # remaining one-shot tokens carry no sticky state.

    # -- replay ------------------------------------------------------------- #
    def check_replayable(self):
        """Raise if the current state cannot be faithfully carried across a boundary."""
        if self.dumps:
            raise RunnerStateError(
                "dump(s) %s are live at a segment boundary; calphy must undump "
                "before a read" % sorted(self.dumps)
            )
        immediate = [n for n, c in self.variables.items() if "$(" in c]
        if immediate:
            raise RunnerStateError(
                "immediate-evaluation variable(s) %s cannot cross a segment "
                "boundary" % immediate
            )

    def build_replay_header(self, segidx, restart_name):
        """Return the replay header lines for segment ``segidx`` (> 0)."""
        lines = []
        if "units" in self._init:
            lines.append(self._init["units"])
        if "atom_style" in self._init:
            lines.append(self._init["atom_style"])
        if self.timestep is not None:
            lines.append(self.timestep)
        lines.append("read_restart %s" % restart_name)
        lines.extend(self.pair_block)
        lines.extend(self.groups.values())
        lines.extend(self.computes.values())
        lines.extend(self.variables.values())
        for entry in self.fixes.values():
            lines.append(_rewrite_replayed_fix(entry["command"], segidx))
            lines.extend(entry["fix_modify"])
        if self.thermo_style is not None:
            lines.append(self.thermo_style)
        if self.thermo is not None:
            lines.append(self.thermo)
        if self.echo is not None:
            lines.append(self.echo)
        return lines


# --------------------------------------------------------------------------- #
# ExecutableRunner
# --------------------------------------------------------------------------- #
def _normalize_cmdargs(cmdargs):
    if cmdargs == "" or cmdargs is None:
        return []
    if isinstance(cmdargs, str):
        return cmdargs.split()
    return list(cmdargs)


class ExecutableRunner:
    """Drives an external LAMMPS binary in restart-continued segments.

    Interface used by calphy: :meth:`command`, :meth:`sync`, :meth:`close`,
    :meth:`rotate_logs`, and the :attr:`logical_commands` property.
    """

    def __init__(self, *, binary, mpi_command, cores, cmdargs, directory,
                 restart_name="calphy.restart", dry_run=False, timeout=None):
        self.binary = binary
        self.mpi_command = mpi_command
        self.cores = cores
        self.cmdargs = _normalize_cmdargs(cmdargs)
        self.directory = directory
        self.restart_name = restart_name
        self.dry_run = dry_run
        self.timeout = timeout

        self.state = SessionState()
        self._buffer = []
        self._logical = []
        self._segidx = -1
        # (segidx, seg_script_basename, seg_log_path) since the last rotate_logs
        self._pending_logs = []

    # -- public interface --------------------------------------------------- #
    @property
    def logical_commands(self):
        """Every command ever passed to :meth:`command`, whitespace-normalized."""
        return list(self._logical)

    def command(self, s):
        """Normalize whitespace, classify the first token, and buffer the command."""
        cmd = " ".join(s.split())
        if not cmd:
            return
        token = cmd.split()[0]
        if token not in KNOWN_TOKENS:
            raise RunnerStateError(
                "unknown command %r (first token %r not in the calphy vocabulary)"
                % (cmd, token)
            )
        self._logical.append(cmd)
        self._buffer.append(cmd)

    def sync(self):
        """Flush the buffer as one segment and execute it (unless dry_run)."""
        if not self._buffer:
            return
        self._segidx += 1
        k = self._segidx
        if k == 0:
            lines = list(self._buffer)
        else:
            self.state.check_replayable()
            lines = self.state.build_replay_header(k, self.restart_name)
            lines.extend(self._buffer)
        lines.append("write_restart %s" % self.restart_name)

        # Commit this segment's commands to the state for the next header.
        for cmd in self._buffer:
            self.state.apply(cmd)
        self._buffer = []

        self._emit_and_run(k, lines)

    def close(self):
        """Final flush (if pending), then drop the rolling restart file."""
        if self._buffer:
            self.sync()
        restart = os.path.join(self.directory, self.restart_name)
        if os.path.exists(restart):
            try:
                os.remove(restart)
            except OSError:
                pass

    def rotate_logs(self, stage_name):
        """Concatenate the segment logs produced since the last rotation.

        Writes ``<stage_name>.log.lammps`` in the sim folder with a separator per
        segment, then resets the pending-log list.  Segment scripts are kept.
        """
        out = os.path.join(self.directory, "%s.log.lammps" % stage_name)
        with open(out, "w") as fout:
            for segidx, seg_script, seg_log in self._pending_logs:
                fout.write("# --- segment %d: %s ---\n" % (segidx, seg_script))
                if os.path.exists(seg_log):
                    with open(seg_log) as fin:
                        fout.write(fin.read())
                fout.write("\n")
        self._pending_logs = []

    # -- internals ---------------------------------------------------------- #
    def _emit_and_run(self, k, lines):
        seg_script = "calphy.seg%d.lmp" % k
        seg_log = "calphy.seg%d.log" % k
        seg_path = os.path.join(self.directory, seg_script)
        log_path = os.path.join(self.directory, seg_log)
        with open(seg_path, "w") as fout:
            fout.write("\n".join(lines) + "\n")
        self._pending_logs.append((k, seg_script, log_path))
        if self.dry_run:
            return
        self._run(seg_path, seg_script, log_path)

    def _run(self, seg_path, seg_script, log_path):
        argv = []
        if self.cores > 1:
            argv += [self.mpi_command, "-np", str(self.cores)]
        argv += [self.binary, "-in", seg_script, "-log",
                 os.path.basename(log_path), "-screen", "none"]
        argv += self.cmdargs
        try:
            proc = subprocess.run(
                argv, cwd=self.directory, capture_output=True, text=True,
                timeout=self.timeout,
            )
        except subprocess.TimeoutExpired as exc:
            raise LammpsExecutionError(
                self._failure_message(
                    seg_path, log_path,
                    "timed out after %s s" % self.timeout,
                    exc.stderr,
                )
            )
        if proc.returncode != 0:
            raise LammpsExecutionError(
                self._failure_message(
                    seg_path, log_path,
                    "exited with return code %d" % proc.returncode,
                    proc.stderr,
                )
            )
        if self._log_reports_error(log_path):
            raise LammpsExecutionError(
                self._failure_message(
                    seg_path, log_path, "LAMMPS reported an ERROR", proc.stderr
                )
            )

    @staticmethod
    def _log_reports_error(log_path):
        if not os.path.exists(log_path):
            return False
        with open(log_path) as fin:
            tail = fin.readlines()[-50:]
        return any("ERROR" in line for line in tail)

    @staticmethod
    def _failure_message(seg_path, log_path, reason, stderr):
        excerpt = ""
        if os.path.exists(log_path):
            with open(log_path) as fin:
                excerpt = "".join(fin.readlines()[-15:])
            source = "log tail (%s)" % log_path
        else:
            excerpt = (stderr or "")[-2000:]
            source = "stderr (log %s missing)" % log_path
        return (
            "LAMMPS segment failed: %s\n"
            "  segment script: %s\n"
            "  segment log:    %s\n"
            "  --- %s ---\n%s"
            % (reason, seg_path, log_path, source, excerpt)
        )


# --------------------------------------------------------------------------- #
# Preflight capability check (Part 4)
# --------------------------------------------------------------------------- #
#: Section headings in `lmp -h` -> the category key we track.
_STYLE_CATEGORIES = ("pair", "fix", "compute")

#: Styles that come from an optional LAMMPS package, for the error message hint.
_STYLE_PACKAGE = {
    ("fix", "ti/spring"): "EXTRA-FIX",
    ("pair", "ufm"): "EXTRA-PAIR",
    ("fix", "atom/swap"): "MC",
    ("fix", "qtb"): "QTB",
    ("pair", "hybrid/scaled"): "EXTRA-PAIR",
}

_STYLE_CACHE = {}


def _parse_styles(text):
    """Parse the ``* Pair/Fix/Compute styles:`` sections of ``lmp -h`` output.

    Captures every whitespace-separated token between a category heading and the
    next ``*`` heading, so it is tolerant of column wrapping and blank lines.
    """
    result = {cat: set() for cat in _STYLE_CATEGORIES}
    current = None
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("*"):
            current = None
            # headings look like "* Pair styles:" or "* Fix styles" (some
            # sections omit the trailing colon), so match tolerantly.
            head = stripped[1:].strip().rstrip(":")
            if head.lower().endswith("styles"):
                word = head.split()[0].lower()
                if word in result:
                    current = word
            continue
        if current is not None and stripped:
            result[current].update(stripped.split())
    return result


def _query_styles(binary):
    try:
        proc = subprocess.run(
            [binary, "-h"], capture_output=True, text=True, timeout=60
        )
    except (OSError, subprocess.TimeoutExpired):
        return {cat: set() for cat in _STYLE_CATEGORIES}
    text = (proc.stdout or "") + "\n" + (proc.stderr or "")
    return _parse_styles(text)


def get_binary_styles(binary):
    """Return ``{"pair": {...}, "fix": {...}, "compute": {...}}`` for ``binary``.

    Parses ``lmp -h``.  Cached per (path, mtime).  If the binary cannot produce a
    style listing (some builds crash on ``-h``), the returned sets are empty --
    callers must treat an all-empty result as "capabilities unknown".
    """
    try:
        mtime = os.path.getmtime(binary)
    except OSError:
        mtime = None
    key = (binary, mtime)
    if key not in _STYLE_CACHE:
        _STYLE_CACHE[key] = _query_styles(binary)
    return _STYLE_CACHE[key]


def required_styles(calc):
    """The pair/fix/compute styles a calculation needs, derived from ``calc``."""
    pair = set(getattr(calc, "_pair_style_names", None) or [])
    fix = {"nve", "nvt", "npt", "langevin", "ave/time", "print", "momentum"}
    compute = {"msd", "temp/com", "pair"}

    mode = calc.mode
    phase = (calc.reference_phase or "").lower()

    if phase == "solid" and mode in ("fe", "ts", "tscale"):
        fix.add("ti/spring")
    if phase == "liquid" or mode == "melting_temperature":
        pair.add("ufm")
        pair.add("hybrid/scaled")
    if mode in ("ts", "tscale") or getattr(calc, "pair_mode", None) == "overlay":
        pair.add("hybrid/scaled")
    if getattr(getattr(calc, "monte_carlo", None), "n_swaps", 0) > 0:
        fix.add("atom/swap")
    if getattr(calc, "_qtb", False):
        fix.add("qtb")

    return {"pair": pair, "fix": fix, "compute": compute}


def preflight(calc, binary):
    """Verify ``binary`` provides every style ``calc`` needs, else raise ValueError.

    Skips (with a logged warning) when ``CALPHY_SKIP_PREFLIGHT=1`` or when the
    binary cannot report its styles.
    """
    if os.environ.get("CALPHY_SKIP_PREFLIGHT") == "1":
        logger.warning(
            "CALPHY_SKIP_PREFLIGHT=1: skipping the LAMMPS capability preflight."
        )
        return

    available = get_binary_styles(binary)
    if not any(available.values()):
        logger.warning(
            "Could not read a style list from `%s -h`; skipping the preflight "
            "capability check. If a required style is missing the run will fail "
            "later with a LAMMPS ERROR.", binary,
        )
        return

    required = required_styles(calc)
    missing = []
    for cat in _STYLE_CATEGORIES:
        for name in sorted(required.get(cat, set())):
            if name not in available.get(cat, set()):
                missing.append((cat, name))

    if missing:
        lines = []
        for cat, name in missing:
            pkg = _STYLE_PACKAGE.get((cat, name), "core / recent LAMMPS")
            lines.append("  - %s style %r  (LAMMPS package: %s)" % (cat, name, pkg))
        raise ValueError(
            "The LAMMPS binary %r is missing styles this calculation needs:\n%s\n"
            "Rebuild LAMMPS with the packages above (e.g. "
            "`-D PKG_EXTRA-FIX=yes`) or use a full binary. Set "
            "CALPHY_SKIP_PREFLIGHT=1 to bypass this check."
            % (binary, "\n".join(lines))
        )
