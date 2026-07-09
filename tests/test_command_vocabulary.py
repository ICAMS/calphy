"""Command-vocabulary contract (Part 2).

Every first token emitted by any calphy driver path must belong to the closed
vocabulary in EXECUTABLE_RUNNER_PLAN.md Appendix A.  The ExecutableRunner (Part 3)
raises RunnerStateError on any command whose first token is not classified, so an
unknown token here would be an unhandled command after the refactor.

The list below is a literal copy of the Appendix A table.  Once Part 3 lands and
`calphy.runner` exposes the classification table, this test should import it and
assert the two agree (reconciliation noted in the plan).
"""
import glob
import os

import pytest

# --- Appendix A: authoritative command vocabulary (first tokens) ------------- #
INIT = {"units", "atom_style", "boundary", "box", "timestep"}
STICKY = {
    "pair_style", "pair_coeff", "mass", "group", "compute", "variable",
    "fix", "fix_modify", "thermo", "thermo_style", "echo", "dump",
}
ONE_SHOT = {
    "run", "velocity", "displace_atoms", "change_box",
    "read_data", "read_restart", "write_data", "write_restart", "print",
}
REMOVAL = {"undump", "unfix", "uncompute"}

VOCABULARY = INIT | STICKY | ONE_SHOT | REMOVAL

GOLDEN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden")
GOLDEN_FILES = sorted(glob.glob(os.path.join(GOLDEN_DIR, "*.txt")))


def test_goldens_exist():
    assert GOLDEN_FILES, "no golden files found -- run the golden-stream tests first"


@pytest.mark.parametrize("path", GOLDEN_FILES, ids=lambda p: os.path.basename(p))
def test_golden_sanity(path):
    """Each golden is non-empty and carries no Python repr leakage."""
    with open(path) as fh:
        text = fh.read()
    assert text.strip(), f"{os.path.basename(path)} is empty"
    assert "object at 0x" not in text, f"{os.path.basename(path)} contains a repr leak"


@pytest.mark.parametrize("path", GOLDEN_FILES, ids=lambda p: os.path.basename(p))
def test_golden_vocabulary(path):
    """Every command's first token is in the Appendix A vocabulary."""
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            token = line.split()[0]
            assert token in VOCABULARY, (
                f"{os.path.basename(path)}:{lineno}: unknown first token "
                f"{token!r} (not in Appendix A vocabulary)\n  line: {line}"
            )
