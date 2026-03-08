"""
Tests for $USER / ${USER} environment-variable expansion in
Calculation.fix_paths (calphy.input).
"""

import os
import getpass
import tempfile
import types

import pytest

from calphy.input import Calculation


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _call_fix_paths(potlist):
    """
    Call Calculation.fix_paths on *potlist* without constructing a full
    Calculation object (fix_paths does not use any instance attributes).
    """
    dummy = types.SimpleNamespace()
    return Calculation.fix_paths(dummy, potlist)


# ---------------------------------------------------------------------------
# Tests for $USER expansion
# ---------------------------------------------------------------------------


class TestUserEnvExpansion:

    def test_dollar_user_expanded_when_file_exists(self):
        """$USER in the potential path should be replaced by the real username
        when the resolved file actually exists on disk."""
        username = getpass.getuser()

        with tempfile.NamedTemporaryFile(suffix=".eam.alloy", delete=False) as f:
            real_path = f.name

        try:
            # Build the $USER variant only if the real path contains the username
            if f"/{username}/" not in real_path:
                pytest.skip(
                    "Temp-file path does not contain the username directory "
                    "component; cannot construct a $USER test path."
                )

            user_path = real_path.replace(f"/{username}/", "/$USER/")
            pair_coeff = f"* * {user_path} Cu"

            result = _call_fix_paths([pair_coeff])

            assert len(result) == 1
            assert "$USER" not in result[0], "Literal $USER should have been expanded"
            assert username in result[0], "Expanded username should appear in result"
            assert (
                real_path in result[0]
            ), "Absolute path to real file should be present"
        finally:
            os.unlink(real_path)

    def test_dollar_user_braces_expanded_when_file_exists(self):
        """${USER} (brace form) should also be expanded."""
        username = getpass.getuser()

        with tempfile.NamedTemporaryFile(suffix=".eam.alloy", delete=False) as f:
            real_path = f.name

        try:
            if f"/{username}/" not in real_path:
                pytest.skip("Temp-file path does not contain the username.")

            user_path = real_path.replace(f"/{username}/", "/${USER}/")
            pair_coeff = f"* * {user_path} Cu"

            result = _call_fix_paths([pair_coeff])

            assert len(result) == 1
            assert "${USER}" not in result[0]
            assert username in result[0]
            assert real_path in result[0]
        finally:
            os.unlink(real_path)

    def test_dollar_user_expanded_even_when_file_missing(self):
        """$USER should be expanded in the output even when the file does not
        exist – calphy must not pass a literal '$USER' string to LAMMPS."""
        username = getpass.getuser()

        pair_coeff = "* * /home/$USER/potentials/Cu.eam Cu"
        result = _call_fix_paths([pair_coeff])

        assert len(result) == 1
        assert "$USER" not in result[0], "Literal $USER must not be passed to LAMMPS"
        assert username in result[0], f"Username '{username}' should appear in result"

    def test_dollar_user_braces_expanded_even_when_file_missing(self):
        """${USER} should be expanded even when the file does not exist."""
        username = getpass.getuser()

        pair_coeff = "* * /home/${USER}/potentials/Cu.eam Cu"
        result = _call_fix_paths([pair_coeff])

        assert len(result) == 1
        assert "${USER}" not in result[0]
        assert username in result[0]

    def test_element_tokens_preserved_after_user_expansion(self):
        """Element tokens (after the potential file) must be preserved
        intact when $USER is expanded."""
        username = getpass.getuser()

        pair_coeff = "* * /home/$USER/potentials/CuZr.eam.fs Cu Zr"
        result = _call_fix_paths([pair_coeff])

        assert len(result) == 1
        tokens = result[0].split()
        # Tokens 0,1 are '* *'; token 2 is the filename; tokens 3+ are elements
        assert tokens[0] == "*"
        assert tokens[1] == "*"
        assert tokens[3] == "Cu"
        assert tokens[4] == "Zr"

    def test_multiple_pair_coeffs_all_expanded(self):
        """$USER should be expanded in every entry when potlist has multiple
        pair_coeff strings."""
        username = getpass.getuser()

        potlist = [
            "* * /home/$USER/pots/A.eam A",
            "* * /home/$USER/pots/B.eam B",
        ]
        result = _call_fix_paths(potlist)

        assert len(result) == 2
        for entry in result:
            assert "$USER" not in entry
            assert username in entry


# ---------------------------------------------------------------------------
# Tests that pre-existing (non-$USER) behaviour is unchanged
# ---------------------------------------------------------------------------


class TestExistingBehaviourPreserved:

    def test_absolute_path_used_when_file_exists(self):
        """A relative path that resolves to an existing file is converted to
        its absolute form (existing behaviour)."""
        with tempfile.NamedTemporaryFile(
            suffix=".eam.alloy", dir=".", delete=False
        ) as f:
            real_path = os.path.abspath(f.name)
            rel_name = os.path.basename(f.name)

        try:
            pair_coeff = f"* * {rel_name} Cu"
            result = _call_fix_paths([pair_coeff])

            assert len(result) == 1
            assert real_path in result[0]
        finally:
            os.unlink(real_path)

    def test_nonexistent_path_without_vars_unchanged(self):
        """A path with no env-var tokens that does not exist should be left
        exactly as supplied (existing behaviour)."""
        pair_coeff = "* * /nonexistent/path/Cu.eam Cu"
        result = _call_fix_paths([pair_coeff])
        assert result == [pair_coeff]

    def test_short_pair_coeff_passed_through_unchanged(self):
        """pair_coeff strings with fewer than 3 tokens must not be touched."""
        pair_coeff = "* *"
        result = _call_fix_paths([pair_coeff])
        assert result == [pair_coeff]

    def test_tilde_home_expanded(self):
        """~ should also be expanded (bonus: expanduser is applied alongside
        expandvars)."""
        home = os.path.expanduser("~")
        pair_coeff = "* * ~/potentials/Cu.eam Cu"
        result = _call_fix_paths([pair_coeff])

        assert len(result) == 1
        assert "~" not in result[0], "Tilde should have been expanded"
        assert home in result[0]
