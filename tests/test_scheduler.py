import os
import re
import sys
import tempfile
import importlib.util
import pytest

import importlib.util, sys

# Import scheduler directly to avoid the pyscal3 dependency in calphy/__init__.py
_spec = importlib.util.spec_from_file_location(
    "calphy.scheduler",
    os.path.join(os.path.dirname(__file__), "..", "calphy", "scheduler.py"),
)
_mod = importlib.util.module_from_spec(_spec)
sys.modules["calphy.scheduler"] = _mod
_spec.loader.exec_module(_mod)
Local = _mod.Local
SLURM = _mod.SLURM
SGE = _mod.SGE


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _read_script(path):
    with open(path) as f:
        return f.read()


def _sbatch_lines(script):
    """Return the list of #SBATCH directive lines (without the '#SBATCH ' prefix)."""
    return [
        line[len("#SBATCH ") :].strip()
        for line in script.splitlines()
        if line.startswith("#SBATCH ")
    ]


# ---------------------------------------------------------------------------
# SLURM – defaults
# ---------------------------------------------------------------------------


class TestSLURMDefaults:
    def test_no_partition_by_default(self):
        s = SLURM({})
        with tempfile.TemporaryDirectory() as d:
            s.queueoptions["directory"] = d
            script_path = os.path.join(d, "job.sub")
            s.write_script(script_path)
            directives = _sbatch_lines(_read_script(script_path))
        assert not any(line.startswith("--partition") for line in directives)

    def test_no_memory_by_default(self):
        """memory defaults to None → no --mem line."""
        s = SLURM({})
        with tempfile.TemporaryDirectory() as d:
            s.queueoptions["directory"] = d
            script_path = os.path.join(d, "job.sub")
            s.write_script(script_path)
            directives = _sbatch_lines(_read_script(script_path))
        assert not any("--mem" in line for line in directives)

    def test_hint_nomultithread_by_default(self):
        s = SLURM({})
        with tempfile.TemporaryDirectory() as d:
            s.queueoptions["directory"] = d
            script_path = os.path.join(d, "job.sub")
            s.write_script(script_path)
            directives = _sbatch_lines(_read_script(script_path))
        assert "--hint=nomultithread" in directives

    def test_no_cpus_per_task_by_default(self):
        s = SLURM({})
        with tempfile.TemporaryDirectory() as d:
            s.queueoptions["directory"] = d
            script_path = os.path.join(d, "job.sub")
            s.write_script(script_path)
            directives = _sbatch_lines(_read_script(script_path))
        assert not any("--cpus-per-task" in line for line in directives)


# ---------------------------------------------------------------------------
# SLURM – target script (matches the user's Raven/GPU setup)
# ---------------------------------------------------------------------------


class TestSLURMTargetScript:
    """
    Reproduce:
        #SBATCH --ntasks=1
        #SBATCH --mem=125000
        #SBATCH --hint=nomultithread
        #SBATCH --cpus-per-task=18
        #SBATCH --constraint="gpu"
        #SBATCH --gres=gpu:a100:1
        source activate calphy-custom-lammps
    (no --partition line)
    """

    def _make_scheduler(self, tmpdir):
        options = {
            "scheduler": "slurm",
            "jobname": "test_job",
            "cores": 1,
            "cpus_per_task": 18,
            "memory": "125000",
            "hint": "nomultithread",
            "options": {"constraint": '"gpu"', "gres": "gpu:a100:1"},
            "commands": ["source activate calphy-custom-lammps"],
            "directory": tmpdir,
        }
        return SLURM(options, cores=1, directory=tmpdir)

    def test_ntasks(self):
        with tempfile.TemporaryDirectory() as d:
            s = self._make_scheduler(d)
            s.write_script(os.path.join(d, "job.sub"))
            directives = _sbatch_lines(_read_script(os.path.join(d, "job.sub")))
        assert "--ntasks=1" in directives

    def test_cpus_per_task(self):
        with tempfile.TemporaryDirectory() as d:
            s = self._make_scheduler(d)
            s.write_script(os.path.join(d, "job.sub"))
            directives = _sbatch_lines(_read_script(os.path.join(d, "job.sub")))
        assert "--cpus-per-task=18" in directives

    def test_mem_total(self):
        """--mem (total) should appear, not --mem-per-cpu."""
        with tempfile.TemporaryDirectory() as d:
            s = self._make_scheduler(d)
            s.write_script(os.path.join(d, "job.sub"))
            script = _read_script(os.path.join(d, "job.sub"))
        assert "--mem=125000" in script
        assert "--mem-per-cpu" not in script

    def test_no_partition(self):
        with tempfile.TemporaryDirectory() as d:
            s = self._make_scheduler(d)
            s.write_script(os.path.join(d, "job.sub"))
            directives = _sbatch_lines(_read_script(os.path.join(d, "job.sub")))
        assert not any("--partition" in line for line in directives)

    def test_constraint_gres(self):
        with tempfile.TemporaryDirectory() as d:
            s = self._make_scheduler(d)
            s.write_script(os.path.join(d, "job.sub"))
            directives = _sbatch_lines(_read_script(os.path.join(d, "job.sub")))
        assert '--constraint="gpu"' in directives
        assert "--gres=gpu:a100:1" in directives

    def test_activate_command(self):
        with tempfile.TemporaryDirectory() as d:
            s = self._make_scheduler(d)
            s.write_script(os.path.join(d, "job.sub"))
            script = _read_script(os.path.join(d, "job.sub"))
        assert "source activate calphy-custom-lammps" in script

    def test_hint_present(self):
        with tempfile.TemporaryDirectory() as d:
            s = self._make_scheduler(d)
            s.write_script(os.path.join(d, "job.sub"))
            directives = _sbatch_lines(_read_script(os.path.join(d, "job.sub")))
        assert "--hint=nomultithread" in directives


# ---------------------------------------------------------------------------
# SLURM – suppressing defaults via None
# ---------------------------------------------------------------------------


class TestSLURMNullOverrides:
    def test_hint_suppressed_with_none(self):
        """Setting hint=None in options dict should remove the --hint line."""
        s = SLURM({"hint": None})
        with tempfile.TemporaryDirectory() as d:
            s.queueoptions["directory"] = d
            s.write_script(os.path.join(d, "job.sub"))
            directives = _sbatch_lines(_read_script(os.path.join(d, "job.sub")))
        assert not any("--hint" in line for line in directives)

    def test_memory_set_explicitly(self):
        s = SLURM({"memory": "64000"})
        with tempfile.TemporaryDirectory() as d:
            s.queueoptions["directory"] = d
            s.write_script(os.path.join(d, "job.sub"))
            directives = _sbatch_lines(_read_script(os.path.join(d, "job.sub")))
        assert "--mem=64000" in directives

    def test_partition_set(self):
        s = SLURM({"queuename": "gpu"})
        with tempfile.TemporaryDirectory() as d:
            s.queueoptions["directory"] = d
            s.write_script(os.path.join(d, "job.sub"))
            directives = _sbatch_lines(_read_script(os.path.join(d, "job.sub")))
        assert "--partition=gpu" in directives


# ---------------------------------------------------------------------------
# SLURM – script structure
# ---------------------------------------------------------------------------


class TestSLURMScriptStructure:
    def test_shebang(self):
        s = SLURM({})
        with tempfile.TemporaryDirectory() as d:
            s.queueoptions["directory"] = d
            s.write_script(os.path.join(d, "job.sub"))
            script = _read_script(os.path.join(d, "job.sub"))
        assert script.startswith("#!/bin/bash")

    def test_chdir_present(self):
        with tempfile.TemporaryDirectory() as d:
            s = SLURM({}, directory=d)
            s.write_script(os.path.join(d, "job.sub"))
            directives = _sbatch_lines(_read_script(os.path.join(d, "job.sub")))
        assert any("--chdir" in line for line in directives)

    def test_walltime(self):
        s = SLURM({"walltime": "12:00:00"})
        with tempfile.TemporaryDirectory() as d:
            s.queueoptions["directory"] = d
            s.write_script(os.path.join(d, "job.sub"))
            directives = _sbatch_lines(_read_script(os.path.join(d, "job.sub")))
        assert "--time=12:00:00" in directives


# ---------------------------------------------------------------------------
# Local scheduler – smoke test
# ---------------------------------------------------------------------------


class TestLocal:
    def test_write_creates_script(self):
        with tempfile.TemporaryDirectory() as d:
            s = Local({}, directory=d)
            s.queueoptions["commands"] = ["echo hello"]
            script_path = os.path.join(d, "local.sub")
            s.write_script(script_path)
            assert os.path.exists(script_path)

    def test_submit_sets_executable(self):
        """submit() must chmod the script before executing it."""
        with tempfile.TemporaryDirectory() as d:
            s = Local({}, directory=d)
            s.maincommand = "true"  # harmless command
            script_path = os.path.join(d, "local.sub")
            s.write_script(script_path)
            s.submit()
            assert os.access(script_path, os.X_OK)

    def test_shebang(self):
        with tempfile.TemporaryDirectory() as d:
            s = Local({}, directory=d)
            script_path = os.path.join(d, "local.sub")
            s.write_script(script_path)
            script = _read_script(script_path)
        assert script.startswith("#!/bin/bash")


# ---------------------------------------------------------------------------
# SGE – smoke test
# ---------------------------------------------------------------------------


class TestSGE:
    def test_write_script(self):
        with tempfile.TemporaryDirectory() as d:
            s = SGE({"queuename": "all.q"}, directory=d)
            script_path = os.path.join(d, "sge.sub")
            s.write_script(script_path)
            script = _read_script(script_path)
        assert "#$ -N" in script
        assert "#$ -l qname=all.q" in script
