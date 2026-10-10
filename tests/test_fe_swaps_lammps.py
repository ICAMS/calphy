"""Solid fe integration with Monte Carlo swaps against the real LAMMPS binary.

Needs the thermoatoms LAMMPS fork: its ``fix atom/swap`` takes an equal-style
variable as swap temperature, which the Frenkel-Ladd path uses for the
Metropolis test at T / (1 - lambda).  Skipped when no ``lmp`` resolves or the
resolved binary is stock LAMMPS.

The system is a degenerate binary: fcc Cu with half of the atoms labelled
type 2, a two-element copy of the Cu EAM in which both elements are Cu, and
equal masses.  Every arrangement of the labels has the same energy at every
lambda, so every swap is accepted and the free energy with swaps must be the
one without swaps plus k_B T ln(1/2), with the same switching work.
"""
import os
import subprocess
import tempfile

import numpy as np
import pytest
import yaml

from calphy.input import Calculation
from calphy.solid import Solid
from calphy.routines import routine_fe
from calphy.runner import resolve_lammps_executable
from calphy.integrators import kb

HERE = os.path.dirname(os.path.abspath(__file__))
POT = os.path.join(HERE, "Cu01.eam.alloy")

try:
    BINARY = resolve_lammps_executable(None)
except ValueError:
    BINARY = None


def _fork_has_variable_swap_temperature(binary):
    """True if ``fix atom/swap`` accepts ``v_name`` as its temperature."""
    if binary is None:
        return False
    with tempfile.TemporaryDirectory() as d:
        script = os.path.join(d, "probe.in")
        with open(script, "w") as fh:
            fh.write(
                "units metal\nlattice fcc 3.615\nregion box block 0 2 0 2 0 2\n"
                "create_box 2 box\ncreate_atoms 1 box\n"
                "set type 1 type/fraction 2 0.5 12345\nmass * 63.546\n"
                "pair_style eam/alloy\npair_coeff * * %s Cu Cu\n"
                "variable t equal 500.0\n"
                "fix s all atom/swap 1 1 1 v_t ke no types 1 2\nrun 0\n" % POT
            )
        proc = subprocess.run(
            [binary, "-in", script, "-log", "none", "-screen", "none"],
            cwd=d, capture_output=True, timeout=120,
        )
        return proc.returncode == 0


pytestmark = [
    pytest.mark.lammps,
    pytest.mark.skipif(BINARY is None, reason="no lmp binary resolves"),
    pytest.mark.skipif(
        not _fork_has_variable_swap_temperature(BINARY),
        reason="lmp is not the thermoatoms fork (no variable swap temperature)",
    ),
]


def write_two_type_cu(folder, n=4, seed=7):
    """fcc Cu n x n x n with half of the atoms labelled type 2, and a two-element
    copy of Cu01.eam.alloy in which both elements are Cu."""
    a = 3.615
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    pos = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                pos += list((basis + [i, j, k]) * a)
    pos = np.array(pos)
    rng = np.random.default_rng(seed)
    types = np.ones(len(pos), int)
    types[rng.choice(len(pos), len(pos) // 2, replace=False)] = 2
    data = os.path.join(folder, "cu2.data")
    with open(data, "w") as fh:
        fh.write("fcc Cu, half of the atoms labelled type 2\n\n%d atoms\n2 atom types\n\n" % len(pos))
        for ax in "xyz":
            fh.write("0.0 %.6f %slo %shi\n" % (n * a, ax, ax))
        fh.write("\nMasses\n\n1 63.546\n2 63.546\n\nAtoms # atomic\n\n")
        for i, (p, t) in enumerate(zip(pos, types)):
            fh.write("%d %d %.6f %.6f %.6f\n" % (i + 1, t, p[0], p[1], p[2]))

    lines = open(POT).read().splitlines()
    nrho, _, nr, _, _ = lines[4].split()
    nrho, nr = int(nrho), int(nr)
    vals = np.array(" ".join(lines[6:]).split(), float)
    F, rho, phi = vals[:nrho], vals[nrho:nrho + nr], vals[nrho + nr:nrho + 2 * nr]
    block = lambda arr: "\n".join("%.16e" % x for x in arr) + "\n"
    pot = os.path.join(folder, "cu2.eam.alloy")
    with open(pot, "w") as fh:
        fh.write("\n".join(lines[:3]) + "\n2 Cu Ni\n" + lines[4] + "\n")
        for _ in range(2):
            fh.write(lines[5] + "\n" + block(F) + block(rho))
        for _ in range(3):
            fh.write(block(phi))
    return data, pot, len(pos)


def make_calc(folder, data, pot, n_swaps):
    mc = {"n_steps": 10, "n_swaps": n_swaps, "use_custom_lammps": True} if n_swaps else {"n_swaps": 0}
    return Calculation(
        element=["Cu", "Ni"], mass=[63.546, 63.546], lattice=data, lattice_constant=0.0,
        mode="fe", reference_phase="solid", temperature=500.0, pressure=0.0,
        pair_style="eam/alloy", pair_coeff="* * %s Cu Ni" % pot,
        lammps_executable=BINARY, monte_carlo=mc,
        md={"timestep": 0.001, "n_small_steps": 500, "n_every_steps": 10,
            "n_repeat_steps": 10, "n_cycles": 20},
        # fixed springs: the reference is then the same in both runs (the
        # measured k would differ between runs and shift W and the reference
        # against each other)
        spring_constants=[2.0, 2.0], tolerance={"pressure": 50.0},
        n_equilibration_steps=1000, n_switching_steps=3000, n_iterations=1,
        queue={"scheduler": "local", "cores": 1},
    )


def run(folder, data, pot, n_swaps):
    calc = make_calc(folder, data, pot, n_swaps)
    simfolder = os.path.join(folder, "swaps%d" % n_swaps)
    os.makedirs(simfolder)
    job = routine_fe(Solid(calculation=calc, simfolder=simfolder))
    return job, simfolder


def test_degenerate_binary_swaps(tmp_path):
    folder = str(tmp_path)
    data, pot, natoms = write_two_type_cu(folder)
    job_s, folder_s = run(folder, data, pot, n_swaps=5)
    job_0, _ = run(folder, data, pot, n_swaps=0)

    res_s, res_0 = job_s.report["results"], job_0.report["results"]

    # bookkeeping: same reference, ideal mixing term only with swaps
    # the volumes converge separately; the reference sees that only through
    # the 1/N ln V term of the CM correction
    assert res_s["reference_system"] == pytest.approx(res_0["reference_system"], abs=2e-4)
    assert res_s["ideal_mixing"] == pytest.approx(kb * 500.0 * np.log(0.5))
    assert res_0["ideal_mixing"] == 0.0
    assert res_s["free_energy"] == pytest.approx(
        res_s["reference_system"] + res_s["ideal_mixing"] + res_s["work"] + res_s["pv"]
    )

    # the labels are degenerate: every swap is accepted, in the configurational
    # equilibration and along both legs
    ceq = np.loadtxt(os.path.join(folder_s, "swap.equilibration.dat"))
    assert np.all(ceq[:, 3] == ceq[:, 2]) and ceq[-1, 2] > 0
    for leg in ("forward", "backward"):
        stats = np.loadtxt(os.path.join(folder_s, "swap.fe.%s_1.dat" % leg))
        assert np.all(stats[:, 3] == stats[:, 2]) and stats[-1, 2] > 0
        # lambda runs 0 -> 1 on the forward leg, 1 -> 0 on the backward one
        lam = stats[:, 1]
        assert (lam[-1] > lam[0]) == (leg == "forward")

    # the switching work is a configurational quantity and does not see the
    # labels: with swaps it is the no-swap work up to the switching noise of
    # these short runs (a few meV/atom)
    assert abs(res_s["work"] - res_0["work"]) < 5e-3
    assert abs(
        (res_s["free_energy"] - res_0["free_energy"]) - res_s["ideal_mixing"]
    ) < 5e-3

    # every atom carried the common mass in the switching run
    conf = open(os.path.join(folder_s, "conf.fe.backward_1.data")).read()
    assert "63.546" in conf
