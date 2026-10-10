"""Degenerate binary for T1: fcc Cu with half of the atoms labelled type 2 and
a two-element copy of tests/Cu01.eam.alloy in which both elements are Cu.

    python make_degenerate_cu.py --repeat 6
    python make_degenerate_cu.py --relabel-only some_fcc.data --fraction 0.3

Every arrangement of the labels has the same energy, so the free energy with
swaps must be the one without swaps plus k_B T sum x ln x.
"""
import argparse
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
POT = os.path.join(HERE, "..", "..", "Cu01.eam.alloy")


def write_fcc(path, n, a, fraction, seed):
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    pos = np.array([(basis + [i, j, k]) * a for i in range(n) for j in range(n) for k in range(n)])
    pos = pos.reshape(-1, 3)
    types = relabel(np.ones(len(pos), int), fraction, seed)
    write_data(path, pos, types, [n * a] * 3)
    return len(pos)


def relabel(types, fraction, seed):
    rng = np.random.default_rng(seed)
    types = np.ones_like(types)
    types[rng.choice(len(types), int(round(fraction * len(types))), replace=False)] = 2
    return types


def write_data(path, pos, types, box):
    with open(path, "w") as fh:
        fh.write("fcc, fraction of atoms labelled type 2\n\n%d atoms\n2 atom types\n\n" % len(pos))
        for ax, l in zip("xyz", box):
            fh.write("0.0 %.6f %slo %shi\n" % (l, ax, ax))
        fh.write("\nMasses\n\n1 63.546\n2 63.546\n\nAtoms # atomic\n\n")
        for i, (p, t) in enumerate(zip(pos, types)):
            fh.write("%d %d %.6f %.6f %.6f\n" % (i + 1, t, p[0], p[1], p[2]))


def write_two_element_potential(src, dst, names=("Cu", "Ni")):
    lines = open(src).read().splitlines()
    assert int(lines[3].split()[0]) == 1, "expects a single element setfl file"
    nrho, _, nr, _, _ = lines[4].split()
    nrho, nr = int(nrho), int(nr)
    vals = np.array(" ".join(lines[6:]).split(), float)
    F, rho, phi = vals[:nrho], vals[nrho:nrho + nr], vals[nrho + nr:nrho + 2 * nr]
    block = lambda arr: "\n".join("%.16e" % x for x in arr) + "\n"
    with open(dst, "w") as fh:
        fh.write("\n".join(lines[:3]) + "\n2 %s %s\n" % names + lines[4] + "\n")
        for _ in range(2):
            fh.write(lines[5] + "\n" + block(F) + block(rho))
        for _ in range(3):
            fh.write(block(phi))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeat", type=int, default=6)
    ap.add_argument("--lattice-constant", type=float, default=3.615)
    ap.add_argument("--fraction", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--relabel-only", metavar="DATA", help="relabel a fraction of the atoms of an existing data file (ASE)")
    ap.add_argument("--out", default=HERE)
    args = ap.parse_args()
    if args.relabel_only:
        from ase.io import read
        at = read(args.relabel_only, format="lammps-data", style="atomic")
        types = relabel(np.ones(len(at), int), args.fraction, args.seed)
        write_data(os.path.join(args.out, "relabelled.data"), at.positions, types, at.cell.lengths())
        print("relabelled.data:", len(at), "atoms,", (types == 2).sum(), "type 2")
        return
    n = write_fcc(os.path.join(args.out, "cu2.data"), args.repeat, args.lattice_constant, args.fraction, args.seed)
    write_two_element_potential(POT, os.path.join(args.out, "cu2.eam.alloy"))
    print("cu2.data: %d atoms; cu2.eam.alloy written" % n)


if __name__ == "__main__":
    main()
