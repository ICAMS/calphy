"""Evaluation of the cluster validation runs, see README.md.

    python analyse.py summary <folder>...
    python analyse.py degenerate <swaps folder> <noswaps folder>
    python analyse.py acceptance <folder> [--bins 20]
    python analyse.py sro <dump> --cutoff 3.0 [--types 1 2]
    python analyse.py rate <folder>...
    python analyse.py compare <FL+swaps folder> <pure fe folder> <composition scaling folder>

A folder is a calphy run folder (contains report.yaml) or a directory that
contains exactly one such folder.
"""
import argparse
import glob
import os
import sys

import numpy as np
import yaml

KB = 8.617333262e-5  # eV/K


def run_folder(path):
    if os.path.exists(os.path.join(path, "report.yaml")):
        return path
    hits = [d for d in glob.glob(os.path.join(path, "*")) if os.path.exists(os.path.join(d, "report.yaml"))]
    if len(hits) != 1:
        sys.exit("%s: expected one run folder with report.yaml, found %d" % (path, len(hits)))
    return hits[0]


def report(path):
    with open(os.path.join(run_folder(path), "report.yaml")) as fh:
        return yaml.safe_load(fh)


def inputs(path):
    files = glob.glob(os.path.join(run_folder(path), "input*.yaml")) + glob.glob(os.path.join(run_folder(path), "*.yaml"))
    for f in files:
        if os.path.basename(f) == "report.yaml":
            continue
        try:
            data = yaml.safe_load(open(f))
        except Exception:
            continue
        if isinstance(data, dict) and "calculations" in data:
            return data["calculations"][0]
    return {}


def summary(args):
    print("%-40s %12s %9s %12s %12s %12s %10s" % ("run", "F", "err", "W", "dissipation", "ideal_mix", "mode"))
    for p in args.folders:
        r = report(p)["results"]
        print("%-40s %12.6f %9.1e %12.6f %12.2e %12.6f %10s" % (
            os.path.basename(run_folder(p).rstrip("/"))[:40], r["free_energy"], r["error"],
            r["work"], r["dissipation"], r.get("ideal_mixing", 0.0),
            "swaps" if r.get("ideal_mixing", 0.0) != 0 else "no swaps"))


def degenerate(args):
    rs, r0 = report(args.swaps)["results"], report(args.noswaps)["results"]
    t = report(args.swaps)["input"]["temperature"]
    conc = [float(x) for x in str(report(args.swaps)["input"]["concentration"]).split()]
    expected = KB * t * sum(x * np.log(x) for x in conc if x > 0)
    dF = rs["free_energy"] - r0["free_energy"]
    err = np.hypot(rs["error"], r0["error"])
    print("F(swaps) - F(noswaps) = %.6f eV/atom" % dF)
    print("k_B T sum x ln x      = %.6f eV/atom (ideal_mixing in report: %.6f)" % (expected, rs["ideal_mixing"]))
    print("residual %.2e eV/atom, combined error %.1e, %.1f sigma" % (dF - expected, err, abs(dF - expected) / err if err > 0 else float("nan")))
    for leg in ("equilibration", "fe.forward_1", "fe.backward_1"):
        f = os.path.join(run_folder(args.swaps), "swap.%s.dat" % leg)
        if os.path.exists(f):
            d = np.loadtxt(f)
            print("%-16s attempted %8.0f accepted %8.0f" % (leg, d[-1, -2], d[-1, -1]))
    ok = abs(dF - expected) < max(5e-4, 2 * err)
    print("PASS" if ok else "FAIL")


def acceptance(args):
    folder = run_folder(args.folder)
    for leg in ("forward", "backward"):
        files = sorted(glob.glob(os.path.join(folder, "swap.fe.%s_*.dat" % leg)))
        if not files:
            continue
        lam_all, acc_all = [], []
        for f in files:
            d = np.loadtxt(f)
            lam = 0.5 * (d[1:, 1] + d[:-1, 1])
            att = np.diff(d[:, 2]); acc = np.diff(d[:, 3])
            lam_all.append(lam); acc_all.append((att, acc))
        lam = np.concatenate(lam_all)
        att = np.concatenate([a for a, _ in acc_all]); acc = np.concatenate([b for _, b in acc_all])
        edges = np.linspace(0, 1, args.bins + 1)
        print("%s leg, %d iterations: acceptance per lambda bin" % (leg, len(files)))
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (lam >= lo) & (lam < hi) if hi < 1 else (lam >= lo) & (lam <= hi)
            a, s = att[m].sum(), acc[m].sum()
            print("  lambda %.2f-%.2f  attempted %8.0f  accepted %8.0f  rate %.4f" % (lo, hi, a, s, s / a if a > 0 else float("nan")))


def sro(args):
    from ase.io import read
    from ase.neighborlist import neighbor_list
    frames = read(args.dump, format="lammps-dump-text", index=":")
    ta, tb = args.types
    print("frame  N_A  N_B  x_B  P_B(A)  alpha_AB   (first shell, cutoff %.2f)" % args.cutoff)
    for k, at in enumerate(frames):
        types = at.get_array("type") if at.has("type") else at.numbers
        i, j = neighbor_list("ij", at, args.cutoff)
        ma = types[i] == ta
        n_b_around_a = np.sum(types[j][ma] == tb)
        p_b_a = n_b_around_a / max(ma.sum(), 1)
        x_b = np.mean(types == tb)
        alpha = 1 - p_b_a / x_b if x_b > 0 else float("nan")
        print("%5d %4d %4d %.3f %.4f %9.4f" % (k, (types == ta).sum(), (types == tb).sum(), x_b, p_b_a, alpha))


def rate(args):
    rows = []
    for p in args.folders:
        r = report(p)
        # switching length from the forward data of iteration 1 (one row per step plus the header)
        fwd = os.path.join(run_folder(p), "forward_1.dat")
        n_switch = sum(1 for l in open(fwd) if not l.startswith("#")) - 1 if os.path.exists(fwd) else 0
        rows.append((n_switch, r["results"]["free_energy"], r["results"]["error"], r["results"]["dissipation"], run_folder(p)))
    rows.sort()
    print("%10s %12s %9s %12s  %s" % ("n_switch", "F", "err", "dissipation", "folder"))
    for n, f, e, q, folder in rows:
        print("%10d %12.6f %9.1e %12.2e  %s" % (n, f, e, q, os.path.basename(folder.rstrip("/"))))
    if len(rows) >= 2:
        d = rows[-1][1] - rows[-2][1]
        e = np.hypot(rows[-1][2], rows[-2][2])
        print("two longest runs differ by %.2e eV/atom (%.1f sigma): %s" % (d, abs(d) / e if e > 0 else float("nan"),
              "PASS" if abs(d) < max(1e-3, 2 * e) else "FAIL"))


def compare(args):
    fl, pure, cs = report(args.fl)["results"], report(args.pure)["results"], report(args.cs)["results"]
    f_route2 = pure["free_energy"] + cs["free_energy"]
    err = np.sqrt(fl["error"] ** 2 + pure["error"] ** 2 + cs["error"] ** 2)
    print("Frenkel-Ladd with swaps:           F = %.6f +- %.1e" % (fl["free_energy"], fl["error"]))
    print("pure fe + composition scaling:     F = %.6f +- %.1e  (pure %.6f, dF_cs %.6f)" % (f_route2, np.hypot(pure["error"], cs["error"]), pure["free_energy"], cs["free_energy"]))
    d = fl["free_energy"] - f_route2
    print("difference %.2e eV/atom, %.1f sigma: %s" % (d, abs(d) / err if err > 0 else float("nan"), "PASS" if abs(d) < 2 * max(err, 5e-4) else "FAIL"))
    print("note: the composition scaling result must include its ideal mixing term the same way (see postprocessing.fix_composition_scaling)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("summary"); s.add_argument("folders", nargs="+"); s.set_defaults(func=summary)
    s = sub.add_parser("degenerate"); s.add_argument("swaps"); s.add_argument("noswaps"); s.set_defaults(func=degenerate)
    s = sub.add_parser("acceptance"); s.add_argument("folder"); s.add_argument("--bins", type=int, default=20); s.set_defaults(func=acceptance)
    s = sub.add_parser("sro"); s.add_argument("dump"); s.add_argument("--cutoff", type=float, required=True); s.add_argument("--types", nargs=2, type=int, default=[1, 2]); s.set_defaults(func=sro)
    s = sub.add_parser("rate"); s.add_argument("folders", nargs="+"); s.set_defaults(func=rate)
    s = sub.add_parser("compare"); s.add_argument("fl"); s.add_argument("pure"); s.add_argument("cs"); s.set_defaults(func=compare)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
