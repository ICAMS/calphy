"""Inputs for T3 (switching rate series) and T4 (cross-check against composition
scaling) for a binary alloy whose potential is not in the repository.

    python make_alloy_inputs.py --elements Cu Ni --masses 63.546 58.693 \
        --potential /path/CuNi.eam.alloy --lattice-file /path/cuni_random.data \
        --temperature 800 --lmp /path/to/fork/lmp --cores 16

--lattice-file is the random solid solution (lammps-data, types 1 and 2 in the
order of --elements).  T4 additionally needs the pure crystal of element 1 of
the same size (--pure-lattice-file); composition scaling then takes it to the
composition of --lattice-file.  Element counts are read from the data files.
"""
import argparse
import os
import re


def counts(data_file):
    txt = open(data_file).read()
    atoms = re.search(r"Atoms.*?\n\n(.*)", txt, re.S).group(1).split("\n")
    types = [int(l.split()[1]) for l in atoms if l.strip() and not l.startswith("Velocities")]
    return {t: types.count(t) for t in sorted(set(types))}


def common(args, lattice, mode):
    return f"""  element: [{args.elements[0]}, {args.elements[1]}]
  mass: [{args.masses[0]}, {args.masses[1]}]
  lattice: {os.path.abspath(lattice)}
  lattice_constant: 0.0
  mode: {mode}
  reference_phase: solid
  temperature: {args.temperature}
  pressure: 0.0
  pair_style: {args.pair_style}
  pair_coeff: '* * {os.path.abspath(args.potential)} {args.elements[0]} {args.elements[1]}'
  lammps_executable: {args.lmp}
  n_print_steps: 1000
  md:
    timestep: 0.001
    n_small_steps: 2000
    n_cycles: 200
  n_equilibration_steps: 25000
  n_iterations: {args.iterations}
  queue:
    scheduler: local
    cores: {args.cores}
"""


MC = """  monte_carlo:
    n_steps: 10
    n_swaps: 10
    use_custom_lammps: true
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--elements", nargs=2, required=True)
    ap.add_argument("--masses", nargs=2, type=float, required=True)
    ap.add_argument("--potential", required=True)
    ap.add_argument("--pair-style", default="eam/alloy")
    ap.add_argument("--lattice-file", required=True, help="random solid solution, lammps-data")
    ap.add_argument("--pure-lattice-file", help="pure crystal of element 1, same cell (T4)")
    ap.add_argument("--temperature", type=float, default=800.0)
    ap.add_argument("--lmp", default="/path/to/thermoatoms/lammps/build/lmp")
    ap.add_argument("--cores", type=int, default=16)
    ap.add_argument("--iterations", type=int, default=3)
    ap.add_argument("--switching-steps", nargs="+", type=int, default=[25000, 50000, 100000, 200000])
    ap.add_argument("--out", default="inputs")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # T3: Frenkel-Ladd with swaps, one input per switching length
    for n in args.switching_steps:
        with open(os.path.join(args.out, "t3_fl_swaps_%d.yaml" % n), "w") as fh:
            fh.write("# T3: Frenkel-Ladd with swaps, n_switching_steps %d\ncalculations:\n-\n" % n)
            fh.write(common(args, args.lattice_file, "fe") + MC + "  n_switching_steps: %d\n" % n)
    with open(os.path.join(args.out, "t3_fl_noswaps.yaml"), "w") as fh:
        fh.write("# T3 reference: the same arrangement without swaps (one frozen configuration)\ncalculations:\n-\n")
        fh.write(common(args, args.lattice_file, "fe") + "  n_switching_steps: %d\n" % max(args.switching_steps))

    # T4: pure crystal by fe, and composition scaling with swaps to the alloy composition
    if args.pure_lattice_file:
        c = counts(args.lattice_file)
        with open(os.path.join(args.out, "t4_pure_fe.yaml"), "w") as fh:
            fh.write("# T4: pure %s crystal, plain fe\ncalculations:\n-\n" % args.elements[0])
            fh.write(common(args, args.pure_lattice_file, "fe")
                     .replace("element: [%s, %s]" % tuple(args.elements), "element: [%s]" % args.elements[0])
                     .replace("mass: [%s, %s]" % tuple(args.masses), "mass: [%s]" % args.masses[0])
                     .replace("%s %s'" % tuple(args.elements), "%s'" % args.elements[0])
                     + "  n_switching_steps: %d\n" % max(args.switching_steps))
        with open(os.path.join(args.out, "t4_composition_scaling.yaml"), "w") as fh:
            fh.write("# T4: composition scaling with swaps from pure %s to the T3 composition\ncalculations:\n-\n" % args.elements[0])
            fh.write(common(args, args.pure_lattice_file, "composition_scaling") + MC)
            fh.write("  composition_scaling:\n    output_chemical_composition:\n      %s: %d\n      %s: %d\n"
                     % (args.elements[0], c.get(1, 0), args.elements[1], c.get(2, 0)))
            fh.write("  n_switching_steps: %d\n" % max(args.switching_steps))
    print("inputs written to", args.out)


if __name__ == "__main__":
    main()
