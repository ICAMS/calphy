"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021-2026 (c) Sarath Menon, Yury Lysogorskiy, Ralf Drautz
Interdisciplinary Centre for Advanced Materials Simulation (ICAMS),
Ruhr University Bochum, 44801 Bochum, Germany

calphy is published and distributed under the Academic Software Licence v1.0 (ASL).
calphy is distributed in the hope that it will be useful for non-commercial academic
research, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the LICENSE file for details.

The ASL permits academic non-commercial use only. Contact
sarath.menon@ruhr-uni-bochum.de to enquire about commercial use rights.

More information about the program can be found in:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz.
"Automated Free Energy Calculation from Atomistic Simulations." Physical Review Materials 5(10), 2021
DOI: 10.1103/PhysRevMaterials.5.103801

For more information contact:
sarath.menon@ruhr-uni-bochum.de
"""

import numpy as np
import yaml
import os
import itertools

from calphy.integrators import *
import calphy.helpers as ph
import calphy.phase as cph


class Alchemy(cph.Phase):
    """
    Class for alchemical transformations

    Parameters
    ----------
    options : dict
        dict of input options

    kernel : int
        the index of the calculation that should be run from
        the list of calculations in the input file

    simfolder : string
        base folder for running calculations

    lmp : LAMMPS object
        The LAMMPS object to use for the calculation
    """

    def __init__(self, calculation=None, simfolder=None, log_to_screen=False, lmp=None):

        # call base class
        super().__init__(
            calculation=calculation, simfolder=simfolder, log_to_screen=log_to_screen, lmp=lmp,
        )

    def _end_states(self):
        """
        The initial and final potential of the alchemical switch.

        The pair lists hold both end states back to back: the first half
        describes the initial potential, the second half the final one. Each
        half is one physical potential made of one component (a plain pair
        style) or several (``pair_mode: overlay``, or the component-wise
        rewrite done by composition scaling). A component is the tuple
        (style with options, style name, pair_coeff).
        """
        styles = self.calc._pair_style_with_options
        names = self.calc._pair_style_names
        coeffs = self.calc.pair_coeff
        if len(styles) % 2 or len(styles) != len(coeffs):
            raise ValueError(
                "alchemical switching needs pair_style and pair_coeff lists of "
                "equal, even length (initial components followed by final "
                "components); got %d pair styles and %d pair coeffs"
                % (len(styles), len(coeffs))
            )
        half = len(styles) // 2
        components = list(zip(styles, names, coeffs))
        return components[:half], components[half:]

    @staticmethod
    def _pure_pair_style_command(components):
        """``pair_style`` for one potential on its own: the plain style for a
        single component, ``hybrid/overlay`` of the components otherwise."""
        if len(components) == 1:
            return "pair_style       %s" % components[0][0]
        return "pair_style       hybrid/overlay %s" % " ".join(
            style for style, _, _ in components
        )

    @staticmethod
    def _pure_pair_coeff_commands(components):
        """``pair_coeff`` lines matching :meth:`_pure_pair_style_command`."""
        if len(components) == 1:
            return ["pair_coeff       %s" % components[0][2]]
        return ph.hybrid_pair_coeff_commands_for(
            [name for _, name, _ in components], [coeff for _, _, coeff in components]
        )

    def _scaled_pair_commands(self, initial, final, initial_scale, final_scale):
        """
        Install ``hybrid/scaled`` with every initial component scaled by
        ``initial_scale`` and every final component by ``final_scale``,
        then define computes ``c1``/``c2`` (or ``c1_1, c1_2, ...`` for
        multi-component ends) holding the two potential energies and the
        per-atom variables ``dU1``/``dU2`` recorded during switching.

        Returns the commands and the compute ids to release afterwards.
        """
        terms = ["%s %s" % (initial_scale, style) for style, _, _ in initial]
        terms += ["%s %s" % (final_scale, style) for style, _, _ in final]
        commands = ["pair_style       hybrid/scaled %s" % " ".join(terms)]

        components = initial + final
        names = [name for _, name, _ in components]
        commands += ph.hybrid_pair_coeff_commands_for(
            names, [coeff for _, _, coeff in components]
        )

        tags = ph.hybrid_component_tags(names)
        compute_ids, energies = [], []
        for label, part, part_tags in (
            ("c1", initial, tags[: len(initial)]),
            ("c2", final, tags[len(initial) :]),
        ):
            ids = (
                [label]
                if len(part) == 1
                else ["%s_%d" % (label, i) for i in range(1, len(part) + 1)]
            )
            compute_commands, energy = ph.hybrid_pair_compute_commands(ids, part_tags)
            commands += compute_commands
            compute_ids += ids
            energies.append(energy if len(ids) == 1 else "(%s)" % energy)

        commands.append("variable        step equal step")
        commands.append("variable        dU1 equal %s/atoms" % energies[0])
        commands.append("variable        dU2 equal %s/atoms" % energies[1])
        return commands, compute_ids

    def _coupling_pair(self, lmp, ramp="0.0", tag=None):
        """Coupling mode: base components (all but last) at constant scale
        1.0, the LAST component scaled by `ramp` (a constant or v_name).
        H(lam) = U_base + lam*U_last — an alchemical coupling ramp of one
        added component over an arbitrary multi-component base. With `tag`,
        defines compute cR<tag> + variable dUc<tag> for the ramped
        component's per-atom energy (unique ids per stage)."""
        ns = self.calc._pair_style_with_options
        names = self.calc._pair_style_names
        terms = " ".join(f"1.0 {s}" for s in ns[:-1]) + f" {ramp} {ns[-1]}"
        lmp.command("pair_style       hybrid/scaled %s" % terms)
        for i, pc in enumerate(self.calc.pair_coeff):
            words = pc.split()
            same = [j for j, n in enumerate(names) if n == names[i]]
            idx = [names[i]] + ([str(same.index(i) + 1)]
                                if len(same) > 1 else [])
            lmp.command("pair_coeff       "
                        + " ".join([*words[:2], *idx, *words[2:]]))
        if tag is not None:
            same = [j for j, n in enumerate(names) if n == names[-1]]
            occ = (" %d" % (same.index(len(names) - 1) + 1)) \
                if len(same) > 1 else ""
            lmp.command("compute         cR%s all pair %s%s"
                        % (tag, names[-1], occ))
            lmp.command("variable        dUc%s equal c_cR%s/atoms"
                        % (tag, tag))
        return lmp

    def _run_integration_coupling(self, iteration=1):
        """Coupling-mode integration: ramp ONLY the last component 0->1
        (forward) and 1->0 (backward) over the fixed base.
        W = int <U_last>_lam dlam;  dF = (W_f + W_b)/2 per iteration."""
        lmp = ph.create_object(self.calc, self.simfolder, lmp=self.lmp)
        conf = os.path.join(self.simfolder, "conf.equilibration.data")
        # style only before read_data (pair_coeff needs the box)
        ns = self.calc._pair_style_with_options
        terms = " ".join(f"1.0 {s}" for s in ns[:-1]) + f" 0.0 {ns[-1]}"
        lmp.command("pair_style       hybrid/scaled %s" % terms)
        lmp = ph.read_data(lmp, conf)
        self._coupling_pair(lmp, ramp="0.0")
        lmp = ph.set_mass(lmp, self.calc)
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)
        lmp.command(
            "velocity          all create %f %d mom yes rot yes dist gaussian"
            % (self.calc._temperature, np.random.randint(1, 10000)))
        if self.calc.npt:
            lmp.command(
                "fix             f1 all npt temp %f %f %f %s %f %f %f"
                % (self.calc._temperature, self.calc._temperature,
                   self.calc.md.thermostat_damping[1], self.iso,
                   self.calc._pressure, self.calc._pressure,
                   self.calc.md.barostat_damping[1]))
        else:
            lmp.command("fix             f1 all nvt temp %f %f %f"
                        % (self.calc._temperature, self.calc._temperature,
                           self.calc.md.thermostat_damping[1]))
        lmp.command("thermo_style    custom step pe")
        lmp.command("thermo          1000")
        lmp.command("run             %d" % self.calc.n_equilibration_steps)

        # forward: 0 -> 1
        lmp.command("variable         clambda equal ramp(0.0,1.0)")
        self._coupling_pair(lmp, ramp="v_clambda", tag="f")
        lmp.command(
            'fix             f2 all print 1 "${dUcf} ${dUcf} ${clambda}" '
            'title "# dU_ramped[eV/atom] dU_ramped[eV/atom] lambda" '
            "screen no file forward_%d.dat" % iteration)
        lmp.command("run             %d" % self.calc._n_switching_steps)
        lmp.command("unfix           f2")
        lmp.command("uncompute       cRf")
        lmp.command("variable        clambda delete")

        # equilibrate at lam = 1 (full base + component)
        self._coupling_pair(lmp, ramp="1.0")
        lmp.command("run             %d" % self.calc.n_equilibration_steps)

        # backward: 1 -> 0
        lmp.command("variable         clambda equal ramp(1.0,0.0)")
        self._coupling_pair(lmp, ramp="v_clambda", tag="b")
        lmp.command(
            'fix             f3 all print 1 "${dUcb} ${dUcb} ${clambda}" '
            'title "# dU_ramped[eV/atom] dU_ramped[eV/atom] lambda" '
            "screen no file backward_%d.dat" % iteration)
        lmp.command("run             %d" % self.calc._n_switching_steps)
        lmp.command("unfix           f3")

        self.lammps_close(lmp=lmp)
        lmp.rotate_logs("integration")

    def run_averaging(self):
        """
        Run averaging routine

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        Run averaging routine using LAMMPS. Starting from the initial lattice two different routines can
        be followed:
        If pressure is specified, MD simulations are run until the pressure converges within the given
        threshold value.
        Fix lattice option is not implemented at present.
        At the end of the run, the averaged box dimensions are calculated.
        """
        lmp = ph.create_object(self.calc, self.simfolder, lmp=self.lmp)

        if self.calc.alchemy_coupling:
            # equilibrate on the BASE system (lam = 0): full base, ramped
            # component at scale 0
            ns = self.calc._pair_style_with_options
            terms = (" ".join(f"1.0 {s}" for s in ns[:-1])
                     + f" 0.0 {ns[-1]}")
            lmp.command("pair_style       hybrid/scaled %s" % terms)
            lmp = ph.create_structure(lmp, self.calc)
            self._coupling_pair(lmp, ramp="0.0")
        else:
            # equilibrate with the initial potential on its own
            initial, _ = self._end_states()
            lmp.command(self._pure_pair_style_command(initial))

            # set up structure
            lmp = ph.create_structure(lmp, self.calc)

            # set up potential
            for command in self._pure_pair_coeff_commands(initial):
                lmp.command(command)
        lmp = ph.set_mass(lmp, self.calc)

        # add some computes
        lmp.command("variable         mvol equal vol")
        lmp.command("variable         mlx equal lx")
        lmp.command("variable         mly equal ly")
        lmp.command("variable         mlz equal lz")
        lmp.command("variable         mpress equal press")
        lmp.command("variable         mpe equal pe/atoms")
        lmp.command("variable         metotal equal etotal/atoms")
        lmp.command("variable         mtemp equal temp")

        # add some computes
        if not self.calc._fix_lattice:
            if self.calc._pressure == 0:
                self.run_zero_pressure_equilibration(lmp)
            else:
                self.run_finite_pressure_equilibration(lmp)

            # equilibration-frame dump (no-op unless
            # n_print_steps_equilibration > 0)
            self.start_equilibration_dump(lmp)

            # this is when the averaging routine starts
            self.run_pressure_convergence(lmp)

        # run if a constrained lattice is used
        else:
            self.start_equilibration_dump(lmp)
            # routine in which lattice constant will not varied, but is set to a given fixed value
            self.run_constrained_pressure_convergence(lmp)

        # check for melting (skip in coupling mode: the base ensemble may
        # legitimately be a liquid — the phase was validated when its
        # baseline free energy was measured)
        self.stop_equilibration_dump(lmp)
        self.dump_current_snapshot(lmp, "traj.equilibration_stage2.dat")
        if not self.calc.alchemy_coupling:
            self.check_if_melted(lmp, "traj.equilibration_stage2.dat")

        # close object and process traj
        lmp = ph.write_data(lmp, "conf.equilibration.data")

        self.lammps_close(lmp=lmp)
        lmp.rotate_logs("averaging")

    def run_integration(self, iteration=1):
        """
        Run integration routine

        Parameters
        ----------
        iteration : int, optional
            iteration number for running independent iterations

        Returns
        -------
        None

        Notes
        -----
        Run the integration routine where the initial and final systems are connected using
        the lambda parameter. See algorithm 4 in publication.
        """
        if self.calc.alchemy_coupling:
            return self._run_integration_coupling(iteration=iteration)

        initial, final = self._end_states()

        # create lammps object
        lmp = ph.create_object(self.calc, self.simfolder, lmp=self.lmp)

        # Adiabatic switching parameters.
        lmp.command("variable        li       equal   1.0")
        lmp.command("variable        lf       equal   0.0")

        # the equilibrated configuration is read in with the initial potential
        lmp.command(self._pure_pair_style_command(initial))
        conf = os.path.join(self.simfolder, "conf.equilibration.data")
        lmp = ph.read_data(lmp, conf)
        for command in self._pure_pair_coeff_commands(initial):
            lmp.command(command)
        lmp = ph.set_mass(lmp, self.calc)

        lmp.command("group g1 type 1")
        lmp.command("group g2 type 2")

        # remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        lmp.command(
            "velocity          all create %f %d mom yes rot yes dist gaussian"
            % (self.calc._temperature, np.random.randint(1, 10000))
        )
        # Integrator & thermostat.
        if self.calc.npt:
            lmp.command(
                "fix             f1 all npt temp %f %f %f %s %f %f %f"
                % (
                    self.calc._temperature,
                    self.calc._temperature,
                    self.calc.md.thermostat_damping[1],
                    self.iso,
                    self.calc._pressure,
                    self.calc._pressure,
                    self.calc.md.barostat_damping[1],
                )
            )
        else:
            lmp.command(
                "fix             f1 all nvt temp %f %f %f"
                % (
                    self.calc._temperature,
                    self.calc._temperature,
                    self.calc.md.thermostat_damping[1],
                )
            )

        lmp.command("thermo_style    custom step pe")
        lmp.command("thermo          1000")
        lmp.command("run             %d" % self.calc.n_equilibration_steps)

        # equilibration run is over

        # ---------------------------------------------------------------
        # FWD cycle: initial potential scaled 1 -> 0, final potential 0 -> 1
        # ---------------------------------------------------------------
        lmp.command("variable         flambda equal ramp(${li},${lf})")
        lmp.command("variable         blambda equal ramp(${lf},${li})")

        commands, compute_ids = self._scaled_pair_commands(
            initial, final, "v_flambda", "v_blambda"
        )
        for command in commands:
            lmp.command(command)

        swap_fixes = self._add_swap_fixes(
            lmp, self.calc.monte_carlo.forward_swap_types, "Forward"
        )

        # Thermo output.
        lmp.command("thermo_style    custom step v_dU1 v_dU2")
        lmp.command("thermo          1000")

        # save the necessary items to a file: first step
        lmp.command(
            'fix             f2 all print 1 "${dU1} ${dU2} ${flambda}" '
            'title "# dU_1[eV/atom] dU_2[eV/atom] lambda" '
            "screen no file forward_%d.dat"
            % iteration
        )
        lmp.command("run             %d" % self.calc._n_switching_steps)

        lmp.command("unfix           f2")
        for compute_id in compute_ids:
            lmp.command("uncompute       %s" % compute_id)
        for swap_fix in swap_fixes:
            lmp.command("unfix %s" % swap_fix)

        # now equilibrate with the final potential on its own
        lmp.command(self._pure_pair_style_command(final))
        for command in self._pure_pair_coeff_commands(final):
            lmp.command(command)

        # Thermo output.
        lmp.command("thermo_style    custom step pe")
        lmp.command("thermo          1000")

        # run eqbrm run
        lmp.command("run             %d" % self.calc.n_equilibration_steps)

        # ---------------------------------------------------------------
        # BKD cycle: initial potential scaled 0 -> 1, final potential 1 -> 0
        # ---------------------------------------------------------------
        lmp.command("variable         flambda equal ramp(${lf},${li})")
        lmp.command("variable         blambda equal ramp(${li},${lf})")

        commands, compute_ids = self._scaled_pair_commands(
            initial, final, "v_flambda", "v_blambda"
        )
        for command in commands:
            lmp.command(command)

        swap_fixes = self._add_swap_fixes(
            lmp, self.calc.monte_carlo.reverse_swap_types, "Reverse"
        )

        # Thermo output.
        lmp.command("thermo_style    custom step v_dU1 v_dU2")
        lmp.command("thermo          1000")

        # save the necessary items to a file: first step
        lmp.command(
            'fix             f2 all print 1 "${dU1} ${dU2} ${flambda}" '
            'title "# dU_1[eV/atom] dU_2[eV/atom] lambda" '
            "screen no file backward_%d.dat"
            % iteration
        )
        lmp.command("run             %d" % self.calc._n_switching_steps)

        lmp.command("unfix           f2")
        for compute_id in compute_ids:
            lmp.command("uncompute       %s" % compute_id)
        for swap_fix in swap_fixes:
            lmp.command("unfix %s" % swap_fix)

        self.lammps_close(lmp=lmp)
        lmp.rotate_logs("integration")

    def _add_swap_fixes(self, lmp, swap_types, pass_name):
        """
        Add ``fix atom/swap`` moves between every pair of ``swap_types`` when
        Monte Carlo swaps are requested; returns the fix ids that were added.
        """
        mc = self.calc.monte_carlo
        if not (mc.n_swaps > 0 and len(swap_types) >= 2):
            return []

        swap_combos = list(itertools.combinations(swap_types, 2))
        self.logger.info(
            f"{pass_name} pass: {mc.n_swaps} swap moves per combo, "
            f"{len(swap_combos)} combinations every {mc.n_steps}"
        )
        for combo in swap_combos:
            self.logger.info(f"  Swapping types: {combo[0]} <-> {combo[1]}")

        fix_ids = []
        for idx, (type1, type2) in enumerate(swap_combos):
            fix_id = "swap%d" % idx
            extra = " noforce yes localE yes" if mc.use_custom_lammps else ""
            lmp.command(
                "fix  %s all atom/swap %d %d %d %f ke no types %s %s%s"
                % (
                    fix_id,
                    mc.n_steps,
                    mc.n_swaps,
                    np.random.randint(1, 10000),
                    self.calc._temperature,
                    type1,
                    type2,
                    extra,
                )
            )
            fix_ids.append(fix_id)
        return fix_ids

    def thermodynamic_integration(self):
        """
        Calculate free energy after integration step

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        Calculates the final work, energy dissipation; In alchemical mode, there is reference system,
        the calculated free energy is the same as the work.
        """
        if self.calc.alchemy_coupling:
            ws, qs = [], []
            for i in range(1, self.calc.n_iterations + 1):
                fwd = np.loadtxt(os.path.join(self.simfolder,
                                              "forward_%d.dat" % i))
                bkd = np.loadtxt(os.path.join(self.simfolder,
                                              "backward_%d.dat" % i))
                wf = np.trapezoid(fwd[:, 0], fwd[:, 2])          # lam 0 -> 1
                wb = -np.trapezoid(bkd[:, 0], bkd[:, 2])         # lam 1 -> 0
                ws.append(0.5 * (wf + wb))
                qs.append(0.5 * (wf - wb))
            self.w = float(np.mean(ws))
            self.qdiss = float(np.mean(qs))
            self.ferr = float(np.std(ws))
            self.fe = self.w
            return

        w, q, qerr = find_w(self.simfolder, self.calc, full=True, solid=False)

        self.w = w
        self.qdiss = q
        self.ferr = qerr
        self.fe = self.w

    def mass_integration(self, ref_mass, target_masses, target_counts):
        mcorsum = integrate_mass(
            ref_mass,
            target_masses,
            target_counts,
            self.calc._temperature,
            self.natoms,
        )
        return mcorsum
