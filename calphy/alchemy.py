"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^2, Ralf Drautz^2
^1: Max Planck Institut für Eisenforschung, Dusseldorf, Germany
^2: Ruhr-University Bochum, Bochum, Germany

calphy is published and distributed under the Academic Software License v1.0 (ASL).
calphy is distributed in the hope that it will be useful for non-commercial academic research,
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
calphy API is published and distributed under the BSD 3-Clause "New" or "Revised" License
See the LICENSE FILE for more details.

More information about the program can be found in:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz.
“Automated Free Energy Calculation from Atomistic Simulations.” Physical Review Materials 5(10), 2021
DOI: 10.1103/PhysRevMaterials.5.103801

For more information contact:
sarath.menon@ruhr-uni-bochum.de/yury.lysogorskiy@icams.rub.de

Notes
-----
- swapping is strictly only performed between types 1 and 2 at the moment; this needs to be refined further
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

    """

    def __init__(self, calculation=None, simfolder=None, log_to_screen=False):

        # call base class
        super().__init__(
            calculation=calculation, simfolder=simfolder, log_to_screen=log_to_screen,
        )

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
        lmp = ph.create_object(self.calc, self.simfolder)
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
        lmp = ph.create_object(self.calc, self.simfolder)

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
            lmp.command(f"pair_style {self.calc._pair_style_with_options[0]}")

            # set up structure
            lmp = ph.create_structure(lmp, self.calc)

            # set up potential
            lmp.command(f"pair_coeff {self.calc.pair_coeff[0]}")
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

        # create lammps object
        lmp = ph.create_object(self.calc, self.simfolder)

        # Adiabatic switching parameters.
        lmp.command("variable        li       equal   1.0")
        lmp.command("variable        lf       equal   0.0")

        lmp.command(f"pair_style {self.calc._pair_style_with_options[0]}")

        # read dump file
        # conf = os.path.join(self.simfolder, "conf.equilibration.dump")
        conf = os.path.join(self.simfolder, "conf.equilibration.data")
        lmp = ph.read_data(lmp, conf)

        # set up hybrid potential
        # here we only need to set one potential
        lmp.command(f"pair_coeff {self.calc.pair_coeff[0]}")
        lmp = ph.set_mass(lmp, self.calc)

        # NEW ADDED
        lmp.command("group g1 type 1")
        lmp.command("group g2 type 2")
        # lmp = ph.set_double_hybrid_potential(lmp, self.options, self.calc._pressureair_style, self.calc._pressureair_coeff)

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
        # FWD cycle
        # ---------------------------------------------------------------
        lmp.command("variable         flambda equal ramp(${li},${lf})")
        lmp.command("variable         blambda equal ramp(${lf},${li})")

        # lmp.command("pair_style       hybrid/scaled v_flambda %s v_blambda ufm 7.5"%self.options["md"]["pair_style"])

        # Compute pair definitions
        if self.calc.pair_style[0] == self.calc.pair_style[1]:
            pc = self.calc.pair_coeff[0]
            pcraw = pc.split()
            pc1 = " ".join(
                [
                    *pcraw[:2],
                    *[
                        self.calc._pair_style_names[0],
                    ],
                    "1",
                    *pcraw[2:],
                ]
            )
            pc = self.calc.pair_coeff[1]
            pcraw = pc.split()
            pc2 = " ".join(
                [
                    *pcraw[:2],
                    *[
                        self.calc._pair_style_names[1],
                    ],
                    "2",
                    *pcraw[2:],
                ]
            )
        else:
            pc = self.calc.pair_coeff[0]
            pcraw = pc.split()
            pc1 = " ".join(
                [
                    *pcraw[:2],
                    *[
                        self.calc._pair_style_names[0],
                    ],
                    *pcraw[2:],
                ]
            )
            pc = self.calc.pair_coeff[1]
            pcraw = pc.split()
            pc2 = " ".join(
                [
                    *pcraw[:2],
                    *[
                        self.calc._pair_style_names[1],
                    ],
                    *pcraw[2:],
                ]
            )

        lmp.command(
            "pair_style       hybrid/scaled v_flambda %s v_blambda %s"
            % (
                self.calc._pair_style_with_options[0],
                self.calc._pair_style_with_options[1],
            )
        )
        lmp.command("pair_coeff       %s" % pc1)
        lmp.command("pair_coeff       %s" % pc2)

        # apply pair force commands
        if self.calc._pair_style_names[0] == self.calc._pair_style_names[1]:
            lmp.command(
                "compute         c1 all pair %s 1" % self.calc._pair_style_names[0]
            )
            lmp.command(
                "compute         c2 all pair %s 2" % self.calc._pair_style_names[1]
            )
        else:
            lmp.command(
                "compute         c1 all pair %s" % self.calc._pair_style_names[0]
            )
            lmp.command(
                "compute         c2 all pair %s" % self.calc._pair_style_names[1]
            )

        # Output variables.
        lmp.command("variable        step equal step")
        lmp.command(
            "variable        dU1 equal c_c1/atoms"
        )  # Driving-force obtained from NEHI procedure.
        lmp.command("variable        dU2 equal c_c2/atoms")

        # add swaps if n_swap is > 0 - forward pass
        if (
            self.calc.monte_carlo.n_swaps > 0
            and len(self.calc.monte_carlo.forward_swap_types) >= 2
        ):
            swap_types = self.calc.monte_carlo.forward_swap_types
            swap_combos = list(itertools.combinations(swap_types, 2))
            self.logger.info(
                f"Forward pass: {self.calc.monte_carlo.n_swaps} swap moves per combo, {len(swap_combos)} combinations every {self.calc.monte_carlo.n_steps}"
            )
            for combo in swap_combos:
                self.logger.info(f"  Swapping types: {combo[0]} <-> {combo[1]}")

            for idx, (type1, type2) in enumerate(swap_combos):
                swap_str = f"{type1} {type2}"
                if self.calc.monte_carlo.use_custom_lammps:
                    lmp.command(
                        "fix  swap%d all atom/swap %d %d %d %f ke no types %s noforce yes localE yes"
                        % (
                            idx,
                            self.calc.monte_carlo.n_steps,
                            self.calc.monte_carlo.n_swaps,
                            np.random.randint(1, 10000),
                            self.calc._temperature,
                            swap_str,
                        )
                    )
                else:
                    lmp.command(
                        "fix  swap%d all atom/swap %d %d %d %f ke no types %s"
                        % (
                            idx,
                            self.calc.monte_carlo.n_steps,
                            self.calc.monte_carlo.n_swaps,
                            np.random.randint(1, 10000),
                            self.calc._temperature,
                            swap_str,
                        )
                    )

            # Use the first swap fix for output tracking
            # lmp.command("variable a equal f_swap0[1]")
            # lmp.command("variable b equal f_swap0[2]")
            # lmp.command(
            #    'fix             swap_print all print 1 "${a} ${b} ${flambda}" screen no file swap.forward_%d.dat'
            #    % iteration
            # )

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

        # now equilibrate at the second potential
        lmp.command("unfix           f2")
        lmp.command("uncompute       c1")
        lmp.command("uncompute       c2")

        # NEW SWAP
        if self.calc.monte_carlo.n_swaps > 0:
            swap_types = self.calc.monte_carlo.forward_swap_types
            swap_combos = list(itertools.combinations(swap_types, 2))
            for idx in range(len(swap_combos)):
                lmp.command(f"unfix swap{idx}")
            # lmp.command("unfix swap_print")

        lmp.command("pair_style      %s" % self.calc._pair_style_with_options[1])
        lmp.command("pair_coeff      %s" % self.calc.pair_coeff[1])

        # Thermo output.
        lmp.command("thermo_style    custom step pe")
        lmp.command("thermo          1000")

        # run eqbrm run
        lmp.command("run             %d" % self.calc.n_equilibration_steps)

        # reverse switching
        lmp.command("variable         flambda equal ramp(${lf},${li})")
        lmp.command("variable         blambda equal ramp(${li},${lf})")

        lmp.command(
            "pair_style       hybrid/scaled v_flambda %s v_blambda %s"
            % (
                self.calc._pair_style_with_options[0],
                self.calc._pair_style_with_options[1],
            )
        )
        lmp.command("pair_coeff       %s" % pc1)
        lmp.command("pair_coeff       %s" % pc2)

        # apply pair force commands
        if self.calc._pair_style_names[0] == self.calc._pair_style_names[1]:
            lmp.command(
                "compute         c1 all pair %s 1" % self.calc._pair_style_names[0]
            )
            lmp.command(
                "compute         c2 all pair %s 2" % self.calc._pair_style_names[1]
            )
        else:
            lmp.command(
                "compute         c1 all pair %s" % self.calc._pair_style_names[0]
            )
            lmp.command(
                "compute         c2 all pair %s" % self.calc._pair_style_names[1]
            )

        # Output variables.
        lmp.command("variable        step equal step")
        lmp.command(
            "variable        dU1 equal c_c1/atoms"
        )  # Driving-force obtained from NEHI procedure.
        lmp.command("variable        dU2 equal c_c2/atoms")

        # add swaps if n_swap is > 0 - reverse pass
        if (
            self.calc.monte_carlo.n_swaps > 0
            and len(self.calc.monte_carlo.reverse_swap_types) >= 2
        ):
            swap_types = self.calc.monte_carlo.reverse_swap_types
            swap_combos = list(itertools.combinations(swap_types, 2))
            self.logger.info(
                f"Reverse pass: {self.calc.monte_carlo.n_swaps} swap moves per combo, {len(swap_combos)} combinations every {self.calc.monte_carlo.n_steps}"
            )
            for combo in swap_combos:
                self.logger.info(f"  Swapping types: {combo[0]} <-> {combo[1]}")

            for idx, (type1, type2) in enumerate(swap_combos):
                swap_str = f"{type1} {type2}"
                if self.calc.monte_carlo.use_custom_lammps:
                    lmp.command(
                        "fix  swap%d all atom/swap %d %d %d %f ke no types %s noforce yes localE yes"
                        % (
                            idx,
                            self.calc.monte_carlo.n_steps,
                            self.calc.monte_carlo.n_swaps,
                            np.random.randint(1, 10000),
                            self.calc._temperature,
                            swap_str,
                        )
                    )
                else:
                    lmp.command(
                        "fix  swap%d all atom/swap %d %d %d %f ke no types %s"
                        % (
                            idx,
                            self.calc.monte_carlo.n_steps,
                            self.calc.monte_carlo.n_swaps,
                            np.random.randint(1, 10000),
                            self.calc._temperature,
                            swap_str,
                        )
                    )

            # Use the first swap fix for output tracking
            # lmp.command("variable a equal f_swap0[1]")
            # lmp.command("variable b equal f_swap0[2]")
            # lmp.command(
            #'fix             swap_print all print 1 "${a} ${b} ${blambda}" screen no file swap.backward_%d.dat'
            #    % iteration
            # )

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

        # now equilibrate at the second potential
        lmp.command("unfix           f2")
        lmp.command("uncompute       c1")
        lmp.command("uncompute       c2")

        if self.calc.monte_carlo.n_swaps > 0:
            swap_types = self.calc.monte_carlo.reverse_swap_types
            swap_combos = list(itertools.combinations(swap_types, 2))
            for idx in range(len(swap_combos)):
                lmp.command(f"unfix swap{idx}")
            # lmp.command("unfix swap_print")

        self.lammps_close(lmp=lmp)
        lmp.rotate_logs("integration")

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
