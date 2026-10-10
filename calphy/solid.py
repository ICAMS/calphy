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

import itertools
import numpy as np
import yaml
import copy
import sys
import os

from calphy.integrators import *
import calphy.helpers as ph
import calphy.phase as cph
from calphy.errors import *


class Solid(cph.Phase):
    """
    Class for free energy calculation with solid as the reference state

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
            calculation=calculation,
            simfolder=simfolder,
            log_to_screen=log_to_screen,
        )

    def run_spring_constant_convergence(self, lmp):
        """ """
        if self.calc._qtb:
            qtb = self.calc.quantum_thermal_bath
            lmp.command("fix              3 all nve")
            lmp.command(
                "fix              3q all qtb temp %f damp %f seed %d f_max %f N_f %d"
                % (
                    self.calc._temperature,
                    qtb.thermostat_damping,
                    np.random.randint(1, 10**8),
                    qtb.f_max,
                    qtb.n_f,
                )
            )
        else:
            lmp.command(
                "fix              3 all nvt temp %f %f %f"
                % (
                    self.calc._temperature,
                    self.calc._temperature,
                    self.calc.md.thermostat_damping[1],
                )
            )

        # apply fix
        lmp = ph.compute_msd(lmp, self.calc)

        if ph.check_if_any_is_none(self.calc.spring_constants):
            # similar averaging routine
            laststd = 0.00
            for i in range(self.calc.md.n_cycles):
                lmp.command("run              %d" % int(self.calc.md.n_small_steps))
                lmp.sync()  # flush before analyse_spring_constants reads msd.dat
                k_mean, k_std = self.analyse_spring_constants(lmp)
                self.logger.info(
                    "At count %d mean k is %f std is %f" % (i + 1, k_mean[0], k_std[0])
                )
                if np.abs(laststd - k_std[0]) < self.calc.tolerance.spring_constant:
                    # now reevaluate spring constants
                    self.assign_spring_constants(k_mean)
                    break
                laststd = k_std[0]

        else:
            if not (len(self.calc.spring_constants) == self.calc.n_elements):
                raise ValueError(
                    "Spring constant input length should be same as number of elements, spring constant length %d, # elements %d"
                    % (len(self.calc.spring_constants), self.calc.n_elements)
                )

            # still run a small NVT cycle
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))
            self.k = self.calc.spring_constants
            self.logger.info("Used user input sprint constants")
            self.logger.info(self.k)

        if self.calc._qtb:
            lmp.command("unfix         3q")
        lmp.command("unfix         3")

    def analyse_spring_constants(self, lmp):
        """
        Analyse spring constant routine
        """
        ncount = int(self.calc.n_equilibration_steps) // int(
            self.calc.md.n_every_steps * self.calc.md.n_repeat_steps
        )

        # now we can check if it converted; read cumulative msd.dat across segments
        msd = lmp.read_timeseries("msd.dat")
        k_mean = []
        k_std = []
        for i in range(self.calc.n_elements):
            quant = msd[:, i + 1][-ncount + 1 :]
            mean_quant = np.round(np.mean(quant), decimals=2)
            std_quant = np.round(np.std(quant), decimals=2)
            if mean_quant == 0:
                self.logger.warning(
                    "MSD for element index %d averaged to ~0; spring-constant "
                    "estimation is unreliable. Consider increasing equilibration "
                    "time or providing spring_constants explicitly." % i
                )
                mean_quant = 1.00
            mean_quant = 3 * kb * self.calc._temperature / mean_quant
            k_mean.append(mean_quant)
            k_std.append(std_quant)
        return k_mean, k_std

    def _add_fe_swap_fixes(self, lmp, temperature, ke, prefix):
        """
        Add ``fix atom/swap`` between every pair of species.

        Parameters
        ----------
        lmp : LAMMPS runner
        temperature : str
            swap temperature argument, a number or an equal-style variable
            reference such as ``v_tswap`` (thermoatoms LAMMPS fork)
        ke : str
            ``yes`` or ``no``, whether velocities are rescaled on a swap
        prefix : str
            fix ids are ``<prefix><index>``

        Returns
        -------
        list of str
            the fix ids that were added
        """
        mc = self.calc.monte_carlo
        pairs = list(itertools.combinations(range(1, self.calc.n_elements + 1), 2))
        self.logger.info(
            "%s: %d swap moves between each of %d species pairs every %d steps",
            prefix, mc.n_swaps, len(pairs), mc.n_steps,
        )
        extra = " noforce yes localE yes" if mc.use_custom_lammps else ""
        fix_ids = []
        for idx, (type1, type2) in enumerate(pairs):
            fix_id = "%s%d" % (prefix, idx)
            lmp.command(
                "fix               %s all atom/swap %d %d %d %s ke %s types %d %d%s"
                % (
                    fix_id,
                    mc.n_steps,
                    mc.n_swaps,
                    np.random.randint(1, 10000),
                    temperature,
                    ke,
                    type1,
                    type2,
                    extra,
                )
            )
            fix_ids.append(fix_id)
        # cumulative attempted and accepted swaps summed over the pairs
        lmp.command(
            "variable          swapatt equal %s"
            % "+".join("f_%s[1]" % fix_id for fix_id in fix_ids)
        )
        lmp.command(
            "variable          swapacc equal %s"
            % "+".join("f_%s[2]" % fix_id for fix_id in fix_ids)
        )
        return fix_ids

    def _fe_swap_statistics_file(self, lmp, filename):
        """
        ``fix print`` of step, lambda and the cumulative attempted and
        accepted swap counts every ``monte_carlo.n_steps`` steps, the
        acceptance diagnostic of a switching leg.
        """
        lmp.command(
            'fix               f5 all print %d "${step} ${lambda} ${swapatt} ${swapacc}" '
            'title "# step lambda swaps_attempted swaps_accepted" screen no file %s'
            % (self.calc.monte_carlo.n_steps, filename)
        )

    def run_configurational_equilibration(self, lmp):
        """
        Equilibrate the arrangement of the species on the lattice sites.

        Runs NVT at the target temperature with the full potential and Monte
        Carlo swaps between all species pairs, in blocks of
        ``md.n_small_steps``, until the block-averaged potential energy per
        atom changes by less than ``tolerance.swap_energy`` between
        successive blocks (at most ``md.n_cycles`` blocks).  This is the
        slow step of a switching run with swaps: the acceptance rate here is
        that of the real alloy at T.  The block averages and the cumulative
        swap counters are written to ``swap.equilibration.dat``.
        """
        temp = self.calc._temperature
        lmp.command(
            "fix               3 all nvt temp %f %f %f"
            % (temp, temp, self.calc.md.thermostat_damping[1])
        )
        swap_fixes = self._add_fe_swap_fixes(lmp, "%f" % temp, "yes", "ceq")

        n_every = int(self.calc.md.n_every_steps)
        n_repeat = int(self.calc.md.n_repeat_steps)
        lmp.command("variable          peatom equal pe/atoms")
        lmp.command(
            "fix               5 all ave/time %d %d %d v_peatom v_swapatt v_swapacc "
            'title2 "# TimeStep pe[eV/atom] swaps_attempted swaps_accepted" '
            "file swap.equilibration.dat"
            % (n_every, n_repeat, n_every * n_repeat)
        )

        ncount = int(self.calc.md.n_small_steps) // (n_every * n_repeat)
        tol = self.calc.tolerance.swap_energy
        last_pe = None
        converged = False
        for i in range(self.calc.md.n_cycles):
            lmp.command("run               %d" % int(self.calc.md.n_small_steps))
            lmp.sync()  # flush before reading swap.equilibration.dat
            data = np.atleast_2d(lmp.read_timeseries("swap.equilibration.dat"))
            pe = float(np.mean(data[-ncount:, 1]))
            attempted = float(data[-1, 2])
            accepted = float(data[-1, 3])
            rate = accepted / attempted if attempted > 0 else 0.0
            self.logger.info(
                "configurational equilibration block %d: pe %.6f eV/atom, "
                "swaps accepted %.0f of %.0f (%.3f)",
                i + 1, pe, accepted, attempted, rate,
            )
            if last_pe is not None and abs(pe - last_pe) < tol:
                converged = True
                break
            last_pe = pe
        if not converged:
            self.logger.warning(
                "configurational equilibration did not converge to %.1e eV/atom "
                "in %d blocks; raise md.n_cycles or md.n_small_steps, or "
                "monte_carlo.n_swaps", tol, self.calc.md.n_cycles,
            )

        lmp.command("unfix             5")
        for fix_id in swap_fixes:
            lmp.command("unfix             %s" % fix_id)
        lmp.command("unfix             3")

    def assign_spring_constants(self, k):
        """
        Here the spring constants are finalised, add added to the class
        """
        # first replace any provided values with user values
        if ph.check_if_any_is_not_none(self.calc.spring_constants):
            spring_constants = copy.copy(self.calc.spring_constants)
            k = ph.replace_nones(spring_constants, k, logger=self.logger)

        # add sanity checks
        k = ph.validate_spring_constants(k, logger=self.logger)

        # now save
        self.k = k

        self.logger.info("finalized sprint constants")
        self.logger.info(self.k)

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
        If `fix_lattice` option is True, then the input structure is used as it is and the corresponding pressure
        is calculated.
        At the end of the run, the averaged box dimensions are calculated.
        """

        lmp = self._open_lammps()

        # set up potential
        lmp = ph.set_pair_style(lmp, self.calc)

        # set up structure
        lmp = ph.create_structure(lmp, self.calc)

        lmp = ph.set_pair_coeff(lmp, self.calc)
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

        # Run if a constrained lattice is not needed
        if not self.calc._fix_lattice:
            if self.calc._pressure == 0:
                self.run_zero_pressure_equilibration(lmp)
            else:
                self.run_finite_pressure_equilibration(lmp)

            # equilibration-frame dump (post warm-up; no-op unless
            # n_print_steps_equilibration > 0)
            self.start_equilibration_dump(lmp)

            # this is when the averaging routine starts
            self.run_pressure_convergence(lmp)

            # dump snapshot and check if melted
            self.dump_current_snapshot(lmp, "traj.equilibration_stage1.dat")
            self.check_if_melted(lmp, "traj.equilibration_stage1.dat")

        # run if a constrained lattice is used
        else:
            self.start_equilibration_dump(lmp)
            # routine in which lattice constant will not varied, but is set to a given fixed value
            self.run_constrained_pressure_convergence(lmp)

        # start MSD calculation routine
        # there two possibilities here - if spring constants are provided, use it. If not, calculate it
        self.run_spring_constant_convergence(lmp)

        # with Monte Carlo swaps, equilibrate the arrangement of the species
        # on the lattice sites at T before the switching starts
        if self.calc._fe_swaps:
            self.run_configurational_equilibration(lmp)

        # check for melting
        self.stop_equilibration_dump(lmp)
        self.dump_current_snapshot(lmp, "traj.equilibration_stage2.dat")
        self.check_if_melted(lmp, "traj.equilibration_stage2.dat")
        lmp = ph.write_data(lmp, "conf.equilibration.data")
        # close object and process traj
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
        lmp = self._open_lammps()

        # set up potential
        lmp = ph.set_pair_style(lmp, self.calc)

        # read in the conf file: the equilibrated configuration, or for later
        # iterations the end of the previous backward leg (see
        # _integration_start_configuration)
        conf = self._integration_start_configuration(iteration)
        lmp = ph.read_data(lmp, conf)

        lmp = ph.set_pair_coeff(lmp, self.calc)
        lmp = ph.set_mass(lmp, self.calc)
        if self.calc._fe_swaps:
            # One common mass for every species.  The switching work is a
            # configurational integral and does not depend on the masses;
            # with equal masses a swap leaves the kinetic energy and the
            # mass weighted centre of mass unchanged, and every arrangement
            # of the species is exactly degenerate at the Einstein end.
            # thermodynamic_integration evaluates the reference with the
            # real masses (Einstein term) and this mass (CM correction).
            lmp = ph.set_uniform_mass(lmp, self.calc)

        # remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        # create groups - each species belong to one group
        for i in range(self.calc.n_elements):
            lmp.command("group  g%d type %d" % (i + 1, i + 1))

        # get counts of each group
        for i in range(self.calc.n_elements):
            lmp.command("variable   count%d equal count(g%d)" % (i + 1, i + 1))

        # initialise everything
        lmp.command("run               0")

        # apply initial fixes
        lmp.command("fix               f1 all nve")

        # apply fix for each spring
        # TODO: Add option to select function
        for i in range(self.calc.n_elements):
            lmp.command(
                "fix               ff%d g%d ti/spring 10.0 100 100 function 2"
                % (i + 1, i + 1)
            )

        # apply temp fix
        if self.calc._qtb:
            qtb = self.calc.quantum_thermal_bath
            lmp.command(
                "fix               f3 all qtb temp %f damp %f seed %d f_max %f N_f %d"
                % (
                    self.calc._temperature,
                    qtb.thermostat_damping,
                    np.random.randint(1, 10**8),
                    qtb.f_max,
                    qtb.n_f,
                )
            )
            # QTB does not consume a base temperature compute, so the temp/com
            # group/correction trick used by langevin does not apply.
            lmp.command("compute           Tcm all temp/com")
        else:
            lmp.command(
                "fix               f3 all langevin %f %f %f %d zero yes"
                % (
                    self.calc._temperature,
                    self.calc._temperature,
                    self.calc.md.thermostat_damping[1],
                    np.random.randint(1, 10000),
                )
            )

            # compute com and apply to fix
            lmp.command("compute           Tcm all temp/com")
            lmp.command("fix_modify        f3 temp Tcm")

        lmp.command("variable          step    equal step")
        lmp.command("variable          dU1      equal pe/atoms")
        for i in range(self.calc.n_elements):
            lmp.command("variable          dU%d      equal f_ff%d" % (i + 2, i + 1))

        lmp.command("variable          lambda  equal f_ff1[1]")

        # add thermo command to force variable evaluation
        lmp.command("thermo_style      custom step pe c_Tcm")
        lmp.command("thermo            10000")

        # Create velocity
        lmp.command(
            "velocity          all create %f %d mom yes rot yes dist gaussian"
            % (self.calc._temperature, np.random.randint(1, 10000))
        )

        # Both equilibration blocks of the cycle are warm starts.  The first
        # one re-thermalises an already equilibrated configuration whose
        # velocities were just regenerated; the second holds the system at the
        # Einstein-crystal end (lambda = 1), a set of damped harmonic
        # oscillators that forgets its initial state within a few thermostat
        # relaxation times.  fix ti/spring drives its lambda schedule off a
        # single t_equil, so the two blocks must have the same length and it
        # is passed to the fix here.
        n_equil = self._warm_start_steps(npt=False)
        self.logger.info(
            "integration iteration %d: equilibration blocks of %d steps "
            "(warm start; n_equilibration_steps = %d)",
            iteration, n_equil, self.calc.n_equilibration_steps,
        )

        # reapply
        for i in range(self.calc.n_elements):
            lmp.command(
                "fix               ff%d g%d ti/spring %f %d %d function 2"
                % (
                    i + 1,
                    i + 1,
                    self.k[i],
                    self.calc._n_switching_steps,
                    n_equil,
                )
            )

        # Monte Carlo swaps of species between lattice sites, active through
        # both equilibration blocks and both switching legs.  The springs
        # belong to the sites (the ti/spring groups are static), so a swap
        # changes only the pair energy, which enters H(lambda) as
        # (1 - lambda) U.  fix atom/swap tests exp(-dU / k T_swap) on the
        # bare pair energy, so T_swap = T / (1 - lambda) gives the Metropolis
        # test of the switched Hamiltonian; the offset keeps the variable
        # finite at lambda = 1, where every swap is accepted and the
        # arrangement is randomised.  Needs the thermoatoms fork, which
        # takes an equal-style variable as swap temperature.
        swap_fixes = []
        if self.calc._fe_swaps:
            lmp.command(
                "variable          tswap   equal %f/(1.0-v_lambda+1.0e-12)"
                % self.calc._temperature
            )
            swap_fixes = self._add_fe_swap_fixes(lmp, "v_tswap", "no", "swap")

        # Equilibriate structure
        lmp.command("run               %d" % n_equil)

        # write out energy
        str1 = 'fix f4 all print 1 "${dU1} '
        str2 = []
        for i in range(self.calc.n_elements):
            str2.append("${dU%d}" % (i + 2))

        str2.append('${lambda}"')
        str2 = " ".join(str2)
        title_cols = (
            ["dU_sys[eV/atom]"]
            + ["dU_ref%d[eV/atom]" % (i + 1) for i in range(self.calc.n_elements)]
            + ["lambda"]
        )
        str3 = ' title "# %s" screen no file forward_%d.dat' % (
            " ".join(title_cols),
            iteration,
        )
        command = str1 + str2 + str3
        lmp.command(command)

        if self.calc.n_print_steps > 0:
            lmp.command(
                "dump              d1 all custom %d traj.fe.forward_%d.dat id type mass x y z fx fy fz"
                % (self.calc.n_print_steps, iteration)
            )

        # turn on swap moves
        # if self.calc.monte_carlo.n_swaps > 0:
        #    self.logger.info(f'{self.calc.monte_carlo.n_swaps} swap moves are performed between 1 and 2 every {self.calc.monte_carlo.n_steps}')
        #    lmp.command("fix  swap all atom/swap %d %d %d %d ke yes types 1 2"%(self.calc.monte_carlo.n_steps,
        #                                                                        self.calc.monte_carlo.n_swaps,
        #                                                                        np.random.randint(1, 10000),
        #                                                                        self.calc._temperature))
        #
        #    lmp.command("variable a equal f_swap[1]")
        #    lmp.command("variable b equal f_swap[2]")
        #    lmp.command("fix             swap2 all print 1 \"${a} ${b}\" screen no file swap.fe.forward_%d.dat"%iteration)

        if swap_fixes:
            self._fe_swap_statistics_file(lmp, "swap.fe.forward_%d.dat" % iteration)

        # Forward switching over ts steps
        lmp.command("run               %d" % self.calc._n_switching_steps)
        lmp.command("unfix             f4")
        if swap_fixes:
            lmp.command("unfix             f5")

        if self.calc.n_print_steps > 0:
            lmp.command("undump           d1")

        # if self.calc.monte_carlo.n_swaps > 0:
        #    lmp.command("unfix swap")
        #    lmp.command("unfix swap2")

        # Equilibrate at the Einstein-crystal end.  fix ti/spring holds
        # lambda at exactly 1 for this whole run, so the interatomic forces
        # enter as (1 - lambda) * f = 0 and the pair style would be evaluated
        # every step only to be multiplied by zero.  Switching the pair
        # compute off leaves the trajectory bitwise identical and makes this
        # block essentially free; it is switched back on before the backward
        # leg, which needs the real forces and energies again.
        # With swaps the pair energies are still needed for the Metropolis
        # test, so the pair compute stays on.
        if not swap_fixes:
            lmp.command("pair_modify       compute no")
        lmp.command("run               %d" % n_equil)
        if not swap_fixes:
            lmp.command("pair_modify       compute yes")

        # write out energy
        str1 = 'fix f4 all print 1 "${dU1} '
        str2 = []
        for i in range(self.calc.n_elements):
            str2.append("${dU%d}" % (i + 2))

        str2.append('${lambda}"')
        str2 = " ".join(str2)
        title_cols = (
            ["dU_sys[eV/atom]"]
            + ["dU_ref%d[eV/atom]" % (i + 1) for i in range(self.calc.n_elements)]
            + ["lambda"]
        )
        str3 = ' title "# %s" screen no file backward_%d.dat' % (
            " ".join(title_cols),
            iteration,
        )
        command = str1 + str2 + str3
        lmp.command(command)

        if self.calc.n_print_steps > 0:
            lmp.command(
                "dump              d1 all custom %d traj.fe.backward_%d.dat id type mass x y z fx fy fz"
                % (self.calc.n_print_steps, iteration)
            )

        # add swaps if n_swap is > 0
        # if self.calc.monte_carlo.n_swaps > 0:
        #    self.logger.info(f'{self.calc.monte_carlo.n_swaps} swap moves are performed between 1 and 2 every {self.calc.monte_carlo.n_steps}')
        #    lmp.command("fix  swap all atom/swap %d %d %d %d ke yes types 2 1"%(self.calc.monte_carlo.n_steps,
        #                                                                        self.calc.monte_carlo.n_swaps,
        #                                                                        np.random.randint(1, 10000),
        #                                                                        self.calc._temperature))
        #
        #    lmp.command("variable a equal f_swap[1]")
        #    lmp.command("variable b equal f_swap[2]")
        #    lmp.command("fix             swap2 all print 1 \"${a} ${b}\" screen no file swap.fe.backward_%d.dat"%iteration)

        if swap_fixes:
            self._fe_swap_statistics_file(lmp, "swap.fe.backward_%d.dat" % iteration)

        # Reverse switching over ts steps
        lmp.command("run               %d" % self.calc._n_switching_steps)
        lmp.command("unfix             f4")
        if swap_fixes:
            lmp.command("unfix             f5")
            for fix_id in swap_fixes:
                lmp.command("unfix             %s" % fix_id)

        if self.calc.n_print_steps > 0:
            lmp.command("undump           d1")

        # if self.calc.monte_carlo.n_swaps > 0:
        #    lmp.command("unfix swap")
        #    lmp.command("unfix swap2")

        # the real system at T again: starting point of the next iteration
        lmp = ph.write_data(lmp, "conf.fe.backward_%d.data" % iteration)

        # close object
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
        Calculates the final work, energy dissipation and free energy by
        matching with Einstein crystal
        """
        use_quantum_reference = self.calc._qtb
        if use_quantum_reference:
            self.logger.info(
                "Using quantum harmonic-oscillator Einstein-crystal reference "
                "(required for self-consistency with QTB sampling)."
            )
        fe, fcm = get_einstein_crystal_fe(
            self.calc, self.vol, self.k,
            return_contributions=True,
            quantum=use_quantum_reference,
        )

        w, q, qerr = find_w(self.simfolder, self.calc, full=True, solid=True)

        self.fref = fe + fcm
        self.feinstein = fe
        self.fcm = fcm
        self.w = w
        self.qdiss = q
        self.ferr = qerr

        # add pressure contribution if required
        if self.calc._pressure != 0:
            p = self.calc._pressure / EV_A3_TO_BAR
            v = self.vol / self.natoms
            self.pv = p * v
        else:
            self.pv = 0

        # calculate final free energy
        self.fe = self.fref + self.w + self.pv
