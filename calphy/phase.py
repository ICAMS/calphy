"""
calphy: a Python library and command line interface for automated free
energy calculations.

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

"""

import numpy as np
import yaml
import copy
import os
import shutil
import itertools

from calphy.integrators import *
import calphy.helpers as ph
from calphy.errors import *
from calphy.input import generate_metadata


class Phase:
    """
    Class for free energy calculation.

    Parameters
    ----------
    input : Calculation class
        input options

    simfolder : string
        base folder for running calculations

    """

    def __init__(self, calculation=None, simfolder=None, log_to_screen=False, lmp=None):

        self.calc = copy.deepcopy(calculation)

        # serialise input
        indict = {"calculations": [self.calc.model_dump()]}
        with open(os.path.join(simfolder, "input_file.yaml"), "w") as fout:
            yaml.safe_dump(indict, fout)

        self.simfolder = simfolder
        self.log_to_screen = log_to_screen
        self.publications = []

        logfile = os.path.join(self.simfolder, "calphy.log")
        self.logger = ph.prepare_log(logfile, screen=log_to_screen)

        if self.calc._pressure is None:
            pressure_string = "None"
        else:
            pressure_string = "%f" % self.calc._pressure

        self.logger.info(
            "Temperature start: %f K, temperature stop: %f K, pressure: %s bar"
            % (self.calc._temperature, self.calc._temperature_stop, pressure_string)
        )

        self.iso = self.calc._pressure_coupling
        self.logger.info("Pressure adjusted in %s" % self.iso)

        self.logger.info("Reference phase is %s" % self.calc.reference_phase)
        if self.calc.reference_phase == "liquid":
            if self.calc.melting_cycle:
                self.logger.info(
                    "Melting cycle will run, this can be turned off using the keyword melting_cycle"
                )
            else:
                self.logger.info("Melting cycle is turned off")

        # now thermostat and barostat damping process
        if self.calc.equilibration_control is None:
            self.logger.info(
                "Thermostat/Barostat combo for equilibration cycle is not explicitely specified"
            )
            self.logger.info(
                "Thermostat/Barostat combo for equilibration cycle can be specified using keyword equilibration_control"
            )
            self.calc.equilibration_control = "nose-hoover"
        else:
            self.logger.info(
                "Equilibration stage is done using %s barostat/thermostat"
                % self.calc.equilibration_control
            )

        if self.calc._qtb:
            qtb = self.calc.quantum_thermal_bath
            self.logger.info(
                "mode=fe-qtb: Dammak quantum thermal bath active for all MD stages"
            )
            self.logger.info("QTB thermostat damping is %f" % qtb.thermostat_damping)
            self.logger.info("QTB barostat damping is %f" % qtb.barostat_damping)
            self.logger.info("QTB f_max is %f THz, N_f is %d" % (qtb.f_max, qtb.n_f))
            self.calc.md.thermostat_damping = [
                qtb.thermostat_damping,
                qtb.thermostat_damping,
            ]
            self.calc.md.barostat_damping = [
                qtb.barostat_damping,
                qtb.barostat_damping,
            ]
            self.logger.info("These values can be tuned by adding in the input file:")
            self.logger.info("quantum_thermal_bath:")
            self.logger.info("   thermostat_damping: <float>  # τ in ps")
            self.logger.info("   barostat_damping: <float>")
            self.logger.info("   f_max: <float>               # THz cutoff")
            self.logger.info("   n_f: <int>                   # spectrum bins")
            if self.calc.equilibration_control == "berendsen":
                self.logger.warning(
                    "mode=fe-qtb overrides equilibration_control=berendsen; "
                    "QTB sampling is used in the equilibration stage too so the "
                    "established volume reflects quantum thermal expansion."
                )
        elif self.calc.equilibration_control == "nose-hoover":
            self.logger.info(
                "Nose-Hoover thermostat damping is %f"
                % self.calc.nose_hoover.thermostat_damping
            )
            self.calc.md.thermostat_damping = [
                self.calc.nose_hoover.thermostat_damping,
                self.calc.md.thermostat_damping,
            ]
            self.logger.info(
                "Nose-Hoover barostat damping is %f"
                % self.calc.nose_hoover.barostat_damping
            )
            self.calc.md.barostat_damping = [
                self.calc.nose_hoover.barostat_damping,
                self.calc.md.barostat_damping,
            ]
            self.logger.info("These values can be tuned by adding in the input file:")
            self.logger.info("nose_hoover:")
            self.logger.info("   thermostat_damping: <float>")
            self.logger.info("   barostat_damping: <float>")
            if self.calc.md.thermostat_damping[0] > 10.0:
                self.logger.warning("Equil. Nose-Hoover thermostat damping is high!")
            if self.calc.md.barostat_damping[0] > 10.0:
                self.logger.warning("Equil. Nose-Hoover barostat damping is high!")
        else:
            self.logger.info(
                "Berendsen thermostat damping is %f"
                % self.calc.berendsen.thermostat_damping
            )
            self.calc.md.thermostat_damping = [
                self.calc.berendsen.thermostat_damping,
                self.calc.md.thermostat_damping,
            ]
            self.logger.info(
                "Berendsen barostat damping is %f"
                % self.calc.berendsen.barostat_damping
            )
            self.calc.md.barostat_damping = [
                self.calc.berendsen.barostat_damping,
                self.calc.md.barostat_damping,
            ]
            self.logger.info("These values can be tuned by adding in the input file:")
            self.logger.info("berendsen:")
            self.logger.info("   thermostat_damping: <float>")
            self.logger.info("   barostat_damping: <float>")
            if self.calc.md.thermostat_damping[0] < 1.0:
                self.logger.warning("Equil. Berendsen thermostat damping is low!")
            if self.calc.md.barostat_damping[0] < 1.0:
                self.logger.warning("Equil. Berendsen barostat damping is high!")

        self.logger.info(
            "Integration stage is done using Nose-Hoover thermostat and barostat when needed"
        )
        self.logger.info(
            "Thermostat damping is %f" % (self.calc.md.thermostat_damping[1])
        )
        self.logger.info("Barostat damping is %f" % (self.calc.md.barostat_damping[1]))

        if self.calc._fix_lattice:
            self.logger.info(
                "Lattice is fixed, pressure convergence criteria is 50*tolerance.pressure; change if needed!"
            )

        self.l = self.calc.lattice
        self.alat = self.calc.lattice_constant
        self.vol = None

        # other properties
        self.cores = self.calc.queue.cores
        self.ncells = np.prod(self.calc.repeat)
        self.natoms = self.calc._natoms
        self.logger.info(
            "%d atoms in %d cells on %d cores" % (self.natoms, self.ncells, self.cores)
        )

        # reference system props; may not be always used
        # TODO : Add option to customize UFM parameters
        self.eps = self.calc._temperature * self.calc.uhlenbeck_ford_model.p * kb

        # Resolve the UFM length scale(s). `sigma` is either a scalar (original
        # single-component reference) or a dict of per-element-pair values for the
        # two-leg reference path. self._ufm_sigma_by_type maps (type_i, type_j) ->
        # sigma; None for the scalar case. self._is_two_leg flags the new path.
        self._ufm_sigma_by_type = self._resolve_ufm_sigmas()
        self._is_two_leg = (
            self.calc.uhlenbeck_ford_model.single_sigma is not None
        )

        if self._ufm_sigma_by_type is None:
            sigma_max = self.calc.uhlenbeck_ford_model.sigma
        else:
            sigma_max = max(self._ufm_sigma_by_type.values())
        # the multi-component (leg-1) UFM cutoff is set by the largest sigma
        self.ufm_cutoff = 5 * sigma_max

        # eps and cutoff for the single-component endpoint (leg 2)
        if self._is_two_leg:
            single_p = self.calc.uhlenbeck_ford_model.single_p
            if single_p is None:
                single_p = self.calc.uhlenbeck_ford_model.p
            self.single_eps = self.calc._temperature * single_p * kb
            self.single_sigma = self.calc.uhlenbeck_ford_model.single_sigma
            self.single_ufm_cutoff = 5 * self.single_sigma

        # properties that will be calculated later
        self.volatom = None
        self.k = None
        self.rho = None

        self.ferr = 0
        self.fref = 0
        self.feinstein = 0
        self.fcm = 0
        self.fideal = 0

        self.w = 0
        self.pv = 0
        self.fe = 0

        # box dimensions that need to be stored
        self.lx = None
        self.ly = None
        self.lz = None

        # now manually tune pair styles
        if self.calc.pair_style is not None:
            self.logger.info("pair_style: %s" % self.calc._pair_style_with_options[0])
            self.logger.info("pair_coeff: %s" % self.calc.pair_coeff[0])

            # log second pair style
            if len(self.calc.pair_style) > 1:
                self.logger.info(
                    "second pair_style: %s" % self.calc._pair_style_with_options[1]
                )
                self.logger.info("second pair_coeff: %s" % self.calc.pair_coeff[1])
        else:
            self.logger.info("pair_style or pair_coeff not provided")
        self._lmp = lmp

    def _resolve_ufm_sigmas(self):
        """
        Resolve the UFM length scale into a per-LAMMPS-type-pair mapping.

        Returns
        -------
        None
            if ``uhlenbeck_ford_model.sigma`` is a scalar (original
            single-component behaviour).
        dict
            mapping (type_i, type_j) with type_i <= type_j (1-indexed LAMMPS
            types) -> sigma, when ``sigma`` is given as a dict of element-pair
            length scales (two-leg path). Keys in the input are
            "<elementA>_<elementB>" and are order-insensitive. Cross terms that
            are not supplied are left out and filled by LAMMPS geometric mixing.
        """
        sigma = self.calc.uhlenbeck_ford_model.sigma
        if not isinstance(sigma, dict):
            return None

        # element symbol -> 1-indexed LAMMPS type (order in `element` list)
        elem_to_type = {el: i + 1 for i, el in enumerate(self.calc.element)}

        resolved = {}
        for key, val in sigma.items():
            parts = key.split("_")
            if len(parts) != 2:
                raise ValueError(
                    "UFM sigma key '%s' must be of the form 'ElementA_ElementB'"
                    % key
                )
            ea, eb = parts
            if ea not in elem_to_type or eb not in elem_to_type:
                raise ValueError(
                    "UFM sigma key '%s' refers to element(s) not in %s"
                    % (key, self.calc.element)
                )
            ti, tj = elem_to_type[ea], elem_to_type[eb]
            if ti > tj:
                ti, tj = tj, ti
            resolved[(ti, tj)] = float(val)
        return resolved

    def ufm_pair_coeff_commands(self, eps, sigma_scalar, sigma_by_type, substyle=""):
        """
        Build the ``pair_coeff`` command(s) for a UFM interaction.

        Parameters
        ----------
        eps : float
            UFM energy scale.
        sigma_scalar : float
            length scale to use when ``sigma_by_type`` is None (single-component).
        sigma_by_type : dict or None
            mapping (type_i, type_j) -> sigma for the multi-component case.
        substyle : str
            hybrid/scaled substyle token to insert after the type pair, e.g.
            "ufm", "ufm 1", "ufm 2". Empty string for a plain (non-hybrid)
            ``pair_style ufm`` where no style keyword is allowed.

        Returns
        -------
        list of str
        """
        tok = (" " + substyle) if substyle else ""
        if sigma_by_type is None:
            return ["pair_coeff       * *%s %f %f" % (tok, eps, sigma_scalar)]
        # Multi-component: LAMMPS requires EVERY declared type pair to be set, even
        # for types not present in the structure. Emit a base "* *" coeff using a
        # default sigma (the largest provided value) so all pairs are defined, then
        # override the explicitly-specified pairs. Later pair_coeff lines override
        # earlier ones in LAMMPS, so order matters: base first, specifics after.
        default_sigma = max(sigma_by_type.values())
        cmds = ["pair_coeff       * *%s %f %f" % (tok, eps, default_sigma)]
        for (ti, tj), sig in sorted(sigma_by_type.items()):
            cmds.append("pair_coeff       %d %d%s %f %f" % (ti, tj, tok, eps, sig))
        return cmds

    def __repr__(self):
        """
        String of the class
        """
        data = self.calc.__repr__()
        return data

    def _from_dict(self, org_dict, indict):
        for key, val in indict.items():
            if isinstance(val, dict):
                if key not in org_dict.keys():
                    org_dict[key] = {}
                self._from_dict(org_dict[key], val)
            else:
                org_dict[key] = val

    def dump_current_snapshot(self, lmp, filename):
        """ """
        lmp.command(
            "dump              2 all custom 1 %s id type mass x y z vx vy vz"
            % (filename)
        )
        lmp.command("run               0")
        lmp.command("undump            2")

    def get_structures(self, stage="fe", direction="forward", n_iteration=1):
        """ """
        species = self.calc.element
        filename = None

        if stage == "fe":
            filename = os.path.join(self.simfolder, "conf.equilibration.data")
        elif stage == "ts":
            if direction == "forward":
                filename = os.path.join(
                    self.simfolder, "traj.ts.forward_%d.dat" % n_iteration
                )
            elif direction == "backward":
                filename = os.path.join(
                    self.simfolder, "traj.ts.backward_%d.dat" % n_iteration
                )

        structures = None
        if filename is not None:
            structures = ph.get_structures(filename, species, index=None)

        return structures

    def check_if_melted(self, lmp, filename):
        """
        Check whether the solid has melted, using a structural solid-fraction
        criterion.

        The solid fraction is computed from the trajectory snapshot in
        ``filename``; if it drops below ``tolerance.solid_fraction`` the run is
        aborted with a MeltedError.  Detection can be turned off by setting
        ``tolerance.solid_fraction: 0``.
        """
        solids = ph.find_solid_fraction(os.path.join(self.simfolder, filename))
        if solids / self.natoms < self.calc.tolerance.solid_fraction:
            self.lammps_close(lmp=lmp)
            # Preserve log file on error
            logfile = os.path.join(self.simfolder, "log.lammps")
            try:
                os.rename(
                    logfile, os.path.join(self.simfolder, "melted_error.log.lammps")
                )
            except OSError as e:
                self.logger.warning(f"Failed to rename log file: {e}")
            raise MeltedError(
                "System melted, increase size or reduce temp!\n Solid detection algorithm only works with BCC/FCC/HCP/SC/DIA. Detection algorithm can be turned off by setting:\n tolerance.solid_fraction: 0"
            )

    def check_if_solidfied(self, lmp, filename):
        """
        Check whether the liquid has solidified, using a structural
        solid-fraction criterion.

        If the solid fraction in ``filename`` exceeds
        ``tolerance.liquid_fraction`` the run is aborted with a SolidifiedError.
        """
        solids = ph.find_solid_fraction(os.path.join(self.simfolder, filename))
        if solids / self.natoms > self.calc.tolerance.liquid_fraction:
            self.lammps_close(lmp=lmp)
            # Preserve log file on error
            logfile = os.path.join(self.simfolder, "log.lammps")
            try:
                os.rename(
                    logfile, os.path.join(self.simfolder, "solidified_error.log.lammps")
                )
            except OSError as e:
                self.logger.warning(f"Failed to rename log file: {e}")
            raise SolidifiedError("System solidified, increase temperature")

    def fix_nose_hoover(
        self,
        lmp,
        temp_start_factor=1.0,
        temp_end_factor=1.0,
        press_start_factor=1.0,
        press_end_factor=1.0,
        stage=0,
        ensemble="npt",
    ):
        """
        Fix Nose-Hoover thermostat and barostat

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        lmp.command(
            "fix              nh1 all npt temp %f %f %f %s %f %f %f"
            % (
                temp_start_factor * self.calc._temperature,
                temp_end_factor * self.calc._temperature,
                self.calc.md.thermostat_damping[stage],
                self.iso,
                press_start_factor * self.calc._pressure,
                press_end_factor * self.calc._pressure,
                self.calc.md.barostat_damping[stage],
            )
        )

    def fix_berendsen(
        self,
        lmp,
        temp_start_factor=1.0,
        temp_end_factor=1.0,
        press_start_factor=1.0,
        press_end_factor=1.0,
        stage=0,
        ensemble="npt",
    ):
        """
        Fix Nose-Hoover thermostat and barostat

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        lmp.command("fix              b1a all nve")
        lmp.command(
            "fix              b1b all temp/berendsen %f %f %f"
            % (
                temp_start_factor * self.calc._temperature,
                temp_end_factor * self.calc._temperature,
                self.calc.md.thermostat_damping[stage],
            )
        )
        lmp.command(
            "fix              b1c all press/berendsen %s %f %f %f"
            % (
                self.iso,
                press_start_factor * self.calc._pressure,
                press_end_factor * self.calc._pressure,
                self.calc.md.barostat_damping[stage],
            )
        )

    def fix_qtb(
        self,
        lmp,
        temp_start_factor=1.0,
        temp_end_factor=1.0,
        press_start_factor=1.0,
        press_end_factor=1.0,
        stage=0,
        ensemble="npt",
    ):
        """
        Apply the Dammak quantum thermal bath as a thermostat.

        QTB only thermostats — it must be paired with an integrator. We use
        fix nph for NPT (pressure-controlled) and fix nve for NVT (canonical).
        Temperature endpoints are read but only the start value is used: LAMMPS
        fix qtb takes a single temperature, so any ramping must be implemented
        externally.
        """
        qtb = self.calc.quantum_thermal_bath
        t_qtb = temp_start_factor * self.calc._temperature
        if abs(temp_end_factor - temp_start_factor) > 1e-12:
            self.logger.warning(
                "fix qtb does not support ramped temperature; using start value %f K"
                % t_qtb
            )
        if ensemble == "npt":
            lmp.command(
                "fix              nh1 all nph %s %f %f %f"
                % (
                    self.iso,
                    press_start_factor * self.calc._pressure,
                    press_end_factor * self.calc._pressure,
                    qtb.barostat_damping,
                )
            )
        else:
            lmp.command("fix              nh1 all nve")
        # qtb damping is in time units (ps in metal units)
        lmp.command(
            "fix              nh1_qtb all qtb temp %f damp %f seed %d f_max %f N_f %d"
            % (
                t_qtb,
                qtb.thermostat_damping,
                np.random.randint(1, 10**8),
                qtb.f_max,
                qtb.n_f,
            )
        )

    def unfix_qtb(self, lmp):
        lmp.command("unfix            nh1_qtb")
        lmp.command("unfix            nh1")

    def unfix_nose_hoover(self, lmp):
        """
        Fix Nose-Hoover thermostat and barostat

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        lmp.command("unfix            nh1")

    def unfix_berendsen(self, lmp):
        """
        Fix Nose-Hoover thermostat and barostat

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        lmp.command("unfix            b1a")
        lmp.command("unfix            b1b")
        lmp.command("unfix            b1c")

    def run_zero_pressure_equilibration(self, lmp):
        """
        Run a zero pressure equilibration

        Parameters
        ----------
        lmp: LAMMPS object

        Returns
        -------
        None

        Notes
        -----
        Each method should close all the fixes. Run a small eqbr routine to achieve zero pressure
        """
        # set velocity
        lmp.command(
            "velocity         all create %f %d"
            % (self.calc._temperature, np.random.randint(1, 10000))
        )

        # apply fixes depending on thermostat/barostat
        if self.calc._qtb:
            self.fix_qtb(lmp, ensemble="npt")
        elif self.calc.equilibration_control == "nose-hoover":
            self.fix_nose_hoover(lmp)
        else:
            self.fix_berendsen(lmp)
        # start thermo logging
        lmp.command("thermo_style     custom step pe press vol etotal temp lx ly lz")
        lmp.command("thermo           10")

        # run MD
        lmp.command("run              %d" % int(self.calc.md.n_small_steps))

        # remove fixes
        if self.calc._qtb:
            self.unfix_qtb(lmp)
        elif self.calc.equilibration_control == "nose-hoover":
            self.unfix_nose_hoover(lmp)
        else:
            self.unfix_berendsen(lmp)

    def run_finite_pressure_equilibration(self, lmp):
        """
        Run a finite pressure equilibration

        Parameters
        ----------
        lmp: LAMMPS object

        Returns
        -------
        None

        Notes
        -----
        Each method should close all the fixes. Run a equilibration routine to reach the given finite pressure.
        The pressure is implemented in one fix, while temperature is gradually ramped.
        The thermostat can work faster than barostat, which means that the structure will melt before the pressure is scaled, this ramping
        can prevent the issue.
        """
        # create velocity
        lmp.command(
            "velocity         all create %f %d"
            % (0.25 * self.calc._temperature, np.random.randint(1, 10000))
        )

        # for QTB / Nose-Hoover thermo/baro combination
        if self.calc._qtb:
            # QTB cannot ramp T; do three cycles at staircase T to mimic the warm-up
            self.fix_qtb(lmp, temp_start_factor=0.5, temp_end_factor=0.5, ensemble="npt")
            lmp.command("thermo_style     custom step pe press vol etotal temp")
            lmp.command("thermo           10")
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))
            self.unfix_qtb(lmp)

            self.fix_qtb(lmp, temp_start_factor=1.0, temp_end_factor=1.0, ensemble="npt")
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))
            self.unfix_qtb(lmp)

            self.fix_qtb(lmp, ensemble="npt")
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))
            self.unfix_qtb(lmp)

        elif self.calc.equilibration_control == "nose-hoover":
            # Cycle 1: 0.25-0.5 temperature, full pressure
            self.fix_nose_hoover(lmp, temp_start_factor=0.25, temp_end_factor=0.5)
            lmp.command("thermo_style     custom step pe press vol etotal temp")
            lmp.command("thermo           10")
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))
            self.unfix_nose_hoover(lmp)

            # Cycle 2: 0.5-1.0 temperature, full pressure
            self.fix_nose_hoover(lmp, temp_start_factor=0.5, temp_end_factor=1.0)
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))
            self.unfix_nose_hoover(lmp)

            # Cycle 3: full temperature, full pressure
            self.fix_nose_hoover(lmp)
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))
            self.unfix_nose_hoover(lmp)

        else:
            # Cycle 1: 0.25-0.5 temperature, full pressure
            self.fix_berendsen(lmp, temp_start_factor=0.25, temp_end_factor=0.5)
            lmp.command("thermo_style     custom step pe press vol etotal temp")
            lmp.command("thermo           10")
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))
            self.unfix_berendsen(lmp)

            # Cycle 2: 0.5-1.0 temperature, full pressure
            self.fix_berendsen(lmp, temp_start_factor=0.5, temp_end_factor=1.0)
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))
            self.unfix_berendsen(lmp)

            # Cycle 3: full temperature, full pressure
            self.fix_berendsen(lmp)
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))
            self.unfix_berendsen(lmp)

    def run_pressure_convergence(self, lmp):
        """
        Run a pressure convergence routine

        Parameters
        ----------
        lmp: LAMMPS object

        Returns
        -------
        None

        Notes
        -----
        Take the equilibrated structure and rigorously check for pressure convergence.
        The cycle is stopped when the average pressure is within the given cutoff of the target pressure.
        """
        if self.calc.script_mode:
            self.run_minimal_pressure_convergence(lmp)
        else:
            self.run_iterative_pressure_convergence(lmp)

    def run_iterative_pressure_convergence(self, lmp):
        """
        Run a pressure convergence routine

        Parameters
        ----------
        lmp: LAMMPS object

        Returns
        -------
        None

        Notes
        -----
        Take the equilibrated structure and rigorously check for pressure convergence.

        Full-length NPT cycles are run throughout.  The first cycle is excluded
        from the running average to discard the initial transient.  After
        ``n_fit_warmup`` cycles (default 5), a linear P-V fit is used to
        predict the equilibrium volume and rescale the box, accelerating
        convergence.  The mean pressure (over all post-transient data) must
        fall within ``tolerance.pressure`` of the target to declare convergence.
        """

        # apply fixes
        if self.calc._qtb:
            self.fix_qtb(lmp, ensemble="npt")
        elif self.calc.equilibration_control == "nose-hoover":
            self.fix_nose_hoover(lmp)
        else:
            self.fix_berendsen(lmp)

        ave_every = int(self.calc.md.n_every_steps)
        ave_repeat = int(self.calc.md.n_repeat_steps)
        ave_freq = ave_every * ave_repeat

        lmp.command(
            "fix              2 all ave/time %d %d %d v_mlx v_mly v_mlz v_mpress v_mpe v_metotal v_mtemp file avg.dat"
            % (ave_every, ave_repeat, ave_freq)
        )

        ncount = int(self.calc.md.n_small_steps) // ave_freq
        target_pressure = self.calc._pressure
        converged = False
        pv_history = []  # [(vol_per_atom, mean_pressure), ...]
        n_fit_warmup = 5  # run this many cycles before attempting linear P-V fit
        n_skip = 1  # drop first cycle(s) from averaging (transient)

        for i in range(int(self.calc.md.n_cycles)):
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))

            file = os.path.join(self.simfolder, "avg.dat")

            lx, ly, lz, ipress = np.loadtxt(file, usecols=(1, 2, 3, 4), unpack=True)
            # Average over all data after dropping the first cycle
            skip_samples = n_skip * ncount
            if len(ipress) <= skip_samples:
                # not enough data yet (still in the skipped transient)
                lxpc = ipress
                lx_avg = lx
                ly_avg = ly
                lz_avg = lz
            else:
                lxpc = ipress[skip_samples:]
                lx_avg = lx[skip_samples:]
                ly_avg = ly[skip_samples:]
                lz_avg = lz[skip_samples:]
            mean = np.mean(lxpc)
            std = np.std(lxpc)
            volatom = np.mean((lx_avg * ly_avg * lz_avg) / self.natoms)
            self.logger.info(
                "At count %d mean pressure is %.2f bar, std %.2f, vol/atom %.4f"
                % (i + 1, mean, std, volatom)
            )

            pv_history.append((volatom, mean))

            if (np.abs(mean - target_pressure)) < self.calc.tolerance.pressure:
                self.logger.info("Pressure within tolerance")
                self.lx = np.round(np.mean(lx_avg), decimals=3)
                self.ly = np.round(np.mean(ly_avg), decimals=3)
                self.lz = np.round(np.mean(lz_avg), decimals=3)
                self.volatom = volatom
                self.vol = self.lx * self.ly * self.lz
                self.rho = self.natoms / (self.lx * self.ly * self.lz)

                self.logger.info(
                    "finalized vol/atom %f at pressure %f" % (self.volatom, mean)
                )
                self.logger.info(
                    "Avg box dimensions x: %f, y: %f, z:%f"
                    % (self.lx, self.ly, self.lz)
                )
                converged = True
                break
            else:
                # After enough warmup cycles, fit P(V) linearly and rescale box
                if len(pv_history) >= n_fit_warmup:
                    scale = self._fit_volume_scale(pv_history, target_pressure)
                    if scale is not None:
                        self.logger.info(
                            "Applying linear P-V fit correction — scale factor %.6f"
                            % scale
                        )
                        lmp.command(
                            "change_box       all x scale %f y scale %f z scale %f remap"
                            % (scale, scale, scale)
                        )

        if not converged:
            self.lammps_close(lmp=lmp)
            # Preserve log file on error
            logfile = os.path.join(self.simfolder, "log.lammps")
            try:
                os.rename(
                    logfile,
                    os.path.join(
                        self.simfolder, "pressure_convergence_error.log.lammps"
                    ),
                )
            except OSError as e:
                self.logger.warning(f"Failed to rename log file: {e}")
            raise ValueError(
                "Pressure did not converge after MD runs, maybe change lattice_constant and try?"
            )

        # unfix thermostat and barostat
        if self.calc._qtb:
            self.unfix_qtb(lmp)
        elif self.calc.equilibration_control == "nose-hoover":
            self.unfix_nose_hoover(lmp)
        else:
            self.unfix_berendsen(lmp)

        lmp.command("unfix            2")

    @staticmethod
    def _fit_volume_scale(pv_history, target_pressure):
        """
        Estimate a box scale factor by fitting a linear P(V) model to the
        accumulated volume/pressure history and predicting the volume at
        ``target_pressure``.

        Parameters
        ----------
        pv_history : list of (vol_per_atom, pressure) tuples
        target_pressure : float
            Target pressure in bar.

        Returns
        -------
        scale : float or None
            Isotropic scale factor to apply to each box dimension, or None if
            the fit is degenerate.  Clamped to [0.96, 1.04] to avoid large jumps.
        """
        vols = np.array([v for v, p in pv_history])
        pres = np.array([p for v, p in pv_history])

        # guard against NaN from early cycles
        mask = np.isfinite(vols) & np.isfinite(pres)
        vols, pres = vols[mask], pres[mask]

        # need at least 2 distinct volumes for a meaningful fit
        if len(vols) < 2 or np.ptp(vols) < 1e-12:
            return None

        # linear fit: P = slope * V + intercept
        slope, intercept = np.polyfit(vols, pres, 1)

        # slope must be negative (pressure decreases as volume increases)
        if slope >= 0:
            return None

        V_target = (target_pressure - intercept) / slope
        V_curr = vols[-1]

        if V_target <= 0:
            return None

        scale = (V_target / V_curr) ** (1.0 / 3.0)
        # clamp to prevent destructive jumps
        scale = max(0.96, min(1.04, scale))
        return scale

    def run_minimal_pressure_convergence(self, lmp):
        """
        Run a pressure convergence routine

        Parameters
        ----------
        lmp: LAMMPS object

        Returns
        -------
        None

        Notes
        -----
        Take the equilibrated structure and rigorously check for pressure convergence.
        The cycle is stopped when the average pressure is within the given cutoff of the target pressure.
        """

        # apply fixes
        if self.calc._qtb:
            self.fix_qtb(lmp, ensemble="npt")
        elif self.calc.equilibration_control == "nose-hoover":
            self.fix_nose_hoover(lmp)
        else:
            self.fix_berendsen(lmp)

        lmp.command(
            "fix              2 all ave/time %d %d %d v_mlx v_mly v_mlz v_mpress v_mpe v_metotal v_mtemp file avg.dat"
            % (
                int(self.calc.md.n_every_steps),
                int(self.calc.md.n_repeat_steps),
                int(self.calc.md.n_every_steps * self.calc.md.n_repeat_steps),
            )
        )

        lmp.command("run              %d" % int(self.calc.md.n_small_steps))
        lmp.command("run              %d" % int(self.calc.n_equilibration_steps))

        # unfix thermostat and barostat
        if self.calc._qtb:
            self.unfix_qtb(lmp)
        elif self.calc.equilibration_control == "nose-hoover":
            self.unfix_nose_hoover(lmp)
        else:
            self.unfix_berendsen(lmp)

        lmp.command("unfix            2")

    def run_constrained_pressure_convergence(self, lmp):
        if self.calc.script_mode:
            self.run_minimal_constrained_pressure_convergence(lmp)
        else:
            self.run_iterative_constrained_pressure_convergence(lmp)

    def run_iterative_constrained_pressure_convergence(self, lmp):
        """ """
        lmp.command(
            "velocity         all create %f %d"
            % (self.calc._temperature, np.random.randint(1, 10000))
        )
        if self.calc._qtb:
            qtb = self.calc.quantum_thermal_bath
            lmp.command("fix              1 all nve")
            lmp.command(
                "fix              1q all qtb temp %f damp %f seed %d f_max %f N_f %d"
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
                "fix              1 all nvt temp %f %f %f"
                % (
                    self.calc._temperature,
                    self.calc._temperature,
                    self.calc.md.thermostat_damping[1],
                )
            )
        lmp.command("thermo_style     custom step pe press vol etotal temp lx ly lz")
        lmp.command("thermo           10")

        # this is when the averaging routine starts
        lmp.command(
            "fix              2 all ave/time %d %d %d v_mlx v_mly v_mlz v_mpress v_mpe v_metotal v_mtemp file avg.dat"
            % (
                int(self.calc.md.n_every_steps),
                int(self.calc.md.n_repeat_steps),
                int(self.calc.md.n_every_steps * self.calc.md.n_repeat_steps),
            )
        )

        lastmean = 100000000
        converged = False
        for i in range(int(self.calc.md.n_cycles)):
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))

            # now we can check if it converted
            mean, std, volatom = self.process_pressure()
            self.logger.info(
                "At count %d mean pressure is %f with %f vol/atom"
                % (i + 1, mean, volatom)
            )

            if (np.abs(mean - lastmean)) < 50 * self.calc.tolerance.pressure:
                # here we actually have to set the pressure
                self.finalise_pressure()
                converged = True
                break

            lastmean = mean

        if self.calc._qtb:
            lmp.command("unfix            1q")
        lmp.command("unfix            1")
        lmp.command("unfix            2")

        if not converged:
            self.lammps_close(lmp=lmp)
            # Preserve log file on error
            logfile = os.path.join(self.simfolder, "log.lammps")
            try:
                os.rename(
                    logfile,
                    os.path.join(
                        self.simfolder, "constrained_pressure_error.log.lammps"
                    ),
                )
            except OSError as e:
                self.logger.warning(f"Failed to rename log file: {e}")
            raise ValueError("pressure did not converge")

    def process_pressure(
        self,
    ):
        if self.calc.script_mode:
            ncount = int(self.calc.n_equilibration_steps) // int(
                self.calc.md.n_every_steps * self.calc.md.n_repeat_steps
            )
        else:
            ncount = int(self.calc.md.n_small_steps) // int(
                self.calc.md.n_every_steps * self.calc.md.n_repeat_steps
            )

        # now we can check if it converted
        file = os.path.join(self.simfolder, "avg.dat")
        # we have to clean the data, so as just the last block is selected
        lx, ly, lz, lxpc = np.loadtxt(file, usecols=(1, 2, 3, 4), unpack=True)
        lx = lx[-ncount + 1 :]
        ly = ly[-ncount + 1 :]
        lz = lz[-ncount + 1 :]
        lxpc = lxpc[-ncount + 1 :]
        mean = np.mean(lxpc)
        std = np.std(lxpc)
        volatom = np.mean((lx * ly * lz) / self.natoms)
        return mean, std, volatom

    def finalise_pressure(
        self,
    ):
        if self.calc.script_mode:
            ncount = int(self.calc.n_equilibration_steps) // int(
                self.calc.md.n_every_steps * self.calc.md.n_repeat_steps
            )
        else:
            ncount = int(self.calc.md.n_small_steps) // int(
                self.calc.md.n_every_steps * self.calc.md.n_repeat_steps
            )

        file = os.path.join(self.simfolder, "avg.dat")
        lx, ly, lz, lxpc = np.loadtxt(file, usecols=(1, 2, 3, 4), unpack=True)
        lx = lx[-ncount + 1 :]
        ly = ly[-ncount + 1 :]
        lz = lz[-ncount + 1 :]
        lxpc = lxpc[-ncount + 1 :]

        mean = np.mean(lxpc)
        std = np.std(lxpc)
        volatom = np.mean((lx * ly * lz) / self.natoms)

        self.calc._pressure = mean
        self.lx = np.round(np.mean(lx), decimals=3)
        self.ly = np.round(np.mean(ly), decimals=3)
        self.lz = np.round(np.mean(lz), decimals=3)
        self.volatom = volatom
        self.vol = self.lx * self.ly * self.lz
        self.rho = self.natoms / (self.lx * self.ly * self.lz)
        self.logger.info("finalized vol/atom %f at pressure %f" % (self.volatom, mean))
        self.logger.info(
            "Avg box dimensions x: %f, y: %f, z:%f" % (self.lx, self.ly, self.lz)
        )

    def run_minimal_constrained_pressure_convergence(self, lmp):
        """ """
        lmp.command(
            "velocity         all create %f %d"
            % (self.calc._temperature, np.random.randint(1, 10000))
        )
        if self.calc._qtb:
            qtb = self.calc.quantum_thermal_bath
            lmp.command("fix              1 all nve")
            lmp.command(
                "fix              1q all qtb temp %f damp %f seed %d f_max %f N_f %d"
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
                "fix              1 all nvt temp %f %f %f"
                % (
                    self.calc._temperature,
                    self.calc._temperature,
                    self.calc.md.thermostat_damping[1],
                )
            )
        lmp.command("thermo_style     custom step pe press vol etotal temp lx ly lz")
        lmp.command("thermo           10")

        # this is when the averaging routine starts
        lmp.command(
            "fix              2 all ave/time %d %d %d v_mlx v_mly v_mlz v_mpress v_mpe v_metotal v_mtemp file avg.dat"
            % (
                int(self.calc.md.n_every_steps),
                int(self.calc.md.n_repeat_steps),
                int(self.calc.md.n_every_steps * self.calc.md.n_repeat_steps),
            )
        )

        lmp.command("run              %d" % int(self.calc.md.n_small_steps))
        lmp.command("run              %d" % int(self.calc.n_equilibration_steps))

        if self.calc._qtb:
            lmp.command("unfix            1q")
        lmp.command("unfix            1")
        lmp.command("unfix            2")

    def submit_report(self, extra_dict=None):
        """
        Submit final report containing results

        Parameters
        ----------
        extra_dict: dict
            extra information to be written out

        Returns
        -------
        None
        """
        report = {}

        # input quantities
        report["input"] = {}
        report["input"]["temperature"] = int(self.calc._temperature)
        report["input"]["pressure"] = float(self.calc._pressure)
        report["input"]["lattice"] = str(self.calc._original_lattice)
        report["input"]["element"] = " ".join(np.array(self.calc.element).astype(str))
        report["input"]["concentration"] = " ".join(
            np.array(
                [val["composition"] for key, val in self.calc._element_dict.items()]
            ).astype(str)
        )

        # average quantities
        report["average"] = {}
        report["average"]["vol_atom"] = float(self.volatom)

        if self.k is not None:
            report["average"]["spring_constant"] = " ".join(
                np.array(self.k).astype(str)
            )
        if self.rho is not None:
            report["average"]["density"] = float(self.rho)

        # results
        report["results"] = {}
        report["results"]["free_energy"] = float(self.fe)
        report["results"]["error"] = float(self.ferr)
        report["results"]["reference_system"] = float(self.fref)
        report["results"]["einstein_crystal"] = float(self.feinstein)
        report["results"]["com_correction"] = float(self.fcm)
        report["results"]["work"] = float(self.w)
        report["results"]["pv"] = float(self.pv)
        report["results"]["unit"] = "eV/atom"

        if extra_dict is not None:
            self._from_dict(report, extra_dict)

        self.report = report

        reportfile = os.path.join(self.simfolder, "report.yaml")
        with open(reportfile, "w") as f:
            yaml.dump(report, f)

        self.logger.info("Report written in %s" % reportfile)

        # now we have to write out the results
        self.logger.info("Please cite the following publications:")
        self.logger.info("- 10.1103/PhysRevMaterials.5.103801")
        self.publications.append("10.1103/PhysRevMaterials.5.103801")

        if self.calc.mode == "fe":
            if self.calc.reference_phase == "solid":
                self.logger.info("- 10.1016/j.commatsci.2015.10.050")
                self.publications.append("10.1016/j.commatsci.2015.10.050")
            else:
                self.logger.info("- 10.1016/j.commatsci.2018.12.029")
                self.logger.info("- 10.1063/1.4967775")
                self.publications.append("10.1016/j.commatsci.2018.12.029")
                self.publications.append("10.1063/1.4967775")

    # ------------------------------------------------------------------
    # Internal helpers for temperature-window block sweeps
    # ------------------------------------------------------------------

    def _run_sweep(
        self,
        lmp,
        lambda_var: str,
        output_file_pattern: str,
        sweep_label: str,
    ) -> None:
        """
        Run a forward or backward sweep as a single continuous LAMMPS run.

        Records ``dU press vol lambda`` at every step to ``output_file_pattern``.

        Parameters
        ----------
        lmp : lammps object
            Active LAMMPS instance.  The pair style and the lambda ramp
            variables (and the integrator fix) must already be defined before
            this method is called.
        lambda_var : str
            Name of the LAMMPS variable to record, e.g. ``"flambda"`` or
            ``"blambda"``.
        output_file_pattern : str
            Name for the output data file, e.g. ``"ts.forward_1.dat"``.
        sweep_label : str
            Human-readable label used in log messages.
        """
        n_sweep = self.calc._n_sweep_steps
        lmp.command(
            'fix               f3 all print 1 "${dU} $(press) $(vol) ${%s}" '
            'screen no file %s' % (lambda_var, output_file_pattern)
        )
        self.logger.info("ts-sweep %s: %d steps", sweep_label, n_sweep)
        lmp.command("run               %d" % n_sweep)
        lmp.command("unfix             f3")

    def _reversible_scaling_forward(self, iteration: int = 1) -> None:
        """
        Perform the forward sweep of a reversible-scaling calculation.

        1. Initial NPT equilibration at T0.
        2. COM-constrained equilibration at T0.
        3. Forward sweep: λ 1 → T0/Tf.
        4. Write ``conf.ts.forward_{iteration}.data`` for the backward sweep.

        Parameters
        ----------
        iteration : int
            Reversible-scaling iteration index.
        """
        solid = self.calc.reference_phase == "solid"
        t0 = self.calc._temperature
        tf = self.calc._temperature_stop
        li = 1.0
        lf = t0 / tf
        pi = self.calc._pressure
        pf = lf * pi

        self.logger.info(
            "forward sweep (iteration %d): T %.1f → %.1f K, "
            "λ %.4f → %.4f, P %.4f → %.4f bar",
            iteration, t0, tf, li, lf, pi, pf,
        )

        lmp = ph.create_object(
            cores=self.cores,
            directory=self.simfolder,
            timestep=self.calc.md.timestep,
            cmdargs=self.calc.md.cmdargs,
            init_commands=self.calc.md.init_commands,
            script_mode=self.calc.script_mode,
            lmp=self._lmp,
        )

        lmp.command("echo              log")
        lmp.command("variable          li equal %f" % li)
        lmp.command("variable          lf equal %f" % lf)

        lmp = ph.set_pair_style(lmp, self.calc)

        conf = os.path.join(self.simfolder, "conf.equilibration.data")
        lmp = ph.read_data(lmp, conf)

        lmp = ph.set_pair_coeff(lmp, self.calc)
        lmp = ph.set_mass(lmp, self.calc)

        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        # ── Initial equilibration ──────────────────────────────────────────
        if self.calc.npt:
            lmp.command(
                "fix               f1 all npt temp %f %f %f %s %f %f %f"
                % (t0, t0, self.calc.md.thermostat_damping[1],
                   self.iso, pi, pi, self.calc.md.barostat_damping[1])
            )
        else:
            lmp.command(
                "fix               f1 all nvt temp %f %f %f"
                % (t0, t0, self.calc.md.thermostat_damping[1])
            )

        self.logger.info("forward sweep (iteration %d): initial equilibration start", iteration)
        lmp.command("run               %d" % self.calc.n_equilibration_steps)
        self.logger.info("forward sweep (iteration %d): initial equilibration done", iteration)

        lmp.command("unfix             f1")

        # ── COM-constrained equilibration ──────────────────────────────────
        lmp.command("variable         xcm equal xcm(all,x)")
        lmp.command("variable         ycm equal xcm(all,y)")
        lmp.command("variable         zcm equal xcm(all,z)")

        if self.calc.npt:
            lmp.command(
                "fix               f1 all npt temp %f %f %f %s %f %f %f "
                "fixedpoint ${xcm} ${ycm} ${zcm}"
                % (t0, t0, self.calc.md.thermostat_damping[1],
                   self.iso, pi, pi, self.calc.md.barostat_damping[1])
            )
        else:
            lmp.command(
                "fix               f1 all nvt temp %f %f %f "
                "fixedpoint ${xcm} ${ycm} ${zcm}"
                % (t0, t0, self.calc.md.thermostat_damping[1])
            )

        lmp.command("compute           tcm all temp/com")
        lmp.command("fix_modify        f1 temp tcm")
        lmp.command("variable          step    equal step")
        lmp.command("variable          dU      equal c_thermo_pe/atoms")
        lmp.command("thermo_style      custom step pe c_tcm press vol")
        lmp.command("thermo            10000")

        lmp.command(
            "velocity          all create %f %d mom yes rot yes dist gaussian"
            % (t0, np.random.randint(1, 10000))
        )
        self.logger.info(
            "forward sweep (iteration %d): COM-constrained equilibration start", iteration
        )
        lmp.command("run               %d" % self.calc.n_equilibration_steps)
        self.logger.info(
            "forward sweep (iteration %d): COM-constrained equilibration done", iteration
        )

        # ----------------------------------------------------------------
        # Lambda schedule for the forward sweep.
        #
        # "linear" (default): lambda = ramp(li, lf) — simple linear
        #   interpolation; LAMMPS ramp() resets automatically each run.
        #
        # "uniform_temperature": T_eq(s) = T0/lambda is linear in step
        #   so every Kelvin bin gets the same number of MD samples.
        #   Requires explicit step0 capture before each sweep.
        # ----------------------------------------------------------------
        lmp.command("variable         T0_rs equal %f" % t0)
        if self.calc.lambda_schedule == "uniform_temperature":
            lmp.command("variable         Nsweep equal %d" % self.calc._n_sweep_steps)
            lmp.command("variable         Tf_rs equal %f" % tf)
            # Capture the step at the START of the sweep so the formula is
            # independent of any prior MD steps (no reset_timestep needed).
            lmp.command("variable         step0 equal $(step)")
            lmp.command(
                "variable         flambda equal "
                "v_T0_rs/(v_T0_rs+(v_Tf_rs-v_T0_rs)*(step-v_step0)/v_Nsweep)"
            )
            lmp.command(
                "variable         blambda equal "
                "v_T0_rs/(v_Tf_rs-(v_Tf_rs-v_T0_rs)*(step-v_step0)/v_Nsweep)"
            )
        else:  # "linear" (default)
            lmp.command("variable         flambda equal ramp(${li},${lf})")
            lmp.command("variable         blambda equal ramp(${lf},${li})")
        lmp.command("variable         fscale equal v_flambda-1.0")
        lmp.command("variable         bscale equal v_blambda-1.0")
        lmp.command("variable         one equal 1.0")
        lmp.command("variable         ftemp equal v_T0_rs/v_flambda")
        lmp.command("variable         btemp equal v_T0_rs/v_blambda")

        lmp.command(ph.scaled_pair_style_command(self.calc, ["v_one", "v_fscale"]))
        for cmd in ph.hybrid_pair_coeff_commands(self.calc, repeat_index=0, total_repeats=2):
            lmp.command(cmd)
        for cmd in ph.hybrid_pair_coeff_commands(self.calc, repeat_index=1, total_repeats=2):
            lmp.command(cmd)

        # ── Optional MC swaps ───────────────────────────────────────────────
        if (
            self.calc.monte_carlo.n_swaps > 0
            and len(self.calc.monte_carlo.forward_swap_types) >= 2
        ):
            swap_types  = self.calc.monte_carlo.forward_swap_types
            swap_combos = list(itertools.combinations(swap_types, 2))
            self.logger.info(
                "forward sweep (iteration %d): %d swap moves/combo, "
                "%d combinations every %d steps",
                iteration, self.calc.monte_carlo.n_swaps,
                len(swap_combos), self.calc.monte_carlo.n_steps,
            )
            for combo in swap_combos:
                self.logger.info("  swapping types %s ↔ %s", combo[0], combo[1])
            for idx, (type1, type2) in enumerate(swap_combos):
                lmp.command(
                    "fix  swap%d all atom/swap %d %d %d ${ftemp} ke yes types %s %s"
                    % (idx, self.calc.monte_carlo.n_steps,
                       self.calc.monte_carlo.n_swaps,
                       np.random.randint(1, 10000), type1, type2)
                )

        if self.calc.n_print_steps > 0:
            lmp.command(
                "dump              d1 all custom %d traj.ts.forward_%d.dat "
                "id type mass x y z vx vy vz"
                % (self.calc.n_print_steps, iteration)
            )

        # ── Forward sweep ───────────────────────────────────────────────────
        self.logger.info("forward sweep (iteration %d): sweep start", iteration)
        try:
            self._run_sweep(
                lmp=lmp,
                lambda_var="flambda",
                output_file_pattern="ts.forward_%d.dat" % iteration,
                sweep_label="forward (iteration %d)" % iteration,
            )
        except Exception:
            # Make sure the LAMMPS instance (especially in interactive
            # pylammpsmpi mode where it is reused for the backward sweep) is
            # cleared and the log file is renamed before the exception
            # propagates so a fresh lmp can be created for the backward sweep.
            try:
                self.lammps_close(lmp=lmp)
            except Exception as _close_exc:
                self.logger.debug(
                    "forward sweep cleanup: lammps_close failed: %s", _close_exc
                )
            logfile = os.path.join(self.simfolder, "log.lammps")
            if os.path.exists(logfile):
                try:
                    os.rename(
                        logfile,
                        os.path.join(
                            self.simfolder, "reversible_scaling_forward.log.lammps"
                        ),
                    )
                except Exception:
                    pass
            raise
        self.logger.info("forward sweep (iteration %d): sweep done", iteration)

        # ── Cleanup swaps / dump ────────────────────────────────────────────
        if self.calc.monte_carlo.n_swaps > 0:
            swap_types  = self.calc.monte_carlo.forward_swap_types
            swap_combos = list(itertools.combinations(swap_types, 2))
            for idx in range(len(swap_combos)):
                lmp.command("unfix swap%d" % idx)

        if self.calc.n_print_steps > 0:
            lmp.command("undump           d1")

        # ── Save forward-sweep end configuration ────────────────────────────
        conf_forward = os.path.join(
            self.simfolder, "conf.ts.forward_%d.data" % iteration
        )
        lmp.command("write_data        %s" % conf_forward)
        self.logger.info(
            "forward sweep (iteration %d): configuration saved to %s",
            iteration, os.path.basename(conf_forward),
        )

        self.lammps_close(lmp=lmp)

        logfile = os.path.join(self.simfolder, "log.lammps")
        if os.path.exists(logfile):
            os.rename(
                logfile,
                os.path.join(
                    self.simfolder, "reversible_scaling_forward.log.lammps"
                ),
            )

    def _reversible_scaling_backward(self, iteration: int = 1) -> None:
        """
        Perform the backward sweep of a reversible-scaling calculation.

        1. Load ``conf.ts.forward_{iteration}.data`` written by the forward
           sweep.
        2. Middle equilibration at Tf.
        3. Backward sweep: λ T0/Tf → 1.

        Parameters
        ----------
        iteration : int
            Reversible-scaling iteration index.
        """
        solid = self.calc.reference_phase == "solid"
        t0 = self.calc._temperature
        tf = self.calc._temperature_stop
        li = 1.0
        lf = t0 / tf
        pi = self.calc._pressure
        pf = lf * pi

        self.logger.info(
            "backward sweep (iteration %d): T %.1f → %.1f K, "
            "λ %.4f → %.4f, P %.4f → %.4f bar",
            iteration, tf, t0, lf, li, pf, pi,
        )

        lmp = ph.create_object(
            cores=self.cores,
            directory=self.simfolder,
            timestep=self.calc.md.timestep,
            cmdargs=self.calc.md.cmdargs,
            init_commands=self.calc.md.init_commands,
            script_mode=False,
            lmp=self._lmp,
        )

        lmp.command("echo              log")
        lmp.command("variable          li equal %f" % li)
        lmp.command("variable          lf equal %f" % lf)

        lmp = ph.set_pair_style(lmp, self.calc)

        conf = os.path.join(
            self.simfolder, "conf.ts.forward_%d.data" % iteration
        )
        lmp = ph.read_data(lmp, conf)

        lmp = ph.set_pair_coeff(lmp, self.calc)
        lmp = ph.set_mass(lmp, self.calc)

        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        # ── Re-install scaled potential at constant λ = lf BEFORE the
        # middle equilibration.  The forward sweep ended with the scaled
        # pair style active at λ = lf, so the snapshot stored in
        # ``conf.ts.forward_<iter>.data`` is in equilibrium with that
        # Hamiltonian (effective temperature Tf, expanded box).  If we
        # equilibrated here under the *unscaled* potential at T0, the
        # thermostat/barostat would re-thermalise to a much colder/denser
        # state, and the first samples of the backward sweep would show a
        # large transient bump in dU as the system re-expanded under the
        # scaled potential.  Using a constant scaling variable (rather
        # than the ramp) keeps λ frozen at lf during this run.
        lmp.command("variable          one equal 1.0")
        lmp.command("variable          bscale_eq equal %f" % (lf - 1.0))
        lmp.command(
            ph.scaled_pair_style_command(self.calc, ["v_one", "v_bscale_eq"])
        )
        for cmd in ph.hybrid_pair_coeff_commands(self.calc, repeat_index=0, total_repeats=2):
            lmp.command(cmd)
        for cmd in ph.hybrid_pair_coeff_commands(self.calc, repeat_index=1, total_repeats=2):
            lmp.command(cmd)

        lmp.command("variable         xcm equal xcm(all,x)")
        lmp.command("variable         ycm equal xcm(all,y)")
        lmp.command("variable         zcm equal xcm(all,z)")

        if self.calc.npt:
            lmp.command(
                "fix               f1 all npt temp %f %f %f %s %f %f %f "
                "fixedpoint ${xcm} ${ycm} ${zcm}"
                % (t0, t0, self.calc.md.thermostat_damping[1],
                   self.iso, pi, pi, self.calc.md.barostat_damping[1])
            )
        else:
            lmp.command(
                "fix               f1 all nvt temp %f %f %f "
                "fixedpoint ${xcm} ${ycm} ${zcm}"
                % (t0, t0, self.calc.md.thermostat_damping[1])
            )

        lmp.command("compute           tcm all temp/com")
        lmp.command("fix_modify        f1 temp tcm")
        lmp.command("variable          step    equal step")
        lmp.command("variable          dU      equal c_thermo_pe/atoms")
        lmp.command("thermo_style      custom step pe c_tcm press vol")
        lmp.command("thermo            10000")

        # ── Middle equilibration at effective Tf (scaled potential, λ=lf) ──
        self.logger.info(
            "backward sweep (iteration %d): middle equilibration start", iteration
        )
        lmp.command("run               %d" % self.calc.n_equilibration_steps)
        self.logger.info(
            "backward sweep (iteration %d): middle equilibration done", iteration
        )

        # Phase-stability check at Tf
        if not self.calc.script_mode:
            self.dump_current_snapshot(lmp, "traj.temp.dat")
            if solid:
                self.check_if_melted(lmp, "traj.temp.dat")
            else:
                self.check_if_solidfied(lmp, "traj.temp.dat")

        # ── Switch from constant-λ scaled potential to ramping scaled
        # potential for the backward sweep.  The scaled potential is
        # already active (set during the constant-lambda middle equil),
        # so no set_potential() call is needed.  We just re-define the
        # lambda variables for the sweep.
        # T0_rs is needed by both schedules for ftemp/btemp.
        lmp.command("variable         T0_rs equal %f" % t0)
        if self.calc.lambda_schedule == "uniform_temperature":
            lmp.command("variable         Nsweep equal %d" % self.calc._n_sweep_steps)
            lmp.command("variable         Tf_rs equal %f" % tf)
            lmp.command("variable         step0 equal $(step)")
            lmp.command(
                "variable         flambda equal "
                "v_T0_rs/(v_T0_rs+(v_Tf_rs-v_T0_rs)*(step-v_step0)/v_Nsweep)"
            )
            lmp.command(
                "variable         blambda equal "
                "v_T0_rs/(v_Tf_rs-(v_Tf_rs-v_T0_rs)*(step-v_step0)/v_Nsweep)"
            )
        else:  # "linear"
            lmp.command("variable         flambda equal ramp(${li},${lf})")
            lmp.command("variable         blambda equal ramp(${lf},${li})")
        lmp.command("variable         fscale equal v_flambda-1.0")
        lmp.command("variable         bscale equal v_blambda-1.0")
        lmp.command("variable         ftemp equal v_T0_rs/v_flambda")
        lmp.command("variable         btemp equal v_T0_rs/v_blambda")

        lmp.command(ph.scaled_pair_style_command(self.calc, ["v_one", "v_bscale"]))
        for cmd in ph.hybrid_pair_coeff_commands(self.calc, repeat_index=0, total_repeats=2):
            lmp.command(cmd)
        for cmd in ph.hybrid_pair_coeff_commands(self.calc, repeat_index=1, total_repeats=2):
            lmp.command(cmd)

        # ── Optional MC swaps ───────────────────────────────────────────────
        if (
            self.calc.monte_carlo.n_swaps > 0
            and len(self.calc.monte_carlo.reverse_swap_types) >= 2
        ):
            swap_types  = self.calc.monte_carlo.reverse_swap_types
            swap_combos = list(itertools.combinations(swap_types, 2))
            self.logger.info(
                "backward sweep (iteration %d): %d swap moves/combo, "
                "%d combinations every %d steps",
                iteration, self.calc.monte_carlo.n_swaps,
                len(swap_combos), self.calc.monte_carlo.n_steps,
            )
            for combo in swap_combos:
                self.logger.info("  swapping types %s ↔ %s", combo[0], combo[1])
            for idx, (type1, type2) in enumerate(swap_combos):
                lmp.command(
                    "fix  swap%d all atom/swap %d %d %d ${btemp} ke yes types %s %s"
                    % (idx, self.calc.monte_carlo.n_steps,
                       self.calc.monte_carlo.n_swaps,
                       np.random.randint(1, 10000), type1, type2)
                )

        if self.calc.n_print_steps > 0:
            lmp.command(
                "dump              d1 all custom %d traj.ts.backward_%d.dat "
                "id type mass x y z vx vy vz"
                % (self.calc.n_print_steps, iteration)
            )

        # ── Backward sweep ──────────────────────────────────────────────────
        self.logger.info("backward sweep (iteration %d): sweep start", iteration)
        self._run_sweep(
            lmp=lmp,
            lambda_var="blambda",
            output_file_pattern="ts.backward_%d.dat" % iteration,
            sweep_label="backward (iteration %d)" % iteration,
        )
        self.logger.info("backward sweep (iteration %d): sweep done", iteration)

        # ── Cleanup swaps / dump ────────────────────────────────────────────
        if self.calc.monte_carlo.n_swaps > 0:
            swap_types  = self.calc.monte_carlo.reverse_swap_types
            swap_combos = list(itertools.combinations(swap_types, 2))
            for idx in range(len(swap_combos)):
                lmp.command("unfix swap%d" % idx)

        if self.calc.n_print_steps > 0:
            lmp.command("undump           d1")

        if self.calc.script_mode:
            file = os.path.join(
                self.simfolder, "reversible_scaling_%d.lmp" % iteration
            )
            lmp.write(file)
            self.logger.info("Please cite the following publications:")
            self.logger.info("- 10.1103/PhysRevLett.83.3973")
            self.publications.append("10.1103/PhysRevLett.83.3973")
            return

        self.lammps_close(lmp=lmp)

        logfile = os.path.join(self.simfolder, "log.lammps")
        if os.path.exists(logfile):
            os.rename(
                logfile,
                os.path.join(
                    self.simfolder, "reversible_scaling_backward.log.lammps"
                ),
            )

    def reversible_scaling(self, iteration=1):
        """
        Perform reversible scaling calculation in NPT.

        Calls :meth:`_reversible_scaling_forward` (initial equilibration +
        forward sweep, saves ``conf.ts.forward_{iteration}.data``) followed by
        :meth:`_reversible_scaling_backward` (middle equilibration at Tf +
        backward sweep).

        Parameters
        ----------
        iteration : int, optional
            Iteration of the calculation.  Default 1.
        """
        self.logger.info("Starting temperature sweep cycle: %d", iteration)

        self._reversible_scaling_forward(iteration=iteration)
        self._reversible_scaling_backward(iteration=iteration)

        self.logger.info("Finished temperature sweep cycle: %d", iteration)
        self.logger.info("Please cite the following publications:")
        self.logger.info("- 10.1103/PhysRevLett.83.3973")
        self.publications.append("10.1103/PhysRevLett.83.3973")

    def integrate_reversible_scaling(self, scale_energy=True, return_values=False):
        """
        Perform integration after reversible scaling

        Parameters
        ----------
        scale_energy : bool, optional
            If True, scale the energy during reversible scaling.

        return_values : bool, optional
            If True, return integrated values

        Returns
        -------
        res : list of lists of shape 1x3
            Only returned if `return_values` is True.
        """
        res, ediss = integrate_rs(
            self.simfolder,
            self.fe,
            self.calc._temperature,
            self.natoms,
            p=self.calc._pressure,
            nsims=self.calc.n_iterations,
            scale_energy=scale_energy,
            return_values=return_values,
        )

        # Cache the max forward/backward energy dissipation on the Phase
        # instance.  Downstream code (e.g. MeltingTemp) reads this as a
        # quality flag: clean reversible sweeps give ediss ~ 1e-4 eV/atom,
        # whereas a hidden phase transition during the sweep produces
        # ediss orders of magnitude larger (the forward and backward
        # integrals can no longer match) and the resulting FE curve is
        # contaminated.
        self.ediss = float(ediss)

        self.logger.info(
            f"Maximum energy dissipation along the temperature scaling part: {ediss} eV/atom"
        )
        if np.abs(ediss) > 1e-4:
            self.logger.warning(
                f"Found max energy dissipation of {ediss} along the temperature scaling path. Please ensure there are no structural changes!"
            )

        if return_values:
            return res

    def scan_temperature_range(self):
        """
        Pre-flight temperature-range scan for a reversible-scaling (ts) run.

        Runs a single fast real-thermostat temperature ramp (T0 -> Tf under
        NPT) and analyses the fluctuation response functions to find the onset
        of a phase transition.  Depending on
        ``phase_transition_detection.mode`` the requested temperature range is
        then left as-is, reduced to the clean sub-range, or the run is aborted:

          'none'  — never called (the caller gates on mode != 'none').
          'adapt' — on detection, reduce ``calc._temperature_stop`` to the
                    detected clean onset (the number of switching steps is left
                    unchanged); a clean scan leaves the range untouched.
          'warn'  — log the detected clean range without modifying anything.
          'stop'  — on detection, raise PhaseTransitionError.

        Unlike the production sweep, this ramp uses a *measured* temperature
        (the thermostat genuinely ramps), so the response functions are the
        plain NPT fluctuation expressions with no lambda reduction — see
        :mod:`calphy.range_scan`.

        Returns
        -------
        None
        """
        from calphy.range_scan import RangeScan, plot_scan
        from calphy.errors import PhaseTransitionError

        td = self.calc.phase_transition_detection

        t0 = float(self.calc._temperature)
        tf = float(self.calc._temperature_stop)
        p0 = self.calc._pressure if self.calc._pressure is not None else 0.0
        n_scan = int(td.prescan_steps)

        self.logger.info("=" * 60)
        self.logger.info("Pre-flight temperature-range scan (mode=%s)", td.mode)
        self.logger.info(
            "pre-scan: requested ts range [%.1f, %.1f] K", t0, tf,
        )
        self.logger.info(
            "pre-scan: diagnostic ramp T %.1f -> %.1f K over %d steps "
            "(onset_fraction=%.2f)",
            t0, tf, n_scan, td.onset_fraction,
        )

        # ── Build the LAMMPS object and load the equilibrated configuration ──
        lmp = ph.create_object(
            cores=self.cores,
            directory=self.simfolder,
            timestep=self.calc.md.timestep,
            cmdargs=self.calc.md.cmdargs,
            init_commands=self.calc.md.init_commands,
            script_mode=False,
            lmp=self._lmp,
        )

        lmp.command("echo              log")
        lmp = ph.set_pair_style(lmp, self.calc)

        conf = os.path.join(self.simfolder, "conf.equilibration.data")
        lmp = ph.read_data(lmp, conf)

        lmp = ph.set_pair_coeff(lmp, self.calc)
        lmp = ph.set_mass(lmp, self.calc)
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        # ── Short equilibration at T0 ───────────────────────────────────────
        lmp.command(
            "fix               1 all npt temp %f %f %f %s %f %f %f"
            % (t0, t0, self.calc.md.thermostat_damping[1],
               self.iso, p0, p0, self.calc.md.barostat_damping[1])
        )
        lmp.command("run               %d" % self.calc.n_equilibration_steps)
        lmp.command("unfix             1")

        # ── Real-thermostat ramp T0 -> Tf, recording every step ─────────────
        pf = (t0 / tf) * p0
        lmp.command("variable          dU      equal pe/atoms")
        lmp.command(
            "fix               f2 all npt temp %f %f %f %s %f %f %f"
            % (t0, tf, self.calc.md.thermostat_damping[1],
               self.iso, p0, pf, self.calc.md.barostat_damping[1])
        )
        scan_file = "prescan.forward.dat"
        lmp.command(
            'fix               fp all print 1 "${dU} $(press) $(vol) $(temp)" '
            'screen no file %s' % scan_file
        )
        self.logger.info("pre-scan: ramp start (%d steps)", n_scan)
        lmp.command("run               %d" % n_scan)
        lmp.command("unfix             fp")
        lmp.command("unfix             f2")

        self.lammps_close(lmp=lmp)
        logfile = os.path.join(self.simfolder, "log.lammps")
        if os.path.exists(logfile):
            try:
                os.rename(logfile, os.path.join(self.simfolder, "prescan.log.lammps"))
            except OSError:
                pass

        # ── Analyse the ramp ────────────────────────────────────────────────
        scan_path = os.path.join(self.simfolder, scan_file)
        try:
            data = np.loadtxt(scan_path, comments="#")
        except Exception as exc:
            self.logger.warning(
                "pre-scan: could not read %s (%s) — skipping range check",
                scan_file, exc,
            )
            return
        if data.ndim < 2 or data.shape[0] < 100:
            self.logger.warning(
                "pre-scan: too few samples in %s — skipping range check", scan_file
            )
            return

        dU, press, vol, temp = data.T
        # Detector calibration (peak_threshold, min_agreement, onset_sigma,
        # onset_level, windows) uses RangeScan's internal defaults; users tune
        # the scan only through mode, prescan_steps and onset_fraction.
        scanner = RangeScan(target_pressure=p0)
        result = scanner.find_clean_range(
            pe=dU, press=press, vol_total=vol, temp=temp,
            natoms=self.natoms, t_start=t0, t_stop=tf,
        )

        # Save the diagnostic signal plot (best-effort; never fatal).
        plot_path = os.path.join(self.simfolder, "prescan_signals.png")
        if plot_scan(
            pe=dU, press=press, vol_total=vol, temp=temp,
            natoms=self.natoms, target_pressure=p0, outpath=plot_path,
            result=result,
        ):
            self.logger.info("pre-scan: signal plot saved to %s",
                             os.path.basename(plot_path))

        if not result.transition_found:
            self.logger.info(
                "pre-scan: RESULT — no phase transition detected over "
                "[%.1f, %.1f] K", t0, tf,
            )
            self.logger.info(
                "pre-scan: clean range = [%.1f, %.1f] K (full requested range); "
                "running ts sweep unchanged", t0, tf,
            )
            self.logger.info("=" * 60)
            return

        # A transition was found.  Apply the fractional safety margin toward T0:
        # the detected onset is the foot of the deviation in a *fast* ramp, but
        # the ts backward sweep equilibrates AT the boundary for many steps, so
        # its practical stability limit is somewhat below the ramp onset (and the
        # onset itself is noisy).  Back the boundary off by a fraction of the
        # super-heated/cooled span:  T_clean = T0 + frac * (T_onset - T0).
        frac = float(td.onset_fraction)
        t_onset = float(result.onset_temperature)
        t_clean = t0 + frac * (t_onset - t0)
        trimmed = abs(tf - t_clean)

        self.logger.warning(
            "pre-scan: RESULT — phase transition DETECTED", )
        self.logger.warning(
            "pre-scan:   triggering signals : %s (confidence %.0f%%)",
            ", ".join(result.triggered_signals), result.confidence * 100,
        )
        self.logger.warning(
            "pre-scan:   onset temperature  : %.1f K  (foot of deviation)",
            t_onset,
        )
        self.logger.warning(
            "pre-scan:   peak / collapse    : %.1f K", result.peak_temperature,
        )
        self.logger.warning(
            "pre-scan:   safety margin      : onset_fraction=%.2f  -> "
            "backed off %.1f K below onset", frac, t_onset - t_clean,
        )
        self.logger.warning(
            "pre-scan:   clean range        : [%.1f, %.1f] K  "
            "(requested [%.1f, %.1f] K, trimmed %.1f K)",
            t0, t_clean, t0, tf, trimmed,
        )

        if td.mode == "warn":
            self.logger.warning(
                "pre-scan (mode='warn'): NOT adapting — ts sweep runs over the "
                "full requested range [%.1f, %.1f] K despite the detected "
                "transition", t0, tf,
            )
            self.logger.info("=" * 60)
            return

        if td.mode == "stop":
            self.logger.info("=" * 60)
            raise PhaseTransitionError(
                "Pre-scan detected a phase transition over [%.1f, %.1f] K "
                "(onset ~ %.1f K, signals: %s).  The clean range is "
                "[%.1f, %.1f] K — re-submit with a corrected temperature range "
                "or set phase_transition_detection.mode: adapt."
                % (t0, tf, t_onset,
                   ", ".join(result.triggered_signals), t0, t_clean)
            )

        # td.mode == "adapt": reduce the upper temperature, keep the same
        # number of switching steps.
        old_tf = self.calc._temperature_stop
        self.calc._temperature_stop = float(t_clean)
        self.logger.warning(
            "pre-scan (mode='adapt'): ADAPTING ts range — T_stop %.1f K -> "
            "%.1f K", old_tf, t_clean,
        )
        self.logger.warning(
            "pre-scan (mode='adapt'): ts sweep will run [%.1f, %.1f] K over %d "
            "switching steps (n_sweep_steps unchanged)",
            t0, t_clean, self.calc._n_sweep_steps,
        )
        self.logger.info("=" * 60)

    def temperature_scaling(self, iteration=1):
        """
        Perform temperature scaling calculation in NPT

        Parameters
        ----------
        iteration : int, optional
            iteration of the calculation. Default 1

        Returns
        -------
        None
        """
        solid = False
        if self.calc.reference_phase == "solid":
            solid = True

        t0 = self.calc._temperature
        tf = self.calc._temperature_stop
        li = 1
        lf = t0 / tf
        p0 = self.calc._pressure
        pf = lf * p0

        # create lammps object
        lmp = ph.create_object(
            cores=self.cores,
            directory=self.simfolder,
            timestep=self.calc.md.timestep,
            cmdargs=self.calc.md.cmdargs,
            init_commands=self.calc.md.init_commands,
            script_mode=self.calc.script_mode,
            lmp=self._lmp,
        )

        lmp.command("echo              log")
        lmp.command("variable          li equal %f" % li)
        lmp.command("variable          lf equal %f" % lf)

        lmp = ph.set_pair_style(lmp, self.calc)

        # read in conf
        # conf = os.path.join(self.simfolder, "conf.dump")
        conf = os.path.join(self.simfolder, "conf.equilibration.data")
        lmp = ph.read_data(lmp, conf)

        # set up potential
        lmp = ph.set_pair_coeff(lmp, self.calc)
        lmp = ph.set_mass(lmp, self.calc)

        # remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        # equilibrate first
        lmp.command(
            "fix               1 all npt temp %f %f %f %s %f %f %f"
            % (
                t0,
                t0,
                self.calc.md.thermostat_damping[1],
                self.iso,
                p0,
                p0,
                self.calc.md.barostat_damping[1],
            )
        )
        lmp.command("run               %d" % self.calc.n_equilibration_steps)
        lmp.command("unfix             1")

        # now scale system to final temp, thereby recording enerfy at every step
        lmp.command("variable          step    equal step")
        lmp.command("variable          dU      equal pe/atoms")
        lmp.command("variable          lambda equal ramp(${li},${lf})")

        lmp.command(
            "fix               f2 all npt temp %f %f %f %s %f %f %f"
            % (
                t0,
                tf,
                self.calc.md.thermostat_damping[1],
                self.iso,
                p0,
                pf,
                self.calc.md.barostat_damping[1],
            )
        )

        self.logger.info(
            "ts-sweep tscale forward (iteration %d): T %.1f → %.1f K, "
            "%d steps",
            iteration, t0, tf, self.calc._n_sweep_steps,
        )
        self._run_sweep(
            lmp=lmp,
            lambda_var="lambda",
            output_file_pattern="ts.forward_%d.dat" % iteration,
            sweep_label="tscale forward (iteration %d)" % iteration,
        )

        lmp.command("unfix             f2")

        lmp.command(
            "fix               1 all npt temp %f %f %f %s %f %f %f"
            % (
                tf,
                tf,
                self.calc.md.thermostat_damping[1],
                self.iso,
                pf,
                pf,
                self.calc.md.barostat_damping[1],
            )
        )
        lmp.command("run               %d" % self.calc.n_equilibration_steps)
        lmp.command("unfix             1")

        # check melting or freezing
        lmp.command(
            "dump              2 all custom 1 traj.temp.dat id type mass x y z vx vy vz"
        )
        lmp.command("run               0")
        lmp.command("undump            2")

        if not self.calc.script_mode:
            self.dump_current_snapshot(lmp, "traj.temp.dat")
            if solid:
                self.check_if_melted(lmp, "traj.temp.dat")
            else:
                self.check_if_solidfied(lmp, "traj.temp.dat")

        # start reverse loop
        lmp.command("variable          lambda equal ramp(${lf},${li})")

        lmp.command(
            "fix               f2 all npt temp %f %f %f %s %f %f %f"
            % (
                t0,
                t0,
                self.calc.md.thermostat_damping[1],
                self.iso,
                p0,
                pf,
                self.calc.md.barostat_damping[1],
            )
        )

        self.logger.info(
            "ts-sweep tscale backward (iteration %d): T %.1f → %.1f K, "
            "%d steps",
            iteration, tf, t0, self.calc._n_sweep_steps,
        )
        self._run_sweep(
            lmp=lmp,
            lambda_var="lambda",
            output_file_pattern="ts.backward_%d.dat" % iteration,
            sweep_label="tscale backward (iteration %d)" % iteration,
        )

        lmp.command("unfix             f2")

        if self.calc.script_mode:
            file = os.path.join(
                self.simfolder, "temperature_scaling_%d.lmp" % iteration
            )
            lmp.write(file)
            return

        self.lammps_close(lmp=lmp)
        # Preserve log file
        logfile = os.path.join(self.simfolder, "log.lammps")
        if os.path.exists(logfile):
            os.rename(
                logfile, os.path.join(self.simfolder, "temperature_scaling.log.lammps")
            )

    def pressure_scaling(self, iteration=1):
        """
        Perform pressure scaling calculation in NPT

        Parameters
        ----------
        iteration : int, optional
            iteration of the calculation. Default 1

        Returns
        -------
        None
        """
        t0 = self.calc._temperature
        li = 1
        lf = self.calc._pressure_stop
        p0 = self.calc._pressure
        pf = self.calc._pressure_stop

        # create lammps object
        lmp = ph.create_object(
            cores=self.cores,
            directory=self.simfolder,
            timestep=self.calc.md.timestep,
            cmdargs=self.calc.md.cmdargs,
            init_commands=self.calc.md.init_commands,
            script_mode=self.calc.script_mode,
            lmp=self._lmp,
        )

        lmp.command("echo              log")
        lmp.command("variable          li equal %f" % li)
        lmp.command("variable          lf equal %f" % lf)
        lmp.command("variable          p0 equal %f" % p0)
        lmp.command("variable          pf equal %f" % pf)

        lmp = ph.set_pair_style(lmp, self.calc)

        # read in conf
        # conf = os.path.join(self.simfolder, "conf.dump")
        conf = os.path.join(self.simfolder, "conf.equilibration.data")
        lmp = ph.read_data(lmp, conf)

        # set up potential
        lmp = ph.set_pair_coeff(lmp, self.calc)
        lmp = ph.set_mass(lmp, self.calc)

        # remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        # equilibrate first
        lmp.command(
            "fix               1 all npt temp %f %f %f %s %f %f %f"
            % (
                t0,
                t0,
                self.calc.md.thermostat_damping[1],
                self.iso,
                p0,
                p0,
                self.calc.md.barostat_damping[1],
            )
        )
        lmp.command("run               %d" % self.calc.n_equilibration_steps)
        lmp.command("unfix             1")

        # now scale system to final temp, thereby recording enerfy at every step
        lmp.command("variable          step    equal step")
        lmp.command("variable          dU      equal pe/atoms")
        lmp.command("variable          lambda equal ramp(${li},${lf})")
        lmp.command("variable          pp equal ramp(${p0},${pf})")

        lmp.command(
            "fix               f2 all npt temp %f %f %f %s %f %f %f"
            % (
                t0,
                t0,
                self.calc.md.thermostat_damping[1],
                self.iso,
                p0,
                pf,
                self.calc.md.barostat_damping[1],
            )
        )
        lmp.command(
            'fix               f3 all print 1 "${dU} ${pp} $(vol) ${lambda}" screen no file ps.forward_%d.dat'
            % iteration
        )
        lmp.command("run               %d" % self.calc._n_sweep_steps)

        lmp.command("unfix             f2")
        lmp.command("unfix             f3")

        lmp.command(
            "fix               1 all npt temp %f %f %f %s %f %f %f"
            % (
                t0,
                t0,
                self.calc.md.thermostat_damping[1],
                self.iso,
                pf,
                pf,
                self.calc.md.barostat_damping[1],
            )
        )
        lmp.command("run               %d" % self.calc.n_equilibration_steps)
        lmp.command("unfix             1")

        # start reverse loop
        lmp.command("variable          lambda equal ramp(${lf},${li})")
        lmp.command("variable          pp equal ramp(${pf},${p0})")

        lmp.command(
            "fix               f2 all npt temp %f %f %f %s %f %f %f"
            % (
                t0,
                t0,
                self.calc.md.thermostat_damping[1],
                self.iso,
                pf,
                p0,
                self.calc.md.barostat_damping[1],
            )
        )
        lmp.command(
            'fix               f3 all print 1 "${dU} ${pp} $(vol) ${lambda}" screen no file ps.backward_%d.dat'
            % iteration
        )
        lmp.command("run               %d" % self.calc._n_sweep_steps)

        # Preserve log file
        logfile = os.path.join(self.simfolder, "log.lammps")
        if os.path.exists(logfile):
            os.rename(
                logfile, os.path.join(self.simfolder, "pressure_scaling.log.lammps")
            )

        self.logger.info("Please cite the following publications:")
        self.logger.info("- 10.1016/j.commatsci.2022.111275")
        self.publications.append("10.1016/j.commatsci.2022.111275")

    def integrate_pressure_scaling(self, return_values=False):
        """
        Perform integration after reversible scaling

        Parameters
        ----------
        scale_energy : bool, optional
            If True, scale the energy during reversible scaling.
        return_values : bool, optional
            If True, return integrated values

        Returns
        -------
        res : list of lists of shape 1x3
            Only returned if `return_values` is True.
        """
        res = integrate_ps(
            self.simfolder,
            self.fe,
            self.natoms,
            self.calc._pressure,
            self.calc._pressure_stop,
            nsims=self.calc.n_iterations,
            return_values=return_values,
        )

        if return_values:
            return res

    def clean_up(self):
        """
        Run a clean up job
        """
        # serialise input configuration
        shutil.copy(
            self.calc.lattice, os.path.join(self.simfolder, "input_configuration.data")
        )

        # write simple metadata
        metadata = generate_metadata()
        metadata["publications"] = self.publications

        with open(os.path.join(self.simfolder, "metadata.yaml"), "w") as fout:
            yaml.safe_dump(metadata, fout)

    def lammps_close(self, lmp):
        if self._lmp is None:
            lmp.close()
        else:
            lmp.clear()
