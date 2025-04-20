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

import pyscal3.traj_process as ptp
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

    def __init__(self, calculation=None, simfolder=None, log_to_screen=False):

        self.calc = copy.deepcopy(calculation)

        # serialise input
        indict = {"calculations": [self.calc.dict()]}
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

        if self.calc._iso:
            self.iso = "iso"
        else:
            self.iso = "aniso"
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

        if self.calc.equilibration_control == "nose-hoover":
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
        self.ufm_cutoff = 5 * self.calc.uhlenbeck_ford_model.sigma

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
            if self.calc.potential_file is not None:
                self.logger.info("potential is being loaded from file instead")

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
        """ """
        solids = ph.find_solid_fraction(os.path.join(self.simfolder, filename))
        if solids / lmp.natoms < self.calc.tolerance.solid_fraction:
            lmp.close()
            raise MeltedError(
                "System melted, increase size or reduce temp!\n Solid detection algorithm only works with BCC/FCC/HCP/SC/DIA. Detection algorithm can be turned off by setting:\n tolerance.solid_fraction: 0"
            )

    def check_if_solidfied(self, lmp, filename):
        """ """
        solids = ph.find_solid_fraction(os.path.join(self.simfolder, filename))
        if solids / lmp.natoms > self.calc.tolerance.liquid_fraction:
            lmp.close()
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
        if self.calc.equilibration_control == "nose-hoover":
            self.fix_nose_hoover(lmp)
        else:
            self.fix_berendsen(lmp)
        # start thermo logging
        lmp.command("thermo_style     custom step pe press vol etotal temp lx ly lz")
        lmp.command("thermo           10")

        # run MD
        lmp.command("run              %d" % int(self.calc.md.n_small_steps))

        # remove fixes
        if self.calc.equilibration_control == "nose-hoover":
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

        # for Nose-Hoover thermo/baro combination
        if self.calc.equilibration_control == "nose-hoover":
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
        The cycle is stopped when the average pressure is within the given cutoff of the target pressure.
        """

        # apply fixes
        if self.calc.equilibration_control == "nose-hoover":
            self.fix_nose_hoover(lmp)
        else:
            self.fix_berendsen(lmp)

        lmp.command(
            "fix              2 all ave/time %d %d %d v_mlx v_mly v_mlz v_mpress file avg.dat"
            % (
                int(self.calc.md.n_every_steps),
                int(self.calc.md.n_repeat_steps),
                int(self.calc.md.n_every_steps * self.calc.md.n_repeat_steps),
            )
        )

        laststd = 0.00
        converged = False

        for i in range(int(self.calc.md.n_cycles)):
            lmp.command("run              %d" % int(self.calc.md.n_small_steps))
            ncount = int(self.calc.md.n_small_steps) // int(
                self.calc.md.n_every_steps * self.calc.md.n_repeat_steps
            )
            # now we can check if it converted
            file = os.path.join(self.simfolder, "avg.dat")
            lx, ly, lz, ipress = np.loadtxt(file, usecols=(1, 2, 3, 4), unpack=True)

            lxpc = ipress
            mean = np.mean(lxpc)
            std = np.std(lxpc)
            volatom = np.mean((lx * ly * lz) / self.natoms)
            self.logger.info(
                "At count %d mean pressure is %f with %f vol/atom"
                % (i + 1, mean, volatom)
            )

            if (np.abs(mean - self.calc._pressure)) < self.calc.tolerance.pressure:

                # process other means
                self.lx = np.round(np.mean(lx[-ncount + 1 :]), decimals=3)
                self.ly = np.round(np.mean(ly[-ncount + 1 :]), decimals=3)
                self.lz = np.round(np.mean(lz[-ncount + 1 :]), decimals=3)
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
            laststd = std

        if not converged:
            lmp.close()
            raise ValueError(
                "Pressure did not converge after MD runs, maybe change lattice_constant and try?"
            )

        # unfix thermostat and barostat
        if self.calc.equilibration_control == "nose-hoover":
            self.unfix_nose_hoover(lmp)
        else:
            self.unfix_berendsen(lmp)

        lmp.command("unfix            2")

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
        if self.calc.equilibration_control == "nose-hoover":
            self.fix_nose_hoover(lmp)
        else:
            self.fix_berendsen(lmp)

        lmp.command(
            "fix              2 all ave/time %d %d %d v_mlx v_mly v_mlz v_mpress file avg.dat"
            % (
                int(self.calc.md.n_every_steps),
                int(self.calc.md.n_repeat_steps),
                int(self.calc.md.n_every_steps * self.calc.md.n_repeat_steps),
            )
        )

        lmp.command("run              %d" % int(self.calc.md.n_small_steps))
        lmp.command("run              %d" % int(self.calc.n_equilibration_steps))

        # unfix thermostat and barostat
        if self.calc.equilibration_control == "nose-hoover":
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
            "fix              2 all ave/time %d %d %d v_mlx v_mly v_mlz v_mpress file avg.dat"
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

        lmp.command("unfix            1")
        lmp.command("unfix            2")

        if not converged:
            lmp.close()
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
        lz = lx[-ncount + 1 :]
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
        lz = lx[-ncount + 1 :]
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
            "fix              2 all ave/time %d %d %d v_mlx v_mly v_mlz v_mpress file avg.dat"
            % (
                int(self.calc.md.n_every_steps),
                int(self.calc.md.n_repeat_steps),
                int(self.calc.md.n_every_steps * self.calc.md.n_repeat_steps),
            )
        )

        lmp.command("run              %d" % int(self.calc.md.n_small_steps))
        lmp.command("run              %d" % int(self.calc.n_equilibration_steps))

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

    def reversible_scaling(self, iteration=1):
        """
        Perform reversible scaling calculation in NPT

        Parameters
        ----------
        iteration : int, optional
            iteration of the calculation. Default 1

        Returns
        -------
        None
        """
        self.logger.info(f"Starting temperature sweep cycle: {iteration}")
        solid = False
        if self.calc.reference_phase == "solid":
            solid = True

        t0 = self.calc._temperature
        tf = self.calc._temperature_stop
        li = 1
        lf = t0 / tf
        pi = self.calc._pressure
        pf = lf * pi

        # create lammps object
        lmp = ph.create_object(
            self.cores,
            self.simfolder,
            self.calc.md.timestep,
            self.calc.md.cmdargs,
            self.calc.md.init_commands,
        )

        lmp.command("echo              log")
        lmp.command("variable          li equal %f" % li)
        lmp.command("variable          lf equal %f" % lf)

        lmp.command(f"pair_style {self.calc._pair_style_with_options[0]}")

        # read in conf file
        # conf = os.path.join(self.simfolder, "conf.equilibration.dump")
        conf = os.path.join(self.simfolder, "conf.equilibration.data")
        lmp = ph.read_data(lmp, conf)

        # set up potential
        lmp.command(f"pair_coeff {self.calc.pair_coeff[0]}")
        lmp = ph.set_mass(lmp, self.calc)

        # remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        # set thermostat and run equilibrium
        if self.calc.npt:
            lmp.command(
                "fix               f1 all npt temp %f %f %f %s %f %f %f"
                % (
                    t0,
                    t0,
                    self.calc.md.thermostat_damping[1],
                    self.iso,
                    pi,
                    pi,
                    self.calc.md.barostat_damping[1],
                )
            )
        else:
            lmp.command(
                "fix               f1 all nvt temp %f %f %f"
                % (t0, t0, self.calc.md.thermostat_damping[1])
            )

        self.logger.info(f"Starting equilibration: {iteration}")
        lmp.command("run               %d" % self.calc.n_equilibration_steps)
        self.logger.info(f"Finished equilibration: {iteration}")

        lmp.command("unfix             f1")

        # now fix com
        lmp.command("variable         xcm equal xcm(all,x)")
        lmp.command("variable         ycm equal xcm(all,y)")
        lmp.command("variable         zcm equal xcm(all,z)")

        if self.calc.npt:
            lmp.command(
                "fix               f1 all npt temp %f %f %f %s %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"
                % (
                    t0,
                    t0,
                    self.calc.md.thermostat_damping[1],
                    self.iso,
                    pi,
                    pi,
                    self.calc.md.barostat_damping[1],
                )
            )
        else:
            lmp.command(
                "fix               f1 all nvt temp %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"
                % (t0, t0, self.calc.md.thermostat_damping[1])
            )

        # compute com and modify fix
        lmp.command("compute           tcm all temp/com")
        lmp.command("fix_modify        f1 temp tcm")

        lmp.command("variable          step    equal step")
        lmp.command("variable          dU      equal c_thermo_pe/atoms")
        lmp.command("thermo_style      custom step pe c_tcm press vol")
        lmp.command("thermo            10000")

        # create velocity and equilibriate
        lmp.command(
            "velocity          all create %f %d mom yes rot yes dist gaussian"
            % (t0, np.random.randint(1, 10000))
        )

        self.logger.info(f"Starting equilibration with constrained com: {iteration}")
        lmp.command("run               %d" % self.calc.n_equilibration_steps)
        self.logger.info(f"Finished equilibration with constrained com: {iteration}")

        lmp.command("variable         flambda equal ramp(${li},${lf})")
        lmp.command("variable         blambda equal ramp(${lf},${li})")
        lmp.command("variable         fscale equal v_flambda-1.0")
        lmp.command("variable         bscale equal v_blambda-1.0")
        lmp.command("variable         one equal 1.0")
        lmp.command(
            f"variable        ftemp equal v_blambda*{self.calc._temperature_stop}"
        )
        lmp.command(
            f"variable        btemp equal v_flambda*{self.calc._temperature_stop}"
        )

        # set up potential
        pc = self.calc.pair_coeff[0]
        pcraw = pc.split()
        pcnew1 = " ".join(
            [
                *pcraw[:2],
                *[
                    self.calc._pair_style_names[0],
                ],
                "1",
                *pcraw[2:],
            ]
        )
        pcnew2 = " ".join(
            [
                *pcraw[:2],
                *[
                    self.calc._pair_style_names[0],
                ],
                "2",
                *pcraw[2:],
            ]
        )

        lmp.command(
            "pair_style       hybrid/scaled v_one %s v_fscale %s"
            % (
                self.calc._pair_style_with_options[0],
                self.calc._pair_style_with_options[0],
            )
        )
        lmp.command("pair_coeff       %s" % pcnew1)
        lmp.command("pair_coeff       %s" % pcnew2)

        lmp.command(
            'fix               f3 all print 1 "${dU} $(press) $(vol) ${flambda}" screen no file ts.forward_%d.dat'
            % iteration
        )

        # add swaps if n_swap is > 0
        if self.calc.monte_carlo.n_swaps > 0:
            self.logger.info(
                f"{self.calc.monte_carlo.n_swaps} swap moves are performed between {self.calc.monte_carlo.swap_types[0]} and {self.calc.monte_carlo.swap_types[1]} every {self.calc.monte_carlo.n_steps}"
            )
            lmp.command(
                "fix  swap all atom/swap %d %d %d ${ftemp} ke yes types %d %d"
                % (
                    self.calc.monte_carlo.n_steps,
                    self.calc.monte_carlo.n_swaps,
                    np.random.randint(1, 10000),
                    self.calc.monte_carlo.swap_types[0],
                    self.calc.monte_carlo.swap_types[1],
                )
            )

            lmp.command("variable a equal f_swap[1]")
            lmp.command("variable b equal f_swap[2]")
            lmp.command(
                'fix             swap2 all print 1 "${a} ${b} ${ftemp}" screen no file swap.rs.forward_%d.dat'
                % iteration
            )

        if self.calc.n_print_steps > 0:
            lmp.command(
                "dump              d1 all custom %d traj.ts.forward_%d.dat id type mass x y z vx vy vz"
                % (self.calc.n_print_steps, iteration)
            )

        self.logger.info(f"Started forward sweep: {iteration}")
        lmp.command("run               %d" % self.calc._n_sweep_steps)
        self.logger.info(f"Finished forward sweep: {iteration}")

        if self.calc.monte_carlo.n_swaps > 0:
            lmp.command("unfix swap")
            lmp.command("unfix swap2")

        # unfix
        lmp.command("unfix             f3")
        # lmp.command("unfix             f1")

        if self.calc.n_print_steps > 0:
            lmp.command("undump           d1")

        # switch potential
        lmp.command("run               %d" % self.calc.n_equilibration_steps)

        # check melting or freezing
        if not self.calc.script_mode:
            self.dump_current_snapshot(lmp, "traj.temp.dat")
            if solid:
                self.check_if_melted(lmp, "traj.temp.dat")
            else:
                self.check_if_solidfied(lmp, "traj.temp.dat")

        lmp = ph.set_potential(lmp, self.calc)

        # reverse scaling
        lmp.command("variable         flambda equal ramp(${li},${lf})")
        lmp.command("variable         blambda equal ramp(${lf},${li})")
        lmp.command("variable         fscale equal v_flambda-1.0")
        lmp.command("variable         bscale equal v_blambda-1.0")
        lmp.command("variable         one equal 1.0")
        lmp.command(
            f"variable        ftemp equal v_blambda*{self.calc._temperature_stop}"
        )
        lmp.command(
            f"variable        btemp equal v_flambda*{self.calc._temperature_stop}"
        )

        lmp.command(
            "pair_style       hybrid/scaled v_one %s v_bscale %s"
            % (
                self.calc._pair_style_with_options[0],
                self.calc._pair_style_with_options[0],
            )
        )
        lmp.command("pair_coeff       %s" % pcnew1)
        lmp.command("pair_coeff       %s" % pcnew2)

        # apply fix and perform switching
        lmp.command(
            'fix               f3 all print 1 "${dU} $(press) $(vol) ${blambda}" screen no file ts.backward_%d.dat'
            % iteration
        )

        if self.calc.n_print_steps > 0:
            lmp.command(
                "dump              d1 all custom %d traj.ts.backward_%d.dat id type mass x y z vx vy vz"
                % (self.calc.n_print_steps, iteration)
            )

        # add swaps if n_swap is > 0
        if self.calc.monte_carlo.n_swaps > 0:
            self.logger.info(
                f"{self.calc.monte_carlo.n_swaps} swap moves are performed between {self.calc.monte_carlo.swap_types[1]} and {self.calc.monte_carlo.swap_types[0]} every {self.calc.monte_carlo.n_steps}"
            )
            lmp.command(
                "fix  swap all atom/swap %d %d %d ${btemp} ke yes types %d %d"
                % (
                    self.calc.monte_carlo.n_steps,
                    self.calc.monte_carlo.n_swaps,
                    np.random.randint(1, 10000),
                    self.calc.monte_carlo.swap_types[1],
                    self.calc.monte_carlo.swap_types[0],
                )
            )

            lmp.command("variable a equal f_swap[1]")
            lmp.command("variable b equal f_swap[2]")
            lmp.command(
                'fix             swap2 all print 1 "${a} ${b} ${btemp}" screen no file swap.rs.backward_%d.dat'
                % iteration
            )

        self.logger.info(f"Started backward sweep: {iteration}")
        lmp.command("run               %d" % self.calc._n_sweep_steps)
        self.logger.info(f"Finished backward sweep: {iteration}")

        if self.calc.monte_carlo.n_swaps > 0:
            lmp.command("unfix swap")
            lmp.command("unfix swap2")

        lmp.command("unfix             f3")

        if self.calc.n_print_steps > 0:
            lmp.command("undump           d1")

        # close the object
        lmp.close()

        self.logger.info("Please cite the following publications:")
        if self.calc.mode == "mts":
            self.logger.info("- 10.1063/1.1420486")
            self.publications.append("10.1063/1.1420486")
        else:
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
        res = integrate_rs(
            self.simfolder,
            self.fe,
            self.calc._temperature,
            self.natoms,
            p=self.calc._pressure,
            nsims=self.calc.n_iterations,
            scale_energy=scale_energy,
            return_values=return_values,
        )

        if return_values:
            return res

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
            self.cores,
            self.simfolder,
            self.calc.md.timestep,
            self.calc.md.cmdargs,
            self.calc.md.init_commands,
        )

        lmp.command("echo              log")
        lmp.command("variable          li equal %f" % li)
        lmp.command("variable          lf equal %f" % lf)

        lmp.command(f"pair_style {self.calc._pair_style_with_options[0]}")

        # read in conf
        # conf = os.path.join(self.simfolder, "conf.equilibration.dump")
        conf = os.path.join(self.simfolder, "conf.equilibration.data")
        lmp = ph.read_data(lmp, conf)

        # set up potential
        lmp.command(f"pair_coeff {self.calc.pair_coeff[0]}")
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
        lmp.command(
            'fix               f3 all print 1 "${dU} $(press) $(vol) ${lambda}" screen no file ts.forward_%d.dat'
            % iteration
        )
        lmp.command("run               %d" % self.calc._n_sweep_steps)

        lmp.command("unfix             f2")
        lmp.command("unfix             f3")

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
                tf,
                t0,
                self.calc.md.thermostat_damping[1],
                self.iso,
                pf,
                p0,
                self.calc.md.barostat_damping[1],
            )
        )
        lmp.command(
            'fix               f3 all print 1 "${dU} $(press) $(vol) ${lambda}" screen no file ts.backward_%d.dat'
            % iteration
        )
        lmp.command("run               %d" % self.calc._n_sweep_steps)

        lmp.close()

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
            self.cores,
            self.simfolder,
            self.calc.md.timestep,
            self.calc.md.cmdargs,
            self.calc.md.init_commands,
        )

        lmp.command("echo              log")
        lmp.command("variable          li equal %f" % li)
        lmp.command("variable          lf equal %f" % lf)
        lmp.command("variable          p0 equal %f" % p0)
        lmp.command("variable          pf equal %f" % pf)

        lmp.command(f"pair_style {self.calc._pair_style_with_options[0]}")

        # read in conf
        # conf = os.path.join(self.simfolder, "conf.dump")
        conf = os.path.join(self.simfolder, "conf.equilibration.data")
        lmp = ph.read_data(lmp, conf)

        # set up potential
        lmp.command(f"pair_coeff {self.calc.pair_coeff[0]}")
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

        lmp.close()

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
