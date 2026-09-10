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
"""Dynamic Clausius-Clapeyron integration (``mode: dcci``).

Traces a solid-liquid coexistence line P_coex(T) from one known coexistence
point in a single nonequilibrium run per direction, following de Koning,
Antonelli and Yip, J. Chem. Phys. 115, 11025 (2001).  Both phases are driven
side by side at the fixed kinetic temperature T0 of the starting point with the
scaled Hamiltonian K + lambda*U; the scaled pressure P_RS is ramped linearly and
lambda follows from the Clausius-Clapeyron condition in its scaled form,

    d lambda / d P_RS = -(v_s - v_l) / (u_s - u_l),

evaluated from block averages of the potential energy u = <pe>/(N lambda) and
volume v = <V>/N per atom of the two cells.  Every point maps back to the real
system through T = T0/lambda and P = P_RS/lambda.  See DCCI_PLAN.md.
"""
import copy
import os
import math

import numpy as np
import yaml

import calphy.helpers as ph
from calphy.input import Calculation
from calphy.integrators import EV_A3_TO_BAR
from calphy.liquid import Liquid
from calphy.solid import Solid
from calphy.routines import enable_phase_detection


# --------------------------------------------------------------------------- #
# Pure numerics (no LAMMPS): kept module-level for unit testing
# --------------------------------------------------------------------------- #
def cce_slope(u_s, u_l, v_s, v_l):
    """
    Slope d lambda / d P_RS of the coexistence line in the scaled formulation.

    Parameters
    ----------
    u_s, u_l : float
        Real (unscaled) potential energy per atom of the solid and liquid cell,
        eV/atom.
    v_s, v_l : float
        Volume per atom of the two cells, A^3/atom.

    Returns
    -------
    float
        ``-(v_s - v_l) / (u_s - u_l)`` in A^3/eV, i.e. per (eV/A^3) of scaled
        pressure.  Along a melting line both differences are negative for a
        normal liquid, so the slope is negative: lambda falls (T rises) as the
        scaled pressure rises.
    """
    du = u_s - u_l
    if du == 0:
        raise ZeroDivisionError("solid and liquid potential energies coincide")
    return -(v_s - v_l) / du


def clausius_clapeyron_slope(u_s, u_l, v_s, v_l, p_real_bar, t_real):
    """
    Real-space slope dP/dT = dH / (T dV) in bar/K, for diagnostics.

    dH = (u_l - u_s) + P (v_l - v_s): the kinetic energy per atom is the same in
    both cells at equal temperature and drops out.
    """
    dv = v_l - v_s
    if dv == 0:
        return float("nan")
    dh = (u_l - u_s) + (p_real_bar / EV_A3_TO_BAR) * dv
    return dh / (t_real * dv) * EV_A3_TO_BAR


def lambda_ramp_command(name, lam_start, lam_end, step_start, n_steps):
    """
    ``variable`` command ramping ``name`` linearly from ``lam_start`` at
    ``step_start`` to ``lam_end`` at ``step_start + n_steps``.

    Written with literal numbers only (no ``$(...)`` immediate evaluation) so
    that the ExecutableRunner can replay it across a segment boundary.
    """
    return "variable          %s equal %.12g+(%.12g)*(step-%d)/%d" % (
        name, lam_start, lam_end - lam_start, step_start, n_steps,
    )


def plan_blocks(n_sweep, n_block):
    """Number of whole blocks covering ``n_sweep`` steps (rounded up)."""
    return max(1, int(math.ceil(n_sweep / float(n_block))))


class SweepState:
    """Where a sweep ended: what the next sweep starts from."""

    def __init__(self, lam, p_rs, n_blocks, direction):
        self.lam = lam
        self.p_rs = p_rs
        self.n_blocks = n_blocks
        self.direction = direction


class DynamicCCI:
    """
    Coordinator of a dynamic Clausius-Clapeyron integration.

    Owns one :class:`~calphy.solid.Solid` and one :class:`~calphy.liquid.Liquid`
    sub-job, equilibrated independently at the starting coexistence point and
    then stepped in lockstep through the coupled sweep.

    Parameters
    ----------
    calculation : Calculation
        The validated ``mode: dcci`` calculation.
    simfolder : str
        Folder for the coupled-sweep output; the cells live in ``solid/`` and
        ``liquid/`` below it.
    log_to_screen : bool
        Also log to the screen.
    """

    CELLS = ("solid", "liquid")

    # Thresholds installed for the two cells unless the user set them.
    DETECTION_SOLID_FRACTION = 0.7
    DETECTION_LIQUID_FRACTION = 0.05

    def __init__(self, calculation=None, simfolder=None, log_to_screen=False):
        self.calc = copy.deepcopy(calculation)
        self.simfolder = simfolder
        self.log_to_screen = log_to_screen
        self.publications = ["10.1063/1.1420486", "10.1103/PhysRevLett.83.3973"]

        if self.calc.md.seed is None:
            self.calc.md.seed = int(np.random.SeedSequence().entropy % (2**31 - 2)) + 1
        np.random.seed(self.calc.md.seed)

        with open(os.path.join(simfolder, "input_file.yaml"), "w") as fout:
            yaml.safe_dump({"calculations": [self.calc.model_dump()]}, fout)

        logfile = os.path.join(self.simfolder, "calphy.log")
        self.logger = ph.prepare_log(logfile, screen=log_to_screen)

        self.t0 = float(self.calc._temperature)
        self.p_start = float(self.calc._pressure)
        self.p_stop = float(self.calc._pressure_stop)
        self.n_block = int(self.calc.dcci.n_block_steps)
        self.n_blocks = plan_blocks(self.calc._n_sweep_steps, self.n_block)
        self.logger.info(
            "Dynamic Clausius-Clapeyron integration from the coexistence point "
            "T = %.2f K, P = %.2f bar towards P = %.2f bar" % (self.t0, self.p_start, self.p_stop)
        )
        self.logger.info(
            "Sweep of %d blocks x %d steps (n_switching_steps %d), corrector %s"
            % (self.n_blocks, self.n_block, self.calc._n_sweep_steps, self.calc.dcci.integrator)
        )
        self.logger.info("Master random seed is %d (md.seed)" % self.calc.md.seed)

        self.jobs = {}
        self.states = {}

    # ------------------------------------------------------------------ cells
    def _cell_calculation(self, phase):
        """
        The sub-calculation of one cell: the parent input with the cell's
        reference phase, the starting pressure and the melt/freeze checks on.
        """
        data = self.calc.model_dump()
        data["reference_phase"] = phase
        data["pressure"] = [self.p_start, self.p_start]
        data["temperature"] = self.t0
        data["kernel"] = 0
        data["inputfile"] = os.path.join(self.simfolder, "input_file.yaml")
        # the structure file has been generated by the parent already
        data["lattice"] = self.calc.lattice
        data["file_format"] = "lammps-data"
        # model_dump fills the defaults in; enable_phase_detection must only
        # act where the user gave no value, so blank those out again
        for key in ("solid_fraction", "liquid_fraction"):
            if key not in self.calc.tolerance.model_fields_set:
                data["tolerance"][key] = None
        enable_phase_detection(
            data, self.logger, "dcci",
            self.DETECTION_SOLID_FRACTION, self.DETECTION_LIQUID_FRACTION,
        )
        return Calculation(**data)

    def prepare_cells(self):
        """Create the ``solid/`` and ``liquid/`` sub-jobs."""
        for phase, cls in (("solid", Solid), ("liquid", Liquid)):
            folder = os.path.join(self.simfolder, phase)
            os.makedirs(folder, exist_ok=True)
            # validation copies the structure file into the working directory:
            # let that be the cell's folder
            cwd = os.getcwd()
            os.chdir(folder)
            try:
                calc = self._cell_calculation(phase)
            finally:
                os.chdir(cwd)
            job = cls(calculation=calc, simfolder=folder)
            for handler in self.logger.handlers:
                job.logger.addHandler(handler)
            self.jobs[phase] = job
        self.logger.info(
            "Cells: solid %s (%d atoms), liquid %s (%d atoms)"
            % (self.jobs["solid"].calc.lattice, self.jobs["solid"].natoms,
               self.jobs["liquid"].calc.lattice, self.jobs["liquid"].natoms)
        )

    def equilibrate_cells(self):
        """
        Equilibrate both cells at (T0, P_start) with the standard averaging
        routines, which converge the pressure and write ``conf.equilibration.data``.
        """
        for phase in self.CELLS:
            self.logger.info("Equilibrating the %s cell at %.2f K, %.2f bar" % (phase, self.t0, self.p_start))
            self.jobs[phase].run_averaging()
            self.logger.info(
                "%s cell: %.4f A^3/atom" % (phase, self.jobs[phase].volatom)
            )

    # ------------------------------------------------------------------ sweep
    def _open_cell(self, phase, lam, p_rs, conf, direction, iteration):
        """
        Open a LAMMPS session for one cell with the lambda-scaled potential at
        constant ``lam`` and a barostat at the constant scaled pressure ``p_rs``,
        and equilibrate it there.
        """
        job = self.jobs[phase]
        lmp = ph.create_object(job.calc, job.simfolder)
        lmp.command("echo              log")
        lmp.command("variable          lam equal %.12g" % lam)
        lmp.command(ph.scaled_pair_style_command(job.calc, ["v_lam"]))
        lmp = ph.read_data(lmp, conf)
        for cmd in ph.hybrid_pair_coeff_commands(job.calc):
            lmp.command(cmd)
        lmp = ph.set_mass(lmp, job.calc)

        if direction == "forward":
            # conf.equilibration.data is an equilibrium sample at (T0, P_start)
            # with fluctuating box: pin it to the averaged box as every sweep does
            lmp = ph.remap_box(lmp, job.lx, job.ly, job.lz)
            lmp.command(
                "velocity          all create %f %d mom yes rot yes dist gaussian"
                % (self.t0, np.random.randint(1, 10000))
            )

        lmp.command(
            "fix               f1 all npt temp %f %f %f %s %f %f %f"
            % (self.t0, self.t0, job.calc.md.thermostat_damping[1],
               job.iso, p_rs, p_rs, job.calc.md.barostat_damping[1])
        )
        lmp.command("thermo_style      custom step pe press vol v_lam")
        lmp.command("thermo            %d" % max(self.n_block, 1))
        self.logger.info(
            "%s sweep (iteration %d), %s cell: equilibration of %d steps at "
            "lambda = %.6f, P_RS = %.2f bar"
            % (direction, iteration, phase, job.calc.n_equilibration_steps, lam, p_rs)
        )
        lmp.command("run               %d" % job.calc.n_equilibration_steps)
        lmp.command("unfix             f1")
        return lmp

    def _arm_sweep(self, phase, lmp, p_rs_start, p_rs_stop, direction, iteration):
        """Install the ramping barostat and the block-average output of one cell."""
        job = self.jobs[phase]
        n_total = self.n_blocks_current * self.n_block
        lmp.command("reset_timestep    0")
        lmp.command(
            "fix               f1 all npt temp %f %f %f %s %f %f %f"
            % (self.t0, self.t0, job.calc.md.thermostat_damping[1],
               job.iso, p_rs_start, p_rs_stop, job.calc.md.barostat_damping[1])
        )
        # real potential energy per atom: LAMMPS reports lambda*U under hybrid/scaled
        lmp.command("variable          u_real equal pe/atoms/v_lam")
        lmp.command("variable          vol_atom equal vol/atoms")
        lmp.command(
            "fix               fav all ave/time 1 %d %d v_u_real v_vol_atom v_lam "
            'title1 "# block averages of the dcci %s sweep, %s cell" '
            'title2 "# step u_real[eV/atom] vol[A^3/atom] lambda" '
            "file %s" % (self.n_block, self.n_block, direction, phase,
                         self._cell_file(direction, iteration))
        )
        if self.calc.n_print_steps > 0:
            lmp.command(
                "dump              d1 all custom %d traj.dcci.%s_%d.dat "
                "id type mass x y z vx vy vz" % (self.calc.n_print_steps, direction, iteration)
            )
        return n_total

    @staticmethod
    def _cell_file(direction, iteration):
        """Raw block averages of one cell, inside the cell folder."""
        return "dcci.blocks.%s_%d.dat" % (direction, iteration)

    def _driver_file(self, direction, iteration):
        return os.path.join(self.simfolder, "dcci.%s_%d.dat" % (direction, iteration))

    def _read_block(self, phase, lmp, direction, iteration, k):
        """Block-k averages (u_real, vol_atom, lambda) of one cell."""
        data = np.atleast_2d(lmp.read_timeseries(self._cell_file(direction, iteration)))
        if data.shape[0] != k + 1:
            raise RuntimeError(
                "dcci %s sweep, %s cell: expected %d block-average rows after block %d, found %d"
                % (direction, phase, k + 1, k, data.shape[0])
            )
        return data[-1, 1], data[-1, 2], data[-1, 3]

    def _check_phases(self, lmps, direction, iteration, tag):
        """Melt/freeze checks on both cells (raise MeltedError / SolidifiedError)."""
        for phase in self.CELLS:
            job = self.jobs[phase]
            snapshot = "traj.dcci.%s_%d.%s.dat" % (direction, iteration, tag)
            job.dump_current_snapshot(lmps[phase], snapshot)
            if phase == "solid":
                job.check_if_melted(lmps[phase], snapshot)
            else:
                job.check_if_solidfied(lmps[phase], snapshot)

    def run_sweep(self, direction="forward", iteration=1):
        """
        Drive both cells through one sweep of the scaled pressure and integrate
        lambda block by block.

        Forward: from (lambda = 1, P_RS = P_start) towards P_RS = P_stop, stopping
        early when the real pressure reaches P_stop (``dcci.stop_at_target_pressure``).
        Backward: from where the forward sweep of the same iteration ended, the
        scaled pressure ramped back to P_start over the same number of blocks.

        Writes ``dcci.<direction>_<iteration>.dat`` (one row per block) and
        ``conf.dcci.<direction>_<iteration>.data`` in each cell folder.

        Returns
        -------
        SweepState
            lambda, scaled pressure and number of blocks at the end.
        """
        if direction == "forward":
            lam = 1.0
            p_rs_start, p_rs_stop = self.p_start, self.p_stop
            self.n_blocks_current = self.n_blocks
            confs = {
                phase: os.path.join(self.jobs[phase].simfolder, "conf.equilibration.data")
                for phase in self.CELLS
            }
        else:
            prev = self.states[(iteration, "forward")]
            lam = prev.lam
            p_rs_start, p_rs_stop = prev.p_rs, self.p_start
            self.n_blocks_current = prev.n_blocks
            confs = {
                phase: os.path.join(
                    self.jobs[phase].simfolder, "conf.dcci.forward_%d.data" % iteration
                )
                for phase in self.CELLS
            }

        n_blocks = self.n_blocks_current
        h_bar = (p_rs_stop - p_rs_start) / n_blocks
        h = h_bar / EV_A3_TO_BAR
        self.logger.info(
            "%s sweep (iteration %d): P_RS %.2f -> %.2f bar in %d blocks of %d steps, "
            "starting at lambda = %.6f (T = %.2f K)"
            % (direction, iteration, p_rs_start, p_rs_stop, n_blocks, self.n_block,
               lam, self.t0 / lam)
        )

        lmps = {
            phase: self._open_cell(phase, lam, p_rs_start, confs[phase], direction, iteration)
            for phase in self.CELLS
        }
        for phase in self.CELLS:
            self._arm_sweep(phase, lmps[phase], p_rs_start, p_rs_stop, direction, iteration)

        rows = []
        f_prev = None
        towards_higher = p_rs_stop > p_rs_start
        for k in range(n_blocks):
            p_rs_k = p_rs_start + k * h_bar
            p_rs_k1 = p_rs_k + h_bar
            # predictor: continue the last measured slope through this block
            lam_pred = lam if f_prev is None else lam + h * f_prev
            for phase in self.CELLS:
                lmps[phase].command(
                    lambda_ramp_command("lam", lam, lam_pred, k * self.n_block, self.n_block)
                )
                lmps[phase].command(
                    "run               %d start 0 stop %d" % (self.n_block, n_blocks * self.n_block)
                )
                lmps[phase].sync()

            u_s, v_s, lam_s = self._read_block("solid", lmps["solid"], direction, iteration, k)
            u_l, v_l, lam_l = self._read_block("liquid", lmps["liquid"], direction, iteration, k)
            lam_mid = 0.5 * (lam_s + lam_l)
            f_k = cce_slope(u_s, u_l, v_s, v_l)
            if self.calc.dcci.integrator == "trapezoid" and f_prev is not None:
                lam_new = lam + h * 0.5 * (f_k + f_prev)
            else:
                lam_new = lam + h * f_k
            if not lam_new > 0:
                raise RuntimeError(
                    "dcci %s sweep: lambda became non-positive (%g) in block %d; the "
                    "cells are no longer at coexistence" % (direction, lam_new, k)
                )

            t_mid = self.t0 / lam_mid
            p_mid = 0.5 * (p_rs_k + p_rs_k1) / lam_mid
            dpdt = clausius_clapeyron_slope(u_s, u_l, v_s, v_l, p_mid, t_mid)
            t_new = self.t0 / lam_new
            p_new = p_rs_k1 / lam_new
            rows.append([k, (k + 1) * self.n_block, lam_new, t_new, p_rs_k1, p_new,
                         u_s, u_l, v_s, v_l, dpdt])
            self.logger.info(
                "%s block %d/%d: T = %.2f K, P = %.1f bar (P_RS %.1f, lambda %.6f), "
                "du = %.4f eV/atom, dv = %.4f A^3/atom, dP/dT = %.2f bar/K"
                % (direction, k + 1, n_blocks, t_new, p_new, p_rs_k1, lam_new,
                   u_l - u_s, v_l - v_s, dpdt)
            )
            lam, f_prev = lam_new, f_k

            n_check = self.calc.dcci.n_check_blocks
            if n_check > 0 and (k + 1) % n_check == 0 and k + 1 < n_blocks:
                self._check_phases(lmps, direction, iteration, "block%d" % (k + 1))

            if (
                direction == "forward"
                and self.calc.dcci.stop_at_target_pressure
                and ((towards_higher and p_new >= self.p_stop)
                     or (not towards_higher and p_new <= self.p_stop))
            ):
                self.logger.info(
                    "forward sweep: real pressure %.1f bar reached the target %.1f bar "
                    "after %d of %d blocks; stopping" % (p_new, self.p_stop, k + 1, n_blocks)
                )
                n_blocks = k + 1
                break

        self._write_driver_file(direction, iteration, rows)

        for phase in self.CELLS:
            lmps[phase].command("unfix             fav")
            if self.calc.n_print_steps > 0:
                lmps[phase].command("undump            d1")
        self._check_phases(lmps, direction, iteration, "end")
        for phase in self.CELLS:
            conf = os.path.join(
                self.jobs[phase].simfolder, "conf.dcci.%s_%d.data" % (direction, iteration)
            )
            lmps[phase].command("write_data        %s" % conf)
            self.jobs[phase].lammps_close(lmp=lmps[phase])
            lmps[phase].rotate_logs("dcci_%s_%d" % (direction, iteration))

        state = SweepState(lam, p_rs_start + n_blocks * h_bar, n_blocks, direction)
        self.states[(iteration, direction)] = state
        self.logger.info(
            "%s sweep (iteration %d) done: T = %.2f K, P = %.1f bar"
            % (direction, iteration, self.t0 / state.lam, state.p_rs / state.lam)
        )
        return state

    def _write_driver_file(self, direction, iteration, rows):
        np.savetxt(
            self._driver_file(direction, iteration),
            np.array(rows),
            fmt=["%d", "%d", "%.10f", "%.6f", "%.6f", "%.6f", "%.8f", "%.8f", "%.6f", "%.6f", "%.6f"],
            header="dynamic Clausius-Clapeyron %s sweep, iteration %d; T0 = %.4f K\n"
            "block step lambda T[K] P_RS[bar] P[bar] u_solid[eV/atom] u_liquid[eV/atom] "
            "v_solid[A^3/atom] v_liquid[A^3/atom] dPdT[bar/K]" % (direction, iteration, self.t0),
        )

    # ------------------------------------------------------------- top level
    def calculate_coexistence_line(self):
        """Run the whole mode: equilibrate both cells, sweep, integrate, report."""
        self.prepare_cells()
        self.equilibrate_cells()
        for i in range(1, self.calc.n_iterations + 1):
            self.run_sweep("forward", iteration=i)
        raise NotImplementedError("mode dcci: backward sweep and integration follow in part 3")
