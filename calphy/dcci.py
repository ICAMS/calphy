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
import shutil
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import yaml

import calphy.helpers as ph
from calphy.input import Calculation, generate_metadata
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


def integrate_dcci(simfolder, t0, p_start, nsims=1, return_values=False):
    """
    Combine the forward and backward sweeps into one coexistence line.

    The forward sweep of every iteration provides the pressure grid (with the
    starting point prepended); the backward sweep of the same iteration is
    interpolated onto it.  The line is the mean of both directions, its error
    half the local hysteresis combined with the standard error over the
    iterations, as ``integrate_rs`` does for the temperature sweep.

    Parameters
    ----------
    simfolder : str
        Folder holding ``dcci.forward_<i>.dat`` and ``dcci.backward_<i>.dat``.
    t0, p_start : float
        The starting coexistence point (K, bar).
    nsims : int
        Number of iterations.
    return_values : bool
        Also return the arrays.

    Returns
    -------
    hysteresis : float
        Mean over iterations of T_backward(P_start) - T0, in K: how far the
        round trip misses the starting temperature.
    (pressure, temperature, error, t_forward, t_backward) : arrays
        Only if ``return_values``; on the grid of the first iteration.
    """
    grid = None
    forwards, backwards, hyst = [], [], []
    for i in range(1, nsims + 1):
        fwd = np.atleast_2d(np.loadtxt(os.path.join(simfolder, "dcci.forward_%d.dat" % i)))
        bwd = np.atleast_2d(np.loadtxt(os.path.join(simfolder, "dcci.backward_%d.dat" % i)))
        # real pressure and temperature, starting point included
        p_f = np.concatenate(([p_start], fwd[:, 5]))
        t_f = np.concatenate(([t0], fwd[:, 3]))
        # the backward sweep starts where the forward one ended
        p_b = np.concatenate(([fwd[-1, 5]], bwd[:, 5]))
        t_b = np.concatenate(([fwd[-1, 3]], bwd[:, 3]))
        order = np.argsort(p_b)
        if grid is None:
            grid = p_f
        t_f_grid = np.interp(grid, p_f, t_f) if p_f[-1] >= p_f[0] else np.interp(grid[::-1], p_f[::-1], t_f[::-1])[::-1]
        t_b_grid = np.interp(grid, p_b[order], t_b[order])
        forwards.append(t_f_grid)
        backwards.append(t_b_grid)
        hyst.append(float(np.interp(p_start, p_b[order], t_b[order]) - t0))

    forwards = np.array(forwards)
    backwards = np.array(backwards)
    means = 0.5 * (forwards + backwards)
    temperature = means.mean(axis=0)
    half_hyst = 0.5 * np.abs(forwards - backwards).mean(axis=0)
    stderr = means.std(axis=0, ddof=1) / np.sqrt(nsims) if nsims > 1 else np.zeros_like(temperature)
    error = np.sqrt(half_hyst**2 + stderr**2)
    hysteresis = float(np.mean(hyst))

    np.savetxt(
        os.path.join(simfolder, "coexistence_line.dat"),
        np.column_stack((grid, temperature, error, forwards.mean(axis=0), backwards.mean(axis=0))),
        header="dynamic Clausius-Clapeyron coexistence line\n"
        "pressure[bar]  temperature[K]  error[K]  T_forward[K]  T_backward[K]",
    )
    if return_values:
        return hysteresis, (grid, temperature, error, forwards.mean(axis=0), backwards.mean(axis=0))
    return hysteresis


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

        # The cells exchange information once per block through this driver;
        # in between they are independent and can run at the same time.
        cores = int(self.calc.queue.cores)
        self.parallel = bool(self.calc.dcci.parallel_cells) and cores >= 2
        if self.parallel:
            self.cell_cores = cores // 2
            self.logger.info(
                "dcci.parallel_cells: the solid and the liquid cell run concurrently "
                "with %d cores each (queue.cores = %d)" % (self.cell_cores, cores)
            )
            if cores % 2:
                self.logger.warning("queue.cores = %d is odd; one core stays idle" % cores)
        else:
            self.cell_cores = cores
            if self.calc.dcci.parallel_cells:
                self.logger.warning(
                    "dcci.parallel_cells needs queue.cores >= 2; running the cells sequentially"
                )
            self.logger.info(
                "The solid and the liquid cell run one after the other, each with %d cores" % cores
            )

        self.jobs = {}
        self.states = {}

    def _for_cells(self, fn):
        """
        Call ``fn(phase)`` for the solid and the liquid cell and return
        ``{phase: result}``; concurrently in two threads when
        ``dcci.parallel_cells`` is on.  An exception in either cell propagates.
        """
        if not self.parallel:
            return {phase: fn(phase) for phase in self.CELLS}
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = {phase: pool.submit(fn, phase) for phase in self.CELLS}
            return {phase: fut.result() for phase, fut in futures.items()}

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
        data["queue"]["cores"] = self.cell_cores
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
        def equilibrate(phase):
            self.logger.info("Equilibrating the %s cell at %.2f K, %.2f bar" % (phase, self.t0, self.p_start))
            self.jobs[phase].run_averaging()
            self.logger.info("%s cell: %.4f A^3/atom" % (phase, self.jobs[phase].volatom))

        self._for_cells(equilibrate)

    # ------------------------------------------------------------------ sweep
    def _open_cell(self, phase, lam, p_rs, conf, direction, iteration, seed):
        """
        Open a LAMMPS session for one cell with the lambda-scaled potential at
        constant ``lam`` and a barostat at the constant scaled pressure ``p_rs``,
        and equilibrate it there.  ``seed`` (drawn by the caller, so that the
        random stream is used from one thread only) seeds the velocities of a
        forward sweep.
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
                % (self.t0, seed)
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

        def check(phase):
            job = self.jobs[phase]
            snapshot = "traj.dcci.%s_%d.%s.dat" % (direction, iteration, tag)
            job.dump_current_snapshot(lmps[phase], snapshot)
            if phase == "solid":
                job.check_if_melted(lmps[phase], snapshot)
            else:
                job.check_if_solidfied(lmps[phase], snapshot)

        self._for_cells(check)

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

        seeds = {phase: int(np.random.randint(1, 10000)) for phase in self.CELLS}
        lmps = self._for_cells(
            lambda phase: self._open_cell(
                phase, lam, p_rs_start, confs[phase], direction, iteration, seeds[phase]
            )
        )
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
            ramp = lambda_ramp_command("lam", lam, lam_pred, k * self.n_block, self.n_block)
            run = "run               %d start 0 stop %d" % (self.n_block, n_blocks * self.n_block)

            def block(phase):
                lmps[phase].command(ramp)
                lmps[phase].command(run)
                lmps[phase].sync()

            self._for_cells(block)

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

        def finish(phase):
            conf = os.path.join(
                self.jobs[phase].simfolder, "conf.dcci.%s_%d.data" % (direction, iteration)
            )
            lmps[phase].command("write_data        %s" % conf)
            self.jobs[phase].lammps_close(lmp=lmps[phase])
            lmps[phase].rotate_logs("dcci_%s_%d" % (direction, iteration))

        self._for_cells(finish)

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

    # ---------------------------------------------------------- integration
    def integrate(self):
        """Combine all sweeps into ``coexistence_line.dat`` and judge the hysteresis."""
        self.hysteresis, values = integrate_dcci(
            self.simfolder, self.t0, self.p_start,
            nsims=self.calc.n_iterations, return_values=True,
        )
        pressure, temperature, error = values[0], values[1], values[2]
        self.hysteresis_high = abs(self.hysteresis) > self.calc.dcci.hysteresis_tolerance
        self.line = (pressure, temperature, error)
        self.logger.info(
            "Coexistence line from %.1f to %.1f bar written to coexistence_line.dat "
            "(%d points)" % (pressure[0], pressure[-1], len(pressure))
        )
        self.logger.info(
            "Round trip misses the starting temperature by %.2f K (tolerance %.2f K)"
            % (self.hysteresis, self.calc.dcci.hysteresis_tolerance)
        )
        if self.hysteresis_high:
            self.logger.warning(
                "The forward and backward sweeps do not close: the integration is "
                "not reversible at this sweep length. Increase n_switching_steps "
                "(or the cell size) before trusting the line."
            )
            self.logger.warning("STATE: coexistence line unreliable, hysteresis too high")
        self.logger.info(
            "STATE: T_coex = %.2f K at %.1f bar (from %.2f K at %.1f bar)"
            % (temperature[-1], pressure[-1], self.t0, self.p_start)
        )

    def submit_report(self):
        """Write ``report.yaml`` at the top level."""
        pressure, temperature, error = self.line
        last = self.states[(self.calc.n_iterations, "forward")]
        report = {
            "input": {
                "temperature": float(self.t0),
                "pressure": float(self.p_start),
                "pressure_stop": float(self.p_stop),
                "lattice": str(self.calc._original_lattice),
                "element": " ".join(np.array(self.calc.element).astype(str)),
                "n_block_steps": int(self.n_block),
                "n_iterations": int(self.calc.n_iterations),
            },
            "average": {
                "vol_atom_solid": float(self.jobs["solid"].volatom),
                "vol_atom_liquid": float(self.jobs["liquid"].volatom),
            },
            "results": {
                "coexistence_line": "coexistence_line.dat",
                "pressure_reached": float(pressure[-1]),
                "temperature_at_pressure_reached": float(temperature[-1]),
                "error_at_pressure_reached": float(error[-1]),
                "n_blocks_used": int(last.n_blocks),
                "hysteresis": float(self.hysteresis),
                "hysteresis_high": bool(self.hysteresis_high),
                "unit": "K, bar",
            },
        }
        self.report = report
        with open(os.path.join(self.simfolder, "report.yaml"), "w") as fout:
            yaml.dump(report, fout)
        self.logger.info("Report written in %s" % os.path.join(self.simfolder, "report.yaml"))
        self.logger.info("Please cite the following publications:")
        for doi in self.publications:
            self.logger.info("- %s" % doi)

    def clean_up(self):
        """Serialise the input configuration and the run metadata."""
        shutil.copy(
            self.calc.lattice, os.path.join(self.simfolder, "input_configuration.data")
        )
        metadata = generate_metadata()
        metadata["publications"] = self.publications
        with open(os.path.join(self.simfolder, "metadata.yaml"), "w") as fout:
            yaml.safe_dump(metadata, fout)

    # ------------------------------------------------------------- top level
    def calculate_coexistence_line(self):
        """Run the whole mode: equilibrate both cells, sweep out and back, integrate, report."""
        self.prepare_cells()
        self.equilibrate_cells()
        for i in range(1, self.calc.n_iterations + 1):
            ts = time.time()
            self.run_sweep("forward", iteration=i)
            self.run_sweep("backward", iteration=i)
            self.logger.info("dcci cycle %d finished in %f s" % (i, time.time() - ts))
        self.integrate()
        self.submit_report()
        self.clean_up()
