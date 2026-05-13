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
from calphy.transition_detector import (
    PhaseTransitionMonitor,
    ThermoRecord,
)


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

        # phase transition monitor (created in _create_monitor)
        self._monitor = None
        self._monitor_fed_rows = 0

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

    def _create_monitor(self, expected_phase):
        """
        Instantiate a PhaseTransitionMonitor for the current calculation.

        Parameters
        ----------
        expected_phase : 'solid' or 'liquid'
        """
        td = self.calc.phase_transition_detection
        if td.mode == "none":
            self._monitor = None
            self._monitor_fed_rows = 0
            self.logger.info("Phase transition detector: disabled")
            return

        self.logger.info(
            "Phase transition detector: enabled for %s phase "
            "(T=%.1f K, P=%.3f bar, peak_threshold=%.1f, min_agreement=%d)",
            expected_phase,
            float(self.calc._temperature),
            float(self.calc._pressure) if self.calc._pressure is not None else 0.0,
            td.peak_threshold,
            td.min_agreement,
        )
        self._monitor = PhaseTransitionMonitor(
            expected_phase=expected_phase,
            target_pressure=float(self.calc._pressure) if self.calc._pressure is not None else 0.0,
            temperature=float(self.calc._temperature),
            baseline_window=td.baseline_window,
            recent_window=td.recent_window,
            min_samples_before_check=td.min_samples_before_check,
            peak_threshold=td.peak_threshold,
            min_signal_agreement=td.min_agreement,
        )
        self._monitor_fed_rows = 0

    def _feed_monitor_from_avg_dat(self, avg_file):
        """
        Read any new rows from avg.dat (columns: step, lx, ly, lz, press,
        pe/atom, etotal/atom, temp) and append them to the monitor buffer.

        avg.dat column layout (1-indexed in LAMMPS output, 0-indexed here):
          0: TimeStep
          1: v_mlx,  2: v_mly,  3: v_mlz
          4: v_mpress
          5: v_mpe      (pe/atoms)
          6: v_metotal  (etotal/atoms)
          7: v_mtemp    (temp)

        Returns the updated number of fed rows so the caller can track
        how many rows have already been processed.
        """
        if self._monitor is None:
            return

        if not os.path.exists(avg_file):
            return

        try:
            data = np.loadtxt(avg_file, usecols=(0, 1, 2, 3, 4, 5, 6, 7), unpack=False)
        except (ValueError, IndexError):
            # File may not yet have enough columns (e.g. first run of older code path)
            return

        if data.ndim == 1:
            data = data.reshape(1, -1)

        n_rows = data.shape[0]
        new_rows = data[self._monitor_fed_rows:]
        self._monitor_fed_rows = n_rows

        if len(new_rows) > 0:
            self.logger.info(
                "Phase transition detector: fed %d new rows from %s "
                "(total buffer size: %d)",
                len(new_rows),
                os.path.basename(avg_file),
                n_rows,
            )

        for row in new_rows:
            step, lx, ly, lz, press, pe, etotal, temp = row
            vol_per_atom = (lx * ly * lz) / self.natoms
            record = ThermoRecord(
                step=step,
                temp=temp,
                pe=pe,
                etotal=etotal,
                vol=vol_per_atom,
                press=press,
            )
            # update() returns a TransitionEvent but we don't raise here;
            # we let evaluate_final() do that at the checkpoint call sites
            self._monitor.update(record, raise_on_transition=False)

    def check_if_melted(self, lmp, filename):
        """
        Check whether the solid has melted, using fluctuation-based detection.

        The monitor buffer is evaluated using the data accumulated during the
        preceding pressure-convergence cycles.  If the buffer is too small for
        reliable detection (fewer than min_samples_before_check samples), the
        check is silently skipped — consistent with the previous behaviour of
        setting tolerance.solid_fraction: 0.
        """
        if self._monitor is None:
            return

        buf_size = len(self._monitor.buffer)
        self.logger.info(
            "Phase transition detector: checking solid — buffer has %d samples",
            buf_size,
        )
        event = self._monitor.evaluate_final(raise_on_transition=False)
        if event is None:
            self.logger.info(
                "Phase transition detector: solid phase stable (%d samples analysed)",
                buf_size,
            )
        else:
            self.logger.warning(
                "Phase transition detector: MELTING detected! "
                "triggered_signals=%s confidence=%.0f%%",
                event.triggered_signals,
                event.confidence * 100,
            )
            self.lammps_close(lmp=lmp)
            logfile = os.path.join(self.simfolder, "log.lammps")
            try:
                os.rename(
                    logfile, os.path.join(self.simfolder, "melted_error.log.lammps")
                )
            except OSError as e:
                self.logger.warning(f"Failed to rename log file: {e}")
            self._monitor._raise_for_event(event)

    def check_if_solidfied(self, lmp, filename):
        """
        Check whether the liquid has solidified, using fluctuation-based detection.

        See check_if_melted for details on behaviour when the buffer is small.
        """
        if self._monitor is None:
            return

        buf_size = len(self._monitor.buffer)
        self.logger.info(
            "Phase transition detector: checking liquid — buffer has %d samples",
            buf_size,
        )
        event = self._monitor.evaluate_final(raise_on_transition=False)
        if event is None:
            self.logger.info(
                "Phase transition detector: liquid phase stable (%d samples analysed)",
                buf_size,
            )
        else:
            self.logger.warning(
                "Phase transition detector: SOLIDIFICATION detected! "
                "triggered_signals=%s confidence=%.0f%%",
                event.triggered_signals,
                event.confidence * 100,
            )
            self.lammps_close(lmp=lmp)
            logfile = os.path.join(self.simfolder, "log.lammps")
            try:
                os.rename(
                    logfile, os.path.join(self.simfolder, "solidified_error.log.lammps")
                )
            except OSError as e:
                self.logger.warning(f"Failed to rename log file: {e}")
            self._monitor._raise_for_event(event)

    def _detect_ts_transitions(self):
        """
        Analyse ts.forward / ts.backward files for phase transitions via
        rolling response-function peaks (Cp, kappa_T, alpha_P).

        Called from integrate_reversible_scaling after the MD run completes.
        Emits a Python UserWarning and logs at WARNING level when a transition
        is detected.  The calculation is not interrupted; the warning is
        informational so the user can inspect output files.
        """
        from calphy.transition_detector import detect_ts_transitions as _dts

        td      = self.calc.phase_transition_detection
        t_start = float(self.calc._temperature)
        t_stop  = float(self.calc._temperature_stop)

        self.logger.info(
            "ts-sweep transition detector: scanning %d iteration(s) "
            "(T %.1f K -> %.1f K, peak_threshold=%.1f, min_agreement=%d)",
            self.calc.n_iterations,
            t_start, t_stop,
            td.peak_threshold,
            td.min_agreement,
        )

        for i in range(1, self.calc.n_iterations + 1):
            for sweep_label, fname in [
                ("forward",  f"ts.forward_{i}.dat"),
                ("backward", f"ts.backward_{i}.dat"),
            ]:
                fpath = os.path.join(self.simfolder, fname)
                if not os.path.exists(fpath):
                    self.logger.debug(
                        "ts-sweep detector: %s not found, skipping", fname
                    )
                    continue
                try:
                    dU, press, vol_total, lam = np.loadtxt(
                        fpath, unpack=True, comments="#"
                    )
                except Exception as exc:
                    self.logger.debug(
                        "ts-sweep detector: could not read %s (%s)", fname, exc
                    )
                    continue

                label = f"{sweep_label} (iteration {i})"
                self.logger.info(
                    "ts-sweep transition detector: analysing %s (%d rows)",
                    label, len(dU),
                )
                events = _dts(
                    dU=dU,
                    press=press,
                    vol_total=vol_total,
                    lam=lam,
                    t_start=t_start,
                    t_stop=t_stop,
                    natoms=self.natoms,
                    peak_threshold=td.peak_threshold,
                    min_signal_agreement=td.min_agreement,
                    sweep_label=label,
                    sweep_mode="ts",
                )

                if not events:
                    self.logger.info(
                        "ts-sweep transition detector: no transition detected in %s",
                        label,
                    )
                for event in events:
                    self.logger.warning(
                        "ts-sweep transition detector: transition in %s at "
                        "T ~ %.1f K. Signals: %s (confidence %.0f%%). "
                        "Verify output files for structural changes.",
                        label,
                        event.temperature,
                        ", ".join(event.triggered_signals),
                        event.confidence * 100,
                    )

    def _plot_ts_response_functions(self):
        """
        Plot Cp, kappa_T, alpha_P vs temperature for every ts sweep file and
        save the figures to the simulation folder.

        Called from integrate_reversible_scaling after the MD run completes.
        A PNG file named ``ts_response_<sweep>_<iteration>.png`` is written for
        each forward/backward sweep.  Detected transition temperatures are
        marked with vertical dashed lines.
        """
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            self.logger.warning(
                "matplotlib not available; skipping response function plots"
            )
            return

        import warnings as _warnings
        from calphy.transition_detector import (
            compute_ts_response_arrays as _compute,
            detect_ts_transitions as _dts,
        )

        td      = self.calc.phase_transition_detection
        t_start = float(self.calc._temperature)
        t_stop  = float(self.calc._temperature_stop)

        for i in range(1, self.calc.n_iterations + 1):
            for sweep_label, fname in [
                ("forward",  f"ts.forward_{i}.dat"),
                ("backward", f"ts.backward_{i}.dat"),
            ]:
                fpath = os.path.join(self.simfolder, fname)
                if not os.path.exists(fpath):
                    continue
                try:
                    dU, press, vol_total, lam = np.loadtxt(
                        fpath, unpack=True, comments="#"
                    )
                except Exception as exc:
                    self.logger.debug(
                        "response plot: could not read %s (%s)", fname, exc
                    )
                    continue

                arrs = _compute(
                    dU, press, vol_total, lam, t_start, t_stop, self.natoms,
                    sweep_mode="ts",
                )
                vs = arrs["valid_start"]
                if vs >= len(arrs["T"]):
                    self.logger.debug(
                        "response plot: %s too short for rolling windows, skipping",
                        fname,
                    )
                    continue

                T      = arrs["T"][vs:]
                Cp     = arrs["Cp"][vs:]
                kappa  = arrs["kappa_T"][vs:]
                alpha  = arrs["alpha_P"][vs:]

                # Detect transitions to mark on the plot
                with _warnings.catch_warnings():
                    _warnings.simplefilter("ignore", UserWarning)
                    events = _dts(
                        dU=dU,
                        press=press,
                        vol_total=vol_total,
                        lam=lam,
                        t_start=t_start,
                        t_stop=t_stop,
                        natoms=self.natoms,
                        peak_threshold=td.peak_threshold,
                        min_signal_agreement=td.min_agreement,
                        sweep_label=sweep_label,
                        sweep_mode="ts",
                    )

                fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
                fig.suptitle(
                    f"Thermodynamic response functions\n"
                    f"{sweep_label} sweep, iteration {i}  "
                    f"(T: {t_start:.0f} K \u2192 {t_stop:.0f} K)",
                    fontsize=11,
                )

                def _modz_series(sig):
                    """Robust mod-Z relative to the global median/MAD."""
                    finite = sig[np.isfinite(sig)]
                    if len(finite) < 10:
                        return np.full_like(sig, np.nan)
                    med = float(np.median(finite))
                    mad = float(np.median(np.abs(finite - med))) * 1.4826
                    if mad < 1e-40:
                        return np.full_like(sig, np.nan)
                    return (sig - med) / mad

                def _add_modz_twin(ax, T_arr, sig, color, threshold):
                    """Overlay mod-Z on a twin y-axis (right side, grey/dashed)."""
                    mz = _modz_series(sig)
                    ax2 = ax.twinx()
                    ax2.plot(T_arr, mz, lw=0.6, color=color, alpha=0.35, ls="--")
                    ax2.axhline(threshold, color="grey", lw=0.8, ls=":", alpha=0.7)
                    ax2.set_ylabel("mod-Z", fontsize=7, color="grey")
                    ax2.tick_params(axis="y", labelsize=6, labelcolor="grey")
                    # Keep the primary axis scale from being dominated by outliers
                    mz_finite = mz[np.isfinite(mz)]
                    if len(mz_finite):
                        ax2.set_ylim(bottom=0,
                                     top=max(threshold * 2, float(np.nanpercentile(mz_finite, 99)) * 1.1))
                    return ax2

                threshold = td.peak_threshold

                axes[0].plot(T, Cp, lw=0.8, color="C0")
                axes[0].set_ylabel("$C_P$ (eV K$^{-1}$ atom$^{-1}$)")
                _add_modz_twin(axes[0], T, Cp, "C0", threshold)

                axes[1].plot(T, kappa, lw=0.8, color="C1")
                axes[1].set_ylabel(r"$\kappa_T$ (eV$^{-1}$)")
                _add_modz_twin(axes[1], T, kappa, "C1", threshold)

                axes[2].plot(T, alpha, lw=0.8, color="C2")
                axes[2].set_ylabel(r"$\alpha_P$ (K$^{-1}$)")
                axes[2].set_xlabel("Temperature (K)")
                _add_modz_twin(axes[2], T, alpha, "C2", threshold)

                for event in events:
                    for ax in axes:
                        ax.axvline(
                            event.temperature,
                            color="red",
                            lw=1.2,
                            ls="--",
                            alpha=0.8,
                            label=f"transition ~{event.temperature:.0f} K",
                        )
                    axes[0].legend(fontsize=8, loc="upper left")

                plt.tight_layout()
                outname = f"ts_response_{sweep_label}_{i}.png"
                outpath = os.path.join(self.simfolder, outname)
                try:
                    fig.savefig(outpath, dpi=150)
                    self.logger.info(
                        "Response function plot saved: %s", outpath
                    )
                except Exception as exc:
                    self.logger.warning(
                        "Could not save response plot %s: %s", outpath, exc
                    )
                finally:
                    plt.close(fig)

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

        Full-length NPT cycles are run throughout.  The first cycle is excluded
        from the running average to discard the initial transient.  After
        ``n_fit_warmup`` cycles (default 5), a linear P-V fit is used to
        predict the equilibrium volume and rescale the box, accelerating
        convergence.  The mean pressure (over all post-transient data) must
        fall within ``tolerance.pressure`` of the target to declare convergence.
        """

        # apply fixes
        if self.calc.equilibration_control == "nose-hoover":
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

            # feed any new avg.dat rows into the phase-transition monitor
            self._feed_monitor_from_avg_dat(file)

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
        if self.calc.equilibration_control == "nose-hoover":
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
        if self.calc.equilibration_control == "nose-hoover":
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

    def _run_split_sweep(
        self,
        lmp,
        lambda_var: str,
        output_file_pattern: str,
        sweep_label: str,
        t0: float,
        tf: float,
        iteration: int,
        sweep_mode: str = "ts",
    ) -> None:
        """
        Run a forward or backward reversible-scaling sweep split into
        temperature-based blocks.

        Uses a **single** ``fix print ... flush yes`` output file that grows
        continuously across all blocks; no per-cycle files are written and no
        post-run merging is needed.  Between consecutive ``run`` commands,
        LAMMPS is idle and the file is fully written to disk (guaranteed by
        ``flush yes``), allowing the incremental transition detector to read the
        accumulated data with the correct full-sweep baseline.

        Parameters
        ----------
        lmp : lammps object
            Active LAMMPS instance.  The pair style and the lambda ramp
            variables must already be defined *before* this method is called.
            The ramp *will* be reset (``reset_timestep 0``) here so that
            ``run N start 0 stop total_steps`` works correctly.
        lambda_var : str
            Name of the LAMMPS variable to record, e.g. ``"flambda"`` or
            ``"blambda"``.
        output_file_pattern : str
            Name for the output data file, e.g. ``"ts.forward_1.dat"``.
        sweep_label : str
            Human-readable label used in log messages.
        t0, tf : float
            Temperature endpoints of the sweep (K) used for the block planner.
        iteration : int
            Current reversible-scaling iteration number.
        sweep_mode : str
            ``'ts'`` (reversible scaling, fixed thermostat) or ``'tscale'``
            (temperature scaling, ramping thermostat).  Used for the correct
            temperature formula in the transition detector.
        """
        from calphy.transition_detector import plan_temperature_blocks as _plan

        td = self.calc.phase_transition_detection
        n_sweep = self.calc._n_sweep_steps
        t_win = td.temperature_window if (td.mode != "none" and td.temperature_window > 0) else 0.0

        if t_win <= 0:
            # ── Single-run path (no block splitting) ──────────────────────
            lmp.command(
                'fix               f3 all print 1 "${dU} $(press) $(vol) ${%s}" '
                'screen no file %s' % (lambda_var, output_file_pattern)
            )
            self.logger.info("ts-sweep %s: single run, %d steps", sweep_label, n_sweep)
            lmp.command("run               %d" % n_sweep)
            lmp.command("unfix             f3")
            return

        # ── Block-splitting path ──────────────────────────────────────────
        blocks = _plan(t0, tf, n_sweep, t_win)
        n_blocks = len(blocks) - 1
        self.logger.info(
            "ts-sweep %s: %d steps split into %d blocks of ~%.0f K "
            "(T %.1f → %.1f K, λ %.4f → %.4f)",
            sweep_label, n_sweep, n_blocks, t_win,
            t0, tf, blocks[0]["lambda"], blocks[-1]["lambda"],
        )
        for k, cp in enumerate(blocks):
            self.logger.info(
                "  checkpoint %d: T=%.2f K  λ=%.6f  step=%d",
                k, cp["temp"], cp["lambda"], cp["step"],
            )

        # Reset timestep so ramp(li,lf) + "run N start 0 stop total" is
        # consistent across all block commands.
        lmp.command("reset_timestep 0")

        # Strategy: run the full sweep to completion, writing a per-block
        # checkpoint at every block boundary so post-hoc recovery (in
        # _post_forward_recovery) can restart the backward sweep from any
        # safe boundary if a transition is detected.  No incremental
        # detection / no aborts inside this method — that decision is made
        # once, post-hoc, on the complete forward dataset.
        fpath = os.path.join(self.simfolder, output_file_pattern)
        # Strip ".dat" so we can build a per-block checkpoint name.
        out_stem = output_file_pattern[:-4] if output_file_pattern.endswith(".dat") else output_file_pattern
        current_step = 0
        for k in range(n_blocks):
            block_steps = blocks[k + 1]["step"] - blocks[k]["step"]
            if block_steps <= 0:
                self.logger.debug(
                    "ts-sweep %s block %d has 0 steps — skipping", sweep_label, k
                )
                continue

            self.logger.info(
                "ts-sweep %s block %d/%d: T %.1f→%.1f K  λ %.4f→%.4f  steps %d–%d",
                sweep_label, k, n_blocks - 1,
                blocks[k]["temp"], blocks[k + 1]["temp"],
                blocks[k]["lambda"], blocks[k + 1]["lambda"],
                current_step, current_step + block_steps,
            )

            # Snapshot the state at the START of this block.  Used by
            # _post_forward_recovery as a candidate restart point if a
            # transition is detected later in the sweep.
            # _post_forward_recovery as a candidate restart point if a
            # transition is detected later in the sweep.
            chk_name = "conf.%s_blk%d.data" % (out_stem, k)
            chk_path = os.path.join(self.simfolder, chk_name)
            lmp.command("write_data %s nocoeff" % chk_path)

            if k == 0:
                lmp.command(
                    'fix f3 all print 1 "${dU} $(press) $(vol) ${%s}" '
                    'screen no file %s' % (lambda_var, output_file_pattern)
                )
            else:
                lmp.command(
                    'fix f3 all print 1 "${dU} $(press) $(vol) ${%s}" '
                    'screen no append %s' % (lambda_var, output_file_pattern)
                )

            lmp.command("run %d start 0 stop %d" % (block_steps, n_sweep))
            current_step += block_steps
            lmp.command("unfix f3")

        # Handle any remaining steps (numerical rounding)
        remaining = n_sweep - current_step
        if remaining > 0:
            self.logger.info(
                "ts-sweep %s remainder block: %d steps", sweep_label, remaining
            )
            lmp.command(
                'fix f3 all print 1 "${dU} $(press) $(vol) ${%s}" '
                'screen no append %s' % (lambda_var, output_file_pattern)
            )
            lmp.command("run %d start 0 stop %d" % (remaining, n_sweep))
            lmp.command("unfix f3")

        self.logger.info("ts-sweep %s: sweep complete (%d steps)", sweep_label, n_sweep)

    def _reversible_scaling_forward(self, iteration: int = 1) -> None:
        """
        Perform the forward sweep of a reversible-scaling calculation.

        1. Initial NPT equilibration at T0 (records thermodynamics for the
           online phase-transition monitor).
        2. COM-constrained equilibration at T0.
        3. Forward sweep: λ 1 → T0/Tf (optionally split into temperature
           windows).
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
            script_mode=False,
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

        # ── Phase-transition monitor variables ────────────────────────────
        lmp.command("variable         rs_mlx equal lx")
        lmp.command("variable         rs_mly equal ly")
        lmp.command("variable         rs_mlz equal lz")
        lmp.command("variable         rs_mpress equal press")
        lmp.command("variable         rs_mpe equal pe/atoms")
        lmp.command("variable         rs_metotal equal etotal/atoms")
        lmp.command("variable         rs_mtemp equal temp")

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
        _rs_ave_every  = int(self.calc.md.n_every_steps)
        _rs_ave_repeat = int(self.calc.md.n_repeat_steps)
        _rs_ave_freq   = _rs_ave_every * _rs_ave_repeat
        lmp.command(
            "fix               rs_ave all ave/time %d %d %d "
            "v_rs_mlx v_rs_mly v_rs_mlz v_rs_mpress v_rs_mpe v_rs_metotal v_rs_mtemp "
            "file rs_equil_avg.dat"
            % (_rs_ave_every, _rs_ave_repeat, _rs_ave_freq)
        )
        lmp.command("run               %d" % self.calc.n_equilibration_steps)
        lmp.command("unfix             rs_ave")
        self.logger.info("forward sweep (iteration %d): initial equilibration done", iteration)

        rs_avg_file = os.path.join(self.simfolder, "rs_equil_avg.dat")
        self._monitor_fed_rows = 0
        self._feed_monitor_from_avg_dat(rs_avg_file)

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

        # ── Define scaling variables ────────────────────────────────────────
        lmp.command("variable         flambda equal ramp(${li},${lf})")
        lmp.command("variable         blambda equal ramp(${lf},${li})")
        lmp.command("variable         fscale equal v_flambda-1.0")
        lmp.command("variable         bscale equal v_blambda-1.0")
        lmp.command("variable         one equal 1.0")
        lmp.command(
            "variable        ftemp equal v_blambda*%f" % self.calc._temperature_stop
        )
        lmp.command(
            "variable        btemp equal v_flambda*%f" % self.calc._temperature_stop
        )

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
            self._run_split_sweep(
                lmp=lmp,
                lambda_var="flambda",
                output_file_pattern="ts.forward_%d.dat" % iteration,
                sweep_label="forward (iteration %d)" % iteration,
                t0=t0,
                tf=tf,
                iteration=iteration,
            )
        except Exception:
            # Make sure the LAMMPS instance (especially in interactive
            # pylammpsmpi mode where it is reused for the backward sweep) is
            # cleared and the log file is renamed before the exception
            # propagates.  This is the only way the recovery path in
            # reversible_scaling() can safely create a fresh lmp for the
            # backward sweep.
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
        2. Middle equilibration at Tf (phase-stability check).
        3. Backward sweep: λ T0/Tf → 1 (optionally split into temperature
           windows).

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

        # ── Middle equilibration at Tf ──────────────────────────────────────
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

        # ── Switch to scaling potential ─────────────────────────────────────
        lmp = ph.set_potential(lmp, self.calc)

        lmp.command("variable         flambda equal ramp(${li},${lf})")
        lmp.command("variable         blambda equal ramp(${lf},${li})")
        lmp.command("variable         fscale equal v_flambda-1.0")
        lmp.command("variable         bscale equal v_blambda-1.0")
        lmp.command("variable         one equal 1.0")
        lmp.command(
            "variable        ftemp equal v_blambda*%f" % self.calc._temperature_stop
        )
        lmp.command(
            "variable        btemp equal v_flambda*%f" % self.calc._temperature_stop
        )

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
        self._run_split_sweep(
            lmp=lmp,
            lambda_var="blambda",
            output_file_pattern="ts.backward_%d.dat" % iteration,
            sweep_label="backward (iteration %d)" % iteration,
            t0=tf,
            tf=t0,
            iteration=iteration,
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

        self.lammps_close(lmp=lmp)

        logfile = os.path.join(self.simfolder, "log.lammps")
        if os.path.exists(logfile):
            os.rename(
                logfile,
                os.path.join(
                    self.simfolder, "reversible_scaling_backward.log.lammps"
                ),
            )

    def _truncate_ts_file(self, fpath: str, n_data_rows: int) -> int:
        """
        Keep only the first ``n_data_rows`` data rows of ``fpath`` (rewriting
        the file in place).  Comment lines starting with '#' are preserved.

        Returns the number of data rows actually retained (clamped to the
        rows currently present in the file).
        """
        if not os.path.exists(fpath):
            return 0
        with open(fpath, "r") as fh:
            lines = fh.readlines()
        kept: list = []
        kept_data = 0
        for line in lines:
            if line.startswith("#") or line.strip() == "":
                kept.append(line)
                continue
            if kept_data < n_data_rows:
                kept.append(line)
                kept_data += 1
            else:
                break
        with open(fpath, "w") as fh:
            fh.writelines(kept)
        return kept_data

    def _post_forward_recovery(self, iteration: int) -> None:
        """
        Post-hoc transition handling after a successful forward sweep.

        Reads the full forward ts file, runs the transition detector on the
        complete dataset, and — if a transition is detected — picks the last
        per-block checkpoint that lies safely *before* the detected
        transition (in the sweep direction), truncates the forward ts file
        to that step, copies the checkpoint to the canonical
        ``conf.ts.forward_<iter>.data`` filename used by the backward sweep,
        and reduces ``calc._temperature_stop`` and ``calc._n_sweep_steps``
        so the subsequent backward sweep covers exactly the same lambda
        range.

        Why post-hoc instead of incremental?  Detecting a transition
        reliably requires the full peak shape (rise + peak + descent +
        post-transition plateau).  Per-block detection produces false
        positives at almost every block boundary because:
          - the rolling window is at the head of the data (startup transient)
          - or the rolling window is at the tail (only the rise of a real
            transition is visible, no descent yet)
          - or smooth physical Cp(T) growth across a single block looks
            like a step jump.
        Running the detector once on the complete dataset removes all of
        these failure modes and makes the implementation drastically
        simpler.

        No-op when ``td.mode == 'none'`` or no transition is found.

        Returns nothing; backward sweep is run unconditionally afterwards.
        """
        from calphy.transition_detector import (
            plan_temperature_blocks as _plan,
            detect_ts_transitions as _dts,
        )

        td = self.calc.phase_transition_detection
        if td.mode == "none":
            return

        out = os.path.join(self.simfolder, "ts.forward_%d.dat" % iteration)
        if not os.path.exists(out):
            self.logger.debug(
                "post-forward recovery: %s missing — skipping", out
            )
            return

        try:
            data = np.loadtxt(out, comments="#")
        except Exception as exc:
            self.logger.debug(
                "post-forward recovery: could not read %s (%s)", out, exc
            )
            return
        if data.ndim < 2 or data.shape[0] < 100:
            return

        dU, press, vol_total, lam = data.T
        t_start = float(self.calc._temperature)
        t_stop  = float(self.calc._temperature_stop)

        events = _dts(
            dU=dU, press=press, vol_total=vol_total, lam=lam,
            t_start=t_start, t_stop=t_stop, natoms=self.natoms,
            peak_threshold=td.peak_threshold,
            min_signal_agreement=td.min_agreement,
            sweep_label="forward (iteration %d) post-hoc" % iteration,
            sweep_mode="ts",
        )
        if not events:
            self.logger.info(
                "post-forward recovery: no phase transition detected in "
                "forward sweep (iteration %d)", iteration,
            )
            return

        ev = events[0]
        T_trans = float(ev.temperature)
        self.logger.warning(
            "post-forward recovery: phase transition detected at T ~ %.1f K "
            "(signals: %s, confidence %.0f%%)",
            T_trans, ", ".join(ev.triggered_signals), ev.confidence * 100,
        )

        if td.mode == "stop":
            from calphy.errors import PhaseTransitionError
            try:
                self._plot_ts_response_functions()
            except Exception as _plot_exc:
                self.logger.debug(
                    "post-forward recovery: plot failed before abort: %s",
                    _plot_exc,
                )
            raise PhaseTransitionError(
                "Phase transition detected in forward sweep at T ~ %.1f K "
                "(signals: %s)." % (T_trans, ", ".join(ev.triggered_signals))
            )

        if td.mode == "warn":
            # Warn-only: keep the full sweep, just log + plot.
            try:
                self._plot_ts_response_functions()
            except Exception as _plot_exc:
                self.logger.debug(
                    "post-forward recovery: plot failed: %s", _plot_exc
                )
            return

        # ── td.mode == 'recover' ──────────────────────────────────────────
        # Re-plan the same block structure used by _run_split_sweep so we
        # know which checkpoint files exist and what their step / T values
        # are.  The block planner is deterministic given (t0, tf, n_sweep,
        # window_K).
        n_sweep = self.calc._n_sweep_steps
        t_win = td.temperature_window if td.temperature_window > 0 else 0.0
        if t_win <= 0:
            self.logger.error(
                "post-forward recovery: cannot recover — temperature_window "
                "is 0 (no per-block checkpoints were written).  Re-run "
                "with a non-zero temperature_window."
            )
            return

        blocks = _plan(t_start, t_stop, n_sweep, t_win)
        n_blocks = len(blocks) - 1

        # Find the latest block whose end-temperature is safely before the
        # detected transition (in the sweep direction).  We use the start
        # of block k as the restart point — i.e. blocks[k]["temp"] strictly
        # before T_trans (heating: <, cooling: >).  Add a one-block safety
        # buffer so we restart at least one full window away from the peak.
        heating = t_stop > t_start
        if heating:
            safe_candidates = [
                k for k in range(n_blocks + 1)
                if blocks[k]["temp"] < T_trans
            ]
        else:
            safe_candidates = [
                k for k in range(n_blocks + 1)
                if blocks[k]["temp"] > T_trans
            ]
        if not safe_candidates:
            self.logger.error(
                "post-forward recovery: no checkpoint lies before the "
                "detected transition (T_trans=%.1f K) — cannot recover.  "
                "Lower the starting temperature and re-run.",
                T_trans,
            )
            return

        # One-block buffer: drop the last candidate (closest to T_trans).
        safe_k = safe_candidates[-1]
        if len(safe_candidates) >= 2:
            safe_k = safe_candidates[-2]

        if safe_k <= 0:
            self.logger.error(
                "post-forward recovery: transition too close to t_start — "
                "no safe sub-range remains.  Lower t_start and re-run."
            )
            return

        safe_T    = blocks[safe_k]["temp"]
        safe_step = blocks[safe_k]["step"]
        chk_name  = "conf.ts.forward_%d_blk%d.data" % (iteration, safe_k)
        chk_path  = os.path.join(self.simfolder, chk_name)

        if not os.path.exists(chk_path):
            self.logger.error(
                "post-forward recovery: expected checkpoint %s not found — "
                "cannot recover", chk_path,
            )
            return

        # 1. Truncate the forward ts file to safe_step rows.
        kept = self._truncate_ts_file(out, safe_step)
        self.logger.warning(
            "post-forward recovery: truncated %s to %d data rows "
            "(safe boundary T=%.1f K, block %d/%d)",
            os.path.basename(out), kept, safe_T, safe_k, n_blocks,
        )

        # 2. Promote the checkpoint to the canonical forward-end configuration
        #    so the backward sweep loads the safe-boundary state instead of
        #    the (post-transition) end state of the full forward sweep.
        target_conf = os.path.join(
            self.simfolder, "conf.ts.forward_%d.data" % iteration
        )
        try:
            shutil.copy(chk_path, target_conf)
            self.logger.warning(
                "post-forward recovery: copied %s -> %s",
                os.path.basename(chk_path), os.path.basename(target_conf),
            )
        except Exception as cp_exc:
            self.logger.error(
                "post-forward recovery: failed to copy checkpoint to %s: %s",
                target_conf, cp_exc,
            )
            return

        # 3. Reduce the calc target so the backward sweep covers the same
        #    (now reduced) lambda range and produces the same number of
        #    rows as the truncated forward file.
        old_tf = self.calc._temperature_stop
        old_n  = self.calc._n_sweep_steps
        self.calc._temperature_stop = float(safe_T)
        self.calc._n_sweep_steps    = int(safe_step)
        self.logger.warning(
            "post-forward recovery: reduced T_stop %.1f K -> %.1f K and "
            "n_sweep_steps %d -> %d for the backward sweep",
            old_tf, safe_T, old_n, safe_step,
        )

        # Disable detection for the backward sweep over the now-safe range
        # so we don't re-trigger and loop.  Plots still produced.
        td.mode = "none"
        self.logger.warning(
            "post-forward recovery: detection disabled (mode='none') for "
            "the subsequent backward sweep"
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

        solid = self.calc.reference_phase == "solid"
        expected_phase = "solid" if solid else "liquid"
        self._create_monitor(expected_phase)

        # Always run the forward sweep to completion (no incremental aborts).
        # Per-block checkpoints are written by _run_split_sweep so that
        # _post_forward_recovery can restart the backward sweep from a safe
        # boundary if a transition is detected post-hoc.
        self._reversible_scaling_forward(iteration=iteration)

        # Post-hoc transition handling on the complete forward dataset.
        # If td.mode == 'recover' and a transition is found, this method
        # truncates the forward ts file, promotes a per-block checkpoint to
        # the canonical conf.ts.forward_N.data, and reduces _temperature_stop
        # and _n_sweep_steps so the backward sweep covers the safe sub-range.
        # If td.mode == 'stop' it raises PhaseTransitionError instead.
        # If td.mode == 'warn' it logs and plots but does not modify anything.
        self._post_forward_recovery(iteration=iteration)

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
        # Post-hoc transition detection on the final merged ts files.
        # When temperature_window > 0, the incremental detector already checked
        # the growing output file after every block, including the last block
        # (which contains the full sweep).  Run post-hoc anyway to ensure the
        # correct full-sweep baseline is used and to update the response-function
        # plots with accurate transition markers.
        td = self.calc.phase_transition_detection
        if td.mode != "none":
            self._detect_ts_transitions()

        # Always plot the response functions so the user can inspect them.
        self._plot_ts_response_functions()

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

        self.logger.info(
            f"Maximum energy dissipation along the temperature scaling part: {ediss} eV/atom"
        )
        if np.abs(ediss) > 1e-4:
            self.logger.warning(
                f"Found max energy dissipation of {ediss} along the temperature scaling path. Please ensure there are no structural changes!"
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

        expected_phase = "solid" if solid else "liquid"
        self._create_monitor(expected_phase)

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
            script_mode=False,
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

        # Variables for the phase-transition monitor (feed from initial equil)
        lmp.command("variable         ts_mlx equal lx")
        lmp.command("variable         ts_mly equal ly")
        lmp.command("variable         ts_mlz equal lz")
        lmp.command("variable         ts_mpress equal press")
        lmp.command("variable         ts_mpe equal pe/atoms")
        lmp.command("variable         ts_metotal equal etotal/atoms")
        lmp.command("variable         ts_mtemp equal temp")

        # equilibrate first
        _ts_ave_every  = int(self.calc.md.n_every_steps)
        _ts_ave_repeat = int(self.calc.md.n_repeat_steps)
        _ts_ave_freq   = _ts_ave_every * _ts_ave_repeat
        lmp.command(
            "fix               ts_ave all ave/time %d %d %d "
            "v_ts_mlx v_ts_mly v_ts_mlz v_ts_mpress v_ts_mpe v_ts_metotal v_ts_mtemp "
            "file tscale_equil_avg.dat"
            % (_ts_ave_every, _ts_ave_repeat, _ts_ave_freq)
        )
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
        lmp.command("unfix             ts_ave")

        # Feed the initial-equilibration data into the monitor so check_if_melted
        # / check_if_solidfied have enough samples at the mid-point check.
        self._monitor_fed_rows = 0
        self._feed_monitor_from_avg_dat(
            os.path.join(self.simfolder, "tscale_equil_avg.dat")
        )

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
        self._run_split_sweep(
            lmp=lmp,
            lambda_var="lambda",
            output_file_pattern="ts.forward_%d.dat" % iteration,
            sweep_label="tscale forward (iteration %d)" % iteration,
            t0=t0,
            tf=tf,
            iteration=iteration,
            sweep_mode="tscale",
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
        self._run_split_sweep(
            lmp=lmp,
            lambda_var="lambda",
            output_file_pattern="ts.backward_%d.dat" % iteration,
            sweep_label="tscale backward (iteration %d)" % iteration,
            t0=tf,
            tf=t0,
            iteration=iteration,
            sweep_mode="tscale",
        )

        lmp.command("unfix             f2")

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
            script_mode=False,
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
