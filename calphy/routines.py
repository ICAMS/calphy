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
"""

import copy
import numpy as np
import os
import time
from mendeleev import element
import yaml

from calphy.input import read_inputfile

# import calphy.queuekernel as cq
from calphy.errors import *
import calphy.helpers as ph

from calphy.liquid import Liquid
from calphy.solid import Solid
from calphy.composition_transformation import CompositionTransformation


class MeltingTemp:
    """
    Class for automated melting temperature calculation.

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
        self.calc = calculation
        self.simfolder = simfolder
        self.log_to_screen = log_to_screen
        self.dtemp = self.calc.melting_temperature.step
        self.maxattempts = self.calc.melting_temperature.attempts
        self.attempts = 0
        self.calculations = []

        self.get_trange()
        self.arg = None

        logfile = os.path.join(os.getcwd(), f"{self.calc.create_identifier()}.log")
        self.logger = ph.prepare_log(logfile, screen=log_to_screen)

    def prepare_calcs(self):
        """
        Prepare calculations list from given object

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        # here, we need to prepare a new calculation
        # protocol, read in, modify, write a output
        # read input again
        calculations = {"calculations": []}

        with open(self.calc.inputfile, "r") as fin:
            data = yaml.safe_load(fin)
        calc = data["calculations"][int(self.calc.kernel)]

        calc["mode"] = "ts"
        calc["temperature"] = [int(self.tmin), int(self.tmax)]
        calc["reference_phase"] = "solid"
        calc["phase_transition_detection"] = {"mode": "recover"}
        calc["folder_prefix"] = "mt%d" % self.attempts
        # Preserve n_iterations from the original melting_temperature calculation
        if "n_iterations" in data["calculations"][int(self.calc.kernel)]:
            calc["n_iterations"] = data["calculations"][int(self.calc.kernel)][
                "n_iterations"
            ]
        calculations["calculations"].append(calc)

        with open(self.calc.inputfile, "r") as fin:
            data = yaml.safe_load(fin)
        calc = data["calculations"][int(self.calc.kernel)]

        calc["mode"] = "ts"
        calc["temperature"] = [int(self.tmax), int(self.tmin)]
        calc["reference_phase"] = "liquid"
        calc["phase_transition_detection"] = {"mode": "recover"}
        calc["folder_prefix"] = "mt%d" % self.attempts
        # Preserve n_iterations from the original melting_temperature calculation
        if "n_iterations" in data["calculations"][int(self.calc.kernel)]:
            calc["n_iterations"] = data["calculations"][int(self.calc.kernel)][
                "n_iterations"
            ]
        calculations["calculations"].append(calc)

        outfile = f"{self.calc.create_identifier()}.{self.attempts}.yaml"
        with open(outfile, "w") as fout:
            yaml.safe_dump(calculations, fout)

        # now read in again, which would allow for checking and so on
        # one could do this smartly, and simply create from here.
        self.calculations = read_inputfile(outfile)

    def get_trange(self):
        """
        Get temperature range for calculations

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        tmin = self.calc._temperature - self.dtemp
        if tmin < 0:
            tmin = 10
        tmax = self.calc._temperature + self.dtemp
        self.tmax = tmax
        self.tmin = tmin

    def run_jobs(self):
        """
        Run calculations

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.prepare_calcs()

        self.soljob = Solid(
            calculation=self.calculations[0],
            simfolder=self.calculations[0].create_folders(),
        )
        self.lqdjob = Liquid(
            calculation=self.calculations[1],
            simfolder=self.calculations[1].create_folders(),
        )

        # Propagate MeltingTemp file handlers to sub-job loggers so that
        # all sub-job output also appears in melting_temperature.log
        for handler in self.logger.handlers:
            self.soljob.logger.addHandler(handler)
            self.lqdjob.logger.addHandler(handler)

        self.logger.info(
            "Free energy of %s and %s phases will be calculated"
            % (self.soljob.calc.lattice, self.lqdjob.calc.lattice)
        )
        self.logger.info("Temperature range of %f-%f" % (self.tmin, self.tmax))
        self.logger.info("STATE: Temperature range of %f-%f K" % (self.tmin, self.tmax))
        self.logger.info("Starting solid fe calculation")

        try:
            self.soljob = routine_fe(self.soljob)
        except (MeltedError, PhaseTransitionError):
            self.logger.info("Solid phase melted")
            return 2

        self.logger.info("Starting solid reversible scaling run")
        for i in range(self.soljob.calc.n_iterations):
            try:
                self.soljob.reversible_scaling(iteration=(i + 1))
            except (MeltedError, PhaseTransitionError):
                self.logger.info("Solid system melted during reversible scaling run")
                return 2

        self.solres = self.soljob.integrate_reversible_scaling(
            scale_energy=True, return_values=True
        )

        self.logger.info("Starting liquid fe calculation")
        try:
            self.lqdjob = routine_fe(self.lqdjob)
        except (SolidifiedError, PhaseTransitionError):
            self.logger.info("Liquid froze")
            return 3

        self.logger.info("Starting liquid reversible scaling calculation")
        for i in range(self.lqdjob.calc.n_iterations):
            try:
                self.lqdjob.reversible_scaling(iteration=(i + 1))
            except (SolidifiedError, PhaseTransitionError):
                self.logger.info("Liquid froze during reversible scaling calculation")
                return 3

        self.lqdres = self.lqdjob.integrate_reversible_scaling(
            scale_energy=True, return_values=True
        )

    # Half-width (K) of the temperature window centred on a predicted Tm
    # used when re-running after extrapolation.
    WINDOW_HALF_WIDTH = 200.0

    def _align_fe_curves(self):
        """
        Return (t_common, sol_f, sol_err, lqd_f, lqd_err) interpolated onto a
        shared temperature grid covering the overlap of the two FE curves.

        Returns None if the two curves have no overlapping temperature range
        (e.g. one was truncated by phase-transition recovery far from the other).
        """
        sol_t, sol_f, sol_err = self.solres
        lqd_t, lqd_f, lqd_err = self.lqdres

        ss = np.argsort(sol_t)
        ls = np.argsort(lqd_t)
        sol_t_s, sol_f_s, sol_err_s = sol_t[ss], sol_f[ss], sol_err[ss]
        lqd_t_s, lqd_f_s, lqd_err_s = lqd_t[ls], lqd_f[ls], lqd_err[ls]

        t_lo = max(sol_t_s[0], lqd_t_s[0])
        t_hi = min(sol_t_s[-1], lqd_t_s[-1])
        if t_hi <= t_lo:
            return None

        t_common = np.arange(t_lo, t_hi + 1.0, 1.0)
        sol_f_c   = np.interp(t_common, sol_t_s, sol_f_s)
        sol_err_c = np.interp(t_common, sol_t_s, sol_err_s)
        lqd_f_c   = np.interp(t_common, lqd_t_s, lqd_f_s)
        lqd_err_c = np.interp(t_common, lqd_t_s, lqd_err_s)
        return t_common, sol_f_c, sol_err_c, lqd_f_c, lqd_err_c

    def _find_intersection_in_overlap(self):
        """
        Look for a sign change of (sol_f - lqd_f) within the overlap of the two
        FE curves.  Returns (tm, tmerr) on success, or None if the curves do
        not cross within their overlap (or have no overlap at all).
        """
        aligned = self._align_fe_curves()
        if aligned is None:
            return None
        t_common, sol_f, sol_err, lqd_f, lqd_err = aligned
        diff = sol_f - lqd_f

        # Look for a sign change.
        sign = np.sign(diff)
        # Treat exact zeros as crossings.
        idx = np.where(sign[:-1] * sign[1:] <= 0)[0]
        if len(idx) == 0:
            return None

        # Use the first sign change.
        i = int(idx[0])
        # Linear interpolation between t_common[i] and t_common[i+1].
        d0, d1 = diff[i], diff[i + 1]
        if d1 == d0:
            tm = t_common[i]
            arg = i
        else:
            frac = d0 / (d0 - d1)
            tm = t_common[i] + frac * (t_common[i + 1] - t_common[i])
            arg = i if frac < 0.5 else i + 1

        # Error estimate: same formula as before, evaluated at the nearest grid point.
        suberr = float(np.sqrt(sol_err[arg] ** 2 + lqd_err[arg] ** 2))
        stride = min(50, arg, len(sol_f) - 1 - arg)
        if stride > 0:
            dt = t_common[arg + stride] - t_common[arg - stride]
            sol_slope = (sol_f[arg + stride] - sol_f[arg - stride]) / dt
            lqd_slope = (lqd_f[arg + stride] - lqd_f[arg - stride]) / dt
            slope_diff = sol_slope - lqd_slope
            tmerr = abs(suberr / slope_diff) if slope_diff != 0 else 0.0
        else:
            tmerr = 0.0
        self.arg = arg
        return float(tm), float(tmerr)

    def _extrapolate_intersection(self):
        """
        Fit a 4th-order polynomial to each of the solid and liquid FE curves
        (or a lower order if not enough points are available) and solve
        poly_sol(T) = poly_lqd(T) for the crossing temperature.

        Returns the predicted Tm (float).  Raises ValueError if no real
        crossing can be identified.
        """
        sol_t, sol_f, _ = self.solres
        lqd_t, lqd_f, _ = self.lqdres
        sol_t = np.asarray(sol_t)
        sol_f = np.asarray(sol_f)
        lqd_t = np.asarray(lqd_t)
        lqd_f = np.asarray(lqd_f)

        if len(sol_t) < 2 or len(lqd_t) < 2:
            raise ValueError(
                "Not enough FE points to fit polynomials for extrapolation "
                "(solid: %d, liquid: %d)" % (len(sol_t), len(lqd_t))
            )

        # Use up to 4th order, but never higher than (n_points - 1).
        sol_order = min(4, len(sol_t) - 1)
        lqd_order = min(4, len(lqd_t) - 1)

        sol_poly = np.polyfit(sol_t, sol_f, sol_order)
        lqd_poly = np.polyfit(lqd_t, lqd_f, lqd_order)

        # Pad both to the same length so we can subtract.
        n = max(len(sol_poly), len(lqd_poly))
        sp = np.concatenate([np.zeros(n - len(sol_poly)), sol_poly])
        lp = np.concatenate([np.zeros(n - len(lqd_poly)), lqd_poly])
        diff_poly = sp - lp

        roots = np.roots(diff_poly)
        # Keep real roots above 0 K.
        real_roots = []
        for r in roots:
            if abs(r.imag) < 1e-6 and r.real > 0:
                real_roots.append(r.real)

        if not real_roots:
            raise ValueError(
                "Polynomial extrapolation produced no positive real "
                "intersection of solid and liquid FE curves"
            )

        # Pick the root closest to the midpoint between the two data ranges
        # (most physically plausible Tm).
        anchor = 0.5 * (sol_t.max() + lqd_t.min())
        if not np.isfinite(anchor):
            anchor = 0.5 * (sol_t.mean() + lqd_t.mean())
        tpred = float(min(real_roots, key=lambda r: abs(r - anchor)))

        self.logger.info(
            "Polynomial extrapolation (orders sol=%d, lqd=%d) predicts Tm = %.1f K "
            "(real roots: %s)",
            sol_order, lqd_order, tpred,
            ", ".join("%.1f" % r for r in sorted(real_roots)),
        )
        self.logger.info("STATE: Predicted Tm from extrapolation: %.1f K", tpred)
        return tpred

    def _shift_window(self, returncode):
        """
        Shift the temperature window when one of the phases failed to even
        produce FE data (averaging crashed because the phase was unstable
        in the chosen window).
        """
        shift = self.WINDOW_HALF_WIDTH
        if returncode == 3:
            # Liquid froze before any FE data was collected: window too cold.
            self.tmin = self.tmin + shift
            self.tmax = self.tmax + shift
            self.logger.info(
                "STATE: Liquid froze, shifting window up by %.0f K -> [%.1f, %.1f]",
                shift, self.tmin, self.tmax,
            )
        elif returncode == 2:
            # Solid melted: window too hot.
            self.tmin = max(10.0, self.tmin - shift)
            self.tmax = max(self.tmin + 1.0, self.tmax - shift)
            self.logger.info(
                "STATE: Solid melted, shifting window down by %.0f K -> [%.1f, %.1f]",
                shift, self.tmin, self.tmax,
            )

    def calculate_tm(self):
        """
        Drive the iterative melting-temperature search.

        Algorithm
        ---------
        1. Run a solid + liquid reversible-scaling pass over the current
           [tmin, tmax] window.  Truncation by phase-transition recovery is
           tolerated — whatever FE curve was produced is used as-is.
        2. If both FE curves exist and their solid/liquid free energies
           cross within their overlapping temperature range, that crossing
           is the melting temperature → done.
        3. Otherwise, fit a 4th-order polynomial to each FE curve and solve
           for the intersection (extrapolation).  Re-centre the window on
           that prediction with a half-width of WINDOW_HALF_WIDTH (200 K)
           and repeat.
        4. If a phase failed to even produce an FE curve (averaging melted
           or froze), shift the window in the appropriate direction by
           WINDOW_HALF_WIDTH and repeat.
        """
        for _ in range(self.maxattempts + 1):
            returncode = self.run_jobs()

            if returncode in (2, 3):
                # No usable FE curve from one of the phases.
                self._shift_window(returncode)
            else:
                # Both FE curves available (possibly truncated).
                hit = self._find_intersection_in_overlap()
                if hit is not None:
                    tm, tmerr = hit
                    self.calc_tm = tm
                    self.tmerr = tmerr
                    self.logger.info(
                        "Found melting temperature = %.2f +/- %.2f K "
                        % (tm, tmerr)
                    )
                    if self.calc._melting_temperature is not None:
                        self.logger.info(
                            "Experimental melting temperature = %.2f K "
                            % (self.calc._melting_temperature)
                        )
                    self.logger.info(
                        "STATE: Tm = %.2f K +/- %.2f K" % (tm, tmerr)
                    )
                    return tm, tmerr

                # No real crossing in overlap → extrapolate, recentre, retry.
                try:
                    tpred = self._extrapolate_intersection()
                except ValueError as exc:
                    self.logger.info(
                        "Extrapolation failed (%s); falling back to window "
                        "midpoint between sol_max and lqd_min.",
                        exc,
                    )
                    sol_t = np.asarray(self.solres[0])
                    lqd_t = np.asarray(self.lqdres[0])
                    tpred = 0.5 * (sol_t.max() + lqd_t.min())

                self.tmin = max(10.0, tpred - self.WINDOW_HALF_WIDTH)
                self.tmax = tpred + self.WINDOW_HALF_WIDTH
                self.logger.info(
                    "Restarting with extrapolated Tm = %.1f K, window "
                    "[%.1f, %.1f] K (+/- %.0f K)",
                    tpred, self.tmin, self.tmax, self.WINDOW_HALF_WIDTH,
                )

            self.attempts += 1
            self.logger.info("Attempt incremented to %d" % self.attempts)
            if self.attempts > self.maxattempts:
                raise ValueError(
                    "Maximum number of melting-temperature attempts (%d) reached"
                    % self.maxattempts
                )

        raise ValueError(
            "Melting-temperature search did not converge within %d attempts"
            % self.maxattempts
        )


def routine_fe(job):
    """
    Perform an FE calculation routine
    """
    ts = time.time()
    job.run_averaging()
    te = time.time() - ts
    job.logger.info("Averaging routine finished in %f s" % te)

    # now run integration loops
    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.run_integration(iteration=(i + 1))
        te = time.time() - ts
        job.logger.info("Integration cycle %d finished in %f s" % (i + 1, te))

    job.thermodynamic_integration()
    job.submit_report()
    job.clean_up()
    return job


def routine_ts(job):
    """
    Perform ts routine
    """
    routine_fe(job)

    # now do rev scale steps
    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.reversible_scaling(iteration=(i + 1))
        te = time.time() - ts
        job.logger.info("TS integration cycle %d finished in %f s" % (i + 1, te))

    job.integrate_reversible_scaling(scale_energy=True)
    job.clean_up()
    return job


def routine_tscale(job):
    """
    Perform tscale routine
    """
    routine_fe(job)

    # now do rev scale steps
    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.temperature_scaling(iteration=(i + 1))
        te = time.time() - ts
        job.logger.info("Temperature scaling cycle %d finished in %f s" % (i + 1, te))

    job.integrate_reversible_scaling(scale_energy=False)
    job.clean_up()
    return job


def routine_pscale(job):
    """
    Perform pscale routine
    """
    routine_fe(job)

    # now do rev scale steps
    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.pressure_scaling(iteration=(i + 1))
        te = time.time() - ts
        job.logger.info("Pressure scaling cycle %d finished in %f s" % (i + 1, te))

    job.integrate_pressure_scaling()
    job.clean_up()
    return job


def routine_alchemy(job):
    """
    Perform an FE calculation routine
    """
    ts = time.time()
    job.run_averaging()
    te = time.time() - ts
    job.logger.info("Averaging routine finished in %f s" % te)

    # now run integration loops
    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.run_integration(iteration=(i + 1))
        te = time.time() - ts
        job.logger.info("Alchemy integration cycle %d finished in %f s" % (i + 1, te))

    job.thermodynamic_integration()
    job.submit_report()
    job.clean_up()
    return job


def routine_composition_scaling(job):
    """
    Perform a compositional scaling routine
    """
    # we set up comp scaling first
    job.logger.info("Calculating composition scaling")
    comp = CompositionTransformation(job.calc)
    forward_swap_types, reverse_swap_types = comp.get_swap_types(
        allow_all_swaps=job.calc.monte_carlo.allow_all_swaps
    )
    job.calc.monte_carlo.forward_swap_types = forward_swap_types
    job.calc.monte_carlo.reverse_swap_types = reverse_swap_types
    job.logger.info(f"Forward swap types: {forward_swap_types}")
    job.logger.info(f"Reverse swap types: {reverse_swap_types}")

    # update pair styles
    res = comp.update_pair_coeff(job.calc.pair_coeff[0])
    job.calc.pair_style.append(job.calc.pair_style[0])
    job.calc._pair_style_with_options.append(job.calc._pair_style_with_options[0])
    job.calc.pair_coeff[0] = res[0]
    job.calc.pair_coeff.append(res[1])
    job.logger.info("Update pair coefficients")
    job.logger.info(f"pair coeff 1: {job.calc.pair_coeff[0]}")
    job.logger.info(f"pair coeff 2: {job.calc.pair_coeff[1]}")
    job.calc._pair_style_names.append(job.calc._pair_style_names[0])
    job.logger.info("Update pair styles")
    job.logger.info(f"pair style 1: {job.calc._pair_style_names[0]}")
    job.logger.info(f"pair style 2: {job.calc._pair_style_names[1]}")

    backup_element = job.calc.element.copy()
    job.calc.element = comp.pair_list_old
    # job.calc._ghost_element_count = len(comp.new_atomtype) - len()

    # write new file out and update lattice
    outfilename = ".".join([job.calc.lattice, "comp", "data"])
    comp.write_structure(outfilename)
    job.calc.lattice = outfilename
    job.logger.info(f"Modified lattice written to {outfilename}")

    # prepare mass change methods
    # update and backup mass
    job.logger.info(f"Original mass: {job.calc.mass}")
    backup_mass = job.calc.mass.copy()
    mass_dict = {key: val for (key, val) in zip(backup_element, backup_mass)}

    target_masses = []
    target_counts = []

    ref_mass_list = []

    for mdict in comp.transformation_list:
        ref_mass_list.append(mass_dict[mdict["primary_element"]])
        target_masses.append(mass_dict[mdict["secondary_element"]])
        target_counts.append(mdict["count"])

    if len(backup_mass) > 2:
        job.logger.warning("Composition scaling is untested for more than 2 elements!")

    if len(np.unique(ref_mass_list)) > 1:
        job.logger.warning("More than one kind of transformation found! Stopping")
        raise RuntimeError("More than one kind of transformation found! Stopping")

    ref_mass = ref_mass_list[0]

    # now replace mass
    job.calc.mass = [ref_mass for x in range(len(job.calc.element))]
    job.logger.info(f"Temporarily replacing mass: {job.calc.mass}")

    # update fict elements if needed
    # job.calc._totalelements = comp.maxtype

    # now start cycle
    ts = time.time()
    job.run_averaging()
    te = time.time() - ts
    job.logger.info("Averaging routine finished in %f s" % te)

    # now run integration loops
    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.run_integration(iteration=(i + 1))
        te = time.time() - ts
        job.logger.info("Alchemy integration cycle %d finished in %f s" % (i + 1, te))

    job.thermodynamic_integration()

    job.logger.info("performing mass rescaling")
    job.logger.info(f"Ref. mass is {ref_mass}")
    job.logger.info(f"Target masses are {target_masses}")

    # read the file
    mcorsum = job.mass_integration(ref_mass, target_masses, target_counts)

    job.fe = job.fe - mcorsum
    job.submit_report(
        extra_dict={
            "results": {
                "mass_correction": float(mcorsum),
                "entropy_contribution": float(comp.entropy_contribution),
            }
        }
    )
    job.clean_up()
    return job
