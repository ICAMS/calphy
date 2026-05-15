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
    Automated melting-temperature search via bracket narrowing on two
    reversible-scaling sweeps.

    Algorithm (one attempt)
    -----------------------
    Two reversible-scaling sweeps are run over the current window
    ``[tmin, tmax]``:

    * **Solid sweep**, ``tmin → tmax``: produces ``G_solid(T)`` on
      ``[tmin, T_sol_max]``.  ``T_sol_max`` equals ``tmax`` if the sweep ran to
      completion, or the safe block-boundary temperature if
      ``phase_transition_detection`` truncated it because the solid lost
      stability.  In either case ``T_sol_max`` is an upper bound on Tm
      (with some superheating margin).
    * **Liquid sweep**, ``tmax → tmin``: produces ``G_liquid(T)`` on
      ``[T_liq_min, tmax]``.  ``T_liq_min`` is similarly a lower bound on
      Tm (with some supercooling margin).

    Three outcomes are then possible:

    1. **Overlap with crossing**: the two FE curves overlap on
       ``[T_liq_min, T_sol_max]`` and ``G_solid - G_liquid`` changes sign
       inside the overlap.  The crossing is the melting temperature →
       done.

    2. **Overlap, no crossing**: one phase is uniformly lower in FE across
       the overlap.  Fit each curve with the standard 6-term CALPHAD
       Gibbs-energy polynomial ``G(T) = a + b·T + c·T·lnT + d·T² + e·T³ +
       f·T⁻¹``, solve ``G_sol(T) = G_liq(T)``, recentre the window on the
       prediction with a half-width of ``melting_temperature.step``, and
       retry.

    3. **No overlap (gap)**: ``T_sol_max < T_liq_min``.  Tm is bracketed
       inside the gap.  Pick the gap midpoint as the next centre with a
       half-width sized to comfortably bracket the gap, and retry.

    A fourth outcome — failure of the *averaging* stage at the window
    endpoints, i.e. the streaming-monitor flipped the phase before the ts
    sweep even started — raises a clear error advising the user to adjust
    ``melting_temperature.guess``.  Window shifting on averaging failure is
    no longer attempted because the bracket-narrowing search above already
    handles every case in which the ts sweeps produce data, including
    badly off-centre initial guesses (which manifest as case 2 or 3 with a
    far extrapolation).

    Parameters
    ----------
    calculation : Calculation
        Parsed ``melting_temperature`` calculation object.
    simfolder : str
        Base folder for running calculations.
    log_to_screen : bool, optional
        Mirror the rotating log file to stdout.
    """

    def __init__(self, calculation=None, simfolder=None, log_to_screen=False):
        self.calc = calculation
        self.simfolder = simfolder
        self.log_to_screen = log_to_screen
        # Half-width (in K) of every search window.  Used for the initial
        # window around `guess`, for case-2 recentring after extrapolation,
        # and as a floor for case-3 gap-bracketing.
        self.dtemp = float(self.calc.melting_temperature.step)
        self.maxattempts = self.calc.melting_temperature.attempts
        self.attempts = 0
        self.calculations = []
        self.attempt_history = []
        # Running list of case-2 CALPHAD-extrapolated Tm predictions.  When
        # the last two agree within :attr:`CASE2_CONVERGENCE_TOLERANCE` K
        # the search is declared converged on extrapolation stability —
        # useful when MD hysteresis is wide enough that no case-1 overlap
        # crossing will ever fire (small systems, sluggish potentials).
        self._case2_tpred_history = []
        # Last CLEAN case-2 prediction (used as the recentring anchor when
        # a subsequent attempt comes back contaminated and its own tpred
        # cannot be trusted).
        self._last_clean_tpred = None
        # Running count of consecutive contaminated attempts.  Reset on a
        # clean attempt; if it reaches MAX_CONSECUTIVE_CONTAMINATED the
        # search aborts because no amount of recentring will fix an
        # under-tuned detector.
        self._consecutive_contaminated = 0
        self.arg = None

        self.get_trange()

        logfile = os.path.join(os.getcwd(), f"{self.calc.create_identifier()}.log")
        self.logger = ph.prepare_log(logfile, screen=log_to_screen)

    # ------------------------------------------------------------------
    # Setup helpers
    # ------------------------------------------------------------------

    def get_trange(self):
        """Set the initial ``[tmin, tmax]`` window centred on the user's guess."""
        tmin = self.calc._temperature - self.dtemp
        if tmin < 0:
            tmin = 10.0
        self.tmin = float(tmin)
        self.tmax = float(self.calc._temperature + self.dtemp)

    def _user_ptd_settings(self):
        """
        Return the user's ``phase_transition_detection`` block from the
        original YAML if any, so it can be propagated into the synthesised
        ts sub-calculations (only ``mode`` is forced to ``recover``).
        """
        with open(self.calc.inputfile, "r") as fin:
            data = yaml.safe_load(fin)
        ci = data["calculations"][int(self.calc.kernel)]
        ptd = {}
        if isinstance(ci.get("phase_transition_detection"), dict):
            ptd.update(ci["phase_transition_detection"])
        if isinstance(data.get("phase_transition_detection"), dict):
            for k, v in data["phase_transition_detection"].items():
                ptd.setdefault(k, v)
        ptd["mode"] = "recover"
        return ptd

    def prepare_calcs(self):
        """
        Build a fresh pair of ts sub-calculations (solid heating sweep,
        liquid cooling sweep) for the current attempt and read them back
        through the input parser so they go through the same validation as
        a user-written input.

        The user's ``phase_transition_detection`` settings are preserved;
        only the ``mode`` is forced to ``recover`` so the sub-sweeps
        truncate at the safe boundary instead of either aborting (``stop``)
        or continuing into a mixed-phase region (``warn``/``none``).
        """
        with open(self.calc.inputfile, "r") as fin:
            data = yaml.safe_load(fin)
        base = data["calculations"][int(self.calc.kernel)]
        ptd  = self._user_ptd_settings()

        def _build(reference_phase, t_start, t_stop):
            calc = copy.deepcopy(base)
            calc["mode"] = "ts"
            calc["temperature"] = [int(round(t_start)), int(round(t_stop))]
            calc["reference_phase"] = reference_phase
            calc["phase_transition_detection"] = copy.deepcopy(ptd)
            calc["folder_prefix"] = "mt%d" % self.attempts
            if "n_iterations" in base:
                calc["n_iterations"] = base["n_iterations"]
            return calc

        calculations = {"calculations": [
            _build("solid",  self.tmin, self.tmax),
            _build("liquid", self.tmax, self.tmin),
        ]}
        outfile = f"{self.calc.create_identifier()}.{self.attempts}.yaml"
        with open(outfile, "w") as fout:
            yaml.safe_dump(calculations, fout)
        self.calculations = read_inputfile(outfile)

    # ------------------------------------------------------------------
    # Per-attempt execution
    # ------------------------------------------------------------------

    # Outcome codes returned by run_jobs.  Truncation by the ts-sweep
    # phase-transition detector is NOT a failure — it produces a usable
    # (truncated) FE curve which the bracket logic then consumes.
    OUTCOME_OK              = "ok"
    OUTCOME_SOLID_AVG_FAIL  = "solid_averaging_failed"
    OUTCOME_LIQUID_AVG_FAIL = "liquid_averaging_failed"

    # Tolerance (K) for declaring two consecutive case-2 CALPHAD
    # extrapolations "converged".  A run that produces predictions
    # 1278.7 and 1305.2 K (|Δ|=26.5 K) does NOT converge at 25 K but
    # WOULD converge at 30 K; tighter than ~15 K is unrealistic given the
    # noise of CALPHAD fits over a narrow overlap.  25 K is a balance:
    # catches genuine convergence while still distinguishing two
    # extrapolations that disagree due to changing window placement.
    CASE2_CONVERGENCE_TOLERANCE = 25.0

    # Maximum forward/backward energy-dissipation (eV/atom) of either ts
    # sweep above which we treat the attempt as CONTAMINATED — i.e. the
    # phase-transition detector silently missed a transition during the
    # sweep and the resulting FE curve is a mix of two phases rather than
    # a clean single-phase G(T).  Clean reversible sweeps come in at
    # ~1e-4 eV/atom; a missed transition is typically 0.01–0.1 eV/atom
    # (two to three orders of magnitude bigger).  1e-3 sits comfortably
    # between those regimes.
    EDISS_CONTAMINATION_THRESHOLD = 1e-3

    # How many consecutive contaminated attempts to tolerate before
    # raising — at the limit, no amount of recentring will fix the
    # underlying problem (detector too insensitive for the system size).
    MAX_CONSECUTIVE_CONTAMINATED = 2

    def run_jobs(self):
        """
        Build and run one solid + liquid ts pair for the current window.

        Returns
        -------
        outcome : str
            One of ``OUTCOME_OK``, ``OUTCOME_SOLID_AVG_FAIL``,
            ``OUTCOME_LIQUID_AVG_FAIL``.  On ``OUTCOME_OK`` both
            ``self.solres`` and ``self.lqdres`` are set to the
            ``(T, F, F_err)`` arrays returned by
            :meth:`Phase.integrate_reversible_scaling`.  These curves may
            be truncated by phase-transition recovery — that's expected
            and is handled by :meth:`calculate_tm`.

        Notes
        -----
        With ``phase_transition_detection.mode = recover`` (forced by
        :meth:`prepare_calcs`) the ts sweeps themselves cannot raise
        :class:`PhaseTransitionError` — they truncate silently.  Only the
        streaming phase-stability monitor used during *averaging* can
        still raise :class:`MeltedError` / :class:`SolidifiedError`, and
        that is the only path through which this method reports failure
        rather than success.
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

        # Mirror sub-job log output into melting_temperature.log
        for handler in self.logger.handlers:
            self.soljob.logger.addHandler(handler)
            self.lqdjob.logger.addHandler(handler)

        self.logger.info(
            "Attempt %d: solid (%s) + liquid (%s), window [%.1f, %.1f] K",
            self.attempts,
            self.soljob.calc.lattice, self.lqdjob.calc.lattice,
            self.tmin, self.tmax,
        )
        self.logger.info("STATE: Temperature range of %f-%f K", self.tmin, self.tmax)

        try:
            self.soljob = routine_fe(self.soljob)
        except MeltedError as exc:
            self.logger.warning(
                "Solid averaging failed at T=%.1f K: %s", self.tmin, exc
            )
            return self.OUTCOME_SOLID_AVG_FAIL

        for i in range(self.soljob.calc.n_iterations):
            self.soljob.reversible_scaling(iteration=(i + 1))
        self.solres = self.soljob.integrate_reversible_scaling(
            scale_energy=True, return_values=True
        )
        # Capture the maximum forward/backward energy dissipation as the
        # data-quality flag for this phase's FE curve (see
        # EDISS_CONTAMINATION_THRESHOLD).
        self.sol_ediss = float(getattr(self.soljob, "ediss", float("nan")))

        try:
            self.lqdjob = routine_fe(self.lqdjob)
        except SolidifiedError as exc:
            self.logger.warning(
                "Liquid averaging failed at T=%.1f K: %s", self.tmax, exc
            )
            return self.OUTCOME_LIQUID_AVG_FAIL

        for i in range(self.lqdjob.calc.n_iterations):
            self.lqdjob.reversible_scaling(iteration=(i + 1))
        self.lqdres = self.lqdjob.integrate_reversible_scaling(
            scale_energy=True, return_values=True
        )
        self.lqd_ediss = float(getattr(self.lqdjob, "ediss", float("nan")))
        return self.OUTCOME_OK

    # ------------------------------------------------------------------
    # FE-curve analysis helpers
    # ------------------------------------------------------------------

    def _align_fe_curves(self):
        """
        Interpolate the solid and liquid FE curves onto a 1 K grid covering
        their overlap.  Returns
        ``(t_common, sol_f, sol_err, lqd_f, lqd_err)`` or ``None`` when
        the curves do not overlap.
        """
        sol_t, sol_f, sol_err = self.solres
        lqd_t, lqd_f, lqd_err = self.lqdres
        ss = np.argsort(sol_t); ls = np.argsort(lqd_t)
        sol_t_s, sol_f_s, sol_err_s = sol_t[ss], sol_f[ss], sol_err[ss]
        lqd_t_s, lqd_f_s, lqd_err_s = lqd_t[ls], lqd_f[ls], lqd_err[ls]

        t_lo = max(sol_t_s[0], lqd_t_s[0])
        t_hi = min(sol_t_s[-1], lqd_t_s[-1])
        if t_hi <= t_lo:
            return None

        t_common = np.arange(t_lo, t_hi + 1.0, 1.0)
        return (
            t_common,
            np.interp(t_common, sol_t_s, sol_f_s),
            np.interp(t_common, sol_t_s, sol_err_s),
            np.interp(t_common, lqd_t_s, lqd_f_s),
            np.interp(t_common, lqd_t_s, lqd_err_s),
        )

    def _find_intersection_in_overlap(self):
        """
        Locate a sign change of ``G_sol(T) − G_liq(T)`` inside the overlap.
        Returns ``(tm, tmerr)`` on success, ``None`` if no crossing.

        When two or more crossings are present (rare; usually means noise
        near a near-tangent intersection) the lowest-T crossing is used
        and the situation is logged.
        """
        aligned = self._align_fe_curves()
        if aligned is None:
            return None
        t_common, sol_f, sol_err, lqd_f, lqd_err = aligned
        diff = sol_f - lqd_f
        sign = np.sign(diff)
        idx = np.where(sign[:-1] * sign[1:] <= 0)[0]
        if len(idx) == 0:
            return None
        if len(idx) > 1:
            self.logger.warning(
                "Multiple sign changes (%d) of (G_sol - G_liq) inside the "
                "overlap; using the lowest-T crossing at T ~ %.1f K. "
                "Crossings at: %s",
                len(idx), float(t_common[int(idx[0])]),
                ", ".join("%.1f" % float(t_common[int(j)]) for j in idx),
            )

        i = int(idx[0])
        d0, d1 = diff[i], diff[i + 1]
        if d1 == d0:
            tm, arg = float(t_common[i]), i
        else:
            frac = d0 / (d0 - d1)
            tm = float(t_common[i] + frac * (t_common[i + 1] - t_common[i]))
            arg = i if frac < 0.5 else i + 1

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
        return tm, float(tmerr)

    def _calphad_extrapolate_tm(self):
        """
        Fit a six-term CALPHAD Gibbs polynomial to each phase and locate
        the lowest-T crossing of ``G_sol(T) = G_liq(T)``.

            G(T) = a + b·T + c·T·ln T + d·T² + e·T³ + f·T⁻¹

        This is the canonical CALPHAD form (SGTE pure-element Gibbs
        polynomial truncated at the standard six terms) and is therefore a
        much better extrapolation basis than a naive polynomial: the
        ``T·ln T`` term encodes the leading entropic contribution
        correctly, so extrapolating modestly outside the data window does
        not produce the wild excursions that ``np.polyfit`` gives.

        Returns
        -------
        tpred : float
            Predicted Tm in K.

        Raises
        ------
        ValueError
            If too few points are available to fit, or if no positive
            crossing of the two fits can be located in a generous search
            window around the data.
        """
        from calphy.phase_diagram import (
            _fit_calphad_poly6 as _fit,
            _eval_calphad_poly6 as _eval,
        )

        sol_t = np.asarray(self.solres[0], dtype=float)
        sol_f = np.asarray(self.solres[1], dtype=float)
        lqd_t = np.asarray(self.lqdres[0], dtype=float)
        lqd_f = np.asarray(self.lqdres[1], dtype=float)

        # poly6 has 6 free parameters; require at least that many points
        # per phase, otherwise the fit is degenerate.
        if len(sol_t) < 6 or len(lqd_t) < 6:
            raise ValueError(
                "Not enough FE points for a CALPHAD poly6 fit "
                "(solid: %d, liquid: %d; need ≥6 each)"
                % (len(sol_t), len(lqd_t))
            )
        if np.any(sol_t <= 0) or np.any(lqd_t <= 0):
            raise ValueError(
                "CALPHAD poly6 requires T > 0 everywhere (T·lnT and 1/T)"
            )

        sol_coef = _fit(sol_t, sol_f)
        lqd_coef = _fit(lqd_t, lqd_f)

        # Search for ΔG = G_sol − G_liq sign change on a fine grid spanning
        # both data ranges with a generous extrapolation margin on each
        # side.  The 6-term form is non-polynomial in T (T·lnT, 1/T) so we
        # scan numerically and refine with brentq, rather than trying to
        # solve symbolically.
        t_lo_data = float(min(sol_t.min(), lqd_t.min()))
        t_hi_data = float(max(sol_t.max(), lqd_t.max()))
        margin = max(0.5 * (t_hi_data - t_lo_data), 200.0)
        t_lo = max(50.0, t_lo_data - margin)
        t_hi = t_hi_data + margin
        T_grid = np.arange(t_lo, t_hi + 1.0, 1.0)

        def _dG(T):
            return _eval(sol_coef, T) - _eval(lqd_coef, T)

        dG_grid = _dG(T_grid)
        sign = np.sign(dG_grid)
        cross_idx = np.where(sign[:-1] * sign[1:] < 0)[0]
        if len(cross_idx) == 0:
            raise ValueError(
                "CALPHAD extrapolation found no G_sol = G_liq crossing in "
                "T ∈ [%.0f, %.0f] K" % (t_lo, t_hi)
            )

        # Refine each bracketing crossing with brentq; pick the root closest
        # to the bracket midpoint of the data (most physically plausible Tm).
        from scipy.optimize import brentq
        anchor = 0.5 * (float(sol_t.max()) + float(lqd_t.min()))
        roots = []
        for j in cross_idx:
            try:
                roots.append(float(brentq(_dG, T_grid[j], T_grid[j + 1])))
            except ValueError:
                continue
        if not roots:
            raise ValueError("CALPHAD extrapolation: brentq found no root")
        tpred = float(min(roots, key=lambda r: abs(r - anchor)))

        self.logger.info(
            "CALPHAD extrapolation: G_sol = G_liq at T = %.1f K "
            "(candidates: %s)",
            tpred, ", ".join("%.1f" % r for r in sorted(roots)),
        )
        return tpred

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------

    def _record_attempt(self, **kw):
        """
        Append one row to :attr:`attempt_history` for the final report.

        Per-phase energy dissipation values (``ediss_solid``,
        ``ediss_liquid``) are pulled automatically from the corresponding
        instance attributes set by :meth:`run_jobs`; explicit kwargs of
        the same name win over the auto-populated values for the
        contaminated-attempt path which wants to emphasise them.
        """
        row = {"attempt": int(self.attempts),
               "tmin": float(self.tmin),
               "tmax": float(self.tmax)}
        sol_e = getattr(self, "sol_ediss", None)
        lqd_e = getattr(self, "lqd_ediss", None)
        if sol_e is not None and np.isfinite(sol_e):
            row["ediss_solid"] = float(sol_e)
        if lqd_e is not None and np.isfinite(lqd_e):
            row["ediss_liquid"] = float(lqd_e)
        row.update(kw)
        self.attempt_history.append(row)

    def _write_report(self, tm, tmerr, status="converged"):
        """
        Write ``melting_temperature_report.yaml`` with the converged Tm,
        its uncertainty, the experimental reference (if known), and the
        full attempt history.  The file lives next to the rotating log so
        downstream tooling can parse it without grepping the log.
        """
        report = {
            "status": status,
            "tm": float(tm) if tm is not None else None,
            "tmerr": float(tmerr) if tmerr is not None else None,
            "tm_experimental": (
                float(self.calc._melting_temperature)
                if self.calc._melting_temperature is not None else None
            ),
            "attempts": int(self.attempts + 1),
            "max_attempts": int(self.maxattempts),
            "step": float(self.dtemp),
            "history": self.attempt_history,
        }
        out = os.path.join(os.getcwd(),
                           f"{self.calc.create_identifier()}.report.yaml")
        with open(out, "w") as fh:
            yaml.safe_dump(report, fh, sort_keys=False)
        self.logger.info("Wrote melting-temperature report: %s", out)

    # ------------------------------------------------------------------
    # Main driver
    # ------------------------------------------------------------------

    def _raise_averaging_failure(self, outcome):
        """Case 4: streaming monitor flipped the phase during averaging."""
        which = "solid" if outcome == self.OUTCOME_SOLID_AVG_FAIL else "liquid"
        msg = (
            "melting_temperature: %s-phase averaging failed in window "
            "[%.1f, %.1f] K — the streaming phase-stability monitor "
            "detected a phase flip at the window endpoint, so no reversible-"
            "scaling data could be collected.  The initial window is too "
            "far from the actual Tm; set "
            "`melting_temperature.guess` closer to the expected melting "
            "point and re-run."
            % (which, self.tmin, self.tmax)
        )
        self._record_attempt(outcome=outcome, case="4: averaging failed")
        self._write_report(tm=None, tmerr=None, status="averaging_failed")
        raise MeltedError(msg) if which == "solid" else SolidifiedError(msg)

    def calculate_tm(self):
        """
        Drive the bracket-narrowing search.  See class docstring for the
        full algorithm.  Returns ``(tm, tmerr)`` on convergence; raises
        :class:`ValueError` if the maximum number of attempts is reached
        without convergence, or :class:`MeltedError` / :class:`SolidifiedError`
        if averaging fails at the window endpoints (case 4).
        """
        for attempt in range(self.maxattempts + 1):
            self.attempts = attempt
            outcome = self.run_jobs()
            if outcome != self.OUTCOME_OK:
                self._raise_averaging_failure(outcome)

            # Truncation points double as one-sided brackets on Tm.
            sol_t = np.asarray(self.solres[0], dtype=float)
            lqd_t = np.asarray(self.lqdres[0], dtype=float)
            T_sol_max = float(sol_t.max())
            T_liq_min = float(lqd_t.min())
            has_overlap = T_liq_min <= T_sol_max

            # --- Contamination guard -------------------------------------
            # Either ts sweep can silently sample a hidden phase transition
            # when the detector is under-tuned for the system size.  The
            # tell-tale signature is forward/backward energy dissipation
            # orders of magnitude above the clean ~1e-4 eV/atom floor.  If
            # either phase tripped that flag, treat the attempt as
            # contaminated: skip case-1/2/3 analysis (the FE curves are
            # mixed-phase, so any tpred from them is unreliable), recentre
            # on the last clean prediction (or the overlap midpoint as a
            # safe fallback) and try again.
            sol_e = float(getattr(self, "sol_ediss", float("nan")))
            lqd_e = float(getattr(self, "lqd_ediss", float("nan")))
            max_e = max(abs(sol_e) if np.isfinite(sol_e) else 0.0,
                        abs(lqd_e) if np.isfinite(lqd_e) else 0.0)
            contaminated = max_e > self.EDISS_CONTAMINATION_THRESHOLD
            if contaminated:
                self._consecutive_contaminated += 1
                self.logger.warning(
                    "Attempt %d CONTAMINATED: max |ediss| = %.4f eV/atom "
                    "(threshold %.0e). Solid=%.4f, liquid=%.4f. The "
                    "phase-transition detector missed a transition during "
                    "one or both sweeps — the resulting FE curve is "
                    "mixed-phase and will give an unreliable Tm.  Skipping "
                    "case analysis for this attempt (consecutive "
                    "contaminated: %d / %d).",
                    attempt, max_e, self.EDISS_CONTAMINATION_THRESHOLD,
                    sol_e, lqd_e,
                    self._consecutive_contaminated,
                    self.MAX_CONSECUTIVE_CONTAMINATED,
                )
                self._record_attempt(
                    T_sol_max=T_sol_max, T_liq_min=T_liq_min,
                    ediss_solid=sol_e, ediss_liquid=lqd_e,
                    contaminated=True,
                    case="contaminated (skipped)",
                )

                if self._consecutive_contaminated >= self.MAX_CONSECUTIVE_CONTAMINATED:
                    self._write_report(
                        tm=None, tmerr=None, status="contaminated_aborted",
                    )
                    raise ValueError(
                        "melting_temperature: %d consecutive contaminated "
                        "attempts (max |ediss| > %.0e eV/atom).  The "
                        "phase-transition detector is too insensitive for "
                        "this system+potential — every sweep ends with the "
                        "system having silently changed phase.  Mitigation "
                        "options:\n"
                        "  • Lower `phase_transition_detection.peak_threshold` "
                        "(default 12; try 6 for ~1000-atom EAM systems).\n"
                        "  • Increase the system size (`repeat: [N,N,N]` "
                        "with bigger N) so variance signals exceed threshold.\n"
                        "  • Raise `n_switching_steps` so the sweep is slower "
                        "and the detector can build statistics."
                        % (self._consecutive_contaminated,
                           self.EDISS_CONTAMINATION_THRESHOLD)
                    )

                # Recentre on the previous clean prediction if any; else
                # fall back to the overlap (or gap) midpoint of THIS
                # attempt.  The midpoint is a noisy estimate but at least
                # has the right order of magnitude even from contaminated
                # data (it depends on window geometry, not FE values).
                if self._last_clean_tpred is not None:
                    tpred = self._last_clean_tpred
                else:
                    tpred = 0.5 * (T_sol_max + T_liq_min) if has_overlap \
                        else 0.5 * (T_sol_max + T_liq_min)
                self.tmin = max(10.0, tpred - self.dtemp)
                self.tmax = tpred + self.dtemp
                self.logger.info(
                    "Retry after contamination: recentring on %.1f K -> "
                    "[%.1f, %.1f] K (+/- %.0f K)",
                    tpred, self.tmin, self.tmax, self.dtemp,
                )
                continue
            else:
                # Clean attempt resets the contamination counter.
                self._consecutive_contaminated = 0
            # -------------------------------------------------------------

            self.logger.info(
                "Attempt %d brackets: T_sol_max=%.1f K, T_liq_min=%.1f K, %s",
                attempt, T_sol_max, T_liq_min,
                "overlap" if has_overlap else "gap (T_sol_max < T_liq_min)",
            )

            if has_overlap:
                hit = self._find_intersection_in_overlap()
                if hit is not None:
                    # Case 1: FE crossing found inside the overlap.
                    tm, tmerr = hit
                    self.calc_tm = tm
                    self.tmerr   = tmerr
                    self._record_attempt(
                        T_sol_max=T_sol_max, T_liq_min=T_liq_min,
                        case="1: crossing in overlap",
                        tm=tm, tmerr=tmerr,
                    )
                    self.logger.info(
                        "Found melting temperature = %.2f +/- %.2f K", tm, tmerr,
                    )
                    if self.calc._melting_temperature is not None:
                        self.logger.info(
                            "Experimental melting temperature = %.2f K",
                            self.calc._melting_temperature,
                        )
                    self.logger.info("STATE: Tm = %.2f K +/- %.2f K", tm, tmerr)
                    self._write_report(tm, tmerr)
                    return tm, tmerr

                # Case 2: overlap present but no crossing → CALPHAD extrapolate.
                try:
                    tpred = self._calphad_extrapolate_tm()
                except ValueError as exc:
                    self.logger.warning(
                        "CALPHAD extrapolation failed (%s); falling back to "
                        "overlap midpoint.", exc,
                    )
                    tpred = 0.5 * (T_liq_min + T_sol_max)
                # Bound the prediction so a runaway fit cannot send us off
                # to the moon: clamp to within one window-half-width beyond
                # the current data range on either side.
                t_lo_data = float(min(sol_t.min(), lqd_t.min()))
                t_hi_data = float(max(sol_t.max(), lqd_t.max()))
                tpred = float(np.clip(tpred,
                                      t_lo_data - self.dtemp,
                                      t_hi_data + self.dtemp))
                self._case2_tpred_history.append(tpred)

                # Stability-based convergence: when two consecutive case-2
                # CALPHAD predictions agree within tolerance the search is
                # done even without a case-1 overlap crossing.  This is the
                # only convergence path available for small systems with
                # MD hysteresis wider than the search window (case 1 is
                # then structurally unreachable because both overlaps sit
                # entirely on one side of Tm).
                if len(self._case2_tpred_history) >= 2:
                    last2 = self._case2_tpred_history[-2:]
                    delta = abs(last2[1] - last2[0])
                    if delta <= self.CASE2_CONVERGENCE_TOLERANCE:
                        tm    = float(np.mean(last2))
                        tmerr = float(np.std(last2))
                        self.calc_tm = tm
                        self.tmerr   = tmerr
                        self._record_attempt(
                            T_sol_max=T_sol_max, T_liq_min=T_liq_min,
                            case="2: CALPHAD extrapolations converged",
                            tpred=tpred,
                            tm=tm, tmerr=tmerr,
                        )
                        self.logger.info(
                            "Case-2 CALPHAD extrapolations converged: "
                            "Tm = %.2f +/- %.2f K (last two predictions: "
                            "%.1f, %.1f K; |Δ| = %.1f K ≤ %.1f K)",
                            tm, tmerr, last2[0], last2[1],
                            delta, self.CASE2_CONVERGENCE_TOLERANCE,
                        )
                        if self.calc._melting_temperature is not None:
                            self.logger.info(
                                "Experimental melting temperature = %.2f K",
                                self.calc._melting_temperature,
                            )
                        self.logger.info("STATE: Tm = %.2f K +/- %.2f K", tm, tmerr)
                        self._write_report(
                            tm, tmerr, status="converged_extrapolation",
                        )
                        return tm, tmerr

                self._record_attempt(
                    T_sol_max=T_sol_max, T_liq_min=T_liq_min,
                    case="2: overlap, no crossing → CALPHAD extrapolation",
                    tpred=tpred,
                )
                self._last_clean_tpred = float(tpred)
                self.tmin = max(10.0, tpred - self.dtemp)
                self.tmax = tpred + self.dtemp
                self.logger.info(
                    "Recentring window on extrapolated Tm = %.1f K -> "
                    "[%.1f, %.1f] K (+/- %.0f K)",
                    tpred, self.tmin, self.tmax, self.dtemp,
                )
            else:
                # Case 3: no overlap.  Tm is bracketed in the gap
                # [T_sol_max, T_liq_min].  Recentre on the midpoint with a
                # half-width sized to comfortably contain the gap.
                gap = T_liq_min - T_sol_max
                tpred = 0.5 * (T_sol_max + T_liq_min)
                # Half-width: never narrower than the initial window.  In
                # systems with significant MD hysteresis (typical Cu EAM at
                # ~1000 atoms shows ~150-200 K of superheating /
                # supercooling) shrinking the window after a case-3 gap
                # turned out to trap subsequent attempts inside the
                # hysteresis band — both phases would then stay metastable
                # past Tm in the overlap, no FE crossing would be visible,
                # and the search would drift on extrapolations.  Keeping
                # the half-width at ``self.dtemp`` (or wider if the gap
                # itself demands it) preserves enough headroom on each
                # side of Tm for the deeper phase to sample cleanly.
                half = max(gap, self.dtemp)
                self._record_attempt(
                    T_sol_max=T_sol_max, T_liq_min=T_liq_min,
                    case="3: gap → midpoint",
                    gap=[float(T_sol_max), float(T_liq_min)],
                    tpred=tpred,
                )
                self._last_clean_tpred = float(tpred)
                self.tmin = max(10.0, tpred - half)
                self.tmax = tpred + half
                self.logger.info(
                    "No overlap: Tm in gap [%.1f, %.1f] K (width %.1f K).  "
                    "Recentring on %.1f K -> [%.1f, %.1f] K (+/- %.1f K)",
                    T_sol_max, T_liq_min, gap, tpred,
                    self.tmin, self.tmax, half,
                )

        # Exhausted attempts without case-1 convergence.
        self._write_report(tm=None, tmerr=None, status="max_attempts_reached")
        raise ValueError(
            "Melting-temperature search did not converge within %d attempts. "
            "See attempt history in the report file."
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
