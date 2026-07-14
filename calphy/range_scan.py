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
"Automated Free Energy Calculation from Atomistic Simulations." Physical Review Materials 5(10), 2021
DOI: 10.1103/PhysRevMaterials.5.103801
"""

"""
Pre-flight temperature-range scan for calphy ``ts`` (reversible-scaling) runs.

Before the production reversible-scaling sweep, calphy can run a single fast
**real-thermostat temperature ramp** (T0 -> Tf under NPT, the same physics as
mode ``tscale``) and watch for the onset of a phase transition.  The clean
sub-range [T0, T_clean] determined here is then used to bound the production
``ts`` sweep, so the sweep never crosses a melting / solid-solid transition.

Why a real ramp (not the reversible-scaling lambda ramp)?
---------------------------------------------------------
During this scan the thermostat genuinely ramps from T0 to Tf, so the
temperature is *measured* directly from the MD (the LAMMPS ``temp`` column).
There is no Hamiltonian scaling: lambda == 1 throughout.  This means the
fluctuation response functions are the plain NPT expressions

    C_P     = Var(H)        / (k_B T**2)        [eV/K per atom]
    kappa_T = Var(V)        / (k_B T <V>)       [1/eV]
    alpha_P = Cov(V, H)     / (k_B T**2 <V>)    [1/K]

with H = pe + P * V the per-atom enthalpy.  No division by lambda is needed
(unlike the reversible-scaling detector, which worked on lambda-reduced
H/lambda and a reconstructed equivalent temperature T0/lambda).  Using the
measured temperature also makes the first-moment slope-break fit clean: H(T)
and V(T) are nearly linear within a single phase, so a low-order polynomial in
the measured T extrapolates well and a transition shows up as a sharp residual.

Two detector families
----------------------
1. Variance-based (second moment): peaks in rolling Cp / kappa_T / alpha_P,
   flagged by a robust modified Z-score against the median+MAD baseline.
2. Slope-break (first moment): normalized residuals of H(T) and V(T) against a
   single-phase EOS polynomial fit, which fire as soon as the *mean* departs
   from the baseline (i.e. at the latent-heat / volume jump), typically earlier
   than the variance peak.

A transition is declared only when at least ``min_agreement`` signals trigger
together.  The reported ``onset_temperature`` (the earliest sign of trouble,
not the peak) is the safe boundary used to bound the production sweep.
"""

import warnings
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

# Boltzmann constant in eV/K (LAMMPS metal units)
KB_EV: float = 8.617333262e-5

# 1 bar * Ang^3 = 6.2415091e-7 eV
BAR_ANG3_TO_EV: float = 6.2415091e-7


# ---------------------------------------------------------------------------
# Result record
# ---------------------------------------------------------------------------

@dataclass
class ScanResult:
    """Outcome of a pre-flight temperature-range scan."""

    transition_found: bool = False
    # Safe upper temperature for the production sweep (K).  Equals the
    # requested t_stop when no transition is found.
    clean_t_stop: float = 0.0
    # Earliest sign of the transition (walk-back onset, K).  Always <= the
    # peak temperature for a heating ramp.  This is the value used for
    # clean_t_stop when a transition is found.
    onset_temperature: float = 0.0
    # Peak temperature of the response functions (K), informational.
    peak_temperature: float = 0.0
    # Names of the signals that triggered (Cp, kappa_T, alpha_P, H_break, V_break).
    triggered_signals: List[str] = field(default_factory=list)
    # Fraction of evaluated signals that triggered.
    confidence: float = 0.0


# ---------------------------------------------------------------------------
# Thermodynamics
# ---------------------------------------------------------------------------

def enthalpy(pe: np.ndarray, vol_atom: np.ndarray, pressure: float) -> np.ndarray:
    """Per-atom enthalpy H = pe + P * V (eV/atom).

    ``pressure`` is the target pressure in bar; ``vol_atom`` is the per-atom
    volume in Ang^3/atom.  No division by lambda — this is a real-temperature
    ramp where lambda == 1.
    """
    return pe + pressure * BAR_ANG3_TO_EV * vol_atom


def _safe_div(num: float, den: float, fallback: float = 0.0) -> float:
    return num / den if abs(den) > 1e-300 else fallback


# ---------------------------------------------------------------------------
# Causal rolling-window helpers
# ---------------------------------------------------------------------------

def _rolling_mean(x: np.ndarray, w: int) -> np.ndarray:
    """Causal rolling mean.

    Returns NaN for the first w-1 positions and for any position whose
    w-sample window contains at least one non-finite value, avoiding silent
    NaN propagation through ``np.cumsum``.
    """
    out = np.full(len(x), np.nan)
    if w <= 0 or w > len(x):
        return out
    nan_mask = ~np.isfinite(x)
    x_fill = np.where(nan_mask, 0.0, x)
    cs  = np.cumsum(x_fill)
    cnn = np.cumsum(nan_mask.astype(np.int32))

    sums = cs[w - 1:] - np.concatenate([[0.0], cs[:len(x) - w]])
    n_nans = cnn[w - 1:] - np.concatenate([[0], cnn[:len(x) - w]])
    valid = n_nans == 0
    out[w - 1:] = np.where(valid, sums / w, np.nan)
    return out


def _rolling_var(x: np.ndarray, w: int) -> np.ndarray:
    """Causal rolling variance (population)."""
    rm  = _rolling_mean(x, w)
    rm2 = _rolling_mean(x ** 2, w)
    return np.maximum(rm2 - rm ** 2, 0.0)


def _rolling_cov(x: np.ndarray, y: np.ndarray, w: int) -> np.ndarray:
    """Causal rolling covariance."""
    return _rolling_mean(x * y, w) - _rolling_mean(x, w) * _rolling_mean(y, w)


# ---------------------------------------------------------------------------
# First-moment slope-break detector
# ---------------------------------------------------------------------------

def _slope_break_signal(
    y: np.ndarray,
    x: np.ndarray,
    valid_start: int,
    baseline_frac: float = 0.15,
    fit_order: int = 1,
    y_raw: Optional[np.ndarray] = None,
) -> Optional[dict]:
    """
    Fit a low-order polynomial y(x) on the early single-phase portion and
    return leverage-aware normalized residuals z = (y - y_fit) / pred_std.

    The slope-break signal is a *first-moment* detector: in a single phase
    y(x) follows a smooth equation of state captured by a low-order polynomial
    in x (here x is the *measured* temperature).  A first-order phase
    transition produces a step or kink in y(x), so |z| grows rapidly.  It fires
    as soon as the mean of y departs from the baseline EOS, typically earlier
    than a fluctuation peak, and is phase-agnostic (solid->solid as well as
    solid->liquid).

    Parameters
    ----------
    y             : signal array (smoothed enthalpy H or volume V), per-sample.
    x             : independent variable — the measured temperature (K).
    valid_start   : index at which rolling arrays become valid (NaN before).
    baseline_frac : fraction of valid samples used to fit the baseline EOS.
    fit_order     : polynomial order for the baseline fit (default 1, linear).
    y_raw         : optional unsmoothed counterpart of ``y``.  When provided,
                    the residual scale sigma is calibrated against the *raw*
                    thermal noise so heavy smoothing of ``y`` does not inflate z
                    (which would cause false positives on clean data).

    Returns
    -------
    dict with keys ``z``, ``sigma``, ``coef``, ``base_end`` — or ``None`` when
    too few finite samples are available to fit.
    """
    n = len(y)
    valid_n = n - valid_start
    if valid_n < 50:
        return None
    base_n = max(int(baseline_frac * valid_n), 50)
    base_end = min(valid_start + base_n, n)

    x_b = x[valid_start:base_end]
    y_b = y[valid_start:base_end]
    mask  = np.isfinite(x_b) & np.isfinite(y_b)
    min_pts = max(fit_order + 5, 10)
    if int(mask.sum()) < min_pts:
        return None

    # Leverage-aware normalization: a naive z = (y - fit) / sigma_baseline
    # overstates significance outside the baseline window because polynomial
    # extrapolation variance grows as dT^(2*fit_order).  Properly,
    #     pred_residual_var(x) = sigma**2 * (1 + h(x))
    # where h(x) is the design-matrix leverage.  This keeps |z|=5 calibrated
    # across the full ramp instead of blowing up on harmless extrapolation.
    x_bm = x_b[mask]
    y_bm = y_b[mask]

    # Center/scale x before building the Vandermonde matrix so X^T X is well
    # conditioned (raw T is 10^2-10^3 K; T^2 overflows for fit_order >= 2).
    x_center = float(np.mean(x_bm))
    x_scale  = float(np.std(x_bm))
    if not np.isfinite(x_scale) or x_scale < 1e-30:
        return None
    u_bm = (x_bm - x_center) / x_scale

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        Xb = np.vander(u_bm, fit_order + 1)
        try:
            coef, *_ = np.linalg.lstsq(Xb, y_bm, rcond=None)
            XtX_inv = np.linalg.inv(Xb.T @ Xb)
        except np.linalg.LinAlgError:
            return None

        resid_b = y_bm - Xb @ coef
        sigma_sq = float(np.sum(resid_b ** 2) /
                         max(len(y_bm) - (fit_order + 1), 1))
        if not np.isfinite(sigma_sq) or sigma_sq < 1e-60:
            return None
        # Calibrate sigma against raw (unsmoothed) thermal noise when available.
        if y_raw is not None:
            raw_b = y_raw[valid_start:base_end]
            raw_fit = Xb @ coef
            raw_resid = raw_b - raw_fit
            rmask = np.isfinite(raw_resid)
            if int(rmask.sum()) >= min_pts:
                sigma_sq_raw = float(np.sum(raw_resid[rmask] ** 2) /
                                      max(int(rmask.sum()) - (fit_order + 1), 1))
                if np.isfinite(sigma_sq_raw) and sigma_sq_raw > 1e-60:
                    sigma_sq = sigma_sq_raw

    sigma = float(np.sqrt(sigma_sq))

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        u_full = (x - x_center) / x_scale
        X_full = np.vander(u_full, fit_order + 1)
        y_fit = X_full @ coef
        h = np.einsum("ij,jk,ik->i", X_full, XtX_inv, X_full)
        h = np.maximum(h, 0.0)
        pred_std = np.sqrt(sigma_sq * (1.0 + h))
        z = (y - y_fit) / pred_std
    bad = ~(np.isfinite(y) & np.isfinite(x))
    z[bad] = np.nan
    if valid_start > 0:
        z[:valid_start] = np.nan

    return {"z": z, "sigma": sigma, "coef": coef, "base_end": base_end}


def _detect_slope_break_event(
    z: np.ndarray,
    T: np.ndarray,
    base_end: int,
    trigger_sigma: float,
    onset_sigma: float,
    persistence: int,
) -> Optional[tuple]:
    """
    Scan a normalized-residual array for the first sustained slope break.

    Returns (peak_idx, T_peak, peak_dev, onset_idx, T_onset) on success, or
    ``None`` when no sustained crossing is found.

    "Sustained" means the trigger threshold is crossed AND the signed residual
    keeps the same sign over all trailing data (mean magnitude >= 0.6 *
    trigger_sigma).  A real first-order transition produces a *permanent*
    one-directional shift of <y>(T) (latent heat / volume jump), so the signed
    residual stays the same sign from onset to the end of the ramp.  A zero-mean
    noise burst random-walks across the burst and then collapses toward zero, so
    its trailing signed mean is near zero and it is correctly rejected.
    """
    n = len(z)
    if base_end >= n:
        return None

    abs_z   = np.abs(z)
    region_abs    = abs_z[base_end:]
    region_signed = z[base_end:]

    candidates = np.where(np.isfinite(region_abs) & (region_abs > trigger_sigma))[0]
    if len(candidates) == 0:
        return None

    sustained_rel = None
    for c in candidates:
        chunk_signed = region_signed[c:]
        finite_mask = np.isfinite(chunk_signed)
        if int(finite_mask.sum()) < persistence:
            continue
        mean_signed = float(np.mean(chunk_signed[finite_mask]))
        if abs(mean_signed) < 0.6 * trigger_sigma:
            continue
        sustained_rel = int(c)
        break
    if sustained_rel is None:
        return None

    peak_idx_rel = sustained_rel
    onset_rel = peak_idx_rel
    for j in range(peak_idx_rel - 1, -1, -1):
        v = region_abs[j]
        if not np.isfinite(v):
            continue
        if v <= onset_sigma:
            break
        onset_rel = j

    peak_idx  = base_end + peak_idx_rel
    onset_idx = base_end + onset_rel
    if not (np.isfinite(T[peak_idx]) and np.isfinite(T[onset_idx])):
        return None
    return (
        peak_idx,
        float(T[peak_idx]),
        float(region_abs[peak_idx_rel]),
        onset_idx,
        float(T[onset_idx]),
    )


# ---------------------------------------------------------------------------
# Response-function arrays (measured temperature, no lambda)
# ---------------------------------------------------------------------------

def compute_response_arrays(
    pe: np.ndarray,
    press: np.ndarray,
    vol_total: np.ndarray,
    temp: np.ndarray,
    natoms: int,
    target_pressure: float,
    window_smooth: int = 500,
    window_fluct: int = 1000,
    slope_break_baseline_frac: float = 0.15,
    slope_break_fit_order: int = 1,
) -> dict:
    """
    Compute rolling response functions from a real-temperature ramp.

    Parameters
    ----------
    pe         : potential energy per atom (eV/atom)
    press      : pressure (bar)
    vol_total  : total volume (Ang^3)
    temp       : measured temperature (K) — used directly, not reconstructed
    natoms     : number of atoms
    target_pressure : target pressure (bar) used to form the enthalpy
    window_smooth   : smoothing window (rows) applied before fluctuations
    window_fluct    : rolling window for variance / covariance

    Returns
    -------
    dict with arrays T, Cp, kappa_T, alpha_P, H_z, V_z and the indices
    valid_start, H_base_end, V_base_end.
    """
    vol_atom = vol_total / natoms
    T = np.asarray(temp, dtype=float)
    H = enthalpy(pe, vol_atom, target_pressure)

    sm_H = _rolling_mean(H,        window_smooth)
    sm_V = _rolling_mean(vol_atom, window_smooth)
    sm_T = _rolling_mean(T,        window_smooth)

    mean_V = _rolling_mean(sm_V, window_fluct)
    # Use the smoothed measured temperature as the representative T per window.
    kBT  = KB_EV * sm_T
    kBT2 = kBT * sm_T

    with np.errstate(divide="ignore", invalid="ignore"):
        Cp    = _rolling_var(sm_H, window_fluct) / kBT2
        kappa = _rolling_var(sm_V, window_fluct) / (kBT * mean_V)
        alpha = np.abs(_rolling_cov(sm_V, sm_H, window_fluct) / (kBT2 * mean_V))

    valid_start = window_smooth + window_fluct - 2

    sb_H = _slope_break_signal(
        sm_H, sm_T, valid_start,
        baseline_frac=slope_break_baseline_frac,
        fit_order=slope_break_fit_order,
        y_raw=H,
    )
    sb_V = _slope_break_signal(
        sm_V, sm_T, valid_start,
        baseline_frac=slope_break_baseline_frac,
        fit_order=slope_break_fit_order,
        y_raw=vol_atom,
    )
    H_z = sb_H["z"] if sb_H is not None else np.full_like(Cp, np.nan)
    V_z = sb_V["z"] if sb_V is not None else np.full_like(Cp, np.nan)
    H_base_end = sb_H["base_end"] if sb_H is not None else valid_start
    V_base_end = sb_V["base_end"] if sb_V is not None else valid_start

    return {
        "T":           sm_T,
        "Cp":          Cp,
        "kappa_T":     kappa,
        "alpha_P":     alpha,
        "valid_start": valid_start,
        "H_z":         H_z,
        "V_z":         V_z,
        "H_base_end":  H_base_end,
        "V_base_end":  V_base_end,
    }


# ---------------------------------------------------------------------------
# The pre-flight scan
# ---------------------------------------------------------------------------

class RangeScan:
    """
    Decide a clean upper temperature for a production reversible-scaling sweep.

    Given the thermodynamic time series of a real-temperature ramp (T0 -> Tf
    under NPT), detect the onset of a phase transition and return the safe
    sub-range [T0, T_clean].  Stateless apart from the configured thresholds —
    call :meth:`find_clean_range` once on the full ramp data.
    """

    def __init__(
        self,
        target_pressure: float = 0.0,
        peak_threshold: float = 12.0,
        min_agreement: int = 2,
        onset_sigma: float = 4.0,
        window_smooth: int = 500,
        window_fluct: int = 1000,
        min_descent_frac: float = 0.3,
        tail_margin_frac: float = 0.05,
        slope_break_sigma: float = 5.0,
        slope_break_baseline_frac: float = 0.15,
        slope_break_fit_order: int = 1,
    ) -> None:
        self.target_pressure = target_pressure
        self.peak_threshold = peak_threshold
        self.min_agreement = min_agreement
        self.onset_sigma = onset_sigma
        self.window_smooth = window_smooth
        self.window_fluct = window_fluct
        self.min_descent_frac = min_descent_frac
        self.tail_margin_frac = tail_margin_frac
        self.slope_break_sigma = slope_break_sigma
        self.slope_break_baseline_frac = slope_break_baseline_frac
        self.slope_break_fit_order = slope_break_fit_order

    def find_clean_range(
        self,
        pe: np.ndarray,
        press: np.ndarray,
        vol_total: np.ndarray,
        temp: np.ndarray,
        natoms: int,
        t_start: float,
        t_stop: float,
    ) -> ScanResult:
        """
        Analyse the ramp time series and return a :class:`ScanResult`.

        Parameters
        ----------
        pe        : potential energy per atom (eV/atom)
        press     : pressure (bar)
        vol_total : total volume (Ang^3)
        temp      : measured temperature (K)
        natoms    : number of atoms
        t_start   : requested ramp start temperature (K)
        t_stop    : requested ramp stop temperature (K)

        Returns
        -------
        ScanResult.  ``transition_found == False`` and ``clean_t_stop ==
        t_stop`` when the ramp is clean over the whole range.
        """
        pe = np.asarray(pe, dtype=float)
        press = np.asarray(press, dtype=float)
        vol_total = np.asarray(vol_total, dtype=float)
        temp = np.asarray(temp, dtype=float)

        no_transition = ScanResult(
            transition_found=False,
            clean_t_stop=float(t_stop),
            onset_temperature=float(t_stop),
            peak_temperature=float(t_stop),
            triggered_signals=[],
            confidence=0.0,
        )

        if len(pe) < self.window_fluct:
            return no_transition

        vol_atom = vol_total / natoms
        T = temp
        H = enthalpy(pe, vol_atom, self.target_pressure)

        sm_H = _rolling_mean(H,        self.window_smooth)
        sm_V = _rolling_mean(vol_atom, self.window_smooth)
        sm_T = _rolling_mean(T,        self.window_smooth)

        mean_V = _rolling_mean(sm_V, self.window_fluct)
        kBT  = KB_EV * sm_T
        kBT2 = kBT * sm_T

        with np.errstate(divide="ignore", invalid="ignore"):
            Cp    = _rolling_var(sm_H, self.window_fluct) / kBT2
            kappa = _rolling_var(sm_V, self.window_fluct) / (kBT * mean_V)
            alpha = np.abs(_rolling_cov(sm_V, sm_H, self.window_fluct) / (kBT2 * mean_V))

        valid_start = self.window_smooth + self.window_fluct - 2
        n = len(Cp)
        if valid_start >= n:
            return no_transition

        def _robust_stats(sig):
            chunk = sig[valid_start:]
            finite = chunk[np.isfinite(chunk)]
            if len(finite) < 10:
                return 0.0, 0.0
            med = float(np.median(finite))
            mad = float(np.median(np.abs(finite - med)))
            return med, 1.4826 * mad

        signals = {
            "Cp":      (Cp,    *_robust_stats(Cp)),
            "kappa_T": (kappa, *_robust_stats(kappa)),
            "alpha_P": (alpha, *_robust_stats(alpha)),
        }

        detected = {}
        for sig_name, (sig, base_med, base_scale) in signals.items():
            if base_scale < 1e-40:
                continue
            threshold_val = base_med + self.peak_threshold * base_scale
            region = sig[valid_start:]
            T_region = sm_T[valid_start:]
            finite_mask = np.isfinite(region)
            if not np.any(finite_mask):
                continue
            masked = np.where(finite_mask, region, -np.inf)
            peak_val = float(np.max(masked))
            if peak_val <= threshold_val:
                continue

            peak_idx_rel = int(np.argmax(masked))
            n_region = len(region)
            margin = max(int(self.tail_margin_frac * n_region), 1)

            # Tail-margin guard: a peak in the last tail_margin_frac of the
            # data has no post-peak descent evidence (ramp barely crossed it).
            if peak_idx_rel >= n_region - margin:
                continue

            # Descent-fraction guard: the curve must actually come back down,
            # distinguishing a real transition (rise then liquid plateau) from
            # smooth monotonic Cp growth of a solid heated toward, but not into,
            # its melting point.
            post = region[peak_idx_rel + 1:]
            post = post[np.isfinite(post)]
            if post.size < max(3, self.window_fluct // 5):
                continue
            tail_n = max(self.window_fluct // 2, post.size // 2)
            tail_n = min(tail_n, post.size)
            trailing_level = float(np.median(post[-tail_n:]))
            excursion = peak_val - base_med
            descent = peak_val - trailing_level
            if excursion > 0 and descent < self.min_descent_frac * excursion:
                continue

            # Onset: walk back from the peak to where the signal last dropped
            # below onset_sigma above the baseline — the start of the rise.
            onset_threshold = base_med + self.onset_sigma * base_scale
            onset_idx_rel = peak_idx_rel
            for j in range(peak_idx_rel - 1, -1, -1):
                v = region[j]
                if not np.isfinite(v):
                    continue
                if v <= onset_threshold:
                    break
                onset_idx_rel = j
            T_onset = float(T_region[onset_idx_rel])
            T_peak = float(T_region[peak_idx_rel])
            mod_z = (peak_val - base_med) / base_scale
            detected[sig_name] = (valid_start + peak_idx_rel, T_peak, mod_z,
                                   valid_start + onset_idx_rel, T_onset)

        # First-moment slope-break detection on H(T) and V(T).
        sb_persistence = max(2 * self.window_smooth, self.window_fluct)
        slope_break_signals = [
            ("H_break", sm_H, H),
            ("V_break", sm_V, vol_atom),
        ]
        for sb_name, sb_y, sb_y_raw in slope_break_signals:
            sb = _slope_break_signal(
                sb_y, sm_T, valid_start,
                baseline_frac=self.slope_break_baseline_frac,
                fit_order=self.slope_break_fit_order,
                y_raw=sb_y_raw,
            )
            if sb is None:
                continue
            ev = _detect_slope_break_event(
                sb["z"], sm_T, sb["base_end"],
                trigger_sigma=self.slope_break_sigma,
                onset_sigma=self.onset_sigma,
                persistence=sb_persistence,
            )
            if ev is None:
                continue
            peak_idx, T_peak_sb, peak_dev, onset_idx, T_onset_sb = ev
            detected[sb_name] = (peak_idx, T_peak_sb, peak_dev, onset_idx, T_onset_sb)

        if len(detected) < self.min_agreement:
            return no_transition

        # Merge: peak temperature is the mean (informational); onset is the
        # earliest across signals (conservatively-early safe boundary).
        T_peak = float(np.mean([v[1] for v in detected.values()]))
        earliest = min(detected.values(), key=lambda v: v[3])
        T_onset = float(earliest[4])
        sig_names = list(detected.keys())
        n_signals_total = len(signals) + len(slope_break_signals)

        warnings.warn(
            f"Pre-scan detected a phase transition at T ~ {T_peak:.1f} K "
            f"(onset ~ {T_onset:.1f} K); signals: {', '.join(sig_names)}.",
            UserWarning,
            stacklevel=2,
        )

        return ScanResult(
            transition_found=True,
            clean_t_stop=float(T_onset),
            onset_temperature=float(T_onset),
            peak_temperature=float(T_peak),
            triggered_signals=sig_names,
            confidence=len(detected) / n_signals_total,
        )


# ---------------------------------------------------------------------------
# Diagnostic plot
# ---------------------------------------------------------------------------

def plot_scan(
    pe: np.ndarray,
    press: np.ndarray,
    vol_total: np.ndarray,
    temp: np.ndarray,
    natoms: int,
    target_pressure: float,
    outpath: str,
    result: Optional[ScanResult] = None,
    window_smooth: int = 500,
    window_fluct: int = 1000,
) -> bool:
    """
    Render the pre-flight scan diagnostic signals to ``outpath`` (a PNG).

    This is a *diagnostic aid only* and is written to be **best-effort / never
    fatal**: any failure (matplotlib missing, a headless-backend problem, an
    unwritable path, degenerate data) is swallowed and reported by returning
    ``False``.  A successful save returns ``True``.  To avoid touching global
    pyplot / backend state inside a compute job it draws through matplotlib's
    object-oriented API on an explicit Agg canvas.

    The figure stacks the five detector signals against the (smoothed, measured)
    temperature, recomputed by :func:`compute_response_arrays` with the *same*
    windows the scan used so the plot shows exactly what the detector saw:

      * variance-based response functions  Cp, kappa_T, alpha_P (a peak flags a
        transition);
      * first-moment slope-break residuals H_z, V_z, in units of sigma (a
        sustained excursion past the +/-slope_break_sigma band flags a break).

    When ``result`` is a transition-positive :class:`ScanResult`, the detected
    onset (the safe boundary) and peak temperatures are marked on every panel,
    and the panels of the triggered signals are drawn in red.

    Parameters
    ----------
    pe, press, vol_total, temp : array-like
        The ramp time series (the same arrays passed to
        :meth:`RangeScan.find_clean_range`).
    natoms : int
        Number of atoms.
    target_pressure : float
        Target pressure (bar) used to form the per-atom enthalpy.
    outpath : str
        Destination PNG path.
    result : ScanResult, optional
        Scan outcome; used only to annotate the onset / peak and colour the
        triggered panels.
    window_smooth, window_fluct : int
        Should match the windows the scan used so the plotted signals are the
        ones the detector actually evaluated.

    Returns
    -------
    bool
        ``True`` if the figure was written, ``False`` on any error.
    """
    try:
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg

        arr = compute_response_arrays(
            np.asarray(pe, dtype=float),
            np.asarray(press, dtype=float),
            np.asarray(vol_total, dtype=float),
            np.asarray(temp, dtype=float),
            natoms,
            target_pressure,
            window_smooth=window_smooth,
            window_fluct=window_fluct,
        )

        T = arr["T"]
        vstart = int(arr["valid_start"])
        if T is None or len(T) <= vstart:
            return False

        Tvalid = T[vstart:]
        finiteT = Tvalid[np.isfinite(Tvalid)]
        if finiteT.size < 2:
            return False
        xlo, xhi = float(np.min(finiteT)), float(np.max(finiteT))

        transition = bool(getattr(result, "transition_found", False)) if result is not None else False
        triggered = set(getattr(result, "triggered_signals", []) or [])

        # (ylabel, signal array, is_slope_break_residual, baseline_end, signal_key)
        panels = [
            (r"$C_P$  [eV/K/atom]",       arr["Cp"],      False, None,              "Cp"),
            (r"$\kappa_T$  [1/eV]",       arr["kappa_T"], False, None,              "kappa_T"),
            (r"$\alpha_P$  [1/K]",        arr["alpha_P"], False, None,              "alpha_P"),
            (r"H slope-break  [$\sigma$]", arr["H_z"],    True,  arr["H_base_end"], "H_break"),
            (r"V slope-break  [$\sigma$]", arr["V_z"],    True,  arr["V_base_end"], "V_break"),
        ]

        fig = Figure(figsize=(8.0, 10.0))
        FigureCanvasAgg(fig)
        axes = fig.subplots(len(panels), 1, sharex=True)

        for i, (ylabel, sig, is_z, base_end, key) in enumerate(panels):
            ax = axes[i]
            fired = key in triggered
            ax.plot(T, sig, lw=0.9, color="tab:red" if fired else "tab:blue")
            ax.set_ylabel(ylabel, fontsize=9)
            ax.tick_params(labelsize=8)
            ax.grid(True, alpha=0.25)

            if is_z:
                # slope-break residual: zero line, the +/-5 sigma trigger band
                # (the RangeScan default slope_break_sigma), and a shaded span
                # over the baseline-EOS fit region.
                ax.axhline(0.0, color="0.6", lw=0.6)
                for s in (-5.0, 5.0):
                    ax.axhline(s, color="0.7", lw=0.6, ls=":")
                if base_end is not None:
                    be = int(base_end)
                    if vstart < be < len(T) and np.isfinite(T[vstart]) and np.isfinite(T[be]):
                        ax.axvspan(float(T[vstart]), float(T[be]),
                                   color="0.85", alpha=0.5, lw=0)

            if transition:
                ax.axvline(float(result.onset_temperature), color="tab:green",
                           lw=1.2, ls="--", label="onset" if i == 0 else None)
                ax.axvline(float(result.peak_temperature), color="tab:purple",
                           lw=1.2, ls="-.", label="peak" if i == 0 else None)

        axes[-1].set_xlabel("temperature  [K]", fontsize=9)
        axes[0].set_xlim(xlo, xhi)

        if transition:
            axes[0].legend(loc="upper left", fontsize=8, framealpha=0.7)
            title = ("pre-scan: transition DETECTED  —  onset %.1f K, peak %.1f K  "
                     "(conf %.0f%%)\nsignals: %s"
                     % (result.onset_temperature, result.peak_temperature,
                        100.0 * result.confidence,
                        ", ".join(result.triggered_signals)))
        else:
            title = "pre-scan: no transition detected over [%.1f, %.1f] K" % (xlo, xhi)
        fig.suptitle(title, fontsize=9)
        fig.subplots_adjust(hspace=0.12, top=0.94)

        fig.savefig(outpath, dpi=120)
        return True
    except Exception:
        return False
