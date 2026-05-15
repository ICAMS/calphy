"""
Fluctuation-based phase transition detector for calphy.

Computes three thermodynamic response functions from rolling-window
fluctuations of lambda-reduced enthalpy H/lambda and volume V:

    C_P     = Var(H/lambda) / (k_B T**2)         [eV/K per atom]
    kappa_T = Var(V)   / (k_B T <V>)             [1/eV]
    alpha_P = Cov(V, H/lambda) / (k_B T**2 <V>)  [1/K]

For equilibration runs (lambda == 1), H/lambda = H, recovering the standard
NPT fluctuation-dissipation expressions.

Two detection modes
-------------------
1. Online (streaming) -- PhaseTransitionMonitor:
   Accumulates ThermoRecord entries one-by-one; compares response functions
   computed over the first ``baseline_window`` samples versus the most recent
   ``recent_window`` samples.  Used by phase.py equilibration checks
   (check_if_melted / check_if_solidfied).

2. Post-hoc (ts sweep) -- detect_ts_transitions:
   Analyses a full ts.forward / ts.backward dataset after the MD run finishes.
   Uses rolling windows to compute response functions and emits a warning when
   peaks exceed ``peak_threshold x baseline`` for at least
   ``min_signal_agreement`` signals simultaneously.
"""

import warnings
from dataclasses import dataclass
from typing import List, Optional, Literal

import numpy as np

# Boltzmann constant in eV/K (LAMMPS metal units)
KB_EV: float = 8.617333262e-5

# 1 bar * Ang^3 = 6.2415091e-7 eV
BAR_ANG3_TO_EV: float = 6.2415091e-7


# ---------------------------------------------------------------------------
# Data record
# ---------------------------------------------------------------------------

@dataclass
class ThermoRecord:
    """One time-averaged sample from the MD run (extensive quantities per atom)."""

    step: float
    temp: float       # temperature (K)
    pe: float         # potential energy per atom (eV/atom)
    etotal: float     # total energy per atom (eV/atom)
    vol: float        # volume per atom (Ang^3/atom)
    press: float      # pressure (bar)
    lam: float = 1.0  # lambda scaling parameter; default 1 for equilibration runs


@dataclass
class TransitionEvent:
    """Returned when a phase transition is detected."""

    sample_index: int
    triggered_signals: List[str]
    confidence: float              # fraction of active signals that triggered
    temperature: float = 0.0      # peak temperature of the response function (K)
    onset_temperature: float = 0.0  # first-rise onset of the transition (K)
                                    # Always <= temperature for heating sweeps.
                                    # This is the safe recovery point — it
                                    # corresponds to where the signal first
                                    # exceeded a low 1-sigma threshold above
                                    # the baseline, i.e. the actual start of
                                    # the transition, not the peak top.


# ---------------------------------------------------------------------------
# Fluctuation buffer
# ---------------------------------------------------------------------------

class FluctuationBuffer:
    """Append-only buffer that stores ThermoRecord entries."""

    def __init__(self) -> None:
        self._steps:   List[float] = []
        self._temp:    List[float] = []
        self._pe:      List[float] = []
        self._etotal:  List[float] = []
        self._vol:     List[float] = []
        self._press:   List[float] = []
        self._lam:     List[float] = []

    def append(self, record: ThermoRecord) -> None:
        self._steps.append(float(record.step))
        self._temp.append(float(record.temp))
        self._pe.append(float(record.pe))
        self._etotal.append(float(record.etotal))
        self._vol.append(float(record.vol))
        self._press.append(float(record.press))
        self._lam.append(float(record.lam))

    def __len__(self) -> int:
        return len(self._steps)

    # ------------------------------------------------------------------
    # Array accessors
    # ------------------------------------------------------------------

    def steps(self)  -> np.ndarray: return np.asarray(self._steps)
    def temp(self)   -> np.ndarray: return np.asarray(self._temp)
    def pe(self)     -> np.ndarray: return np.asarray(self._pe)
    def etotal(self) -> np.ndarray: return np.asarray(self._etotal)
    def vol(self)    -> np.ndarray: return np.asarray(self._vol)
    def press(self)  -> np.ndarray: return np.asarray(self._press)
    def lam(self)    -> np.ndarray: return np.asarray(self._lam)

    def enthalpy(self, target_pressure: float) -> np.ndarray:
        """H per atom (eV/atom) = pe + P_target * BAR_ANG3_TO_EV * vol."""
        return self.pe() + target_pressure * BAR_ANG3_TO_EV * self.vol()

    def enthalpy_reduced(self, target_pressure: float) -> np.ndarray:
        """H/lambda per atom -- removes Hamiltonian scaling trend.  Equal to H when lambda=1."""
        lam = self.lam()
        lam_safe = np.where(np.abs(lam) > 1e-30, lam, 1.0)
        return self.enthalpy(target_pressure) / lam_safe

    def slice(self, start: int, end: Optional[int] = None) -> "FluctuationBuffer":
        """Return a new FluctuationBuffer containing only samples [start:end]."""
        sub = FluctuationBuffer()
        sl = slice(start, end)
        for s, t, p, e, v, pr, l in zip(
            self._steps[sl], self._temp[sl], self._pe[sl],
            self._etotal[sl], self._vol[sl], self._press[sl], self._lam[sl],
        ):
            sub.append(ThermoRecord(step=s, temp=t, pe=p, etotal=e, vol=v, press=pr, lam=l))
        return sub


# ---------------------------------------------------------------------------
# Window-level response functions
# ---------------------------------------------------------------------------

def _safe_div(num: float, den: float, fallback: float = 0.0) -> float:
    return num / den if abs(den) > 1e-300 else fallback


def _window_response_functions(
    Hl: np.ndarray,
    vol: np.ndarray,
    T: float,
) -> dict:
    """
    Compute Cp, kappa_T, alpha_P for a fixed window of data.

    Parameters
    ----------
    Hl  : H/lambda per atom (eV/atom)
    vol : volume per atom (Ang^3/atom)
    T   : representative temperature (K)
    """
    T = max(T, 1.0)
    kBT  = KB_EV * T
    kBT2 = kBT * T

    n = len(Hl)
    var_Hl  = float(np.var(Hl,  ddof=1)) if n > 1 else 0.0
    var_V   = float(np.var(vol, ddof=1)) if n > 1 else 0.0
    mean_V  = float(np.mean(vol))
    mean_Hl = float(np.mean(Hl))
    cov_VH  = float(np.mean(vol * Hl) - mean_V * mean_Hl)

    Cp      = _safe_div(var_Hl, kBT2)
    kappa_T = _safe_div(var_V, kBT * mean_V)
    alpha_P = abs(_safe_div(cov_VH, kBT2 * mean_V))

    return {"Cp": Cp, "kappa_T": kappa_T, "alpha_P": alpha_P}


# ---------------------------------------------------------------------------
# Rolling helpers for ts-sweep post-hoc detection
# ---------------------------------------------------------------------------

def _rolling_mean(x: np.ndarray, w: int) -> np.ndarray:
    """Causal rolling mean.

    Returns NaN for the first w-1 positions and for any position whose
    w-sample window contains at least one NaN value.  This avoids silent
    NaN propagation through ``np.cumsum`` when the input has a NaN prefix
    (e.g. after dividing smoothed H by smoothed lambda).
    """
    out = np.full(len(x), np.nan)
    nan_mask = ~np.isfinite(x)
    x_fill = np.where(nan_mask, 0.0, x)  # treat NaN as 0 in sum
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
# First-moment slope-break detector (phase-agnostic)
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
    Fit a low-order polynomial y(x) on the early single-phase portion of the
    valid region and return normalized residuals z = (y - y_fit) / sigma_resid
    over the full data.

    The slope-break signal is a *first-moment* phase-transition detector: in a
    single-phase regime y(x) follows a smooth equation of state captured by
    a low-order polynomial in x.  Any phase transition (solid->solid,
    solid->liquid) produces a step or kink in y(x), causing |z| to grow
    rapidly.

    Unlike fluctuation-based detectors (Cp, kappa_T, alpha_P), this signal
    does not require a variance peak to accumulate — it fires as soon as the
    *mean* of y departs from the baseline EOS.  This catches transitions
    earlier and is phase-agnostic: it relies only on the assumption that
    y(x) is smooth within a single phase, which is true for solid->solid as
    well as solid->liquid transitions.

    Parameters
    ----------
    y             : signal array (e.g. enthalpy H or volume V), per-sample.
    x             : independent variable, typically the *equivalent
                    temperature* T(λ).  Fitting in T rather than λ keeps the
                    baseline polynomial well-behaved over the extrapolated
                    range — V(T) and H(T) are nearly linear in T for a
                    solid, whereas V(λ) and H/λ(λ) are hyperbolic and a
                    polynomial extrapolation breaks down.
    valid_start   : index at which the smoothed/rolling-window arrays become
                    valid (NaN before this).
    baseline_frac : fraction of valid samples (starting at valid_start) used
                    to fit the baseline polynomial.  Should be small enough
                    that the baseline lies entirely in the initial phase
                    (default 0.15 ≈ 15%).
    fit_order     : polynomial order for the baseline fit.  Default 1
                    (linear).  Higher orders increase the false-positive
                    rate because polynomial extrapolation error grows as
                    ΔT^order: with a baseline that covers only ~15% of the
                    sweep range, a quadratic fit can easily produce
                    multi-σ "deviations" purely from extrapolation runaway.
                    Real first-order transitions produce slope changes
                    that linear fits catch cleanly; use order ≥2 only on
                    very wide sweeps where anharmonic curvature dominates.
    y_raw         : optional unsmoothed counterpart of ``y``.  When provided,
                    ``sigma`` (the residual scale used to normalize z) is
                    computed from raw-data residuals at the baseline.  This
                    is essential when ``y`` is heavily smoothed: smoothing
                    reduces the per-sample residual std by a factor of √W,
                    which would inflate z and cause false positives on
                    clean data.  Calibrating σ against the *raw* thermal
                    noise restores the intended "5σ" interpretation.

    Returns
    -------
    dict with keys
      ``z``        — normalized residual array, length len(y); NaN before
                     valid_start or wherever inputs are non-finite.
      ``sigma``    — baseline residual standard deviation (used to normalize).
      ``coef``     — polynomial coefficients (np.polyfit order; highest first).
      ``base_end`` — exclusive end index of the baseline window.
    or ``None`` if too few finite samples are available to fit.
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

    # ------------------------------------------------------------------
    # Leverage-aware normalization.
    #
    # A naive z = (y - fit) / sigma_baseline overstates the significance
    # of points outside the baseline window because polynomial
    # extrapolation variance grows as ΔT^(2 * fit_order).  Properly,
    #
    #     pred_residual_var(x) = sigma² * (1 + h(x))
    #
    # where h(x) = v(x)^T (X_b^T X_b)^{-1} v(x) is the design-matrix
    # leverage at point x and v(x) is the Vandermonde row.  This grows
    # polynomially with extrapolation distance, so |z|=5 stays
    # calibrated everywhere on the sweep instead of blowing up on
    # harmless extrapolated noise.
    # ------------------------------------------------------------------
    x_bm = x_b[mask]
    y_bm = y_b[mask]

    # Center and scale x before building the Vandermonde matrix.  Raw T
    # values are 10²–10³ K, so T² is 10⁴–10⁶ and X^T X overflows / is
    # ill-conditioned for fit_order ≥ 2 (numpy emits overflow + invalid
    # warnings on matmul, and the leverage term loses precision).  Working
    # in u = (x - μ) / σ gives well-conditioned columns of order unity;
    # predictions are identical under the change of variables.
    x_center = float(np.mean(x_bm))
    x_scale  = float(np.std(x_bm))
    if not np.isfinite(x_scale) or x_scale < 1e-30:
        return None
    u_bm = (x_bm - x_center) / x_scale

    # Build the baseline Vandermonde matrix (highest power first, matching
    # np.polyfit's coefficient convention) and solve in least-squares.
    # Some BLAS backends emit spurious overflow/invalid warnings during
    # well-conditioned matmuls on Apple Silicon; suppress while doing the
    # provably-finite baseline math.
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        Xb = np.vander(u_bm, fit_order + 1)
        try:
            coef, *_ = np.linalg.lstsq(Xb, y_bm, rcond=None)
            XtX_inv = np.linalg.inv(Xb.T @ Xb)
        except np.linalg.LinAlgError:
            return None

        # Smoothed-data residual variance on the baseline.  Biased low by
        # autocorrelation when ``y`` is heavily smoothed, but the leverage
        # term compensates: the test statistic (y - y_fit) is itself made
        # of smoothed (autocorrelated) data, so the smoothed-baseline σ is
        # the matching scale.  When ``y_raw`` is supplied we override σ²
        # with the raw thermal noise estimate, which gives the strict-sense
        # per-sample fluctuation scale.
        resid_b = y_bm - Xb @ coef
        sigma_sq = float(np.sum(resid_b ** 2) /
                         max(len(y_bm) - (fit_order + 1), 1))
        if not np.isfinite(sigma_sq) or sigma_sq < 1e-60:
            return None
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

    # Predicted values and leverage h(x) at every sample, in the centered/
    # scaled u-coordinate.  Predictions are identical to the unscaled basis.
    # Suppress overflow / invalid warnings from NaN rows of ``x`` (positions
    # before valid_start propagate through np.vander); those positions are
    # masked to NaN at the end anyway.
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

    return {
        "z": z,
        "sigma": sigma,
        "coef": coef,
        "base_end": base_end,
    }


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

    Returns (peak_idx, T_peak, peak_dev, onset_idx, T_onset) on success,
    where peak_idx/onset_idx are absolute indices into the full data array.
    Returns None when no sustained crossing is found.

    "Sustained" means the trigger threshold is crossed AND, over the next
    ``persistence`` finite samples:
      (a) the mean of |z| stays at least 0.6 * trigger_sigma, AND
      (b) ≥80 % of the signed-z values share the sign of the trigger sample.

    Condition (b) is the key discriminator for zero-mean noise: a real
    first-order phase transition produces a one-directional shift of <y>(T)
    (latent heat positive on melting, volume positive on melting, etc.), so
    the signed residual stays the same sign post-onset.  A noise excursion
    that briefly crosses the trigger threshold but oscillates around zero
    afterwards (sample-mean drift in a noisy tail region) has near-50 %
    sign agreement and is correctly rejected.
    """
    n = len(z)
    if base_end >= n:
        return None

    abs_z   = np.abs(z)
    region_abs    = abs_z[base_end:]
    region_signed = z[base_end:]
    n_region = len(region_abs)

    candidates = np.where(np.isfinite(region_abs) & (region_abs > trigger_sigma))[0]
    if len(candidates) == 0:
        return None

    sustained_rel = None
    for c in candidates:
        # Use ALL trailing data from the trigger onward (not a fixed
        # ``persistence`` window) when checking signed-mean persistence.
        # A real first-order transition produces a *permanent* shift of
        # <y>(T): the latent-heat or volume jump persists from onset to
        # the end of the sweep, so the signed mean over the entire
        # trailing region remains ~trigger_sigma in magnitude.  A
        # zero-mean noise burst (localized variance increase, no mean
        # shift) random-walks across the burst region and then returns
        # to baseline noise afterwards, so the trailing signed mean
        # collapses toward zero once the burst ends.  The ``persistence``
        # argument now functions as the *minimum number of finite
        # trailing samples* required to accept the trigger — i.e. a
        # head-of-tail margin that ensures we have enough data to make
        # the call.
        chunk_signed = region_signed[c:]
        finite_mask = np.isfinite(chunk_signed)
        if int(finite_mask.sum()) < persistence:
            continue
        chunk_signed_f = chunk_signed[finite_mask]
        mean_signed = float(np.mean(chunk_signed_f))
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
# Temperature-block planner
# ---------------------------------------------------------------------------

def plan_temperature_blocks(
    t0: float,
    tf: float,
    total_steps: int,
    window_K: float,
    lambda_schedule: str = "linear",
) -> List[dict]:
    """
    Partition a reversible-scaling sweep into temperature-based blocks.

    Each block covers ``window_K`` Kelvin of the sweep.  The sweep direction
    (heating vs. cooling) is inferred from ``t0`` and ``tf``.

    The step-to-temperature mapping depends on ``lambda_schedule``:

    * ``"linear"`` (default): LAMMPS ramps λ linearly in step,
      so λ(s) = λ₀ + (λ_f − λ₀) × s / total_steps, and temperature is
      non-linear in step:
      ``s(T) = (t0/T − 1) / (t0/tf − 1) × total_steps``

    * ``"uniform_temperature"``: temperature is ramped linearly in step,
      so ``s(T) = (T − t0) / (tf − t0) × total_steps``

    Parameters
    ----------
    t0 : float
        Starting temperature of the sweep (K).  This is always
        ``self.calc._temperature`` (the reference temperature).
    tf : float
        Ending temperature of the sweep (K).
    total_steps : int
        Total number of MD steps in the sweep (``_n_sweep_steps``).
    window_K : float
        Desired temperature width of each block (K).  Must be > 0.
    lambda_schedule : str, optional
        Lambda schedule used by LAMMPS.  ``"linear"`` (default) or
        ``"uniform_temperature"``.  Must match the schedule passed to the
        LAMMPS run so checkpoint step numbers are consistent with the actual
        simulation.

    Returns
    -------
    list of dict
        Each entry has keys ``"temp"`` (K), ``"lambda"`` (dimensionless),
        and ``"step"`` (int, clamped to [0, total_steps]).  The list includes
        both the starting and ending checkpoints, so there are
        ``len(result) - 1`` blocks.  The final entry always has
        ``"step" == total_steps`` and ``"temp" == tf``.

    Examples
    --------
    Linear lambda schedule (default):

    >>> plan_temperature_blocks(1200, 2000, 100000, 400)
    [{'temp': 1200, 'lambda': 1.0, 'step': 0},
     {'temp': 1600, 'lambda': 0.75, 'step': 62500},
     {'temp': 2000, 'lambda': 0.6,  'step': 100000}]

    Uniform-temperature schedule:

    >>> plan_temperature_blocks(1200, 2000, 100000, 400, lambda_schedule='uniform_temperature')
    [{'temp': 1200, 'lambda': 1.0, 'step': 0},
     {'temp': 1600, 'lambda': 0.75, 'step': 40000},
     {'temp': 2000, 'lambda': 0.6,  'step': 100000}]
    """
    if window_K <= 0:
        raise ValueError("window_K must be > 0")
    if total_steps <= 0:
        raise ValueError("total_steps must be > 0")

    # Generate temperature checkpoints
    checkpoints: List[float] = []
    if tf > t0:
        # Heating: t0 → tf
        cur = t0
        while cur < tf:
            checkpoints.append(cur)
            cur += window_K
    else:
        # Cooling: t0 → tf
        cur = t0
        while cur > tf:
            checkpoints.append(cur)
            cur -= window_K

    # Always include the final temperature exactly
    if not checkpoints or checkpoints[-1] != tf:
        checkpoints.append(tf)

    # Pre-compute lambda endpoints for the linear-lambda formula.
    lam_start = 1.0          # t0 / t0
    lam_end   = t0 / tf      # lambda at T = tf

    result = []
    for temp in checkpoints:
        lam = t0 / temp
        if lambda_schedule == "uniform_temperature":
            # T is linear in step: s(T) = (T - t0) / (tf - t0) * total_steps
            step = int((temp - t0) / (tf - t0) * total_steps)
        else:
            # "linear" (default): lambda is linear in step.
            # lambda(s) = lam_start + (lam_end - lam_start) * s / total_steps
            # Inverting: s(T) = (lam - lam_start) / (lam_end - lam_start) * total_steps
            step = int((lam - lam_start) / (lam_end - lam_start) * total_steps)
        step = max(0, min(step, total_steps))
        result.append({"temp": temp, "lambda": lam, "step": step})

    # Guarantee the last checkpoint lands exactly on total_steps
    result[-1]["step"] = total_steps
    return result


# ---------------------------------------------------------------------------
# Post-hoc ts-sweep transition detector
# ---------------------------------------------------------------------------

def detect_ts_transitions(
    dU: np.ndarray,
    press: np.ndarray,
    vol_total: np.ndarray,
    lam: np.ndarray,
    t_start: float,
    t_stop: float,
    natoms: int,
    window_smooth: int = 500,
    window_fluct: int = 1000,
    peak_threshold: float = 8.0,
    baseline_frac: float = 0.2,
    min_signal_agreement: int = 2,
    min_descent_frac: float = 0.3,
    tail_margin_frac: float = 0.05,
    slope_break_sigma: float = 5.0,
    onset_sigma: float = 4.0,
    slope_break_baseline_frac: float = 0.15,
    slope_break_fit_order: int = 1,
    sweep_label: str = "forward",
    sweep_mode: str = "ts",
) -> List[TransitionEvent]:
    """
    Detect phase transitions in a reversible-scaling ts sweep.

    Parameters
    ----------
    dU        : pe/atom array (eV/atom), column 0 of ts file
    press     : pressure array (bar), column 1 of ts file
    vol_total : total volume array (Ang^3), column 2 of ts file
    lam       : lambda array, column 3 of ts file
    t_start   : starting temperature (K)
    t_stop    : ending temperature (K)
    natoms    : number of atoms
    window_smooth        : rows to smooth before computing fluctuations
    window_fluct         : rolling window for variance / covariance
    peak_threshold       : flag when modified Z-score of the peak exceeds this
                           value, where mod_z = (peak - median) / (1.4826 * MAD)
                           computed over all valid data accumulated so far.
                           Default 8.0 cleanly separates real transitions
                           (mod_z > 25) from natural Cp growth in a single
                           solid window (mod_z ~ 2-7).
    baseline_frac        : (kept for backward compatibility; unused with the
                           robust median+MAD criterion)
    min_signal_agreement : signals (Cp, kappa_T, alpha_P) that must peak simultaneously
    min_descent_frac     : after the peak, the signal must descend by at least
                           this fraction of the excursion (peak - baseline median)
                           for the peak to count.  Suppresses false positives in
                           which the response function rises monotonically toward
                           the end of the analysed window without ever turning
                           over (typical of a solid heated toward, but not into,
                           its melting transition).
    tail_margin_frac     : the peak must sit at least ``tail_margin_frac`` of
                           the valid samples away from *both* the beginning and
                           end of the data.  Prevents two classes of false
                           positive:
                           - "max is the last sample" (tail): only a rise is
                             observed, no descent evidence yet (solid near Tm).
                           - "max is the first sample" (head): rolling-window
                             startup transient at the very beginning of a sweep
                             (e.g. liquid backward sweep starting from a cold
                             equilibration) inflates variance temporarily.
    slope_break_sigma          : trigger threshold (in baseline-σ) for the
                                 first-moment slope-break signals on H/lambda
                                 and V.  A break is declared at the first
                                 sample whose normalized residual exceeds this
                                 value (with persistence).  These signals fire
                                 earlier than the variance-based Cp/κ/α peaks
                                 because they detect deviation of the *mean*
                                 from the single-phase equation-of-state rather
                                 than waiting for a fluctuation peak.  Default
                                 5.0.
    onset_sigma                : how far above the single-phase baseline noise
                                 a signal must be to count as "the transition
                                 has started."  Applied uniformly as a walk-back
                                 threshold for ALL signals: variance-based
                                 (Cp, kappa_T, alpha_P) use onset_sigma × MAD
                                 above the baseline median; slope-break signals
                                 (H_break, V_break) use onset_sigma × fit-σ.
                                 Lower values give earlier (more conservative)
                                 onsets; higher values place the onset closer to
                                 where the signal is unambiguous.  Default 4.0.
    slope_break_baseline_frac  : fraction of valid samples (starting at the
                                 first valid row) used to fit the single-phase
                                 EOS polynomial baseline.  Default 0.15.
    slope_break_fit_order      : polynomial order for the baseline fit in T.
                                 Default 1 (linear).  Higher orders increase
                                 the false-positive rate via extrapolation
                                 runaway from a narrow baseline window;
                                 first-order transitions are picked up
                                 cleanly by a linear fit.
    sweep_label          : label for warning messages
    sweep_mode           : ``'ts'`` (reversible scaling, fixed thermostat at
                           t_start, equivalent T = t_start / lambda) or
                           ``'tscale'`` (temperature scaling, ramping thermostat,
                           actual T = t_start + t_stop * (1 - lambda))

    Returns
    -------
    List of TransitionEvent.  Empty if no transition detected.
    A UserWarning is also emitted for each detected event.

    Notes
    -----
    For ``sweep_mode='ts'`` (reversible scaling) the thermostat is fixed at
    *t_start* throughout the sweep; lambda scales the Hamiltonian.  The
    equivalent free-energy reference temperature is::

        T = t_start / lambda

    For ``sweep_mode='tscale'`` (temperature scaling) the thermostat ramps
    from *t_start* to *t_stop* and the actual simulation temperature is::

        T = t_start + t_stop * (1 - lambda)
    """
    if len(dU) < window_fluct:
        return []

    vol_atom = vol_total / natoms

    # Equivalent / actual temperature depending on sweep mode
    if sweep_mode == "ts":
        # Fixed thermostat at t_start; equivalent free-energy T = t_start / lambda
        lam_safe_T = np.where(np.abs(lam) > 1e-30, lam, np.nan)
        T = t_start / lam_safe_T
    else:
        # Ramping thermostat: T = t_start + t_stop * (1 - lambda)
        T = t_start + t_stop * (1.0 - lam)

    H = dU + press * vol_atom * BAR_ANG3_TO_EV  # eV/atom

    # Smooth to reduce thermal noise before computing fluctuations
    sm_H   = _rolling_mean(H,        window_smooth)
    sm_lam = _rolling_mean(lam,      window_smooth)
    sm_V   = _rolling_mean(vol_atom, window_smooth)

    # H/lambda -- removes Hamiltonian-scaling trend
    lam_safe = np.where(np.abs(sm_lam) > 1e-30, sm_lam, np.nan)
    Hl = sm_H / lam_safe

    # Rolling response functions
    mean_V = _rolling_mean(sm_V, window_fluct)
    kBT  = KB_EV * T
    kBT2 = kBT * T

    with np.errstate(divide="ignore", invalid="ignore"):
        Cp    = _rolling_var(Hl,   window_fluct) / kBT2
        kappa = _rolling_var(sm_V, window_fluct) / (kBT * mean_V)
        alpha = np.abs(_rolling_cov(sm_V, Hl, window_fluct) / (kBT2 * mean_V))

    # Baseline statistics: median + MAD over ALL valid data accumulated so far.
    # Using a global robust statistic (rather than the cold-start "first 20%")
    # is essential for incremental block-by-block detection: response functions
    # of a solid grow ~8x naturally over a 200 K temperature window
    # (equilibrium physics), so a "first 20% baseline" makes everything later
    # look like an outlier. Median + MAD is robust against the peak itself
    # influencing the baseline.
    valid_start = window_smooth + window_fluct - 2
    n = len(Cp)
    if valid_start >= n:
        return []

    def _robust_stats(sig):
        chunk = sig[valid_start:]
        finite = chunk[np.isfinite(chunk)]
        if len(finite) < 10:
            return 0.0, 0.0
        med = float(np.median(finite))
        mad = float(np.median(np.abs(finite - med)))
        # 1.4826 * MAD = consistent estimate of std for normal data
        return med, 1.4826 * mad

    base_Cp_med,    base_Cp_scale    = _robust_stats(Cp)
    base_kappa_med, base_kappa_scale = _robust_stats(kappa)
    base_alpha_med, base_alpha_scale = _robust_stats(alpha)

    signals = {
        "Cp":      (Cp,    base_Cp_med,    base_Cp_scale),
        "kappa_T": (kappa, base_kappa_med, base_kappa_scale),
        "alpha_P": (alpha, base_alpha_med, base_alpha_scale),
    }

    # Modified Z-score outlier test: flag when (peak - median) / scale > peak_threshold.
    # This robust criterion correctly distinguishes:
    #   - smooth monotonic growth of Cp through a solid window (mod_z ~ 2)
    #   - sharp transition spike (mod_z > 25)
    # while being insensitive to the magnitude of the peak (no near-zero median
    # divisors, no first-20% baseline that becomes obsolete as data accumulates).
    #
    # Three additional shape guards are required to suppress false positives in
    # incremental block-by-block detection:
    #   (a) tail-margin: the peak must not sit in the last `tail_margin_frac`
    #       of valid samples — this fires only when T_stop is very close to the
    #       transition (sweep barely crosses the peak, no post-peak data for
    #       the descent estimate).  The descent-fraction guard handles the
    #       general case; tail-margin is a hard cutoff for that edge case.
    #   (b) descent-fraction: after the peak the signal must have dropped by
    #       at least `min_descent_frac` of the excursion (peak - baseline).
    #       This is the primary discriminator: smooth monotonic Cp growth in a
    #       solid (which never turns over) fails; a real first-order transition
    #       (Cp rises sharply then settles to a higher liquid plateau) passes.
    #   Note: a head-margin guard (symmetric to tail-margin) is NOT applied.
    #   In post-hoc detection the rolling-window warmup (valid_start rows) has
    #   already excluded startup transients; an additional head cutoff is
    #   redundant and would suppress real early transitions.
    detected = {}
    for sig_name, (sig, base_med, base_scale) in signals.items():
        if base_scale < 1e-40:
            continue
        threshold_val = base_med + peak_threshold * base_scale
        region = sig[valid_start:]
        T_region = T[valid_start:]
        finite_mask = np.isfinite(region)
        if not np.any(finite_mask):
            continue
        masked = np.where(finite_mask, region, -np.inf)
        peak_val = float(np.max(masked))
        if peak_val <= threshold_val:
            continue

        peak_idx_rel = int(np.argmax(masked))
        n_region = len(region)
        margin = max(int(tail_margin_frac * n_region), 1)

        # Tail-margin guard: peak in the last tail_margin_frac of the data.
        # In post-hoc detection on the full sweep, this fires only when T_stop
        # is very close to the transition — the sweep barely crosses the peak
        # and there is insufficient post-peak data for a reliable descent
        # estimate.  The descent-fraction check below handles the general case;
        # this guard is a hard cutoff for that edge case.
        if peak_idx_rel >= n_region - margin:
            continue

        # Descent-fraction guard: the curve must have actually come back
        #     down by the end of the data, not just dipped briefly due to
        #     noise.  Compare the peak to the *median* of the trailing
        #     samples after the peak (a robust estimate of the post-peak
        #     plateau) rather than the noisy minimum.
        post = region[peak_idx_rel + 1:]
        post = post[np.isfinite(post)]
        if post.size < max(3, window_fluct // 5):
            continue
        tail_n = max(window_fluct // 2, post.size // 2)
        tail_n = min(tail_n, post.size)
        trailing_level = float(np.median(post[-tail_n:]))
        excursion = peak_val - base_med
        descent = peak_val - trailing_level
        if excursion > 0 and descent < min_descent_frac * excursion:
            continue

        # Find the onset: scan backward from the peak to the first index
        # where the signal dropped back below 1-sigma above the baseline.
        # This marks where the transition *started* rising, which corresponds
        # far more closely to the actual Tm than the peak top does.
        onset_threshold = base_med + onset_sigma * base_scale
        onset_idx_rel = peak_idx_rel
        for j in range(peak_idx_rel - 1, -1, -1):
            v = region[j]
            if not np.isfinite(v):
                continue
            if v <= onset_threshold:
                break
            onset_idx_rel = j
        T_onset = float(T_region[onset_idx_rel])

        T_trans = float(T_region[peak_idx_rel])
        mod_z = (peak_val - base_med) / base_scale
        detected[sig_name] = (valid_start + peak_idx_rel, T_trans, mod_z,
                               valid_start + onset_idx_rel, T_onset)

    # ------------------------------------------------------------------
    # First-moment slope-break detection (phase-agnostic).
    #
    # Detects deviation of <H/lambda>(lambda) and <V>(lambda) from the
    # single-phase equation-of-state fit through the early baseline.  Fires
    # earlier than the variance peaks and applies to any first-order
    # transition (solid->solid, solid->liquid) because it relies only on the
    # smoothness of y(lam) within a single phase.
    # ------------------------------------------------------------------
    # Persistence: must exceed the autocorrelation length of the smoothed
    # signal (~2 * window_smooth) so an autocorrelated noise excursion that
    # briefly crosses the trigger threshold cannot fake a sustained shift.
    sb_persistence = max(2 * window_smooth, window_fluct)
    # First-moment signals: fit y(T) (equivalent temperature) rather than
    # y(lambda) — single-phase EOS curves H(T), V(T) are nearly linear in T,
    # while the corresponding y(lambda) functions are hyperbolic in lambda and
    # cannot be reliably extrapolated by a low-order polynomial.  Use the
    # raw smoothed enthalpy sm_H (not sm_H/lambda) because the latent-heat
    # jump at a transition is intrinsic to H itself; dividing by lambda only
    # adds the trivial Hamiltonian scaling back into the signal.
    slope_break_signals = [
        ("H_break", sm_H,  H),
        ("V_break", sm_V,  vol_atom),
    ]
    for sb_name, sb_y, sb_y_raw in slope_break_signals:
        sb = _slope_break_signal(
            sb_y, T, valid_start,
            baseline_frac=slope_break_baseline_frac,
            fit_order=slope_break_fit_order,
            y_raw=sb_y_raw,
        )
        if sb is None:
            continue
        ev = _detect_slope_break_event(
            sb["z"], T, sb["base_end"],
            trigger_sigma=slope_break_sigma,
            onset_sigma=onset_sigma,
            persistence=sb_persistence,
        )
        if ev is None:
            continue
        peak_idx, T_peak_sb, peak_dev, onset_idx, T_onset_sb = ev
        detected[sb_name] = (peak_idx, T_peak_sb, peak_dev, onset_idx, T_onset_sb)

    if len(detected) < min_signal_agreement:
        return []

    # Merge across triggered signals.  Use the MINIMUM onset (earliest sign
    # of trouble) rather than the mean — the user-visible failure mode is
    # "didn't stop in time", so the conservatively-early recovery point is
    # preferred when signals disagree.  Peak temperature stays as the mean
    # for informational purposes.
    T_trans   = float(np.mean([v[1] for v in detected.values()]))
    peak_idx  = int(np.mean([v[0] for v in detected.values()]))
    earliest  = min(detected.values(), key=lambda v: v[3])
    onset_idx = int(earliest[3])
    T_onset   = float(earliest[4])
    sig_names = list(detected.keys())
    # Cp/kappa_T/alpha_P report modified Z-score of the variance peak;
    # H_break/V_break report normalized residual sigma of the slope-break.
    def _fmt_metric(name, val):
        if name in ("H_break", "V_break"):
            return f"{name}=|z| {val:.1f}σ"
        return f"{name}=modZ {val:.1f}"
    ratio_str = ", ".join(_fmt_metric(k, v[2]) for k, v in detected.items())

    n_signals_total = len(signals) + len(slope_break_signals)
    event = TransitionEvent(
        sample_index=peak_idx,
        triggered_signals=sig_names,
        confidence=len(detected) / n_signals_total,
        temperature=T_trans,
        onset_temperature=T_onset,
    )

    warnings.warn(
        f"Phase transition detected in {sweep_label} sweep at T ~ {T_trans:.1f} K. "
        f"Response function peaks above baseline: {ratio_str}. "
        f"Verify output files for structural changes.",
        UserWarning,
        stacklevel=2,
    )

    return [event]


# ---------------------------------------------------------------------------
# Response-function arrays for post-run plotting
# ---------------------------------------------------------------------------

def compute_ts_response_arrays(
    dU: np.ndarray,
    press: np.ndarray,
    vol_total: np.ndarray,
    lam: np.ndarray,
    t_start: float,
    t_stop: float,
    natoms: int,
    window_smooth: int = 500,
    window_fluct: int = 1000,
    sweep_mode: str = "ts",
    slope_break_baseline_frac: float = 0.15,
    slope_break_fit_order: int = 1,
) -> dict:
    """
    Compute rolling response functions from a ts sweep file.

    Uses the same algorithm as :func:`detect_ts_transitions` but returns the
    raw arrays (T, Cp, kappa_T, alpha_P) for plotting rather than running peak
    detection.

    Parameters
    ----------
    dU, press, vol_total, lam : arrays as read from ts.forward_N.dat
    t_start, t_stop           : temperature endpoints (K)
    natoms                    : number of atoms
    window_smooth             : smoothing window (rows)
    window_fluct              : fluctuation window (rows)
    sweep_mode                : ``'ts'`` or ``'tscale'`` (see
                                :func:`detect_ts_transitions`)

    Returns
    -------
    dict with keys:
      T           – temperature array (K), shape (N,)
      Cp          – heat capacity, eV/K/atom, shape (N,)
      kappa_T     – isothermal compressibility, 1/eV, shape (N,)
      alpha_P     – thermal expansion coeff, 1/K, shape (N,)
      valid_start – index from which the rolling windows are fully populated
    """
    vol_atom = vol_total / natoms
    if sweep_mode == "ts":
        lam_safe_T = np.where(np.abs(lam) > 1e-30, lam, np.nan)
        T = t_start / lam_safe_T
    else:
        T = t_start + t_stop * (1.0 - lam)

    H = dU + press * vol_atom * BAR_ANG3_TO_EV

    sm_H   = _rolling_mean(H,        window_smooth)
    sm_lam = _rolling_mean(lam,      window_smooth)
    sm_V   = _rolling_mean(vol_atom, window_smooth)

    lam_safe = np.where(np.abs(sm_lam) > 1e-30, sm_lam, np.nan)
    Hl = sm_H / lam_safe

    mean_V = _rolling_mean(sm_V, window_fluct)
    kBT  = KB_EV * T
    kBT2 = kBT * T

    with np.errstate(divide="ignore", invalid="ignore"):
        Cp    = _rolling_var(Hl,   window_fluct) / kBT2
        kappa = _rolling_var(sm_V, window_fluct) / (kBT * mean_V)
        alpha = np.abs(_rolling_cov(sm_V, Hl, window_fluct) / (kBT2 * mean_V))

    valid_start = window_smooth + window_fluct - 2

    # First-moment slope-break residuals (normalized by baseline sigma).
    # Fit y(T) on raw enthalpy and volume — see _slope_break_signal docstring.
    sb_H = _slope_break_signal(
        sm_H, T, valid_start,
        baseline_frac=slope_break_baseline_frac,
        fit_order=slope_break_fit_order,
        y_raw=H,
    )
    sb_V = _slope_break_signal(
        sm_V, T, valid_start,
        baseline_frac=slope_break_baseline_frac,
        fit_order=slope_break_fit_order,
        y_raw=vol_atom,
    )
    H_z      = sb_H["z"]        if sb_H is not None else np.full_like(Cp, np.nan)
    V_z      = sb_V["z"]        if sb_V is not None else np.full_like(Cp, np.nan)
    H_base_end = sb_H["base_end"] if sb_H is not None else valid_start
    V_base_end = sb_V["base_end"] if sb_V is not None else valid_start

    return {
        "T":           T,
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
# Online (streaming) monitor -- used during equilibration
# ---------------------------------------------------------------------------

class PhaseTransitionMonitor:
    """
    Streaming monitor that accumulates ThermoRecord entries and detects phase
    transitions using thermodynamic response functions (Cp, kappa_T, alpha_P).

    Detects transitions by comparing response functions computed over the first
    ``baseline_window`` samples versus the most recent ``recent_window``
    samples.  A transition is flagged when at least ``min_signal_agreement``
    signals exceed ``peak_threshold x baseline``.

    Parameters
    ----------
    expected_phase : 'solid' or 'liquid'
    target_pressure : float
        Target pressure in bar (used for H = pe + P_target * V * BAR_ANG3_TO_EV).
    temperature : float
        Target simulation temperature in K.
    baseline_window : int
        Number of early samples to use as baseline.
    recent_window : int
        Number of latest samples to compare against baseline.
    min_samples_before_check : int
        Minimum buffer size before any check is attempted.
    peak_threshold : float
        Flag when recent signal > peak_threshold x baseline signal.
    min_signal_agreement : int
        Number of signals that must trigger simultaneously.
    """

    def __init__(
        self,
        expected_phase,
        target_pressure: float = 0.0,
        temperature: float = 300.0,
        baseline_window: int = 50,
        recent_window: int = 50,
        min_samples_before_check: int = 100,
        peak_threshold: float = 5.0,
        min_signal_agreement: int = 2,
    ) -> None:
        from calphy.errors import MeltedError, SolidifiedError
        self._MeltedError     = MeltedError
        self._SolidifiedError = SolidifiedError

        self.expected_phase           = expected_phase
        self.target_pressure          = target_pressure
        self.temperature              = temperature
        self.baseline_window          = baseline_window
        self.recent_window            = recent_window
        self.min_samples_before_check = min_samples_before_check
        self.peak_threshold           = peak_threshold
        self.min_signal_agreement     = min_signal_agreement

        self.buffer = FluctuationBuffer()

    # ------------------------------------------------------------------
    def update(
        self,
        record,
        raise_on_transition: bool = True,
    ):
        """Append one record; full evaluation happens at evaluate_final."""
        self.buffer.append(record)
        return None

    # ------------------------------------------------------------------
    def evaluate_final(
        self,
        raise_on_transition: bool = True,
    ):
        """
        Run the detector on the full accumulated buffer.

        Returns None (without raising) when the buffer has too few samples.
        """
        n = len(self.buffer)
        if n < self.min_samples_before_check:
            return None

        base_buf   = self.buffer.slice(0, self.baseline_window)
        recent_buf = self.buffer.slice(n - self.recent_window, n)

        base_Hl  = base_buf.enthalpy_reduced(self.target_pressure)
        base_vol = base_buf.vol()
        base_T   = max(float(np.mean(base_buf.temp())), 1.0)

        rec_Hl  = recent_buf.enthalpy_reduced(self.target_pressure)
        rec_vol = recent_buf.vol()
        rec_T   = max(float(np.mean(recent_buf.temp())), 1.0)

        base_stats   = _window_response_functions(base_Hl, base_vol, base_T)
        recent_stats = _window_response_functions(rec_Hl,  rec_vol,  rec_T)

        triggered = []
        for sig in ("Cp", "kappa_T", "alpha_P"):
            b = base_stats[sig]
            r = recent_stats[sig]
            if b > 1e-30 and _safe_div(r, b) > self.peak_threshold:
                triggered.append(sig)

        if len(triggered) < self.min_signal_agreement:
            return None

        event = TransitionEvent(
            sample_index=n,
            triggered_signals=triggered,
            confidence=len(triggered) / 3.0,
            temperature=rec_T,
        )
        if raise_on_transition:
            self._raise_for_event(event)
        return event

    # ------------------------------------------------------------------
    def reset(self) -> None:
        """Clear accumulated data (e.g. between reversible-scaling iterations)."""
        self.buffer = FluctuationBuffer()

    # ------------------------------------------------------------------
    def _raise_for_event(self, event) -> None:
        signals_str = ", ".join(event.triggered_signals)
        if self.expected_phase == "solid":
            raise self._MeltedError(
                "Phase transition detected (solid -> liquid) at sample "
                f"{event.sample_index} (T ~ {event.temperature:.1f} K). "
                f"Triggered signals: {signals_str} "
                f"(confidence {event.confidence:.0%}). "
                "System likely melted -- consider reducing temperature or "
                "increasing system size.\n"
                "Detection can be tuned via the phase_transition_detection input "
                "block or disabled with phase_transition_detection.mode: none."
            )
        else:
            raise self._SolidifiedError(
                "Phase transition detected (liquid -> solid) at sample "
                f"{event.sample_index} (T ~ {event.temperature:.1f} K). "
                f"Triggered signals: {signals_str} "
                f"(confidence {event.confidence:.0%}). "
                "System likely solidified -- consider increasing temperature.\n"
                "Detection can be tuned via the phase_transition_detection input "
                "block or disabled with phase_transition_detection.mode: none."
            )
