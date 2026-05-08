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
    confidence: float         # fraction of active signals that triggered
    temperature: float = 0.0  # estimated transition temperature (K)


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
# Temperature-block planner
# ---------------------------------------------------------------------------

def plan_temperature_blocks(
    t0: float,
    tf: float,
    total_steps: int,
    window_K: float,
) -> List[dict]:
    """
    Partition a reversible-scaling sweep into temperature-based blocks.

    Each block covers ``window_K`` Kelvin of the sweep.  The sweep direction
    (heating vs. cooling) is inferred from ``t0`` and ``tf``.

    The mapping from temperature to step number uses the relationship
    ``λ = T0 / T`` and the linear ramp ``λ(step) = li + (lf - li) * step / total_steps``,
    where ``li = 1`` and ``lf = T0 / Tf``.

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
    >>> plan_temperature_blocks(1200, 2000, 100000, 400)
    [{'temp': 1200, 'lambda': 1.0, 'step': 0},
     {'temp': 1600, 'lambda': 0.75, 'step': 33333},
     {'temp': 2000, 'lambda': 0.6,  'step': 100000}]
    """
    if window_K <= 0:
        raise ValueError("window_K must be > 0")
    if total_steps <= 0:
        raise ValueError("total_steps must be > 0")

    li = 1.0
    lf = t0 / tf

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

    result = []
    for temp in checkpoints:
        lam = t0 / temp
        # λ(step) = li + (lf - li) * step / total_steps  →  solve for step
        step = int((lam - li) * total_steps / (lf - li))
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
    peak_threshold: float = 5.0,
    baseline_frac: float = 0.2,
    min_signal_agreement: int = 2,
    sweep_label: str = "forward",
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
    peak_threshold       : flag when peak > peak_threshold x baseline median
    baseline_frac        : fraction of valid data to use as baseline
    min_signal_agreement : signals (Cp, kappa_T, alpha_P) that must peak simultaneously
    sweep_label          : label for warning messages

    Returns
    -------
    List of TransitionEvent.  Empty if no transition detected.
    A UserWarning is also emitted for each detected event.

    Notes
    -----
    The actual simulation temperature during a reversible-scaling sweep is::

        T = t_start + t_stop * (1 - lambda)

    This is consistent with the LAMMPS thermostat variable
    ``blambda * t_stop`` (forward sweep) or ``flambda * t_stop`` (backward).
    """
    if len(dU) < window_fluct:
        return []

    vol_atom = vol_total / natoms

    # Actual simulation temperature: T = t_start + t_stop * (1 - lambda)
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

    # Baseline: median over first baseline_frac of valid data.
    # Hl has NaN for the first window_smooth-1 positions (from the sm_lam
    # rolling mean), so rolling_var(Hl) is only valid from position
    # window_smooth + window_fluct - 2 onwards.
    valid_start = window_smooth + window_fluct - 2
    n = len(Cp)
    if valid_start >= n:
        return []

    n_base = max(1, int(baseline_frac * (n - valid_start)))

    def _base_median(sig):
        chunk = sig[valid_start: valid_start + n_base]
        finite = chunk[np.isfinite(chunk)]
        return float(np.median(finite)) if len(finite) > 0 else 0.0

    base_Cp    = _base_median(Cp)
    base_kappa = _base_median(kappa)
    base_alpha = _base_median(alpha)

    signals = {
        "Cp":      (Cp,    base_Cp),
        "kappa_T": (kappa, base_kappa),
        "alpha_P": (alpha, base_alpha),
    }

    # Find peak of each signal; flag if it exceeds threshold x baseline
    detected = {}
    for sig_name, (sig, base) in signals.items():
        if base <= 1e-30:
            continue
        region = sig[valid_start:]
        T_region = T[valid_start:]
        finite_mask = np.isfinite(region)
        if not np.any(finite_mask):
            continue
        masked = np.where(finite_mask, region, -np.inf)
        peak_val = float(np.max(masked))
        if peak_val > peak_threshold * base:
            peak_idx_rel = int(np.argmax(masked))
            T_trans = float(T_region[peak_idx_rel])
            detected[sig_name] = (valid_start + peak_idx_rel, T_trans, peak_val / base)

    if len(detected) < min_signal_agreement:
        return []

    T_trans  = float(np.mean([v[1] for v in detected.values()]))
    peak_idx = int(np.mean([v[0] for v in detected.values()]))
    sig_names = list(detected.keys())
    ratio_str = ", ".join(f"{k}={v[2]:.1f}x" for k, v in detected.items())

    event = TransitionEvent(
        sample_index=peak_idx,
        triggered_signals=sig_names,
        confidence=len(detected) / len(signals),
        temperature=T_trans,
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

    return {
        "T":           T,
        "Cp":          Cp,
        "kappa_T":     kappa,
        "alpha_P":     alpha,
        "valid_start": valid_start,
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
                "Detection can be tuned via the transition_detector input "
                "block or disabled with transition_detector.enabled: false."
            )
        else:
            raise self._SolidifiedError(
                "Phase transition detected (liquid -> solid) at sample "
                f"{event.sample_index} (T ~ {event.temperature:.1f} K). "
                f"Triggered signals: {signals_str} "
                f"(confidence {event.confidence:.0%}). "
                "System likely solidified -- consider increasing temperature.\n"
                "Detection can be tuned via the transition_detector input "
                "block or disabled with transition_detector.enabled: false."
            )
