"""
Fluctuation-based phase transition detector for calphy.

Uses thermodynamic response functions computed from ensemble fluctuations to
detect phase transitions without structural analysis (Steinhardt parameters,
solid fraction, etc.) and without being limited to specific crystal structures.

Derived quantities (fluctuation–dissipation theorem):
  - Cv   (NVT-style):  Var(U) / (kB T²)
  - Cp   (NPT-style):  Var(H) / (kB T²),  H = U + P_target·V
  - alpha_P           : Cov(V, H) / (kB T² <V>)
  - kappa_T           : Var(V) / (kB T <V>)
  - plain variances   : var_pe, var_H, var_vol, var_etotal
  - mean shifts       : mean_H, mean_pe, mean_vol

Detection strategy
------------------
Two non-overlapping sliding windows over the accumulated record buffer:
  - baseline window  : early samples (post burn-in)
  - recent window    : latest ``recent_window`` samples

A transition is flagged when *at least* ``min_agreement`` out of the active
signals simultaneously exceed their respective thresholds:
  1. Variance-ratio test     : Var_recent / Var_baseline  > variance_ratio_threshold
  2. CUSUM mean-shift test   : |mean_recent - mean_baseline| / std_baseline
                                 > mean_shift_threshold
  3. Response-function peak  : |value_recent| / |value_baseline|
                                 > cv_peak_threshold  (Cv, Cp, alpha, kappa_T)
"""

import warnings
from dataclasses import dataclass, field
from typing import List, Optional, Literal

import numpy as np

# Boltzmann constant in eV/K (LAMMPS metal units)
KB_EV: float = 8.617333262e-5

# 1 bar·Å³ = 6.2415091e-7 eV  (pressure × volume → energy, metal units)
BAR_ANG3_TO_EV: float = 6.2415091e-7


# ---------------------------------------------------------------------------
# Data record
# ---------------------------------------------------------------------------

@dataclass
class ThermoRecord:
    """One time-averaged sample from the MD run (all extensive quantities per atom)."""
    step: float
    temp: float       # temperature (K)
    pe: float         # potential energy per atom (eV/atom)
    etotal: float     # total energy per atom (eV/atom)
    vol: float        # volume per atom (Å³/atom)
    press: float      # pressure (bar)


@dataclass
class TransitionEvent:
    """Returned by PhaseTransitionMonitor when a transition is detected."""
    sample_index: int
    triggered_signals: List[str]
    confidence: float   # fraction of active signals that triggered


# ---------------------------------------------------------------------------
# Fluctuation buffer
# ---------------------------------------------------------------------------

class FluctuationBuffer:
    """
    Append-only buffer that stores ThermoRecord entries.

    All arrays grow monotonically; slicing returns a lightweight view object
    so that window computations avoid unnecessary copies.
    """

    def __init__(self) -> None:
        self._steps:   List[float] = []
        self._temp:    List[float] = []
        self._pe:      List[float] = []
        self._etotal:  List[float] = []
        self._vol:     List[float] = []
        self._press:   List[float] = []

    # ------------------------------------------------------------------
    def append(self, record: ThermoRecord) -> None:
        self._steps.append(float(record.step))
        self._temp.append(float(record.temp))
        self._pe.append(float(record.pe))
        self._etotal.append(float(record.etotal))
        self._vol.append(float(record.vol))
        self._press.append(float(record.press))

    def __len__(self) -> int:
        return len(self._steps)

    # ------------------------------------------------------------------
    # Array accessors (always return numpy arrays)
    # ------------------------------------------------------------------

    def steps(self) -> np.ndarray:
        return np.asarray(self._steps)

    def temp(self) -> np.ndarray:
        return np.asarray(self._temp)

    def pe(self) -> np.ndarray:
        return np.asarray(self._pe)

    def etotal(self) -> np.ndarray:
        return np.asarray(self._etotal)

    def vol(self) -> np.ndarray:
        return np.asarray(self._vol)

    def press(self) -> np.ndarray:
        return np.asarray(self._press)

    def enthalpy(self, target_pressure: float) -> np.ndarray:
        """
        H per atom (eV/atom) = pe + P_target * vol.

        Uses the *target* pressure (a constant in NPT) rather than the
        fluctuating instantaneous pressure to avoid pressure-noise coupling.
        """
        return self.pe() + target_pressure * BAR_ANG3_TO_EV * self.vol()

    # ------------------------------------------------------------------
    def slice(self, start: int, end: Optional[int] = None) -> "FluctuationBuffer":
        """Return a new FluctuationBuffer containing only samples [start:end]."""
        sub = FluctuationBuffer()
        sl = slice(start, end)
        for s, t, p, e, v, pr in zip(
            self._steps[sl],
            self._temp[sl],
            self._pe[sl],
            self._etotal[sl],
            self._vol[sl],
            self._press[sl],
        ):
            sub.append(ThermoRecord(step=s, temp=t, pe=p, etotal=e, vol=v, press=pr))
        return sub


# ---------------------------------------------------------------------------
# Response function computation
# ---------------------------------------------------------------------------

def _safe_div(numerator: float, denominator: float, fallback: float = 0.0) -> float:
    if abs(denominator) < 1e-300:
        return fallback
    return numerator / denominator


def compute_response_functions(
    buf: "FluctuationBuffer",
    target_pressure: float,
    temperature: float,
) -> Optional[dict]:
    """
    Compute ensemble averages and fluctuation-based response functions.

    Parameters
    ----------
    buf : FluctuationBuffer
    target_pressure : float
        Target pressure in bar (used for H = pe + P_target * vol).
    temperature : float
        Representative temperature in K (used for kB T² denominator).

    Returns
    -------
    dict with keys:
        var_pe, var_H, var_vol, var_etotal,
        cv, cp, alpha, kappa_T,
        mean_pe, mean_H, mean_vol, mean_temp
    or None if fewer than 2 samples are available.
    """
    if len(buf) < 2:
        return None

    T = max(temperature, 1.0)
    kT2 = KB_EV * T ** 2
    kT  = KB_EV * T

    pe     = buf.pe()
    H      = buf.enthalpy(target_pressure)
    vol    = buf.vol()
    etotal = buf.etotal()

    var_pe     = float(np.var(pe,     ddof=1))
    var_H      = float(np.var(H,      ddof=1))
    var_vol    = float(np.var(vol,    ddof=1))
    var_etotal = float(np.var(etotal, ddof=1))

    mean_vol  = float(np.mean(vol))
    mean_H    = float(np.mean(H))
    mean_pe   = float(np.mean(pe))
    mean_temp = float(np.mean(buf.temp()))

    cv       = _safe_div(var_pe,  kT2)
    cp       = _safe_div(var_H,   kT2)
    cov_VH   = float(np.mean(vol * H) - mean_vol * mean_H)
    alpha    = _safe_div(cov_VH,  kT2 * mean_vol)
    kappa_T  = _safe_div(var_vol, kT  * mean_vol) if kT > 0 else 0.0

    return dict(
        var_pe=var_pe,
        var_H=var_H,
        var_vol=var_vol,
        var_etotal=var_etotal,
        cv=cv,
        cp=cp,
        alpha=alpha,
        kappa_T=kappa_T,
        mean_pe=mean_pe,
        mean_H=mean_H,
        mean_vol=mean_vol,
        mean_temp=mean_temp,
    )


# ---------------------------------------------------------------------------
# Change-point detector
# ---------------------------------------------------------------------------

class OnlineChangePointDetector:
    """
    Detects phase transitions from a FluctuationBuffer using three tests:

    1. **Variance-ratio**  :  Var_recent / Var_baseline  > variance_ratio_threshold
    2. **CUSUM mean-shift**:  |mean_recent - mean_baseline| / std_baseline
                                > mean_shift_threshold
    3. **Response-peak**   :  |stat_recent| / |stat_baseline|  > cv_peak_threshold
                              (applied to Cv, Cp, alpha, kappa_T)

    A transition is reported when ``min_agreement`` or more tests fire
    simultaneously.

    Parameters
    ----------
    active_signals : list[str]
        Signals to evaluate.  Any subset of ALL_SIGNALS.
    baseline_window : int
        Number of samples in the baseline window (taken from early data).
    recent_window : int
        Number of samples in the recent window (taken from latest data).
    min_samples_before_check : int
        Minimum total samples before any check is attempted.
    variance_ratio_threshold : float
        Threshold for the variance-ratio test.
    mean_shift_threshold : float
        Threshold (in units of baseline σ) for the mean-shift test.
    cv_peak_threshold : float
        Threshold ratio for response-function peak test.
    min_agreement : int
        Minimum number of signals that must trigger simultaneously.
    target_pressure : float
        Target pressure in bar (passed to enthalpy calculation).
    temperature : float
        Simulation temperature in K.
    """

    VARIANCE_SIGNALS  = frozenset({"var_pe", "var_H", "var_vol", "var_etotal"})
    RESPONSE_SIGNALS  = frozenset({"cv", "cp", "alpha", "kappa_T"})
    MEAN_SHIFT_SIGNALS = frozenset({"mean_H", "mean_pe", "mean_vol"})

    ALL_SIGNALS: List[str] = [
        "var_pe", "var_H", "var_vol", "var_etotal",
        "cv", "cp", "alpha", "kappa_T",
        "mean_H", "mean_pe", "mean_vol",
    ]

    def __init__(
        self,
        active_signals: List[str],
        baseline_window: int = 50,
        recent_window: int = 50,
        min_samples_before_check: int = 100,
        variance_ratio_threshold: float = 5.0,
        mean_shift_threshold: float = 5.0,
        cv_peak_threshold: float = 5.0,
        min_agreement: int = 2,
        target_pressure: float = 0.0,
        temperature: float = 300.0,
    ) -> None:
        self.active_signals             = list(active_signals)
        self.baseline_window            = baseline_window
        self.recent_window              = recent_window
        self.min_samples_before_check   = min_samples_before_check
        self.variance_ratio_threshold   = variance_ratio_threshold
        self.mean_shift_threshold       = mean_shift_threshold
        self.cv_peak_threshold          = cv_peak_threshold
        self.min_agreement              = min_agreement
        self.target_pressure            = target_pressure
        self.temperature                = temperature

    # ------------------------------------------------------------------
    def check(self, buf: "FluctuationBuffer") -> Optional[TransitionEvent]:
        """
        Inspect current buffer for a phase transition.

        Returns a TransitionEvent if one is detected, otherwise None.
        """
        n = len(buf)
        if n < self.min_samples_before_check:
            return None

        # baseline: skip first 10 % of data as burn-in, take next baseline_window samples
        burn_in   = max(0, n // 10)
        base_start = burn_in
        base_end   = base_start + self.baseline_window

        # ensure windows do not overlap
        if base_end > n - self.recent_window:
            return None

        base_buf   = buf.slice(base_start, base_end)
        recent_buf = buf.slice(n - self.recent_window, n)

        T_base   = max(float(np.mean(base_buf.temp())),   1.0)
        T_recent = max(float(np.mean(recent_buf.temp())), 1.0)

        base_stats   = compute_response_functions(base_buf,   self.target_pressure, T_base)
        recent_stats = compute_response_functions(recent_buf, self.target_pressure, T_recent)

        if base_stats is None or recent_stats is None:
            return None

        triggered: List[str] = []

        for sig in self.active_signals:
            if sig in self.VARIANCE_SIGNALS:
                b_val = base_stats.get(sig, 0.0)
                r_val = recent_stats.get(sig, 0.0)
                if b_val > 1e-30 and _safe_div(r_val, b_val) > self.variance_ratio_threshold:
                    triggered.append(sig)

            elif sig in self.RESPONSE_SIGNALS:
                b_val = abs(base_stats.get(sig, 0.0))
                r_val = abs(recent_stats.get(sig, 0.0))
                if b_val > 1e-30 and _safe_div(r_val, b_val) > self.cv_peak_threshold:
                    triggered.append(sig)

            elif sig in self.MEAN_SHIFT_SIGNALS:
                # retrieve the raw array for baseline and recent windows
                raw_key = sig[len("mean_"):]   # "H", "pe", or "vol"
                if raw_key == "H":
                    b_arr = base_buf.enthalpy(self.target_pressure)
                    r_arr = recent_buf.enthalpy(self.target_pressure)
                elif raw_key == "pe":
                    b_arr = base_buf.pe()
                    r_arr = recent_buf.pe()
                elif raw_key == "vol":
                    b_arr = base_buf.vol()
                    r_arr = recent_buf.vol()
                else:
                    continue

                b_mean = float(np.mean(b_arr))
                b_std  = float(np.std(b_arr, ddof=1))
                r_mean = float(np.mean(r_arr))
                if b_std > 1e-30:
                    shift = abs(r_mean - b_mean) / b_std
                    if shift > self.mean_shift_threshold:
                        triggered.append(sig)

        if not self.active_signals:
            return None

        if len(triggered) >= self.min_agreement:
            confidence = len(triggered) / len(self.active_signals)
            return TransitionEvent(
                sample_index=n,
                triggered_signals=triggered,
                confidence=confidence,
            )

        return None


# ---------------------------------------------------------------------------
# High-level monitor
# ---------------------------------------------------------------------------

class PhaseTransitionMonitor:
    """
    Streaming monitor that accumulates ThermoRecords and checks for phase
    transitions using fluctuation-based statistics.

    Designed as a drop-in replacement for the pyscal3 ``find_solid_fraction``
    snapshot checks in calphy.  The same ``MeltedError`` / ``SolidifiedError``
    exceptions are raised so downstream code (routines.py) is unchanged.

    Parameters
    ----------
    expected_phase : 'solid' or 'liquid'
        The phase that should remain stable.  A detected transition raises
        ``MeltedError`` for 'solid' or ``SolidifiedError`` for 'liquid'.
    target_pressure : float
        Target simulation pressure in bar.
    temperature : float
        Target simulation temperature in K.
    active_signals : list[str] or None
        Signals to include; defaults to all available signals.
    **detector_kwargs
        Extra keyword arguments forwarded to OnlineChangePointDetector.
    """

    ALL_SIGNALS: List[str] = OnlineChangePointDetector.ALL_SIGNALS

    def __init__(
        self,
        expected_phase: Literal["solid", "liquid"],
        target_pressure: float = 0.0,
        temperature: float = 300.0,
        active_signals: Optional[List[str]] = None,
        **detector_kwargs,
    ) -> None:
        # local import to avoid circular dependency
        from calphy.errors import MeltedError, SolidifiedError
        self._MeltedError    = MeltedError
        self._SolidifiedError = SolidifiedError

        self.expected_phase = expected_phase
        self.buffer = FluctuationBuffer()

        if active_signals is None:
            active_signals = list(self.ALL_SIGNALS)

        self.detector = OnlineChangePointDetector(
            active_signals=active_signals,
            target_pressure=target_pressure,
            temperature=temperature,
            **detector_kwargs,
        )

    # ------------------------------------------------------------------
    def update(
        self,
        record: ThermoRecord,
        raise_on_transition: bool = True,
    ) -> Optional[TransitionEvent]:
        """
        Append one record and check for a transition.

        Parameters
        ----------
        record : ThermoRecord
        raise_on_transition : bool
            When True (default), immediately raise MeltedError or
            SolidifiedError.  When False, return the event for the
            caller to handle.

        Returns
        -------
        TransitionEvent or None
        """
        self.buffer.append(record)
        event = self.detector.check(self.buffer)
        if event is None:
            return None
        if raise_on_transition:
            self._raise_for_event(event)
        return event

    # ------------------------------------------------------------------
    def evaluate_final(
        self, raise_on_transition: bool = True
    ) -> Optional[TransitionEvent]:
        """
        Run the detector once on the full accumulated buffer.

        Used as a drop-in for snapshot-based ``check_if_melted`` /
        ``check_if_solidfied`` calls at equilibration checkpoints.
        Returns None (no error raised) if the buffer is too small for
        detection — this preserves the previous behaviour of skipping
        the check when ``tolerance.solid_fraction == 0``.
        """
        event = self.detector.check(self.buffer)
        if event is None:
            return None
        if raise_on_transition:
            self._raise_for_event(event)
        return event

    # ------------------------------------------------------------------
    def reset(self) -> None:
        """Clear accumulated data (e.g. between reversible-scaling iterations)."""
        self.buffer = FluctuationBuffer()

    # ------------------------------------------------------------------
    def _raise_for_event(self, event: TransitionEvent) -> None:
        signals_str = ", ".join(event.triggered_signals)
        if self.expected_phase == "solid":
            raise self._MeltedError(
                f"Phase transition detected (solid → liquid) at sample "
                f"{event.sample_index}. "
                f"Triggered signals: {signals_str} "
                f"(confidence {event.confidence:.0%}). "
                f"System likely melted — consider reducing temperature or "
                f"increasing system size.\n"
                f"Detection can be tuned via the transition_detector input "
                f"block or disabled with transition_detector.enabled: false."
            )
        else:
            raise self._SolidifiedError(
                f"Phase transition detected (liquid → solid) at sample "
                f"{event.sample_index}. "
                f"Triggered signals: {signals_str} "
                f"(confidence {event.confidence:.0%}). "
                f"System likely solidified — consider increasing temperature.\n"
                f"Detection can be tuned via the transition_detector input "
                f"block or disabled with transition_detector.enabled: false."
            )
