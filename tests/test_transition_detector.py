"""
Unit tests for calphy.transition_detector

Tests cover:
  - FluctuationBuffer: append, slice, enthalpy_reduced (lambda-aware)
  - _window_response_functions: correctness of Cp, kappa_T, alpha_P
  - PhaseTransitionMonitor: stable baseline gives no event; step-jump triggers
  - detect_ts_transitions: synthetic ramp data with injected peak
"""

import numpy as np
import pytest
import warnings

from calphy.transition_detector import (
    FluctuationBuffer,
    ThermoRecord,
    TransitionEvent,
    PhaseTransitionMonitor,
    detect_ts_transitions,
    plan_temperature_blocks,
    compute_ts_response_arrays,
    KB_EV,
    BAR_ANG3_TO_EV,
    _window_response_functions,
    _rolling_mean,
    _rolling_var,
    _rolling_cov,
    _slope_break_signal,
    _detect_slope_break_event,
)
from calphy.errors import MeltedError, SolidifiedError


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_record(step, temp=1000.0, pe=-4.0, etotal=-3.5, vol=12.0, press=0.0, lam=1.0):
    return ThermoRecord(step=step, temp=temp, pe=pe, etotal=etotal, vol=vol, press=press, lam=lam)


def _fill_buffer(n, temp=1000.0, pe=-4.0, etotal=-3.5, vol=12.0, press=0.0, lam=1.0,
                 noise=0.0, rng=None):
    """Fill a FluctuationBuffer with n records, optionally with Gaussian noise."""
    if rng is None:
        rng = np.random.default_rng(42)
    buf = FluctuationBuffer()
    for i in range(n):
        buf.append(ThermoRecord(
            step=float(i),
            temp=temp   + noise * rng.standard_normal(),
            pe=pe       + noise * rng.standard_normal(),
            etotal=etotal,
            vol=vol     + noise * rng.standard_normal(),
            press=press + noise * rng.standard_normal(),
            lam=lam,
        ))
    return buf


def _make_step_jump_buffer(n_baseline=120, n_jump=120, jump_vol=2.0, noise=0.01):
    """Buffer where volume jumps after n_baseline records (simulates melting)."""
    rng = np.random.default_rng(0)
    buf = FluctuationBuffer()
    for i in range(n_baseline):
        buf.append(ThermoRecord(
            step=float(i), temp=1000.0,
            pe=-4.0 + noise * rng.standard_normal(),
            etotal=-3.5,
            vol=12.0 + noise * rng.standard_normal(),
            press=0.0,
        ))
    for i in range(n_jump):
        buf.append(ThermoRecord(
            step=float(n_baseline + i), temp=1000.0,
            pe=-3.0 + 0.5 * rng.standard_normal(),
            etotal=-2.5,
            vol=12.0 + jump_vol + 0.5 * rng.standard_normal(),
            press=0.0,
        ))
    return buf


# ---------------------------------------------------------------------------
# FluctuationBuffer tests
# ---------------------------------------------------------------------------

class TestFluctuationBuffer:
    def test_empty(self):
        buf = FluctuationBuffer()
        assert len(buf) == 0

    def test_append_and_len(self):
        buf = FluctuationBuffer()
        for i in range(5):
            buf.append(_make_record(i))
        assert len(buf) == 5

    def test_arrays_correct(self):
        buf = FluctuationBuffer()
        buf.append(ThermoRecord(step=1, temp=500, pe=-3.0, etotal=-2.5, vol=10.0, press=1.0, lam=0.8))
        buf.append(ThermoRecord(step=2, temp=600, pe=-4.0, etotal=-3.5, vol=11.0, press=2.0, lam=0.6))
        np.testing.assert_array_equal(buf.steps(),  [1, 2])
        np.testing.assert_array_equal(buf.temp(),   [500, 600])
        np.testing.assert_array_equal(buf.pe(),     [-3.0, -4.0])
        np.testing.assert_array_equal(buf.vol(),    [10.0, 11.0])
        np.testing.assert_array_equal(buf.lam(),    [0.8, 0.6])

    def test_enthalpy_uses_target_pressure(self):
        buf = FluctuationBuffer()
        target_p = 1000.0  # bar
        buf.append(ThermoRecord(step=0, temp=1000, pe=-4.0, etotal=-3.5, vol=12.0, press=0.0))
        H = buf.enthalpy(target_p)
        expected = -4.0 + target_p * BAR_ANG3_TO_EV * 12.0
        np.testing.assert_allclose(H, [expected], rtol=1e-10)

    def test_enthalpy_reduced_lam_one(self):
        """When lam=1, enthalpy_reduced should equal enthalpy."""
        buf = FluctuationBuffer()
        buf.append(ThermoRecord(step=0, temp=1000, pe=-4.0, etotal=-3.5, vol=12.0, press=0.0, lam=1.0))
        np.testing.assert_allclose(
            buf.enthalpy_reduced(0.0),
            buf.enthalpy(0.0),
            rtol=1e-12,
        )

    def test_enthalpy_reduced_halved_by_lam(self):
        """H/0.5 should be 2*H."""
        buf = FluctuationBuffer()
        buf.append(ThermoRecord(step=0, temp=1000, pe=-4.0, etotal=-3.5, vol=12.0, press=0.0, lam=0.5))
        np.testing.assert_allclose(
            buf.enthalpy_reduced(0.0),
            buf.enthalpy(0.0) / 0.5,
            rtol=1e-12,
        )

    def test_slice(self):
        buf = _fill_buffer(10)
        sub = buf.slice(2, 5)
        assert len(sub) == 3
        np.testing.assert_array_equal(sub.steps(), buf.steps()[2:5])

    def test_slice_no_end(self):
        buf = _fill_buffer(10)
        sub = buf.slice(3)
        assert len(sub) == 7

    def test_lam_default_one(self):
        """Records created without lam should default to 1.0."""
        buf = FluctuationBuffer()
        buf.append(ThermoRecord(step=0, temp=300, pe=-4.0, etotal=-3.5, vol=12.0, press=0.0))
        np.testing.assert_array_equal(buf.lam(), [1.0])


# ---------------------------------------------------------------------------
# _window_response_functions tests
# ---------------------------------------------------------------------------

class TestWindowResponseFunctions:
    def test_zero_variance_identical_records(self):
        Hl  = np.full(20, -4.0)
        vol = np.full(20, 12.0)
        rf = _window_response_functions(Hl, vol, 1000.0)
        assert rf["Cp"]      == pytest.approx(0.0, abs=1e-30)
        assert rf["kappa_T"] == pytest.approx(0.0, abs=1e-30)
        assert rf["alpha_P"] == pytest.approx(0.0, abs=1e-30)

    def test_Cp_formula(self):
        """Cp = Var(H/lambda) / (kB T^2)."""
        rng = np.random.default_rng(7)
        T = 1000.0
        Hl  = -4.0 + 0.1 * rng.standard_normal(200)
        vol = np.full(200, 12.0)
        rf = _window_response_functions(Hl, vol, T)
        expected_Cp = np.var(Hl, ddof=1) / (KB_EV * T**2)
        assert rf["Cp"] == pytest.approx(expected_Cp, rel=1e-6)

    def test_kappa_formula(self):
        """kappa_T = Var(V) / (kB T <V>)."""
        rng = np.random.default_rng(8)
        T = 1000.0
        Hl  = np.full(100, -4.0)
        vol = 12.0 + 0.1 * rng.standard_normal(100)
        rf = _window_response_functions(Hl, vol, T)
        kBT = KB_EV * T
        expected_kappa = np.var(vol, ddof=1) / (kBT * np.mean(vol))
        assert rf["kappa_T"] == pytest.approx(expected_kappa, rel=1e-6)

    def test_all_non_negative(self):
        rng = np.random.default_rng(42)
        Hl  = -4.0 + 0.05 * rng.standard_normal(50)
        vol = 12.0 + 0.05 * rng.standard_normal(50)
        rf = _window_response_functions(Hl, vol, 1000.0)
        assert rf["Cp"]      >= 0.0
        assert rf["kappa_T"] >= 0.0
        assert rf["alpha_P"] >= 0.0


# ---------------------------------------------------------------------------
# Rolling helpers
# ---------------------------------------------------------------------------

class TestRollingHelpers:
    def test_rolling_mean_constant(self):
        x = np.ones(20) * 3.0
        rm = _rolling_mean(x, 5)
        np.testing.assert_allclose(rm[4:], 3.0, rtol=1e-12)
        assert np.all(np.isnan(rm[:4]))

    def test_rolling_mean_values(self):
        x = np.arange(10, dtype=float)
        rm = _rolling_mean(x, 3)
        # rm[2] should be mean(0,1,2) = 1
        assert rm[2] == pytest.approx(1.0)
        assert rm[4] == pytest.approx(3.0)

    def test_rolling_var_constant(self):
        x = np.ones(20) * 5.0
        rv = _rolling_var(x, 4)
        np.testing.assert_allclose(rv[3:], 0.0, atol=1e-25)

    def test_rolling_cov_uncorrelated(self):
        rng = np.random.default_rng(1)
        x = rng.standard_normal(2000)
        y = rng.standard_normal(2000)
        cov = _rolling_cov(x, y, 200)
        # rolling covariance of uncorrelated series should be close to 0
        assert np.nanmean(np.abs(cov[200:])) < 0.15


# ---------------------------------------------------------------------------
# PhaseTransitionMonitor tests
# ---------------------------------------------------------------------------

class TestPhaseTransitionMonitor:
    def _make_monitor(self, expected_phase="solid", **kwargs):
        defaults = dict(
            target_pressure=0.0,
            temperature=1000.0,
            baseline_window=50,
            recent_window=50,
            min_samples_before_check=100,
            peak_threshold=5.0,
            min_signal_agreement=2,
        )
        defaults.update(kwargs)
        return PhaseTransitionMonitor(expected_phase=expected_phase, **defaults)

    def test_returns_none_when_buffer_too_small(self):
        monitor = self._make_monitor(min_samples_before_check=200)
        for i in range(50):
            monitor.update(_make_record(i), raise_on_transition=False)
        assert monitor.evaluate_final(raise_on_transition=False) is None

    def test_no_false_positive_stationary_noise(self):
        """Stationary low-noise data should not trigger a false positive."""
        monitor = self._make_monitor(min_signal_agreement=3)
        buf = _fill_buffer(300, noise=0.001)
        for i in range(len(buf)):
            monitor.buffer.append(
                ThermoRecord(
                    step=buf.steps()[i], temp=buf.temp()[i], pe=buf.pe()[i],
                    etotal=buf.etotal()[i], vol=buf.vol()[i], press=buf.press()[i],
                )
            )
        event = monitor.evaluate_final(raise_on_transition=False)
        assert event is None

    def test_detects_step_jump(self):
        """Large step-jump in volume + energy should trigger detection."""
        monitor = self._make_monitor(
            peak_threshold=2.0,
            min_signal_agreement=1,
        )
        buf = _make_step_jump_buffer(n_baseline=120, n_jump=120, jump_vol=3.0)
        for i in range(len(buf)):
            monitor.buffer.append(
                ThermoRecord(
                    step=buf.steps()[i], temp=buf.temp()[i], pe=buf.pe()[i],
                    etotal=buf.etotal()[i], vol=buf.vol()[i], press=buf.press()[i],
                )
            )
        event = monitor.evaluate_final(raise_on_transition=False)
        assert event is not None
        assert isinstance(event, TransitionEvent)
        assert len(event.triggered_signals) >= 1

    def test_raises_melted_error_for_solid(self):
        monitor = self._make_monitor(
            expected_phase="solid",
            peak_threshold=2.0,
            min_signal_agreement=1,
        )
        buf = _make_step_jump_buffer(n_baseline=120, n_jump=120, jump_vol=3.0)
        for i in range(len(buf)):
            monitor.buffer.append(
                ThermoRecord(
                    step=buf.steps()[i], temp=buf.temp()[i], pe=buf.pe()[i],
                    etotal=buf.etotal()[i], vol=buf.vol()[i], press=buf.press()[i],
                )
            )
        event = monitor.evaluate_final(raise_on_transition=False)
        if event is not None:
            with pytest.raises(MeltedError):
                monitor._raise_for_event(event)

    def test_raises_solidified_error_for_liquid(self):
        monitor = self._make_monitor(
            expected_phase="liquid",
            peak_threshold=2.0,
            min_signal_agreement=1,
        )
        buf = _make_step_jump_buffer(n_baseline=120, n_jump=120, jump_vol=3.0)
        for i in range(len(buf)):
            monitor.buffer.append(
                ThermoRecord(
                    step=buf.steps()[i], temp=buf.temp()[i], pe=buf.pe()[i],
                    etotal=buf.etotal()[i], vol=buf.vol()[i], press=buf.press()[i],
                )
            )
        event = monitor.evaluate_final(raise_on_transition=False)
        if event is not None:
            with pytest.raises(SolidifiedError):
                monitor._raise_for_event(event)

    def test_confidence_between_0_and_1(self):
        monitor = self._make_monitor(
            peak_threshold=2.0,
            min_signal_agreement=1,
        )
        buf = _make_step_jump_buffer(n_baseline=120, n_jump=120, jump_vol=3.0)
        for i in range(len(buf)):
            monitor.buffer.append(
                ThermoRecord(
                    step=buf.steps()[i], temp=buf.temp()[i], pe=buf.pe()[i],
                    etotal=buf.etotal()[i], vol=buf.vol()[i], press=buf.press()[i],
                )
            )
        event = monitor.evaluate_final(raise_on_transition=False)
        if event is not None:
            assert 0.0 < event.confidence <= 1.0

    def test_reset_clears_buffer(self):
        monitor = self._make_monitor()
        for i in range(10):
            monitor.update(_make_record(i), raise_on_transition=False)
        assert len(monitor.buffer) == 10
        monitor.reset()
        assert len(monitor.buffer) == 0

    def test_temperature_in_event(self):
        """TransitionEvent should include the estimated temperature."""
        monitor = self._make_monitor(
            peak_threshold=2.0,
            min_signal_agreement=1,
        )
        buf = _make_step_jump_buffer(n_baseline=120, n_jump=120, jump_vol=3.0)
        for i in range(len(buf)):
            monitor.buffer.append(
                ThermoRecord(
                    step=buf.steps()[i], temp=buf.temp()[i], pe=buf.pe()[i],
                    etotal=buf.etotal()[i], vol=buf.vol()[i], press=buf.press()[i],
                )
            )
        event = monitor.evaluate_final(raise_on_transition=False)
        if event is not None:
            assert event.temperature > 0.0


# ---------------------------------------------------------------------------
# detect_ts_transitions tests
# ---------------------------------------------------------------------------

class TestDetectTsTransitions:
    def _make_ts_data(self, n=5000, T_start=1200.0, T_stop=2500.0, natoms=500,
                       peak_at=0.5, peak_width=200, noise=0.01, rng_seed=42):
        """
        Synthetic ts data with an artificial variance spike at a given position.
        lambda goes from 1 -> T_start/T_stop (forward sweep convention).
        """
        rng = np.random.default_rng(rng_seed)
        lf = T_start / T_stop
        lam = np.linspace(1.0, lf, n)

        # Base pe/atom and vol/atom (per atom already)
        pe   = np.full(n, -4.0) + noise * rng.standard_normal(n)
        press = np.zeros(n)
        vol_atom = np.full(n, 12.0) + noise * rng.standard_normal(n)

        # Inject a large variance spike around peak_at fraction of the array
        idx_peak = int(peak_at * n)
        half_w   = peak_width // 2
        lo, hi   = max(0, idx_peak - half_w), min(n, idx_peak + half_w)
        pe[lo:hi]       += 0.5 * rng.standard_normal(hi - lo)
        vol_atom[lo:hi] += 0.5 * rng.standard_normal(hi - lo)

        vol_total = vol_atom * natoms
        return pe, press, vol_total, lam

    def test_detects_synthetic_peak(self):
        """Injected variance spike should be detected with default settings."""
        pe, press, vol_total, lam = self._make_ts_data(n=5000, peak_at=0.5, peak_width=300)
        events = detect_ts_transitions(
            dU=pe, press=press, vol_total=vol_total, lam=lam,
            t_start=1200.0, t_stop=2500.0, natoms=500,
            window_smooth=100,
            window_fluct=200,
            peak_threshold=3.0,
            min_signal_agreement=1,
            sweep_label="test_forward",
        )
        assert len(events) >= 1

    def test_no_detection_clean_data(self):
        """Clean (no injected peak) data should not trigger a detection."""
        rng = np.random.default_rng(99)
        n = 5000
        lf = 1200.0 / 2500.0
        lam = np.linspace(1.0, lf, n)
        pe  = -4.0 + 0.005 * rng.standard_normal(n)
        press = np.zeros(n)
        vol_total = (12.0 + 0.005 * rng.standard_normal(n)) * 500

        events = detect_ts_transitions(
            dU=pe, press=press, vol_total=vol_total, lam=lam,
            t_start=1200.0, t_stop=2500.0, natoms=500,
            window_smooth=100, window_fluct=200,
            peak_threshold=20.0,
            min_signal_agreement=2,
        )
        assert len(events) == 0

    def test_returns_empty_for_short_data(self):
        """Data shorter than window_fluct should return empty list."""
        n = 50
        lam = np.linspace(1.0, 0.48, n)
        pe  = np.zeros(n)
        press = np.zeros(n)
        vol_total = np.ones(n) * 6000.0

        events = detect_ts_transitions(
            dU=pe, press=press, vol_total=vol_total, lam=lam,
            t_start=1200.0, t_stop=2500.0, natoms=500,
            window_fluct=1000,
        )
        assert events == []

    def test_emits_userwarning(self):
        """detect_ts_transitions should emit a UserWarning when a transition is found."""
        pe, press, vol_total, lam = self._make_ts_data(n=5000, peak_at=0.5, peak_width=400)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            events = detect_ts_transitions(
                dU=pe, press=press, vol_total=vol_total, lam=lam,
                t_start=1200.0, t_stop=2500.0, natoms=500,
                window_smooth=100, window_fluct=200,
                peak_threshold=3.0,
                min_signal_agreement=1,
                sweep_label="test",
            )
            if events:
                assert any(issubclass(wi.category, UserWarning) for wi in w)

    def test_event_temperature_in_range(self):
        """Detected transition temperature should be within [t_start, t_stop]."""
        pe, press, vol_total, lam = self._make_ts_data(n=5000, peak_at=0.5, peak_width=400)
        events = detect_ts_transitions(
            dU=pe, press=press, vol_total=vol_total, lam=lam,
            t_start=1200.0, t_stop=2500.0, natoms=500,
            window_smooth=100, window_fluct=200,
            peak_threshold=3.0,
            min_signal_agreement=1,
        )
        for ev in events:
            assert 1100.0 <= ev.temperature <= 2600.0

    def test_detects_peak_with_descent(self):
        """
        A real first-order-style transition (sharp rise + descent back to a
        higher liquid plateau) must still be detected after the shape guards.
        """
        rng = np.random.default_rng(11)
        n = 6000
        natoms = 500
        T_start, T_stop = 1200.0, 2500.0
        lf = T_start / T_stop
        lam = np.linspace(1.0, lf, n)

        # Baseline solid fluctuations
        pe = -4.0 + 0.005 * rng.standard_normal(n)
        vol_atom = 12.0 + 0.005 * rng.standard_normal(n)

        # Localised spike at index ~ 0.5 (well inside the data) that decays
        # back to a slightly higher plateau (liquid).
        idx_peak = int(0.5 * n)
        spike_w = 250
        lo, hi = idx_peak - spike_w // 2, idx_peak + spike_w // 2
        pe[lo:hi]       += 0.4 * rng.standard_normal(hi - lo)
        vol_atom[lo:hi] += 0.4 * rng.standard_normal(hi - lo)

        press = np.zeros(n)
        vol_total = vol_atom * natoms
        events = detect_ts_transitions(
            dU=pe, press=press, vol_total=vol_total, lam=lam,
            t_start=T_start, t_stop=T_stop, natoms=natoms,
            window_smooth=100, window_fluct=200,
            peak_threshold=5.0,
            min_signal_agreement=1,
            sweep_label="real-transition",
        )
        assert len(events) >= 1

    def test_no_detection_when_peak_at_end(self):
        """
        Even if the modified-Z criterion is satisfied, a peak that sits in
        the trailing tail of the data (no observed descent) must be
        suppressed by the tail-margin guard.

        This test exercises the variance-based tail-margin guard
        specifically: the injected noise burst has zero mean shift but
        elevated variance, so a real transition is not present.  The
        first-moment slope-break detector is disabled (slope_break_sigma
        set high) because random-walk drift in the smoothed mean of an
        injected noise burst can produce a one-directional excursion that
        is genuinely indistinguishable from a small mean shift on the
        slope-break signal alone — by design, slope-break responds to
        sustained mean shifts wherever they occur.  The variance-based
        Cp/kappa_T/alpha_P peak detector with tail-margin remains the
        correct guard for "is there a real peak in the response function?"
        which this test asks.
        """
        rng = np.random.default_rng(3)
        n = 6000
        natoms = 500
        T_start, T_stop = 1200.0, 2500.0
        lf = T_start / T_stop
        lam = np.linspace(1.0, lf, n)

        pe = -4.0 + 0.003 * rng.standard_normal(n)
        vol_atom = 12.0 + 0.003 * rng.standard_normal(n)

        # Inject a spike right at the end of the array — no data after it.
        spike_w = 400
        pe[-spike_w:]       += 0.4 * rng.standard_normal(spike_w)
        vol_atom[-spike_w:] += 0.4 * rng.standard_normal(spike_w)

        press = np.zeros(n)
        vol_total = vol_atom * natoms
        events = detect_ts_transitions(
            dU=pe, press=press, vol_total=vol_total, lam=lam,
            t_start=T_start, t_stop=T_stop, natoms=natoms,
            window_smooth=100, window_fluct=200,
            peak_threshold=5.0,
            min_signal_agreement=1,
            slope_break_sigma=1e6,  # disable first-moment detector for this test
            sweep_label="tail-spike",
        )
        assert events == []


# ---------------------------------------------------------------------------
# Slope-break (first-moment) detector
# ---------------------------------------------------------------------------

class TestSlopeBreak:
    """Tests for the leverage-aware H/V slope-break detector."""

    def _build_sweep(self, n=8000, T_start=1000.0, T_stop=2000.0, natoms=500,
                    noise=0.005, rng_seed=11):
        rng = np.random.default_rng(rng_seed)
        lf = T_start / T_stop
        lam = np.linspace(1.0, lf, n)
        pe = -4.0 + noise * rng.standard_normal(n)
        press = np.zeros(n)
        vol_atom = 12.0 + 0.5 * noise * rng.standard_normal(n)
        vol_total = vol_atom * natoms
        return pe, press, vol_total, lam, rng

    def test_clean_baseline_no_detection(self):
        """Stationary noise gives small leverage-aware z -> no slope break."""
        pe, press, vol_total, lam, _ = self._build_sweep(noise=0.005)
        events = detect_ts_transitions(
            dU=pe, press=press, vol_total=vol_total, lam=lam,
            t_start=1000.0, t_stop=2000.0, natoms=500,
            window_smooth=100, window_fluct=200,
            peak_threshold=1e6,            # disable variance peaks
            min_signal_agreement=2,
            sweep_label="clean",
        )
        assert events == []

    def test_mean_shift_triggers(self):
        """A sustained mean shift in H mid-sweep should fire H_break."""
        n = 8000
        natoms = 500
        T_start, T_stop = 1000.0, 2000.0
        lf = T_start / T_stop
        lam = np.linspace(1.0, lf, n)
        rng = np.random.default_rng(7)
        pe = -4.0 + 0.005 * rng.standard_normal(n)
        # Step up pe by 0.05 eV/atom (mimicking a latent-heat jump) at
        # 60% through the sweep — well inside the data window so onset
        # is far from the tail margin.
        jump_at = int(0.6 * n)
        pe[jump_at:] += 0.05
        press = np.zeros(n)
        vol_total = (12.0 + 0.0025 * rng.standard_normal(n)) * natoms

        events = detect_ts_transitions(
            dU=pe, press=press, vol_total=vol_total, lam=lam,
            t_start=T_start, t_stop=T_stop, natoms=natoms,
            window_smooth=100, window_fluct=200,
            peak_threshold=1e6,            # disable variance peaks
            min_signal_agreement=1,        # only need H_break to fire
            sweep_label="h-step",
        )
        assert len(events) == 1
        assert "H_break" in events[0].triggered_signals

    def test_zero_mean_noise_burst_rejected(self):
        """
        Variance burst with zero mean shift (sign-flipping residuals)
        must NOT fire slope-break thanks to the sign-agreement guard.
        """
        n = 8000
        natoms = 500
        T_start, T_stop = 1000.0, 2000.0
        lf = T_start / T_stop
        lam = np.linspace(1.0, lf, n)
        rng = np.random.default_rng(3)
        pe = -4.0 + 0.003 * rng.standard_normal(n)
        vol_atom = 12.0 + 0.003 * rng.standard_normal(n)
        # Inject a *strictly* zero-mean noise burst in the middle (not the
        # tail).  We subtract the realised sample mean from each injection;
        # otherwise finite-sample drift of a std=0.3 noise series over
        # 1600 samples creates a real ~2.5σ_raw sustained shift inside
        # the burst, which the detector would (correctly) classify as a
        # phase change.  The test's intent is to verify that *purely
        # variance-only* bursts — no mean shift at all — are rejected.
        burst_lo, burst_hi = int(0.5 * n), int(0.7 * n)
        b = burst_hi - burst_lo
        pe_burst  = 0.3 * rng.standard_normal(b)
        pe_burst -= pe_burst.mean()
        vol_burst = 0.3 * rng.standard_normal(b)
        vol_burst -= vol_burst.mean()
        pe[burst_lo:burst_hi]       += pe_burst
        vol_atom[burst_lo:burst_hi] += vol_burst
        press = np.zeros(n)
        vol_total = vol_atom * natoms

        events = detect_ts_transitions(
            dU=pe, press=press, vol_total=vol_total, lam=lam,
            t_start=T_start, t_stop=T_stop, natoms=natoms,
            window_smooth=100, window_fluct=200,
            peak_threshold=1e6,            # disable variance peaks
            min_signal_agreement=1,
            sweep_label="zero-mean-burst",
        )
        # No detection from slope-break alone — sign disagreement filters
        # the burst out even though smoothed residuals briefly cross 5σ.
        assert events == []

    def test_compute_response_arrays_returns_residuals(self):
        """compute_ts_response_arrays must surface H_z, V_z for plotting."""
        pe, press, vol_total, lam, _ = self._build_sweep()
        arrs = compute_ts_response_arrays(
            dU=pe, press=press, vol_total=vol_total, lam=lam,
            t_start=1000.0, t_stop=2000.0, natoms=500,
            window_smooth=100, window_fluct=200,
        )
        for key in ("T", "Cp", "kappa_T", "alpha_P", "valid_start",
                    "H_z", "V_z"):
            assert key in arrs, f"missing key {key} in response arrays"
        assert arrs["H_z"].shape == arrs["T"].shape
        assert arrs["V_z"].shape == arrs["T"].shape

    def test_leverage_grows_with_extrapolation(self):
        """
        Sanity check: pred_std (encoded in z) grows away from baseline so
        the SAME absolute residual magnitude gives smaller |z| at the far
        extrapolation end than just outside the baseline window.

        We inject a flat baseline (y=0 + small noise) for the fit, then
        a constant non-zero residual ``bump`` that the linear fit cannot
        absorb (it would tilt the line so badly that baseline residuals
        explode).  At each point outside the baseline the raw residual
        is the same; leverage-aware normalization shrinks |z| as we move
        further out.
        """
        rng = np.random.default_rng(0)
        n = 5000
        T = np.linspace(1000.0, 2000.0, n)
        valid_start = 100
        # Flat baseline with thermal-scale noise.  Tiny noise lets the
        # fit converge tightly to y=0 over the baseline, so the constant
        # bump applied later is NOT absorbed by intercept drift.
        noise = 1e-3 * rng.standard_normal(n)
        y = noise.copy()
        baseline_frac = 0.25
        base_end_expected = valid_start + int(baseline_frac * (n - valid_start))
        # Apply a constant residual ONLY after the baseline window so the
        # fit (over the baseline) is unaffected.
        bump = 0.5
        y[base_end_expected:] += bump

        sb = _slope_break_signal(
            y, T, valid_start, baseline_frac=baseline_frac, fit_order=1,
        )
        assert sb is not None
        z = sb["z"]
        base_end = sb["base_end"]
        # Compare |z| at three points moving away from baseline.
        z_just_after = float(np.abs(z[base_end + 5]))
        z_mid        = float(np.abs(z[(base_end + n) // 2]))
        z_far_end    = float(np.abs(z[-10]))
        # All three are positive (constant residual = bump),
        # and they decrease monotonically with extrapolation distance.
        assert z_just_after > z_mid > z_far_end, (
            "leverage-aware |z| should shrink as extrapolation distance grows: "
            f"just-after={z_just_after:.2f}, mid={z_mid:.2f}, far={z_far_end:.2f}"
        )


# ---------------------------------------------------------------------------
# plan_temperature_blocks
# ---------------------------------------------------------------------------

class TestPlanTemperatureBlocks:
    """Tests for plan_temperature_blocks() block planner."""

    def test_heating_exact_multiple(self):
        """Heating with window that divides range exactly."""
        blocks = plan_temperature_blocks(1200.0, 2000.0, 100_000, 400.0)
        temps = [b["temp"] for b in blocks]
        assert temps == [1200.0, 1600.0, 2000.0]
        # First step must be 0, last must equal total_steps
        assert blocks[0]["step"] == 0
        assert blocks[-1]["step"] == 100_000

    def test_cooling_exact_multiple(self):
        """Cooling with window that divides range exactly."""
        blocks = plan_temperature_blocks(2000.0, 1200.0, 100_000, 400.0)
        temps = [b["temp"] for b in blocks]
        assert temps == [2000.0, 1600.0, 1200.0]
        assert blocks[0]["step"] == 0
        assert blocks[-1]["step"] == 100_000

    def test_endpoint_always_included(self):
        """Final checkpoint must be exactly tf regardless of window alignment."""
        blocks = plan_temperature_blocks(1200.0, 2000.0, 50_000, 300.0)
        assert blocks[-1]["temp"] == 2000.0
        assert blocks[-1]["step"] == 50_000

    def test_lambda_correct(self):
        """λ values must equal t0/T."""
        t0 = 1200.0
        blocks = plan_temperature_blocks(t0, 2400.0, 100_000, 600.0)
        for b in blocks:
            assert abs(b["lambda"] - t0 / b["temp"]) < 1e-10

    def test_steps_monotone_heating(self):
        """Step numbers must be non-decreasing for a heating sweep."""
        blocks = plan_temperature_blocks(1000.0, 3000.0, 200_000, 250.0)
        steps = [b["step"] for b in blocks]
        assert steps == sorted(steps)

    def test_steps_monotone_cooling(self):
        """Step numbers must be non-decreasing for a cooling sweep."""
        blocks = plan_temperature_blocks(3000.0, 1000.0, 200_000, 250.0)
        steps = [b["step"] for b in blocks]
        assert steps == sorted(steps)

    def test_single_block_large_window(self):
        """Window larger than total range → two checkpoints (start + end)."""
        blocks = plan_temperature_blocks(1200.0, 2000.0, 50_000, 5000.0)
        assert len(blocks) == 2
        assert blocks[0]["temp"] == 1200.0
        assert blocks[-1]["temp"] == 2000.0
        assert blocks[-1]["step"] == 50_000

    def test_invalid_window_raises(self):
        with pytest.raises(ValueError):
            plan_temperature_blocks(1200.0, 2000.0, 50_000, 0.0)

    def test_invalid_steps_raises(self):
        with pytest.raises(ValueError):
            plan_temperature_blocks(1200.0, 2000.0, 0, 100.0)

    def test_step_bounds_clamped(self):
        """No step value should fall outside [0, total_steps]."""
        blocks = plan_temperature_blocks(1500.0, 2500.0, 75_000, 123.7)
        for b in blocks:
            assert 0 <= b["step"] <= 75_000

    def test_linear_schedule_step_matches_lammps(self):
        """
        With the default 'linear' lambda schedule, the step number for each
        checkpoint must correspond to the LAMMPS step at which the linear-ramp
        lambda equals t0/T.  Verify with example 2 parameters: t0=500 K,
        tf=1500 K, n=25000, window=50 K.
        """
        t0, tf, n = 500.0, 1500.0, 25_000
        blocks = plan_temperature_blocks(t0, tf, n, 50.0)  # linear (default)
        lam_start = 1.0
        lam_end = t0 / tf
        for b in blocks:
            # Lambda that LAMMPS produces at this step under ramp(li, lf)
            lam_lammps = lam_start + (lam_end - lam_start) * b["step"] / n
            T_from_lammps = t0 / lam_lammps
            # Must agree with the checkpoint temperature to within 1 K
            assert abs(T_from_lammps - b["temp"]) < 1.0, (
                f"At step {b['step']}: planner T={b['temp']} K, "
                f"LAMMPS T={T_from_lammps:.1f} K (lambda={lam_lammps:.4f})"
            )

    def test_uniform_temperature_schedule_step_matches_linear_T(self):
        """
        With 'uniform_temperature' schedule, the step formula is
        s(T) = (T - t0) / (tf - t0) * total_steps (T linear in step).
        """
        t0, tf, n = 500.0, 1500.0, 25_000
        blocks = plan_temperature_blocks(t0, tf, n, 50.0,
                                         lambda_schedule="uniform_temperature")
        for b in blocks:
            expected_step = int((b["temp"] - t0) / (tf - t0) * n)
            assert b["step"] == min(expected_step, n), (
                f"T={b['temp']}: expected step {expected_step}, got {b['step']}"
            )
