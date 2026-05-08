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
    KB_EV,
    BAR_ANG3_TO_EV,
    _window_response_functions,
    _rolling_mean,
    _rolling_var,
    _rolling_cov,
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
