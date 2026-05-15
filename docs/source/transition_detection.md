# Phase-Transition Detection in Reversible-Scaling Sweeps

## 1. Why we need it

A reversible-scaling (RS) calculation drags a system continuously from a
reference temperature $T_0$ to a target temperature $T_f$ along a single
molecular-dynamics trajectory.  The free energy at every intermediate
temperature is recovered by integrating the per-atom potential energy
$U(\lambda)$ along the scaling parameter $\lambda$.

If a **first-order phase transition** (melting, freezing, allotropic change)
occurs anywhere on the path $[T_0, T_f]$, the reversible-scaling identity is
broken: the system has crossed a coexistence line and the integral no longer
gives a meaningful free energy.  The simulation will still finish, but the
results are silently wrong.

The job of the transition detector is to notice this crossing and either warn,
stop, or automatically recover the calculation, so the user ends up with a
free-energy curve that stays inside one phase.

---

## 2. Two families of signals

Five signals are computed from the same ts data file.  They fall into two
families that target **different statistical moments** of the same physical
event:

| Family | Signals | Moment | Quantity that "lights up" | Lag vs. transition |
|--------|---------|--------|---------------------------|--------------------|
| **Variance-based** (second moment) | $C_p$, $\kappa_T$, $\alpha_P$ | Var / Cov of fluctuations | Fluctuations diverge in the two-phase region | Lags — peak appears $\sim w_f/2$ rows *after* the structural change |
| **Slope-break** (first moment) | $H_{\text{break}}$, $V_{\text{break}}$ | Mean trajectory $\langle y\rangle(T)$ | Latent-heat / volume jump perturbs the single-phase equation of state | Fires near onset — detects mean shift directly |

Both families are run on the same data and the detector requires at least
`min_agreement` (default 2) of the five to fire simultaneously, with the
recovery point set by the **earliest onset** across them.

---

## 3. What is recorded

At every MD step LAMMPS writes one row to a file `ts.forward_<i>.dat`
(or `ts.backward_<i>.dat`) with four columns:

| Column | Symbol     | Meaning                                |
|--------|------------|----------------------------------------|
| 1      | $u$        | Potential energy per atom (eV/atom)    |
| 2      | $P$        | Instantaneous pressure (bar)           |
| 3      | $V$        | Total simulation-box volume ($\mathrm{\AA}^3$) |
| 4      | $\lambda$  | Current value of the scaling parameter |

The equivalent thermodynamic temperature on the path is

$$
T(\lambda) \;=\; \frac{T_0}{\lambda}
\qquad (\text{reversible scaling, fixed thermostat}),
$$
or
$$
T(\lambda) \;=\; T_{\text{start}} + (T_{\text{stop}} - T_{\text{start}})\,(1-\lambda)
\qquad (\text{temperature scaling}).
$$

Per-atom enthalpy is reconstructed as

$$
H \;=\; u + P V / N,
$$

with $N$ the number of atoms.  Per-atom volume is $v = V/N$.

---

## 4. Smoothing pipeline

Raw per-step data are very noisy.  Before any detector runs, two
NaN-aware causal rolling-window operations are applied:

1. **Smoothing window** ($w_s$, default 500 rows): rolling mean of the raw
   columns, suppressing high-frequency thermal noise.

2. **Fluctuation window** ($w_f$, default 1000 rows): rolling variance and
   covariance of the smoothed columns, producing $C_p(t)$, $\kappa_T(t)$,
   $\alpha_P(t)$.

Rows up to index $w_s + w_f - 2$ have an incompletely populated window and
are excluded from detection.

The slope-break signals share the smoothing step but skip the second
window — they operate on $\langle y \rangle(T)$ directly.

---

## 5. Variance-based signals (second moment)

These detect the **divergence of fluctuations** at a first-order
transition.  Near coexistence, both enthalpy and volume can independently
sample either phase within one rolling window, so $\sigma_H^2$ and
$\sigma_V^2$ blow up.

### 5.1 Heat capacity at constant pressure

$$
C_p(T) \;=\; \frac{N\,\sigma_{H/\lambda}^2(T)}{k_B\,T^{2}}
$$

where $\sigma_{H/\lambda}^2$ is the rolling variance of the lambda-reduced
enthalpy per atom $H/\lambda$.  Using $H/\lambda$ rather than $H$ removes
the smooth Hamiltonian-scaling trend in reversible scaling and recovers
the standard equilibrium $C_p$ formula.

### 5.2 Isothermal compressibility

$$
\kappa_T(T) \;=\; \frac{\sigma_V^2(T)}{\langle V\rangle\, k_B\, T},
$$

with $\sigma_V^2$ the rolling variance of the box volume.

### 5.3 Isobaric thermal-expansion coefficient

$$
\alpha_P(T) \;=\; \left|\frac{\mathrm{Cov}(V,\,H/\lambda)(T)}{\langle V \rangle\, k_B\, T^2}\right|.
$$

### 5.4 Modified Z-score criterion

Inside a single-phase window the heat capacity grows naturally by a factor
of $\sim 8$ over a 200 K interval (equilibrium physics).  A naive
mean-plus-3$\sigma$ criterion would fire on this growth.  Instead we use
the **modified Z-score** built from the median and **median absolute
deviation (MAD)**, both of which are robust to the transition peak
itself:

$$
\mathrm{med} \;=\; \mathrm{median}\bigl(X_{\text{valid}}\bigr),
\quad
\mathrm{MAD} \;=\; \mathrm{median}\bigl(|X_{\text{valid}} - \mathrm{med}|\bigr),
$$
$$
\mathrm{modZ}(t) \;=\; \frac{X(t) - \mathrm{med}}{1.4826\;\mathrm{MAD}}.
$$

The baseline ($\mathrm{med}$, $\mathrm{MAD}$) is computed over **all valid
rows of the complete dataset**, not just the first 20 %.  With a global
robust baseline, equilibrium growth contributes $\mathrm{modZ} \lesssim 7$
while a true transition spike gives $\mathrm{modZ} > 20$.

A variance-based signal is flagged when

$$
\max_{t}\,\mathrm{modZ}(t) \;>\; \tau \quad (\text{default }\tau = 12).
$$

The threshold $\tau$ is the user-facing parameter `peak_threshold`.

### 5.5 Shape guards (suppressing false positives)

Two additional shape tests are applied after the modified-Z threshold:

| Guard | Condition rejected | Rationale |
|---|---|---|
| **Tail-margin** | Peak in the last 5 % of valid rows | Only the rising edge is visible; post-peak descent cannot be estimated reliably |
| **Descent-fraction** | Post-peak median has not dropped by ≥ 30 % of the excursion (peak − baseline) | Primary discriminator: monotonic solid $C_p$ growth (which never turns over) fails this; a real spike followed by a liquid plateau passes |

A head-margin guard is *not* applied: the rolling-window warmup already
excludes the first $w_s + w_f - 2 \approx 1500$ rows.

All three guards must be satisfied for a variance signal to count toward
`min_agreement`.

---

## 6. Slope-break signals (first moment)

The variance peaks are a **lagging** indicator: by the time the rolling
window holds enough mixed-phase samples for its variance to peak, the
structural change has already happened.  A first-order transition also
produces a much earlier, sharper signature in the **mean** trajectory:
the latent-heat jump perturbs $\langle H \rangle(T)$ away from the
single-phase equation of state, and the density jump perturbs
$\langle V \rangle(T)$ in the same way.  Two slope-break signals exploit
this directly.

### 6.1 Baseline equation-of-state fit

For each first-moment signal $y \in \{H, V\}$ a low-order polynomial is
fit through the first `slope_break_baseline_frac` (default 15 %) of the
valid samples — the early portion of the sweep that is guaranteed to be
single-phase:

$$
\hat y(T) \;=\; \sum_{k=0}^{p} c_k\,\left(\frac{T - T_c}{T_s}\right)^{k},
\qquad p = \mathrm{slope\_break\_fit\_order}\;(\text{default }1).
$$

$T_c$ and $T_s$ are an internal centering and scaling that improve the
condition number of the least-squares system.  A linear fit ($p=1$) is the
default because higher orders extrapolate poorly from a narrow baseline
window and inflate the false-positive rate.

### 6.2 Leverage-aware normalized residual

The fit residual at every $T$ is normalized by an effective noise scale:

$$
z(T) \;=\; \frac{y(T) - \hat y(T)}{\sigma_{\text{raw}}\,\sqrt{1 + h(T)}},
$$

where

* $\sigma_{\text{raw}}$ is the standard deviation of the residuals
  computed on the **raw** (unsmoothed) baseline data.  Calibrating on
  smoothed data would underestimate $\sigma$ by $\sqrt{w_s}$ and inflate
  $z$ across the board.
* $h(T) = x^\top\,(X^\top X)^{-1}\,x$ is the **leverage** of the design
  vector at $T$.  Inside the baseline window $h \in [0, 1]$; far outside
  it grows unbounded.  The $\sqrt{1+h}$ factor widens the noise envelope
  where the polynomial is being extrapolated, preventing a spurious blow-up
  of $|z|$ purely from extrapolation error.

The result is a $z$-statistic that is approximately $\mathcal N(0,1)$
inside a single phase regardless of how far we extrapolate, and which
spikes sharply when the system steps off the single-phase EOS.

### 6.3 Trigger condition

A slope-break signal fires when **all three** of the following hold:

1. **Threshold crossing**: $|z(T^\star)| > \mathrm{slope\_break\_sigma}$
   (default 5).

2. **Persistence**: the trigger sample $T^\star$ has at least
   $\max(2 w_s,\, w_f)$ finite samples after it.  This is longer than the
   autocorrelation length of the smoothed signal, so an autocorrelated
   noise excursion that briefly crosses the threshold cannot fake a
   sustained shift.

3. **Sign coherence**: the **signed** mean of $z$ over all samples from
   $T^\star$ onward satisfies
   $$
   \bigl|\overline{z}_{\,T \ge T^\star}\bigr| \;\geq\; 0.6\,\mathrm{slope\_break\_sigma}.
   $$
   This is the key discriminator against zero-mean noise: a real
   first-order transition produces a **one-directional** shift of
   $\langle y\rangle(T)$ (latent heat is intrinsically positive on
   heating, the volume jump on melting is positive, etc.), so the signed
   mean past the trigger stays $\sim$ trigger-sigma in magnitude.  A
   localized noise burst that briefly crosses the threshold but oscillates
   around zero afterward fails this test and is correctly rejected.

### 6.4 Why both $H$ and $V$

Some transitions barely move $\langle H \rangle$ but produce a large
$\Delta V$ (e.g. some allotropic transitions); others reverse this
ordering (e.g. liquid-liquid).  Running both signals independently
guarantees that whichever moment carries the transition signature will
fire.  Conversely, a glitch that affects only one of $H$ or $V$ will not
clear `min_agreement` on its own.

### 6.5 Why first-moment detection is phase-agnostic

The variance peaks are at heart melting/freezing detectors: they require
enough latent-heat-driven sampling of two phases inside one rolling window
to inflate the variance.  The slope-break signals require only that the
single-phase $H(T)$ and $V(T)$ curves are smooth (which they are for any
single-phase equation of state) and detect any sustained departure.  This
makes them effective for solid → solid transitions as well as
solid → liquid.

---

## 7. Onset temperature (recovery point)

The peak of a response function lags the structural change; the trigger
sample of a slope-break is somewhere along the rising edge.  Neither is
the temperature at which the transition *started* — and the recovery
checkpoint must lie *before* the start, not at the middle or the top.

For every flagged signal we therefore back-scan from the peak / trigger to
find the **onset**: the earliest sample on the rising edge that still
exceeds the baseline by `onset_sigma` "noise widths".  The walk-back
threshold is uniform across signal families:

| Signal family | Onset definition |
|---------------|------------------|
| Variance: $C_p$, $\kappa_T$, $\alpha_P$ | walk back from peak to the first sample where $X \le \mathrm{med} + \mathrm{onset\_sigma} \cdot 1.4826\,\mathrm{MAD}$ |
| Slope-break: $H_{\text{break}}$, $V_{\text{break}}$ | walk back from trigger to the first sample where $\lvert z \rvert \le \mathrm{onset\_sigma}$ |

The single parameter `onset_sigma` (default 4.0) controls all five.

Across signals that fire, the detector then uses the **minimum onset
temperature** as the recovery point — the most conservatively-early
candidate.  This matches the user-visible failure mode: "didn't stop in
time" is what hurts; "stopped a little too early" only costs you a few
extra K of useful range.

```
sweep direction →   T0 ─── onset ─── trigger ─── peak ─── (post-transition)
                            ↑
                    recovery checkpoint is the last block boundary
                    strictly before this temperature
```

The `TransitionEvent` object carries both temperatures:

| Attribute            | Meaning                                         |
|----------------------|-------------------------------------------------|
| `temperature`        | Temperature at the response-function peak (K)   |
| `onset_temperature`  | Temperature at the back-scanned onset (K)       |

Recovery uses `onset_temperature`; diagnostics and plots show both, plus
the block-boundary recovery cut (which is conservative by another half a
`temperature_window` on average).

---

## 8. Post-hoc detection on the full forward sweep

### 8.1 Why post-hoc, not incremental

Reliable detection of a first-order transition requires the full shape of
the response-function peak — **rise, maximum, and descent** back toward a
post-transition plateau.  Checking after every temperature block
(incremental detection) sees only a partial rise at each block boundary
and produces systematic false positives:

* The rolling-window head transient fires at the very start of any sweep.
* Monotonic solid $C_p$ growth near $T_m$ looks like a transition peak
  when only the rise has been observed (tail-margin false positive).
* The descent-fraction guard can only be applied once the curve has
  turned over.

Running the detector once on the **complete forward dataset** removes all
of these failure modes because the full peak shape, including the
descent, is always visible.  The slope-break signals share this benefit:
their persistence and sign-coherence tests both require a long trailing
window, which is only available after the full sweep is on disk.

### 8.2 Per-block checkpoints

While the forward sweep runs to completion, LAMMPS atomistic
configurations are saved at each temperature-block boundary:

```
conf.ts.forward_<iter>_blk0.data   (at T0)
conf.ts.forward_<iter>_blk1.data   (at T0 + ΔT)
...
conf.ts.forward_<iter>_blkN.data   (at Tf)
```

These checkpoints are the restart points for recovery.  No MD time is
wasted saving them — they are written as LAMMPS `write_data` calls between
contiguous `run` commands.

### 8.3 Sequence of operations

```
forward sweep to completion  →  post-hoc detector  →  backward sweep
      (writes per-block checkpoints)       (may modify T_stop)
```

1. `_reversible_scaling_forward` runs the full sweep and writes all
   checkpoints.
2. `_post_forward_recovery` reads the complete `ts.forward_<iter>.dat`,
   runs `detect_ts_transitions`, and acts according to `mode` (see §10).
3. `_reversible_scaling_backward` runs the backward sweep over the
   (possibly reduced) temperature range determined in step 2.

---

## 9. Recovery procedure (`mode: recover`)

When the detector fires in `recover` mode:

1. **Preserve the full forward data** for plotting only: copy
   `ts.forward_<iter>.dat` → `ts.forward_<iter>_full.dat`.

2. **Locate the safe block**: find the last checkpoint boundary $T_k$ that
   lies strictly *before* `onset_temperature` in the sweep direction.

3. **Truncate the forward ts file**: remove all rows beyond the step count
   corresponding to $T_k$, so the forward dataset covers only the
   single-phase region $[T_0, T_k]$.

4. **Promote the checkpoint**: copy
   `conf.ts.forward_<iter>_blk{k}.data` → `conf.ts.forward_<iter>.data`.
   The backward sweep will load this file for its middle equilibration.

5. **Reduce the sweep range**: set `_temperature_stop = T_k` and
   `_n_sweep_steps` to the step count at block $k$.  The backward sweep
   then covers the same lambda range as the truncated forward file, which
   is required by `integrate_rs` (it pairs `forward[i]` with
   `backward[N-1-i]`).

6. **Disable detection for the backward sweep**: set `td.mode = 'none'`
   so the now-safe backward sweep does not re-trigger and loop.

The result is a consistent forward / backward pair over $[T_0, T_k]$ from
which a valid free-energy curve is obtained.  The
`ts.forward_<iter>_full.dat` copy preserved in step 1 is used by the
response-function plot so the user can see the full temperature range
including the detected transition, while the calculation itself uses only
the safely truncated data.

---

## 10. Behaviour by mode

```yaml
phase_transition_detection:
  mode: recover             # none (default) | warn | recover | stop
  temperature_window: 50.0  # block width ΔT in K
  peak_threshold: 12.0      # modified-Z threshold τ (variance signals)
  min_agreement: 2          # signals that must trigger simultaneously
  onset_sigma: 4.0          # walk-back threshold for onset temperature
```

| `mode`    | What happens when a transition is detected |
|-----------|--------------------------------------------|
| `none`    | Detection disabled entirely. Sweep always completes. *(default)* |
| `warn`    | Warning logged and response-function plots generated; sweep continues unchanged. |
| `recover` | Forward ts file truncated, backward sweep re-run over safe sub-range (see §9). |
| `stop`    | `PhaseTransitionError` raised; simulation stops. Lower `temperature_stop` and re-run. |

---

## 11. Tunable parameters at a glance

| Parameter            | Default  | Affects | Meaning                                           |
|----------------------|---------:|---------|---------------------------------------------------|
| `mode`               | `'none'` | All     | Action on detection (see §10)                     |
| `temperature_window` | `50.0`   | Recovery| Block width $\Delta T$ in K                       |
| `peak_threshold`     | `12.0`   | Variance only | Modified-Z trigger $\tau$ for $C_p$, $\kappa_T$, $\alpha_P$ |
| `min_agreement`      | `2`      | All     | Number of signals (of 5) that must trigger        |
| `onset_sigma`        | `4.0`    | All     | Onset walk-back threshold (noise widths)          |

Two further hyperparameters control the slope-break family and are not
currently exposed in the input file (they have function-level defaults
that have been validated on Cu EAM and remain stable):

| Internal parameter           | Default | Meaning                                 |
|------------------------------|--------:|-----------------------------------------|
| `slope_break_sigma`          | `5.0`   | Trigger threshold on $|z|$ (sigmas)     |
| `slope_break_baseline_frac`  | `0.15`  | Fraction of early data used for EOS fit |
| `slope_break_fit_order`      | `1`     | Polynomial order for the baseline fit   |

The defaults have been validated against Cu EAM solid and liquid sweeps
spanning the melting transition (~1357 K).  In normal use no tuning should
be required.  The most useful knob in difficult cases is `onset_sigma` —
lower it to recover earlier (more conservatively), raise it to keep more
of the sweep at the cost of cutting closer to the transition.

---

## 12. Reading the diagnostic plot

In `warn` and `recover` modes a response-function plot
`ts_response_<sweep>_<iter>.png` is written into the run folder, with one
panel per signal family:

1. $C_p(T)$ with a twin axis showing $\mathrm{modZ}$ normalized by
   `peak_threshold` (so the trigger line sits at $y=1$).
2. $\kappa_T(T)$ with the same modZ twin.
3. $\alpha_P(T)$ with the same modZ twin.
4. $\lvert z_H(T) \rvert$ and $\lvert z_V(T) \rvert$ — the slope-break
   residuals — with a dashed line at `slope_break_sigma`.
5. (Forward sweep only, if available) $\langle q_{12}\rangle$ Steinhardt
   order parameter from per-block checkpoints — a structural sanity
   check, not wired into detection.

Vertical lines overlay each detected event:

* **Red dashed** at `onset_temperature` — the recovery-relevant
  temperature.
* **Dark-red dotted** at `temperature` (peak / trigger).
* **Orange solid** at the block-boundary recovery cut (if recovery was
  applied).  The forward calculation uses only data up to this line; the
  curves to its right come from the preserved `_full` copy and are shown
  for context only.
