# p95 Tail Root-Cause Analysis (Measured 2026-09-25, Sprint 61 investigation)

> **Finding**: The urban-canyon $p_{95}$ tail is **not** false integer fixes and **not**
> outage divergence, which is what `TIER1_ROADMAP.md` §2 currently claims. It is a
> **kinematic motion-model defect**: the filter declares a velocity state it never
> observes, so its constant-velocity propagation contributes no motion and the
> estimate falls behind the vehicle whenever measurement updates are gated.

---

## 1. Measured baseline (Smooth PPK, all epochs)

`cargo run --release --bin eval_f9p_rover -- all`

| Dataset | Fix% | $p_{50}$ | $p_{95}$ | $h_{RMS}$ | Tier-1 target $p_{95}$ |
|---|---:|---:|---:|---:|---:|
| Tokyo Odaiba | 25.0% | 1.308 m | 4.787 m | 2.796 m | 0.15 m |
| Hong Kong TST1 (Patch) | 3.8% | 2.488 m | 14.013 m | 5.354 m | 3.00 m |
| Hong Kong TST1 (Survey) | 6.4% | 1.424 m | 13.739 m | 6.366 m | 1.80 m |
| Tokyo Shinjuku | 7.1% | 1.608 m | 22.377 m | 7.363 m | 1.50 m |
| Hong Kong Whampoa (Patch) | 10.4% | 2.264 m | 25.295 m | 10.097 m | 3.50 m |
| Hong Kong Whampoa (Survey) | 13.4% | 2.045 m | 20.140 m | 9.067 m | 2.00 m |

The $p_{95}/p_{50}$ ratio is **7–14×**. A tail that wide is structural, not noise.

---

## 2. Tail decomposition (Whampoa Patch, 1534 epochs)

| Population | $p_{50}$ | $p_{90}$ | $p_{95}$ | $p_{99}$ | max |
|---|---:|---:|---:|---:|---:|
| All epochs | 2.377 | 21.052 | 27.283 | 36.970 | 48.734 |
| Float epochs only (1370) | 2.666 | 22.475 | 28.793 | 37.262 | 48.734 |

Threshold census: **>5 m in 29.7% of epochs, >10 m in 20.1%, >20 m in 10.9%.**

The bad epochs are **not scattered outliers**. They form long contiguous runs:

| Run length (epochs) | Peak error | All float? |
|---:|---:|:---:|
| 98 | 48.7 m | yes |
| 79 | 35.2 m | yes |
| 58 | 37.7 m | yes |
| 33 | 37.1 m | yes |
| 25 | 25.1 m | yes |

40 such runs cover 455 of 1534 epochs. At 1 Hz that is a run of **up to 98 seconds**
of sustained multi-metre error — categorically different from the isolated
"wrong fix locked in for dozens of epochs" the roadmap describes.

---

## 3. Ruling out the two documented root causes

**False fixes — ruled out.** Every one of the worst runs reports `all_fixed = false`.
The AR gating machinery (Sprints 40–42, 54–59) is working as designed and is *refusing*
to fix. Fixed-subset accuracy is healthy: $p_{50}$ 1.08–1.36 m, $p_{95}$ 2.5–2.8 m.
The tail lives entirely in the float population.

**Outage divergence — ruled out.** Counting satellites per epoch directly in
`datasets/urbannav/hk_whampoa/rover_f9p.obs` (RINEX 3.03):

```
minimum satellites in any epoch ....... 9
longest run with < 10 satellites ...... 1 epoch   (i.e. none)
modal satellites per epoch ............ 15–20
```

There are no outages. The filter always has geometry. The roadmap's §2 item 4
("complete blockages caused unbounded drift") does not describe this dataset.

---

## 4. Actual root cause: the velocity state is never observed

`predict.rs` implements a constant-velocity kinematic model:

```rust
state.pos_ecef += state.vel_ecef * dt;   // predict.rs:34
f[(i, i + 3)] = dt;                      // predict.rs:41
```

But every write to `vel_ecef` in the module is one of:

- `vel_ecef: Vector3::zeros()` — initialisation (`state/mod.rs:64`)
- `state.vel_ecef * dt` — the propagation itself (`predict.rs:34`)
- `state.vel_ecef = Vector3::new(1.0, 2.0, 3.0)` — unit-test fixture (`predict.rs:164`)
- `velocity_ecef: Some(self.state.vel_ecef)` — output serialisation (`mod.rs:320`)

**Nothing ever writes a measured velocity.** The velocity block stays exactly zero for
the whole run, so `pos += vel*dt` is a no-op and the declared motion model is inert.
Position can advance *only* through measurement updates.

This predicts a specific, falsifiable signature — and it is what the data shows:

| Quantity | Value |
|---|---:|
| Truth speed over the dataset (median) | 2.523 m/epoch |
| Estimated displacement during the 98-epoch run | 69.8 m (0.713 m/s) |
| Estimated displacement during the 33-epoch run | 89.8 m (2.722 m/s) |

The estimate is **lagging, not frozen**: it moves, but at ~28% of true speed. Whenever
the NLOS/multipath gating correctly suppresses position updates, the filter cannot
coast on motion and the shortfall accumulates without bound — producing exactly the
observed 20–50 m excursions in multi-epoch runs.

---

## 5. Negative result: naive DD Doppler velocity destabilises the filter

The obvious remedy is to observe velocity from Doppler, which needs no ambiguity state
and whose Jacobian is the geometry vector the pseudorange row already uses:

$$\frac{d(r_{sat}-r_{ref})}{dt} = -e_{sat}\cdot v + e_{ref}\cdot v = d_{geom}\cdot v$$

This was implemented in full (DD Doppler observable plumbed through formation, a
velocity row appended to the IEKF system, Joseph covariance, 4 unit tests). It was
**measured and reverted** rather than shipped.

| Doppler sign | assumed $\sigma$ (m/s) | Whampoa Survey $p_{95}$ |
|:---:|---:|---:|
| −1 | 0.1 | 1500 m+ |
| +1 | 0.1 | 3933 m |
| −1 | 1.0 | 554 m |
| +1 | 1.0 | 1747 m |
| −1 | 2.0 | 192 m |
| *baseline (no Doppler)* | — | **20.14 m** |

Every configuration is 10–200× worse than doing nothing. Inflating $\sigma$ walks the
update back toward zero influence but never below baseline.

**Why it fails.** The modelled range rate omits the satellite-geometry terms, which do
not cancel in a double difference:

$$\big(e_{sat}^{base}-e_{sat}^{rov}\big)\cdot v_{sat} \;+\; \big(e_{ref}^{rov}-e_{ref}^{base}\big)\cdot v_{ref}$$

With a ~1 km baseline and ~20,000 km orbital radius the LOS directions differ by
$\sim 5\times10^{-5}$ rad, and $v_{sat}\approx 3.8$ km/s, so this residual is
**~0.2 m/s — the same order as the rover's own motion**. Gravitational/relativistic
Doppler and the ionospheric range rate add further unmodelled terms. Critically this
contamination is **correlated across satellites** (shared satellite-motion terms), so
it cannot be absorbed as independent white noise in $R$; inflating $\sigma$ only
discards the signal while the filter still couples the corrupted rows into position
through $P_{pos,vel}$.

A correct implementation must therefore add, at minimum: satellite-geometry
(relativity) Doppler correction, Sagnac correction in the rotating frame,
ionospheric range-rate handling, and a correlation-aware noise model. Do not re-attempt
the naive version.

---

## 6. What a correct fix requires

1. **Sagnac + relativity-corrected DD Doppler** as a first-class observation, with the
   satellite-geometry residual modelled rather than absorbed into $R$.
2. An initial velocity covariance in `RtkState` consistent with the declared CV model —
   currently the velocity block is presumably held at an arbitrary variance that no
   measurement ever informs.
3. Re-validate the fixed-subset accuracy as a regression gate: the fixed population
   (currently $p_{95}$ 2.5–2.8 m) must not degrade, since Doppler rows couple into
   position through the cross-covariance.

The lower-risk intermediate step — deriving velocity from successive *validated* float
position fixes and feeding it as a pseudo-measurement — introduces no new physics and
would make the existing CV model functional, but it cannot recover velocity while the
position track itself is lagging.

---

## 7. Benchmark hygiene note

`eval_f9p_rover.rs` sets `lever_arm: None` for both Whampoa datasets while the other four
carry explicit lever arms. Whampoa is therefore scored against the un-offset IMU/SPAN
reference. Its 2-pass calibration converges to `[0.9362, -0.0924, -0.6509]` against a
measured mean body-frame offset of `[0.9282, -0.2094, -0.6462]` — independent
confirmation of a real ~1.16 m mounting offset. This inflates the Whampoa $p_{50}$ floor
and should be resolved before Whampoa numbers are quoted against Tier-1 targets.

`eval_f9p_rover.rs` is at 476 LOC against the <500 limit, so tail diagnostics have no
room in place. They require extracting a shared `bench` harness module first.

---

## 8. Follow-up investigation (Sprint 61, rounds 2–3)

### 8.1 Which measurement survives as error grows

Instrumenting the IEKF to count accepted rows per epoch and bucketing by actual
horizontal error (Whampoa Survey, 1535 epochs):

| Horizontal error | Epochs | Code rows/ep | Phase rows/ep | Pairs offered |
|---|---:|---:|---:|---:|
| < 2 m | 791 | 22.77 | 15.78 | 26 |
| 2–5 m | 286 | 18.56 | 9.01 | 15 |
| 5–20 m | 388 | 16.04 | 6.76 | 15 |
| > 20 m | 70 | 17.20 | 5.54 | 16 |

This resolves the question left open in §4. Measurements are **not** being gated
away: 16–17 pseudorange rows are still accepted at the worst epochs. But carrier
collapses by 65% (15.78 → 5.54). The position is therefore being driven largely by
**raw pseudorange**, where NLOS bias is the dominant error term.

This reframes the problem: it is not a "filter cannot coast" problem, it is a
**code-bias** problem. Adding a motion model would not remove a biased measurement.

### 8.2 Negative result: sticky CMC protection is a no-op

§8.1 implies the CMC multipath detector should be the natural defence — but its
protection was gated behind carrier availability:

```rust
if let Some(cp) = dd_cp { /* downweight */ } else { base_pr_var }   // old
```

Losing carrier phase is itself a *symptom* of NLOS, so the pairs most likely to be
contaminated received no protection at all. This was implemented using the
tracker's remembered (stateful) estimate, decayed per epoch
(`CMC_STICKY_DECAY = 0.7`, self-clearing via `CmcTrack::update`).

**Measured result: no effect.** Across all 6 datasets every row was bit-identical
to baseline except TST1 Survey Fwd, which regressed slightly (p50 1.446 → 1.455 m).
Reverted.

Why inert: `apply_cmc_deweight` only sees `None` when the *band* carries no phase at
all, whereas the low phase-row count in §8.1 is dominated by pairs whose ambiguity
state is not yet allocated (`append_dd_phase_row` returns early when
`get_amb_idx` is `None`) — those pairs still take the `Some` branch and are
unaffected. The `multipath_m` estimate is also only non-zero after
`CMC_WARMUP_EPOCHS` of arc plus a deviation beyond `CMC_MULTIPATH_THRESHOLD_M`.

### 8.3 `robust_inflate` is a deliberate guard — do not weaken it

`robust_inflate` (`update/robust.rs`) inflates R in proportion to the *squared*
innovation, which caps any observation's normalised weight at a constant. A
lagging filter therefore recovers slowly: a 20 m innovation cannot snap the
solution back, because the guard refuses to let one observation dominate.

This is **intentional** and is the mechanism that guarantees zero false fixes, the
project's first-ranked requirement. Trading it for faster tail recovery would
re-open exactly the failure mode that produced the original 7.6–23.5 m tail.
Left untouched.

### 8.4 Concrete defect found: 2-pass lever-arm calibration is unstable

`eval_f9p_rover.rs` declares `lever_arm: None` for both Whampoa datasets and relies
on `execute_calibrated_post_process` to recover the mounting offset. Running that
calibration in **both** passes (instead of PPK only) exposed that the estimator does
not converge consistently:

| Dataset | Fwd-pass arm | PPK-pass arm | Fwd p50 | Fwd p95 |
|---|---|---|---:|---:|
| Whampoa Patch | `[0.9282, -0.2094, -0.6462]` | `[0.9362, -0.0924, -0.6509]` | 2.377 → 2.205 m | 27.283 → 27.723 m |
| Whampoa Survey | `[0.0004, -0.8163, -0.4414]` | `[-0.3110, -0.4597, 0.0840]` | 1.873 → **2.407 m** | 19.302 → **20.236 m** |

For Whampoa Patch the two passes agree closely and $p_{50}$ improves. For Whampoa
Survey they disagree materially (magnitudes 0.93 m vs 0.56 m, different directions)
and every Fwd metric gets worse — the estimator is fitting noise, not a physical
offset. Reverted; an unstable estimate must not be used to score results.

**Actionable follow-up:** make `CalibrationConvergenceCriteria` actually reject a
disagreement between passes (e.g. require the recovered arm to be stable across
independent segments before it is applied), and only then use the recovered arm for
both passes so the two passes measure the same thing. Until then the Whampoa Fwd
numbers carry an uncorrected mounting offset and must not be quoted against
Tier-1 targets.
