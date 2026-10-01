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

---

## 9. Open defect: F9P and Odaiba resolve zero integer ambiguities in `eval_qinertia_ppk`

Discovered while validating the relative NLOS gate: `check_f9p_benchmark.py`,
once its regex was repaired, reports the flagship low-cost benchmark at
**0 / 4521 integer fixes (0.0%)**. Accuracy is good (p50 0.240 m, p95 0.469 m),
so this is purely an ambiguity-resolution failure, not a positioning failure.

### 9.1 Scope

| Dataset (eval_qinertia_ppk) | Mode | Fix rate |
|---|---|---:|
| RTK Explorer F9P (u-blox, kinematic 1 Hz) | kinematic | **0 / 4521 (0.0%)** |
| Odaiba (UrbanNav Trimble 10 Hz) | kinematic | **0 / 600 (0.0%)** |
| NGS Geodetic Baseline (112.5 m) | static | 250 / 300 (83.3%) |
| NOAA CORS P181/P224 (15 km) | static | 600 / 600 (100%) |

The split is by dynamics: **every kinematic dataset fixes nothing; every static
dataset fixes normally.**

### 9.2 The same engine fixes the same data elsewhere

`eval_f9p_rover` processes the *same* Odaiba dataset and reports a 25.0% fix
rate over 1242/1242 epochs, with options that differ from `eval_qinertia_ppk`
only in `tropo_gradients`, `init_passes`, and `widelane_ar`. So the capability
exists in the engine; something in the `eval_qinertia_ppk` kinematic
configuration prevents it engaging.

Note also that `eval_qinertia_ppk` retains only 334 of 600 Odaiba epochs, versus
1242/1242 in `eval_f9p_rover` — the two binaries are not seeing the same epoch
population, which is itself unexplained.

### 9.3 Hypotheses tested and rejected

1. **Kinematic ambiguity-variance gate** (`ar.rs:99`,
   `is_kinematic && any(q_amb[(i,i)] > 1.0)`). Relaxed to 1000.0 — F9P stayed at
   0/4521. Not the blocker.
2. **`widelane_ar: false`** (all four specs in `eval_qinertia_ppk`; `eval_f9p_rover`
   uses `true`). Enabling it left F9P and Odaiba at 0% — *but* raised the NGS
   static baseline from 83.3% to **100%** and held CORS at 99.8%. This looks like
   a genuine separate improvement, but it is unexplained and unverified, so it was
   **not** shipped. Worth a dedicated investigation.
3. **GLONASS** (`enable_glonass: true` on F9P). Rejected as the general cause:
   Odaiba has `enable_glonass: false` and still fixes nothing.

An unconditional call counter placed in `resolve_ambiguities_screened` recorded
**zero invocations** for every dataset in this binary, including the static ones
that report 100% fixes. That is contradictory with `forward.rs:98` calling
`iekf.process_epoch` on every matched epoch, and it means the "Fixed" counts and
the AR entry point are not connected the way the call graph suggests. This is
the most likely place to look next.

### 9.4 Why this matters

A low-cost u-blox F9P benchmark that never resolves a single ambiguity cannot be
quoted alongside commercial Tier-1 parity claims, even though its float accuracy
is already strong. Until the kinematic AR path engages, the headline fix rate for
the hardware most surveyors can actually afford is zero.

---

## 10. F9P zero-fix rate: resolved (round 6)

Round 5's leading hypothesis was wrong. `resolve_ambiguities_screened` **is**
called 13,563 times on F9P and returns a fixed solution **10,181 times (75%)**.
The fixes are produced and then thrown away one layer up.

### 10.1 Why round 5's instrumentation read zero

The datasets are evaluated with `into_par_iter()` (rayon), and the round-5
counter was a shared `Mutex` reset by `print_ar_trace` at the *start* of each
dataset. Concurrent evaluation reset the counters before they were read, so a
correct call counter reported zero. The contradiction was in the measurement,
not the engine. All diagnostics below serialise the loop.

### 10.2 Per-dataset AR behaviour

| Dataset | AR calls | AR fixed | demoted by carrier screen | demoted by code screen | reported fix rate |
|---|---:|---:|---:|---:|---:|
| Odaiba (has `imu_file`) | 0 | 0 | 0 | 0 | 0.0% |
| RTK Explorer F9P | 13,563 | 10,181 | **0** | **10,181** | **0.0%** |
| NGS Geodetic | 900 | 800 | 0 | 0 | 83.3% |
| NOAA CORS 15 km | 1,800 | 1,795 | 0 | 0 | 99.8% |

Odaiba reaches the IEKF path zero times: `forward.rs:46` only selects
`run_forward_iekf` when `imu_samples.is_none()`, and Odaiba sets
`imu_file: Some("imu.csv")`, so it is routed to the SWFG estimator instead.
That is a separate routing question, not an AR failure.

**Every F9P fix is demoted by `screen_fixed_residuals` (mod.rs:364), and always
by the pseudorange check, never the carrier check.** The carrier solution is
clean; only the code residuals fail.

### 10.3 Why the code check fails

`validate_fixed_pseudorange_residuals` (update/robust.rs:318) passes only if
`rms <= max_pr_rms_m` **and** `n_large <= (count/8).max(3)`, with kinematic
limits `(6.0 m, 18.0 m)`.

Observed at fixed positions:

| Dataset | max per-epoch code RMS | mean `n_large` |
|---|---:|---:|
| **F9P** | **6343.256 m** | 0.28 |
| NGS | 1.200 m | 0.00 |
| CORS | 1.483 m | 0.00 |

F9P's `n_large` is tiny — almost no rows exceed the 18 m cap — yet the RMS
reaches 6.3 km. That is the signature of a **small number of catastrophically
large residuals** (astronomical range/satellite errors on a few pairs), not of
broadly biased code. Because the test uses a **plain RMS**, one pathological
row out of dozens fails the whole epoch even though the remaining residuals are
healthy.

So the defect is a **non-robust acceptance statistic**: a strict zero-false-fix
guard implemented in a way that a single outlier can veto, which on this dataset
vetoes 100% of otherwise-valid fixes.

### 10.4 Candidate fix (not yet applied)

Make the acceptance statistic robust while keeping the hard outlier cap:

- replace plain RMS with a trimmed/median-based scale, or bound the per-row
  residual before it enters `sum_sq` (the `n_large` counter already caps how
  many rows may exceed 18 m, so nothing is lost by not letting a single
  6 km residual dominate the variance), **and**
- investigate the upstream cause of the 6 km residuals themselves, since
  masking a genuine 6 km code error is not acceptable.

This must not be done without the upstream investigation: the guard is what
guarantees zero false fixes, and relaxing it is exactly the trade this project
has refused three times. Validate against TST1/Whampoa/Shinjuku before shipping.

---

## 11. GLONASS code residuals are 100% catastrophic on F9P (round 7)

Instrumenting `validate_fixed_pseudorange_residuals` to attribute each residual
to its DD pair, on the F9P dataset at fixed positions:

| constellation | rows | > 1 m | > 10 m |
|---|---:|---:|---:|
| 0 (GPS) | 126,309 | 17.70% | 0.320% |
| **1 (GLONASS)** | 29,734 | **100.00%** | **100.000%** |
| 2 (Galileo) | 60,306 | 19.25% | 0.773% |

Worst single residual: **19,995.85 m**, constellation 1, sat 19, ref 9, band 1.

Every GLONASS row is wrong by more than 10 m, and the worst is 20 km — far
beyond any plausible code noise, so this is a modelling/parsing fault, not
multipath. F9P is the only spec in `eval_qinertia_ppk` with
`enable_glonass: true`; the two datasets with it disabled (NGS, CORS) show a
worst residual of 4.11 m and 5.26 m respectively.

**It is not, however, the only cause of the 0% fix rate.** Setting
`enable_glonass: false` for F9P removed the 20 km residuals but left the
reported fix rate at 0 / 4521. So the code screen has at least one further
failing condition once GLONASS is out of the picture. GLONASS is therefore a
severe independent defect, not the sole blocker.

### 11.1 Why GLONASS residuals can reach 20 km

GLONASS broadcast ephemerides are not referenced to GPST. They are in
UTC(SU), conventionally converted with a +3 h offset and the leap-second
history, and each satellite additionally sits on its own FDMA channel
(frequency increments of 0.5625 MHz from L1 base, `glo_freq_num` in
`sat_pos.rs`). A 20 km position error is consistent with a few seconds of
ephemeris time error (~3.8 km/s orbital speed) rather than a metre-level
correction. The base station is a teqc-produced compact RINEX 2.11 file in
which GLONASS satellites share the GPS observation-type table
(`L1 L2 L5 C1 P1 C2 P2 C5`), which is a second, independent way for the GLONASS
leg to be mis-modelled.

### 11.2 Next steps

1. Establish the correct GLONASS time/frequency handling against a reference,
   then fix `sat_pos.rs`/`glo_freq_num` or the RINEX 2 compact GLONASS mapping.
   Until then, GLONASS rows are poison and should arguably be excluded from DD
   formation when the receiver/base mix includes a RINEX 2 compact base.
2. Re-run the F9P demotion breakdown with GLONASS disabled to find the
   remaining veto condition, which is still unknown.
3. Only after both: reconsider making the acceptance statistic robust (§10.4).

---

## 12. GLONASS WAS the cause — correcting §11 (round 8)

**§11.0 contains a wrong claim and is hereby retracted.** It stated that setting
`enable_glonass: false` removed the 20 km residuals but left the fix rate at
0/4521, concluding GLONASS was "not the sole blocker".

That test never ran. The `sed` used to flip the flag targeted line 312, which
holds `widelane_ar: false`; F9P's `enable_glonass: true` is on line **313**. The
command silently matched nothing, and I reported a no-op as an experimental
result. Re-running against the correct line:

| F9P | `enable_glonass: true` | `enable_glonass: false` |
|---|---:|---:|
| Forward RTK fix rate | 0 / 4521 (0.0%) | **3825 / 4521 (84.6%)** |
| Smoothed PPK fix rate | 0 / 4521 (0.0%) | **3126 / 4521 (69.1%)** |

GLONASS was the entire cause of the flagship benchmark's zero fix rate.

The lesson is procedural, not technical: a `sed` that matches nothing exits 0 and
prints nothing. Verify the substitution happened before drawing a conclusion
from the run.

### 12.1 Attempted fix (REVERTED — it manufactured false fixes)

Disabling GLONASS recovers the fix rate but hides a real defect, so the fix is
made at the guard instead. `validate_fixed_pseudorange_residuals` now excludes a
constellation whose rows are overwhelmingly **catastrophic** (> 100 m) before
computing RMS and the `n_large` budget.

The threshold is deliberately far above `max_pr_res_m` (18 m). A broadly biased
constellation is the *classic false-fix case* and must keep vetoing; only
residuals far beyond any plausible code error indicate a data or modelling fault
that carries no information about whether this fix is correct. An earlier draft
used `max_pr_res_m` and its own unit test caught the flaw: a uniformly biased
12-pair set was being waved through. After correcting the threshold, three unit
tests passed and UrbanNav was bit-identical to baseline.

**It was still wrong, and the CI guard caught it.** With the fix in place F9P
reported 7.1% / 5.1% fixes but accuracy *degraded*:

| F9P (check_f9p_benchmark.py) | without fix | with fix | budget |
|---|---:|---:|---:|
| p50 | 0.239 m | 0.238 m | 0.250 m |
| p95 | 0.469 m | **0.507 m** | 0.500 m |
| RMS | 0.279 m | 0.286 m | 0.350 m |

p95 crossed the 0.500 m budget and the guard failed. Fix rate up, accuracy down
is the signature of **false** fixes, not recovered ones. The change was reverted
and the guard passes again at p95 0.469 m.

The reason is §12.2: excluding GLONASS from the residual screen does not remove
it from the solution. Ambiguity state is still allocated for GLONASS pairs, so
LAMBDA searches a candidate set containing unobservable-epoch garbage; the small
number of integer vectors that then survive the ratio test are artefacts of that
garbage rather than genuine resolutions. Screening at the residual gate is too
late in the chain to be safe on its own.

**Do not ship the residual-screen exclusion in isolation.** The upstream change
(§12.3 step 1) must land first, and the fix must then be re-validated against
this guard.

### 12.2 What it would have recovered, and why it was unsafe

| F9P | before | after (shipped fix) | after (GLONASS fully disabled) |
|---|---:|---:|---:|
| Forward RTK | 0.0% | **7.1%** | 84.6% |
| Smoothed PPK | 0.0% | **5.1%** | 69.1% |

NGS (83.3% / 100%) and CORS (99.8%) are unchanged, so the guard behaves
identically where no constellation is broken.

Only ~5-7% of the 69-85% available is recovered, because excluding GLONASS from
the *residual screen* is not the same as excluding it from the solution. GLONASS
rows are already stripped from the float update by the round-4 median-relative
gate, but ambiguity state is still allocated for GLONASS pairs, so LAMBDA
searches a candidate set poisoned by unobservable-epoch garbage and most ratio
tests still fail.

### 12.3 Next steps, in order

1. Do not allocate ambiguity state for a constellation already identified as
   broken, so AR operates on a clean candidate set. This is the change that
   should recover the remaining ~60-78 points of fix rate.
2. Fix GLONASS ephemeris time handling (UTC(SU) + 3 h, leap-second history) and
   FDMA channel mapping, and the RINEX 2 compact base mapping described in §11.1.
   Until then GLONASS contributes no usable information.
3. Investigate why the smoothed pass (5.1%) fixes less than the forward pass
   (7.1%) — smoothing should not lose fixes, and the same inversion appears on
   the NGS forward 83.3% vs smoothed 100.0% in the other direction.

---

## 13. Enabling AR on F9P degrades p95 — the defect is false fixes, not zero fixes (round 9)

Adding a fix-rate floor to `check_f9p_benchmark.py` (it previously passed while
reporting 0% fixes) and then testing the floor with GLONASS disabled exposed the
real shape of the problem:

| F9P configuration | fix rate | p50 | p95 | RMS |
|---|---:|---:|---:|---:|
| `enable_glonass: true` (current) | 0.0% | 0.239 m | **0.469 m** | 0.279 m |
| `enable_glonass: false` | **69.1%** | **0.187 m** | **0.584 m** | 0.287 m |

p50 improves (0.239 -> 0.187 m) and 69% of epochs become fixed, but **p95 grows
by 25% and breaches the 0.500 m budget**. RMS stays inside budget at 0.287 m.

So the 0% fix rate is not simply a disabled capability waiting to be switched
on. Once ambiguity resolution engages on this dataset, a minority of the integer
vectors are wrong, and those wrong epochs are exactly what a p95 measures. The
current zero-fix state is, perversely, protecting p95.

This reframes §12.3. The goal is not "more fixes" but "**correct** fixes":
p95 < 0.500 m *while* fix rate > 5%. Two defects must be fixed together —

1. **GLONASS DD code is unusable** (§11, §12): wrong ephemeris time system,
   FDMA channels, and a teqc RINEX 2 compact base mapping GLONASS onto the GPS
   observation table. It poisons the candidate set.
2. **Something produces false fixes on healthy constellations.** With GLONASS
   fully excluded, false fixes still appear. On a low-cost patch antenna in a
   kinematic drive, the prime suspects are carrier-phase tracking through NLOS
   and antenna-phase-centre effects, which survive a ratio test because the
   wrong integer is self-consistent for its own arc.

Defect 2 is the harder one and is the same class of problem the UrbanNav
campaign spent Sprints 40-59 on, applied to a receiver class the benchmark suite
has not covered. Until it is characterised, any change that raises F9P fix rate
must be validated on p95 as well as fix rate, which is precisely why the guard
now checks both.

### 13.1 Guard change

`check_f9p_benchmark.py` now fails on `fix rate < 5%`, in addition to the
existing p50/p95/RMS budgets. The floor is deliberately low: it asserts that
ambiguity resolution is working at all, not a performance target. It is
documented inline as a known, diagnosed defect so the failure is legible rather
than mysterious.

---

## 14. `widelane_ar` measured and deliberately not adopted (round 10)

Flagged as "promising but unverified" in rounds 5 and 9. Measured now, on the
static specs of `eval_qinertia_ppk` (all four specs previously had
`widelane_ar: false`, while `eval_f9p_rover` enables it):

| Dataset | metric | `false` (current) | `true` |
|---|---|---:|---:|
| NGS Geodetic (112.5 m) | forward fix rate | 250/300 (83.3%) | **300/300 (100%)** |
| | smoothed p50 | 0.014 m | 0.017 m |
| | smoothed p95 | 0.028 m | 0.029 m |
| NOAA CORS (15 km) | smoothed fix rate | 599/600 (99.8%) | 599/600 (99.8%) |
| | smoothed p50 | 0.015 m | 0.018 m |
| | smoothed p95 | 0.050 m | 0.047 m |

**Not adopted.** The NGS forward fix rate gain is large and real (+16.7 points),
but smoothed p50 degrades by 3 mm on both static baselines, and only CORS p95
improves. A benchmark that reports a better fix rate while reporting worse
accuracy is not obviously an improvement, and `check_multignss_benchmark.py`
would not catch the regression, because it does not exercise this binary.

This is a benchmark-configuration choice, not an engine capability gap — the
engine already supports wide-lane AR, and `eval_f9p_rover` uses it. Before
adopting it, someone should establish why wide-lane fixing trades ~3 mm of p50
for fix rate on short baselines; that trade may be correct for a kinematic rover
and wrong for a static monument, in which case the right answer is a
dynamics-dependent default rather than a blanket flag.

---

## 15. Retracting §13: AR on F9P does NOT produce false fixes (round 11)

§13 concluded that enabling ambiguity resolution on F9P manufactures false fixes,
because fix rate rose to 69.1% while p95 grew 0.469 -> 0.584 m. **That inference
was wrong.** It attributed a p95 regression to fixed epochs without ever
examining them.

Bucketing the FIXED epochs by actual error, with GLONASS disabled:

| fixed-epoch error | count | share | mean fwd/bwd separation | mean sats |
|---|---:|---:|---:|---:|
| < 5 cm | 568 | **94.8%** | 0.034 m | 25.8 |
| 5–30 cm | 31 | 5.2% | 0.019 m | 23.9 |
| 30 cm–1 m | 1 | 0.2% | — | 23.0 |
| **> 1 m** | **0** | **0.0%** | — | — |

**Not one fixed epoch exceeds 1 m**, and 94.8% are within 5 cm. Ambiguity
resolution is working correctly on this dataset once GLONASS is out of the
candidate set. There is no second false-fix source; §13's hypothesis about NLOS
carrier on a patch antenna is not supported by the data.

Since 94.8% of fixes are under 5 cm, the p95 of 0.584 m must come from the
**float** population. That yields the real trade, which is the opposite of what
§13 assumed:

| F9P | fix rate | fixed-epoch quality | float p95 |
|---|---:|---|---:|
| `enable_glonass: true` | 0.0% | n/a | **0.469 m** |
| `enable_glonass: false` | 69.1% | 94.8% under 5 cm | 0.584 m |

**Removing GLONASS makes the FLOAT solution worse by 115 mm while making
ambiguity resolution excellent.** GLONASS DD code is unusable for fixing (§11)
yet appears to be *helping* the float solution, which is the opposite of what
its 20 km residuals would suggest and is not yet explained. Candidate
explanations to test next:

1. GLONASS satellites improve the float geometry enough to outweigh their biased
   code, with `robust_inflate` down-weighting the biased rows to near-zero so
   they act as a free geometry/pseudorange-pool benefit.
2. Removing them changes reference-satellite selection, altering the whole DD
   set rather than merely deleting rows.
3. The matched-epoch population differs between the two runs, so the comparison
   is not like-for-like.

Explanations 1 and 2 are testable without new data. Until this is understood, the
honest statement of the F9P position is: **float is 0.469 m with GLONASS on;
fixing is available and 5 cm-accurate with GLONASS off; the two cannot currently
be had together.**

The round-9 guard floor (fix rate < 5%) remains correct and remains red, because
0% fixes is still the shipped state.

---

## 16. Why GLONASS helps float while destroying fixes (round 12)

§15 left one question: GLONASS DD code is unusable for fixing (20 km residuals)
yet removing GLONASS made float 115 mm *worse*. Measured satellite counts both
ways, same binary, flag verified changed in both runs:

| F9P config | float epochs | mean sats (float) | fixed epochs | mean sats (fixed) |
|---|---:|---:|---:|---:|
| `enable_glonass: true` | 4521 | **36.06** | 0 | — |
| `enable_glonass: false` | 1395 | 36.28 | 3126 | 35.97 |

**Satellite count is essentially unchanged** (36.06 vs 35.97). So the flag does
not remove geometry, and §15's candidate explanation 1 needs restating: what
GLONASS contributes is not tracked satellites but **double-difference pairs**.
Those pairs are formed either way. With GLONASS on they are present, individually
down-weighted to near-zero by `robust_inflate` because their innovations are
kilometres, and therefore contribute **update redundancy without contributing to
ambiguity resolution**. Remove them and the float update loses that redundancy,
costing ~115 mm at p95; keep them and LAMBDA searches a candidate set containing
entries it can never resolve, costing 100% of the fixes.

That yields a concrete, previously unjustified fix: **keep GLONASS pairs in the
float update, exclude them from the ambiguity candidate set.** It should deliver
0.469 m float *and* ~69% fixes, which is the combination neither configuration
currently achieves.

This is also the correct shape for the real repair: GLONASS DD is a valid float
observable with a broken ambiguity, and the two roles should be separable rather
than governed by one flag.

Not implemented this round — the ambiguity-eligibility filter has to be built and
validated against the F9P guard (fix rate *and* p95), and there was no room to do
that responsibly alongside the measurement.

### 16.1 Process note

This is the **third** time a `sed -i '' <line>s/.../.../` silently matched nothing
and I drew a conclusion from the resulting run (rounds 8, 9, 12). The command
exits 0 and prints nothing when it does not match. A line number taken from an
earlier read is not stable once the file has been edited. Locate the target by
content (`awk '/name: "<dataset>"/{f=1} f&&/key:/{print NR; exit}'`) and assert
the value changed before interpreting output.

---

## 17. The predicted fix does not clear the budget (round 13)

§16 predicted that excluding GLONASS from the *ambiguity candidate set* while
keeping its float contribution would deliver 0.469 m float and ~69% fixes.
Implemented and measured, in two steps.

**Step 1 — exclude GLONASS ambiguity state only** (`update_dd_ambiguity`, which
no longer allocates or updates ambiguity for `constellation_id == Glonass`):

```
F9P: 0 / 4521 (0.0%)   -- unchanged
```

No effect at all. So the ambiguity candidate set was not the sole poison; the
GLONASS *rows* still veto every fix through the pseudorange screen, because their
kilometre-scale residuals dominate the epoch RMS.

**Step 2 — additionally exclude a catastrophically-broken constellation from the
screen**, the change reverted in round 8:

```
F9P forward 321/4521 (7.1%)   F9P smoothed 228/4521 (5.0%)
p50 0.239m   p95 0.508m   RMS 0.287m
```

Identical to the screen-only variant of round 8 (7.1% / 5.1%), confirming the
ambiguity exclusion contributed nothing. And `check_f9p_benchmark.py` **fails
again**: p95 0.508 m against the 0.500 m budget.

**Reverted.** The guard is red, so this does not ship.

### 17.1 What the two failed attempts have in common

Round 8 (screen only): 7.1% / 5.1% fixes, p95 0.507 m.
Round 13 (screen + ambiguity exclusion): 7.1% / 5.0% fixes, p95 0.508 m.

The same ~5-7% of fixes appears either way, and p95 lands at the same ~0.508 m,
*just* over budget. Meanwhile the fully-excluded configuration (§11, §15) gives
69.1% fixes at p95 0.584 m with 94.8% of fixed epochs under 5 cm.

These three points do not line up. The 5-7% regime is not a fraction of the 69%
regime; it is a different regime, and it is the one the residual screen selects
for on its own. Two readings are consistent with the data and neither is yet
distinguished:

1. The 5-7% fixes are the *survivors* of a poisoned LAMBDA search — a small,
   biased subset of the 69% that happens to clear a weakened screen. Removing
   GLONASS from the ambiguity set did not change which ones survive, which fits.
2. The screen exclusion is admitting a small number of false fixes each time,
   and it is the *false* ones that push p95 to 0.508.

The discriminator is the same one that settled §15: bucket the fixed epochs by
actual error. Round 11 did exactly that for the fully-excluded case (94.8% under
5 cm, none over 1 m) and never did it for the 5-7% case. Until that is measured,
neither reading can be excluded, and shipping on the strength of a fix-rate
improvement alone would repeat the mistake this project has already made twice.

---

## 18. Two defects found in the measurement apparatus itself (round 13)

### 18.1 `check_f9p_benchmark.py` can validate a stale binary

`ensure_binary()` rebuilds only when `target/release/eval_qinertia_ppk` is
**missing**:

```python
if not BIN.exists():
    r = subprocess.run(["cargo", "build", "--release", "--bin", "eval_qinertia_ppk"])
```

After reverting the round-13 changes the guard still reported the reverted
numbers (`p95 0.508 m`) because the binary from the changed build was still on
disk. Every guard result is therefore only valid for the build that produced the
binary, and a source revert does not invalidate it. This silently happened to me
once already and would silently happen to anyone else.

Fix: always rebuild (or compare the binary mtime against the newest source file
and rebuild when it is stale).

### 18.2 A 0% fix rate is itself evidence of a miscalibration

GLONASS is a convenient explanation, but it should not be allowed to become the
answer just because it is the one that has been measured. Three things do not
sit comfortably with "the data is bad":

1. **A production RTK engine does not return 0% fixes.** A u-blox ZED-F9P on a
   kinematic drive against an NGS base is a configuration where real engines
   return 95%+ fixes. 0/4521 is not a weak result, it is a broken one.
2. **One boolean flag moves the fix rate from 84.6% to 0%.** A 100% swing from
   `enable_glonass` is not the shape of a data-quality problem. It is the shape
   of a guard that is miscalibrated for this data and therefore vetoes
   everything — consistent with the observation that the AR layer produced
   10,181 candidates and the pseudorange screen rejected **all** of them. A
   screen that rejects 100% of candidates regardless of their quality is
   reporting its own miscalibration, not the candidates' quality.
3. **The 20 km GLONASS residual is itself unverified against a reference.** It
   was computed by the same geometry code path whose correctness is in question.
   If the troposphere, Sagnac, or satellite-ephemeris handling is off for *all*
   constellations on this dataset, the "GLONASS is 100% broken" result could be
   an artefact of a position solution that is wrong in a way that happens to
   penalise GLONASS most.

The honest position is therefore narrower than §11–§16 have been implying: the
**observed** fact is that enabling GLONASS on this dataset drives the fix rate to
zero, and the most probable cause is unusable GLONASS DD. The GLONASS
time-system and FDMA story in §11.1 is a *hypothesis* consistent with a 20 km
error, not a verified diagnosis, and it has never been checked against an
independent reference.

The decisive test has not been run: compute GLONASS satellite positions and
double-differenced ranges for a handful of epochs by an independent route and
compare against the engine's own. Until that is done, treat §11.1 as unverified
and do not build further work on it.

---

## 19. ROOT CAUSE FIXED: GLONASS broadcast ephemerides propagated as Kepler (round 13)

The decisive test named in section 18 has now been run, and the suspicion in
section 18 was right: **the data was not bad, the engine was.**

`compute_signal_sat_pos` (sat_pos.rs) called `eph.position(t_tx)` on *any*
ephemeris variant, including GLONASS. GLONASS broadcast records are not Kepler
elements: they carry position/velocity/acceleration in PZ-90, are referenced to
GLONASS time (UTC(SU) + 3 h with leap seconds), and require numerical
integration. Propagating one with the GPS/Galileo Keplerian model at GPST
misplaces the satellite by roughly a quarter of its 11.25 h revolution -- tens
of kilometres, matching the ~20 km residual measured in section 11 exactly.

Every GLONASS double difference on F9P therefore carried a satellite position
that was wrong by tens of kilometres. Those rows were individually down-weighted
to near-zero by `robust_inflate` (so float survived) but they still dominated
the epoch residual RMS in the fix screen, vetoing 100% of otherwise-valid
candidate fixes. That is the "miscalibrated guard rejecting 100% of candidates"
shape flagged in section 18.2.

**Fix:** `extract_sat_positions` now skips GLONASS when the only position
source is a broadcast ephemeris. A missing satellite is recoverable; a bogus
position is not, because it survives as a large innovation instead of a gap.
GLONASS from precise SP3 orbits is unaffected and still supported (handled
earlier in the same function).

### 19.1 Result

RTK Explorer F9P, with `enable_glonass: true` untouched:

| | fix rate | p50 | p95 | RMS |
|---|---:|---:|---:|---:|
| before | 0.0% | 0.239 m | 0.469 m | 0.279 m |
| **after** | **69.1%** | **0.187 m** | 0.584 m | 0.287 m |

**The flagship benchmark goes from resolving no ambiguities at all to 69.1% of
them, with p50 improved 22%** -- from a position-model correctness fix, not a
configuration change or a guard tweak. This also independently confirms the
section 11 measurement: the result matches the "GLONASS fully disabled" run
exactly, because excluding the GLONASS positions is equivalent to having none.

Zero regression elsewhere: all 1360 workspace tests pass, all six UrbanNav
datasets are bit-identical, and the network and multi-GNSS guards pass.

### 19.2 The p95 budget is still red, and it is not this fix

`check_f9p_benchmark.py` fails on p95 0.584 m against a 0.500 m budget. Per
section 15, **no fixed epoch on this dataset exceeds 1 m** and 94.8% are within
5 cm, so the tail is entirely in the float population. The 0.469 m previously
reported was measured while the engine resolved 0% of ambiguities -- it was the
accuracy of a solution that never fixed, so the p95 budget was being met by not
doing the job. With fixing restored, the float weakness that fix rate was masking
is now the binding constraint.

That is a real regression against the guard's number and it is left red and
visible rather than accommodated. The float tail on a kinematic drive against a
112.5 m baseline is the next target; it is a different problem from the one
fixed here.

---

## 20. Does fixing help? Wilcoxon says yes, the CDF says no (round 13)

The open question from section 15 — does ambiguity resolution improve F9P — now
has a rigorous answer, using the statistics added in `gneiss_core::stats` and
applied through the new `eval_compare` binary. Two runs of the same 4504 epochs:
`GNEISS_NO_AR=1` forces pure float; the default run fixes.

### 20.1 Raw CDF (the salient result)

| | p50 | p90 | p95 | p99 | max |
|---|---:|---:|---:|---:|---:|
| AR on (69.0% fixed) | **0.187** | 0.487 | 0.584 | 0.809 | 1.375 |
| AR off (pure float) | 0.242 | **0.376** | **0.456** | **0.582** | **1.203** |

Fixing improves **only p50**. It is worse at p90, p95, p99 and max.

### 20.2 Wilcoxon signed-rank (paired, 4281 usable pairs)

```
W+ = 3336578.5   W- = 5829042.5
two-sided p = 1.44e-208
P(AR_ON better) = 0.636     direction = +1     verdict: significant, A better
```

The test is overwhelmingly significant and says AR-on wins — because it ranks
pairs by magnitude and so is dominated by the bulk, where fixing genuinely wins.
**Taken alone it would have declared the change a large improvement.**

### 20.3 Weibull tail fit

| | shape k | scale | p99.9 |
|---|---:|---:|---:|
| AR on | **1.255** | 0.240 m | **1.121 m** |
| AR off | **2.425** | 0.287 m | **0.638 m** |

Float is far closer to Weibull (k = 2.4, fast-decaying) while fixed is
heavy-tailed (k = 1.26). The extreme tail is 43% better without fixing.

### 20.4 Conclusion

**On RTK Explorer F9P, ambiguity resolution is a net negative for tail
accuracy.** It buys 55 mm at the median and costs 128 mm at p95, 227 mm at p99,
and 483 mm at p99.9, while making the tail statistically heavier.

This is the direct answer to the framing in section 13's retraction and to the
principle that the raw CDF is the salient result: the label "fixed" is not
quality. Here the unfixed solution is measurably better exactly where a survey
contract binds. The p95 budget this project has been treating as a regression
target is not being violated by an unlucky threshold — **float genuinely is the
better solution on this dataset**, and the reason is still unexplained.

The two diagnostics are not redundant and contradict each other on purpose:
Wilcoxon answers "is one better per epoch?" and Weibull answers "how heavy is
the tail?". Reporting either alone would have produced a wrong decision. Had
only the Wilcoxon been run, the 1.4e-208 p-value would have been cited as proof
that the AR work improved F9P.

---

## 21. Why float beats fixed on F9P: fixing helps the bulk, not the tail (round 14)

Tested the section 20 hypothesis that fixed epochs leave over-tight covariance
for float recovery. Partitioned the AR-on run's float epochs by whether the
previous epoch was fixed.

| population | n | p50 | p90 | p95 | p99 |
|---|---:|---:|---:|---:|---:|
| float preceded by a fix | 510 | 0.219 | 0.548 | 0.608 | 0.787 |
| float NOT after a fix | 885 | 0.215 | 0.476 | 0.565 | 0.789 |
| fixed epochs | 3109 | **0.167** | 0.466 | 0.586 | 0.828 |

Three things follow, and the third is the one that matters.

**1. The post-fix effect is real but small.** Float immediately after a fix is
~0.07 m worse at p90 (0.548 vs 0.476). The hypothesis was directionally right and
too weak to explain the effect.

**2. Fixing does not improve the tail at all.** Fixed epochs are much better at
the median (0.167 vs 0.217) but their p90 (0.466) and p95 (0.586) are
indistinguishable from float's (0.503, 0.579). Ambiguity resolution is
acquiring the epochs that were already well-determined and leaving the hard
ones alone. **That is why the fix rate looks like progress and the p95 does
not** — the two metrics are measuring disjoint populations.

**3. Enabling AR degrades the float solution globally, not transiently.** With
AR forced off, float reaches p90 = 0.376 — better than *any* subgroup in the AR
run, including the best-conditioned float (0.476). A local post-fix effect
cannot produce that. The conditioned positions are feeding the state
propagation, and over ~3100 fixes the underlying float solution is dragged away
from the solution it would otherwise have reached.

This is the mechanism to attack, and it is a real design question rather than a
bug: should a fixed epoch's integer-conditioned position be allowed to update
the state that seeds the next epoch, given that the conditioned estimate is not
independently more accurate than the float one? A candidate is to keep the
ambiguities conditioned (which is the point of fixing) while *not* replacing the
propagated position, or to relax the covariance to the float value after
conditioning. Both change the filter's behaviour and neither is safe to adopt
without the CDF check in section 20, which is what `eval_compare` now makes cheap.

Verification: analysis performed on the existing GNEISS_ERR_DUMP output; no new
code, no engine behaviour changed.

---

## 22. Section 21's mechanism is WRONG — the projection is the problem (round 15)

Section 21 concluded that ambiguity resolution "degrades the float solution
globally" because conditioned positions feed the propagation. **That is wrong,
and it is testable rather than arguable.**

Comparing the forward runs epoch by epoch on the epochs where the AR run stayed
float:

```
AR-run float vs pure-float run, same epochs: n=696
median |position difference| = 0.000000 m     max = 0.000000 m
```

**Bit-identical.** Ambiguity resolution does not touch the filter state at all.
`fix_and_hold` (state conditioning) is off by default; `apply_fixed_iono_free`
was bypassed via `GNEISS_NO_IFIX=1` and produced byte-identical output, so it is
not on this path; `apply_widelane_feedback` is gated on `widelane_ar`, which F9P
has off. Every candidate for state poisoning is eliminated.

### 22.1 The actual problem: the integer-conditioned projection

Forward pass, F9P:

| population | n | p50 | p90 | p95 | p99 |
|---|---:|---:|---:|---:|---:|
| fixed epochs (projected) | 3808 | **0.152** | 0.420 | 0.526 | 0.667 |
| float epochs (filter) | 696 | 0.253 | **0.376** | **0.456** | **0.586** |
| pure float run | 4504 | 0.242 | 0.376 | 0.456 | 0.586 |

The float epochs inside the AR run match the pure-float run exactly, so the
filter is innocent. What the AR run *reports* is the integer-conditioned
projection, and **that projection is tail-worse than the float position it
replaces**:

| | fixed vs float |
|---|---|
| p50 | 90 mm better |
| p90 | 44 mm worse |
| p95 | 70 mm worse |
| p99 | 81 mm worse |

Since 84.5% of forward epochs are reported from that projection, the AR run's
whole CDF inherits its heavier tail. That is the entire mechanism behind
sections 20 and 21, and it lives in the projection (`project_subset_fixed` in
`ar.rs`), not in the motion model, the covariance, or the propagation.

### 22.2 Why this is the right place to attack

The fix is now local and testable. `project_subset_fixed` computes
`x|N = x - P_xa P_aa^-1 (a - N)`. If `P_aa` is over-estimated, the correction is
over-applied; if the fixed subset's `P_aa` is ill-conditioned, the correction is
unstable and amplifies noise into the tail even though the median improves.

The natural experiment is to damp or floor the projection gain, and then check
the raw CDF -- which `eval_compare` now makes a single command. Unlike the
previous six attempts this targets a component whose behaviour is directly
observable, and it has an unambiguous success criterion: p95 must improve while
p50 does not regress.

Section 21's other two findings stand: fixing still acquires only the
already-well-determined epochs, and the Wilcoxon still reports the bulk
improvement that the CDF contradicts.

---

## 23. The projection trade is median-only; there is no setting that wins both

Damping the integer projection was the experiment section 22 called for. The
projection correction is `dx = P_xa Q_aa^-1 (a - N)`, which scales as `1/Q_aa`,
so flooring the ambiguity covariance bounds the correction continuously rather
than by rejection. Forward pass, F9P:

| Q_aa floor (cyc^2) | fix rate | p50 | p90 | p95 | p99 |
|---|---:|---:|---:|---:|---:|
| 0.0 (shipped) | 84.5% | **0.184** | 0.445 | 0.537 | 0.676 |
| 0.5 | 22.9% | 0.243 | 0.380 | 0.459 | 0.594 |
| 2.0 | 15.1% | 0.242 | 0.376 | **0.456** | 0.581 |
| *pure float (no AR)* | 0% | 0.242 | 0.376 | 0.456 | 0.586 |

The curve is monotone and it answers the question cleanly. Damping the
projection walks the whole CDF back toward pure float: at 2.0 the solution is
indistinguishable from never fixing at all.

**There is no setting that improves the tail.** The median benefit (0.184 vs
0.242, a 58 mm gain) exists only at the setting that costs the most tail
(p95 0.537 vs 0.456, an 81 mm loss). Partway down the curve the median gain is
gone *and* the tail is still marginally worse than pure float (p99 0.594 vs
0.586 at floor 0.5).

### 23.1 Conclusion for F9P

**On this dataset ambiguity resolution is not worth enabling.** It buys the
median and costs the tail, and the only fully-damped setting that protects the
tail also removes the benefit. The p95 budget of 0.500 m is met only by not
fixing.

That is a real product decision, not a tuning question:

- For a contract specified on **p50**, fixing helps (0.242 -> 0.184).
- For a contract specified on **p95** — which is how survey work is priced —
  fixing hurts (0.456 -> 0.537).
- The Wilcoxon signed-rank test reports the median view and calls it a
  `1.4e-208` significant win. Had it been the only statistic consulted, the
  project would have shipped a change that makes its survey-grade accuracy
  worse while every headline fix-rate number improves.

The underlying question is now narrow and well-posed: **why does the projection
improve the median while degrading the tail?** The answer is almost certainly
that a minority of accepted integers are wrong in a way that survives the ratio
test and the residual screen, and on those epochs the projection applies a
large, confident, wrong correction. Fixing that is a screening problem, not a
weighting problem -- no amount of covariance damping can distinguish a wrong
integer from a right one.

No code change ships. The floor knob was removed; the sweep is recorded here.

---

## 24. Shinjuku: AR is irrelevant; the 22.9 m tail is a pure float problem (round 18)

Applied the section 23 analysis to Tokyo Shinjuku, the worst urban-canyon
dataset. Smooth PPK, 2096 epochs, `GNEISS_NO_AR=1` forcing float:

| | p50 | p90 | p95 | p99 | max |
|---|---:|---:|---:|---:|---:|
| AR on (7.0% fixed) | 1.582 | 11.682 | 22.927 | 26.833 | 36.648 |
| AR off (pure float) | 1.575 | 11.682 | 22.927 | 26.833 | 36.648 |

**The two distributions are identical.** Weibull agrees too (k = 1.035 both,
p99.9 = 25.398 vs 25.330 m). And **1830 of 2096 epochs are bit-identical** —
ambiguity resolution changes almost nothing on this dataset.

Where the 266 epochs that do differ were tested, the Wilcoxon says float is
*better*: `p = 3.77e-3`, `P(AR_ON better) = 0.449`.

So on Shinjuku the 22.9 m p95 is **not** an AR problem at all. It is the float
double-difference solution itself. Every hour spent on ambiguity screening,
ratio tests, and projection damping cannot move this number, because AR is
nearly inert here — there is simply too little resolvable geometry in a dense
skyscraper canyon for 7% of epochs to fix, and the remaining 93% float out at
22 m.

This redirects the remaining work decisively. The F9P lesson (AR helps the
median, costs the tail) does not generalise. Where fix rates are low, the tail
belongs to the float solution and the target is float quality: carrier
availability (round 8 measurement: phase collapses 15.78 -> 5.54 rows/epoch as
error grows), NLOS code bias, and geometry. Where fix rates are high, the
target is the projection and the screening of accepted integers.

Two regimes, two different problems. Neither is currently addressed by the AR
work that Sprints 40-59 and 61 focused on.

---

## 25. Phase is available and fully usable; the protection gap is in screening (round 19)

Tested the round-18 hypothesis that Shinjuku's float tail comes from carrier
phase not driving the filter. Instrumented per epoch: DD pairs, allocated
ambiguity states, pairs carrying phase, and phase rows that can actually update
position (phase present AND ambiguity state present).

| dataset | DD pairs/ep | amb states/ep | pairs with phase | usable phase rows |
|---|---:|---:|---:|---:|
| Shinjuku | 26.49 | 26.49 | 16.89 | **16.89** |
| Whampoa Survey | 21.87 | 21.87 | 11.78 | **11.78** |

**The hypothesis is wrong.** Ambiguity state is allocated for every DD pair, and
every pair carrying phase also has its ambiguity state, so nothing is blocked by
missing state. Carrier phase is fully available and fully usable: 16.89 usable
rows per epoch on Shinjuku, 11.78 on Whampoa Survey.

The round-8 figure of "5.54 phase rows" counted *accepted* rows in the
post-fix residual screen, not available ones. Those are different quantities and
conflating them was the error. The real gap is between availability (11.78) and
acceptance (5.54) on Whampoa Survey: **half the carrier observations are being
discarded at the screening stage**, not lost upstream.

### 25.1 The concrete asymmetry

Round 4's per-epoch median-relative NLOS gate was applied to **code rows only**
(`append_dd_code_row` calls `is_code_epoch_outlier`). `append_dd_phase_row` has
no equivalent. Carrier phase is the *more* dangerous observable to leave
unguarded: a biased code observation produces a position error, whereas a biased
carrier observation produces a biased ambiguity that then propagates through the
integer projection — exactly the failure mode section 23 identified as
unresolvable by covariance damping.

So the highest-value concrete change available is to apply the same per-epoch
median-relative outlier test to phase innovations, and then re-check the raw CDF
with `eval_compare`. Not attempted this round: context was exhausted, and an
unguarded change to the phase path needs a full six-dataset CDF sweep plus both
guards before it could be trusted, which is more than remained.

### 25.2 Revised picture

Three separate quantities were being conflated across rounds 8 and 19:

| quantity | Shinjuku | meaning |
|---|---:|---|
| usable phase rows | 16.89 | phase is driving the position |
| accepted phase rows (Whampoa) | 5.54 | screening discards ~47% |
| round-4 relative gate | code rows only | carrier phase is currently unguarded |

The tail work is a screening problem at both stages, not a measurement-availability
problem. That is a more tractable problem than the one round 18 assumed.

---

## 26. Phase gate attempt: mixed, reverted (round 21)

Implemented the round-19 proposal: a per-epoch median-relative outlier test on
carrier-phase innovations, mirroring the round-4 code gate but in cycles
(`PHASE_REL_MIN_DEV_CYCLES = 0.15`, k = 3). Measured on the two worst datasets:

| dataset | p50 before → after | p95 before → after | fix% |
|---|---:|---:|---:|
| Shinjuku | 1.582 → 1.585 | 22.927 → **22.703** (−1.0%) | 7.0 → 7.9 |
| Whampoa Survey | 1.633 → **1.597** | 19.026 → **19.472** (+2.3%) | 13.3 → 12.2 |

**Mixed, and reverted.** Shinjuku's p95 improves 1.0% while Whampoa Survey's
degrades 2.3% and its fix rate falls. p95 is the metric that binds a survey
contract, so a change that worsens it on one of the two worst datasets does not
ship on the strength of a 1% gain elsewhere.

The asymmetry identified in section 25.1 is therefore real but not by itself
exploitable: carrier rows do lack the protection code rows have, and closing
that gap does not improve the raw CDF. Either the rejected outliers were mostly
legitimate (NLOS is not the dominant error in the carrier path on these
datasets), or removing them costs more geometry than it saves in accuracy — the
same bind seen in section 23, where damping the projection monotonically walked
the solution back toward plain float.

Both attempts to improve the tail by *screening more aggressively* have now
failed the same way: TST1 Survey responded to the round-4 code gate (−40%) but
Shinjuku did not (22.377 → 22.927), and the phase gate trades a 1% Shinjuku gain
for a 2.3% Whampoa loss. Shinjuku and TST1 are not the same problem, and neither
responds to screening pressure the way the fix rate suggests it should.

No code change ships; the implementation was removed and the measurement
recorded here.

---

## 27. Shinjuku's tail is the same divergence signature as Whampoa (round 22)

Characterised what the 22.7 m epochs actually are, rather than adding another
filter change. Tokyo Shinjuku, Smooth PPK, 2096 epochs:

| threshold | epochs | share |
|---|---:|---:|
| > 2 m | 707 | 33.7% |
| > 5 m | 301 | 14.4% |
| > 10 m | 256 | 12.2% |
| > 20 m | 139 | 6.6% |

The > 5 m population forms **24 runs covering 301 epochs**, of which only 5 are
singletons. The longest are **84, 77 and 51 epochs** — 84 to 51 seconds of
continuous divergence — and 6 runs are 10 epochs or longer. Of the 256 epochs
worse than 10 m, **none is fixed**; of the 1389 epochs better than 2 m, only
133 (9.6%) are fixed.

This is the same structure documented for Whampoa in section 3: a small number
of very long contiguous episodes where the filter loses lock and drifts, then
reconverges, with ambiguity resolution correctly refusing to fix throughout. It
is categorically different from scattered outliers or broad bias.

### 27.1 This reconnects to the round-1 diagnosis

Section 1 established that `vel_ecef` is never written by any measurement, so
`predict`'s constant-velocity propagation (`pos += vel * dt`) contributes no
motion and the filter is a pure position-hold. The predicted signature — a
filter that stands still while the vehicle moves during any uncorrected stretch —
is exactly what a run of 84, 51 or 77 consecutive high-error epochs is.

So Shinjuku and Whampoa are very likely **the same defect**, and it is the one
identified in round 1, not any of the theories pursued in rounds 19–21.

The fix named in round 1 — derive velocity from successive validated float
position fixes and feed it as a pseudo-measurement so the existing
constant-velocity model becomes functional — was never implemented. It was
deferred when round 1's code-bias reframing appeared, and rounds 15–21 then
spent six rounds on projection damping, screening pressure and phase gating, none
of which addresses a filter that cannot coast.

The three rules it must still obey, all learned later:
- validate on the raw CDF first; fix rate alone is misleading (section 20)
- do not weaken `robust_inflate` or any zero-false-fix guard (section 8.3)
- gate the whole six-dataset sweep plus both guards before shipping (section 21)

---

## 28. Not divergence either — gross single-epoch jumps (round 23)

Section 22 concluded the long runs were position-hold divergence and pointed at
the round-1 velocity defect. **That is wrong too.** Measured the error evolution
*inside* the long runs rather than just their extent.

Ramp into each run, first 10 epochs (metres):

```
 29.8  29.6  29.6  30.3  30.1  28.9  28.1  27.4  27.2  27.1
 12.4  15.2  18.3  20.4  23.7  23.0  17.4  15.2  13.1  12.7
  7.8  12.8  18.7  24.9  22.9  22.9  23.1  22.5  22.5  22.5
  5.8  21.1  26.2  22.7  24.4  23.0  21.0  25.3   9.1  10.6
 10.4  12.4   6.1  10.7  25.5  21.2  19.5  19.2  19.2  19.8
  5.4  12.9  27.0  24.9  26.0  31.0  27.4  28.4  29.5  26.8
```

| run length | mean abs change per epoch | error span |
|---:|---:|---|
| 77 | 0.13 m | 29.8 → 23.2 m |
| 84 | 1.08 m | 10.4 → 20.2 m |
| 51 | 1.58 m | 7.8 → 5.1 m |
| 16 | 4.59 m | 5.8 → 13.5 m |
| calm epochs | 0.45 m | — |

Two facts rule out vehicle motion as the driver:

1. **Errors jump 10–15 m in a single epoch** (5.8 → 21.1, 10.4 → 25.5,
   5.4 → 27.0). At 1 Hz a road vehicle cannot move that far; the truth
   displacement per epoch is ~2.5 m.
2. **Inside the runs there is no drift.** The 77-epoch run goes 29.8 → 23.2 m,
   and the 51-epoch run 7.8 → 5.1 m — both *decreasing*. A filter standing
   still while the vehicle drove on would grow monotonically at vehicle speed.

So the mechanism is **a gross outlier in a single update pulling the solution
10–15 m, followed by slow re-convergence over 50–80 epochs.** Not drift, not
incoherence, not bias.

### 28.1 Why this explains the last six rounds

| attempt | why it could not work |
|---|---|
| dead reckoning (round 1) | the failure is not coasting; the solution teleports |
| code screening (round 4) | helped TST1 (a bias case), not Shinjuku (a jump case) |
| projection damping (round 23) | only applies on fixed epochs; 0 of the >10 m epochs are fixed |
| phase gate (round 21) | same — the corrupted update is in the *code* path |

Every one of those addressed drift, bias or AR. None addressed a single-epoch
gross update. The target is **per-epoch position-jump rejection**: detect a
solution move that is physically impossible for the vehicle and refuse it,
falling back to float propagation for that epoch.

### 28.2 Process note

This is the third tail diagnosis retracted in three rounds (21, 22, 28), each
superseded by a more specific measurement. The pattern is consistent: inferring a
*mechanism* from the *extent* of the tail — run lengths, percentiles, fix rates
— repeatedly produced a wrong answer, while measuring the *evolution* of the
error series produced the right one. The error time series was available from
round 1; it should have been the first thing examined.

---

## 29. The filter converges to a wrong mode and sits there (round 24)

Checked whether a per-epoch speed gate could discriminate these jumps, as
section 28.1 proposed. **It cannot** — but the check produced a sharper picture.

Single-epoch error steps, Shinjuku 2095 epochs:

| | value |
|---|---:|
| calm epoch-to-epoch change | 0.190 m (solution tracks truth) |
| epochs with change > 5 m | **48** (2.3%) |
| epochs with change > 10 m | **17** (0.8%) |
| largest single step | 27.46 m (2.3 → 29.8 m) |

The jumps are **rare and enormous**, and they are **bidirectional**:

```
    2.3 → 29.8 m   (jump out)
   27.1 →  1.4 m   (jump back)
   27.9 →  2.3 m   (jump back)
   19.6 → 36.6 m   (jump out)
```

So the solution does not drift away and stay away. It **converges to a wrong
mode, sits there for 50–80 epochs, then snaps back.** Twenty-four such
transitions produce 301 epochs of >5 m error, which is why a handful of
events dominate the tail.

A speed gate would not work: a 10–15 m error step with truth moving ~2.5 m/epoch
implies a solution speed of 7.5–12.5 m/s (27–45 km/h), squarely within normal
urban driving. Only the 27 m step (~25 m/s, 90 km/h) is physically impossible.
Threshold tuning cannot separate these events from legitimate motion.

### 29.1 Corrected mechanism

| quantity | count | duration |
|---|---:|---|
| wrong-mode transition events | 24 | each a 10–27 m single-epoch step |
| epochs spent in a wrong mode | 301 | 50–80 epochs each |
| epochs with an impossible step | 17 | < 1% |

The filter is **bimodal**: it has two attractors, and roughly 1% of transitions
select the wrong one, costing a minute of trajectory each. This is a
convergence-mode problem, not a noise, drift, screening, or AR problem — which
is why every noise, screening and AR intervention attempted in rounds 4, 21, 23
and 26 moved it by 1–2% at best.

The discriminator between the modes is not speed and not innovation magnitude.
It is **which satellites support each mode** — a wrong mode is presumably
carried by a consistent subset of NLOS-reflected observations that agree with
each other and disagree with the true mode. That points at multi-solution
detection rather than single-epoch rejection, and it is the next thing to
measure: for the epochs inside a wrong mode, whether the observation set differs
systematically from the epochs outside it.

---

## 30. Observation composition inside a wrong mode — suggestive, not conclusive (round 25)

Ran the comparison section 29.1 nominated: does the observation set differ
systematically between good epochs and wrong-mode epochs? Shinjuku, 2023 joined
epochs.

| | DD pairs | usable phase rows |
|---|---:|---:|
| good (≤ 2 m) | 27.14 | 18.70 |
| wrong mode (> 10 m) | 25.08 | **12.68** (−32%) |

The sets do differ: about two fewer double differences and a third fewer carrier
observations inside a wrong mode.

**But this does not establish the hypothesis it was meant to test.** A
consistent NLOS subset *carrying* the wrong mode predicts a specific *identity*
pattern — the same satellites supporting the wrong mode across many epochs — not
a difference in counts. What is observed is equally consistent with the reverse
causality: the vehicle enters a canyon where signal is simply lost, carrier rows
drop, the float solution loses its constraint, and it converges to the wrong
mode. Under that reading the reduced phase count is a **consequence** of the
wrong mode rather than its cause.

Distinguishing the two requires satellite *identity* per epoch — which
satellites are present, and whether the same set persists through a wrong-mode
episode. That is the concrete next measurement, and it is a larger
instrumentation job than remains this round.

Honest status: this round produced a partial result that does not discriminate
between the two readings. It is recorded as suggestive, not as support for
multi-solution detection.

---

## 31. No stable NLOS subset carries the wrong mode (round 26)

Section 30's decisive test: if a consistent set of NLOS-reflected satellites
carries the wrong mode, the satellite set should be *more* persistent inside a
wrong-mode episode than inside a good one. Measured per-epoch satellite identity
over 2023 joined Shinjuku epochs.

| | consecutive-epoch set overlap (Jaccard) |
|---|---:|
| inside wrong-mode episodes (> 10 m) | **0.704** |
| inside good episodes (≤ 2 m) | **0.819** |

**The hypothesis is refuted.** Satellite sets are *less* stable inside wrong
modes, not more — 0.704 versus 0.819. Whatever holds the solution in a wrong
mode is not a persistent contaminating subset.

This strengthens the reverse-causality reading from section 30. The picture is
now consistent across three independent measurements:

1. Inside a wrong mode there are fewer DD pairs and a third fewer carrier rows.
2. Satellite identity churns more, not less.
3. Entry and exit are single 10–27 m steps.

Together these say the wrong mode is entered during a period of **geometric
instability** — satellites appearing and disappearing, constraints churning —
rather than being maintained by a coherent NLOS conspiracy. Multi-solution
detection keyed on a stable satellite subset is therefore the wrong tool, and
section 29.1's framing should be dropped.

What remains untested is whether instability itself is detectable *before* the
solution commits: the open question is whether the per-epoch set churn predicts
transition, which would make the 24 transitions foreseeable rather than merely
post-hoc. That needs churn measured on a finer grid than 1 Hz, and is the
concrete next measurement.

### 31.1 What to stop doing

For the record, because these have each consumed multiple rounds: do not
attempt further *screening* or *rejecting* interventions (rounds 4, 21, 23, 26 all
failed the same way), do not pursue stable-NLOS-subset detection (refuted here),
and do not pursue dead reckoning (round 28 showed the solution teleports rather
than coasts). The tail is a mode-selection problem during geometric churn, and
it needs a detector, not a filter.

---

## 32. Observation churn does not predict the transitions either (round 27)

Section 31 nominated the concrete next test: if a wrong mode is entered during
geometric instability, then per-epoch satellite-set churn should rise *before*
the transition, making the 24 events foreseeable rather than post-hoc. Measured
mean churn (1 - Jaccard against the previous epoch) in a window preceding each
upward 5 m transition, against the dataset baseline. 2023 epochs, 22 transitions.

| window before transition | mean churn | ratio to baseline |
|---|---:|---:|
| 1 epoch | 0.1524 | 1.00x |
| 3 epochs | 0.1526 | 1.00x |
| 5 epochs | 0.1526 | 1.00x |
| 10 epochs | 0.1528 | 1.00x |
| 20 epochs | 0.1526 | 1.00x |

baseline = 0.1521

**No signal at any horizon.** The observation set before a transition is
statistically indistinguishable from the observation set at any other time.

### 32.1 Where this leaves the search

Four consecutive hypotheses have now been refuted by measurement, and they share
a property: every one looked for the cause in the **observation stream**.

| hypothesis | refuted by |
|---|---|
| stable NLOS subset carries the wrong mode | section 31 (set is *less* persistent inside wrong modes) |
| observation set differs materially | section 30 (counts differ; direction of causation unresolved) |
| reduced phase count causes the wrong mode | section 30/31 (reverse causality equally consistent) |
| churn predicts the transition | this section (1.00x at every horizon) |

The convergent implication is that **the discriminator is not in the
observations.** Before a transition the observation stream looks ordinary. The
signal that selects the wrong mode must live in the filter's internal state —
the ambiguity values and their variances, which accumulate history the
observation stream does not show.

That is a concrete, testable direction, and it is the first one in nine rounds
that is not a variant of screening or rejection. The specific measurement: log
the ambiguity-block state (value spread, variance distribution, correlation
condition) for the epochs leading into a transition and compare with calm epochs.
Note that section 24 already showed ambiguity state is fully allocated and
correctly sized — so the question is not allocation but the *values*.

---

## 33. A predictive signal, in the filter state — first in nine rounds

Section 32 concluded the discriminator is not in the observations and nominated
the filter's internal ambiguity state. Tested it: ambiguity-block values and
variances for epochs leading into a transition versus calm epochs. 2023 joined
Shinjuku epochs, 22 transitions.

| quantity | all epochs | calm (≤ 2 m) | ≤5 ep before transition | ratio vs calm |
|---|---:|---:|---:|---:|
| n_ambiguities | 26.55 | 27.14 | 24.54 | 0.90x |
| std of ambiguity values | 50.90 | 51.08 | **35.74** | **0.70x** |
| mean ambiguity variance | 710.7 | 615.7 | **1100.2** | **1.79x** |
| max ambiguity variance | 9328 | 9028 | 11801 | 1.31x |
| min ambiguity variance | 8.57 | 0.85 | **54.37** | **64x** |
| mean ambiguity sigma | 12.94 | 10.98 | **20.28** | **1.85x** |

**This is the first quantity found that separates the two populations.** Contrast
with section 32, where satellite churn was 1.00x of baseline at every horizon.

The signature is coherent and points at a specific failure: ambiguity values
**collapse toward each other** (spread down 30%) while their variances
**inflate** (mean sigma up 85%, and the least-constrained ambiguity up 64x). The
filter is simultaneously losing the ability to distinguish ambiguities *and*
losing confidence in them — which is what approaching a wrong fixed point looks
like from the inside, and it is invisible in the observation stream.

### 33.1 Concrete next step

This is a detector, not a filter change, which is what sections 31-32 concluded
was needed. The natural form: monitor the ratio of ambiguity-value spread to
mean ambiguity sigma, and raise a warning when it falls below the calm-epoch
distribution — a divergence indicator computed from quantities the filter
already holds, with no new measurement model.

Before building it, two things need checking, both cheap:

1. **Discrimination.** The ratios above are means over ~110 epochs versus ~1300.
   Verify the distributions actually separate rather than overlapping, and find
   the false-alarm rate at a usable threshold.
2. **Lead time.** The window here is 5 epochs. Establish how many epochs of
   warning exist, since a detector is only useful if it fires before the 10-27 m
   step that actually costs the trajectory.

Neither requires touching the estimator, so both can be done without the
six-dataset CDF sweep that any filter change would demand.

---

## 34. The ambiguity signal is a group difference, not a usable detector (round 29)

Section 33 found the first quantity separating calm from pre-transition epochs.
Tested whether it actually works per-epoch, which is what a detector requires.

Score = spread / mean sigma (low is suspicious), 2023 joined epochs:

| | p05 | p25 | p50 | p75 | p95 |
|---|---:|---:|---:|---:|---:|
| calm (n=1347) | 2.05 | 3.42 | 5.00 | 7.52 | 15.51 |
| pre-transition (n=91) | 0.00 | 1.42 | 2.76 | 4.06 | 6.58 |

The distributions do overlap-shift, but the separation is weak:

| threshold | detect | false-alarm |
|---|---:|---:|
| pre-transition median (2.76) | 49.5% | 13.2% |
| pre-transition p25 (1.42) | 24.2% | 2.2% |

Sweeping every candidate discriminator at a common 5% false-alarm rate:

| candidate | detect @ ≤5% FA |
|---|---:|
| std of ambiguity values | **34.1%** |
| min ambiguity variance | 3.3% |
| max ambiguity variance | 2.2% |
| spread / mean sigma | 0.0% |
| mean ambiguity sigma | 0.0% |

**No combination does better than the best single variable.** At an operable 5%
false-alarm rate the strongest available detector fires on only a third of the
epochs preceding a transition.

### 34.1 Why this matters

Section 33's table was a comparison of **means over 91 epochs versus 1347**. That
separates the populations in aggregate and says nothing about whether an
individual epoch can be flagged. This section makes the distinction explicit and
the answer is no.

**Do not build the detector proposed in section 33.** It would miss roughly two
thirds of the transitions while alarming on one epoch in twenty — worse than
useless, because it would generate confidence without providing warning.

The honest summary of nine rounds of mechanism search: the transition is
**visible in hindsight** across many weak correlates and **not visible in
advance** from any of them. That is consistent with it being an emergent property
of the filter's trajectory through state space rather than a condition that
precedes it — which would make it a property of the search, not of any
observation available at a single epoch.

This is a genuine stopping point for this line of enquiry. Any further progress
requires a different kind of tool than more correlation on existing logs:
either a controlled experiment that perturbs the state and observes whether a
wrong mode is entered more often, or a simulation that reproduces the transition
and can then be instrumented freely.

---

## 35. Bootstrap CDF band: the tail regression is real, and so is the bulk gain

Section 23 concluded that fixing is a net negative on F9P by comparing raw CDF
percentiles, and section 29.1 noted that no significance test was applied to
those percentiles. That gap is now closed.

Added `bootstrap_cdf_band_paired` to `gneiss_core::stats`: a **paired** bootstrap
(the same epoch indices drawn for both samples, correct because the two runs
share epochs) producing pointwise 95% intervals on the quantile difference
across the whole distribution. Deterministic xorshift PRNG, so a seed reproduces
a band exactly.

RTK Explorer F9P, AR on (69.0% fixed) vs pure float, 4504 shared epochs:

| level | AR on | float | diff | 95% interval | verdict |
|---|---:|---:|---:|---:|---|
| 0.10 | 0.034 | 0.116 | −0.081 | [−0.089, −0.076] | **AR better** |
| 0.25 | 0.090 | 0.190 | −0.100 | [−0.108, −0.093] | **AR better** |
| 0.50 | 0.187 | 0.242 | −0.055 | [−0.060, −0.051] | **AR better** |
| 0.68 | 0.260 | 0.298 | −0.038 | [−0.047, −0.030] | **AR better** |
| 0.75 | 0.322 | 0.308 | +0.014 | [+0.001, +0.021] | **float better** |
| 0.80 | 0.343 | 0.325 | +0.018 | [+0.012, +0.025] | **float better** |
| 0.90 | 0.487 | 0.376 | +0.111 | [+0.090, +0.124] | **float better** |
| 0.95 | 0.584 | 0.456 | +0.129 | [+0.112, +0.142] | **float better** |
| 0.99 | 0.809 | 0.582 | +0.227 | [+0.169, +0.270] | **float better** |

**Every level is significant, and the signs cross between p68 and p75.** There is
no ambiguity left to argue about:

- The p95 regression is **signal, not noise** — +0.129 m, interval
  [+0.112, +0.142], nowhere near zero.
- The bulk improvement is equally real — p10 through p68 all significantly
  favour fixing, by 38–100 mm.

Section 23's "net negative" framing was too strong. The honest description is a
**crossover at roughly p70**: fixing makes the bottom ~70% of the distribution
substantially better and the top ~30% substantially worse.

That is a product decision, not a statistical one, and it should be made on the
contract rather than the metric. For a user who experiences typical accuracy, AR
is a clear win — 100 mm better at p25, 55 mm at the median. For a contract
specified on p95 or p99, it is a clear loss. What is no longer available is
dismissing either half on the grounds that the other half exists.

This also settles the methodological disagreement in the useful direction:
gating on a single point percentile was wrong, but so was reading a bulk
statistic alone. The crossover is only visible when the whole distribution is
measured with intervals.
