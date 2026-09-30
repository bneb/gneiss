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
