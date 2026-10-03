# Shinjuku p95 Tail — Root Cause Analysis

**Status:** mechanism characterised; fix not yet designed.
**Measured at** `55cb240` (+ `eval_f9p_rover` error-dump change, numerics unchanged).
**Reproduce:** `GNEISS_ERR_DUMP=/tmp/sj ./target/release/eval_f9p_rover shinjuku`

## Headline

| subset | n | p50 | p68 | p95 | max |
|---|---|---|---|---|---|
| **fixed** (`quality==1`) | 159 | 0.821 | 1.185 | **2.024** | 3.644 |
| **float** (`quality==2`) | 1937 | 1.628 | 2.193 | **22.783** | 39.290 |
| all | 2096 | 1.545 | 2.043 | 22.737 | 39.290 |

Fixed solutions are already Tier-1. **The entire objective gap is the float tail** —
16.1% of float epochs exceed 5 m. Every epoch >5 m is labelled float; zero are labelled fixed.

## Excursion morphology

21 contiguous runs; 312 epochs total. Median run length 3, but the tail mass is in
three long runs (90, 90, 52 epochs) that alone account for 232 of 312 bad epochs.

Run 1 (epochs 760–849, tow 282035–282124), with truth speed from `reference.csv`:

| idx | tow | err (m) | truth speed (m/s) |
|---|---|---|---|
| 755 | 282030 | 1.44 | 7.46 |
| 760 | 282035 | **8.11** | 7.70 |
| 765 | 282040 | 29.08 | 4.43 |
| 770–835 | 282045–282110 | ~23 | **0.00 (stopped)** |
| 850 | 282125 | 4.82 | 2.52 |

Three facts define the mechanism:

1. **Onset is discontinuous** — 0.50 m → 8.11 m in a single epoch at 7.7 m/s. This is
   not gradual filter divergence.
2. **The error is frozen while stationary.** For 70 s the rover does not move (0.00 m/s)
   and the error moves only 0.6 m (23.4 → 22.8). A float filter with full GNSS at a
   standstill collapses 23 m → <1 m in seconds. Ours does not move at all.
3. **Recovery is instantaneous on motion resume**, not gradual.

## Ruled out

- **A hard-coded 90-epoch window.** Runs 1 and 2 are both exactly 90 epochs, which
  suggests a fixed constant. Grepped `crates/gneiss-rtk/src/` — no 90-epoch/90-second
  convergence, reset, or lock constant exists in the pipeline.
- **Dead reckoning.** The error is largest while the rover is *stopped*, and stationary
  dead reckoning cannot accumulate 31 m.
- **Measurement-count / gross-error screening.** `screen_gross_pr_errors`
  (`estimators/rtk_iekf/screen.rs:28`) caps removals at 3 per epoch and de-weights rather
  than drops when carrier phase is present; it cannot freeze the solution.
- **A single dropped satellite.** 14 of 21 runs occur at ≥0.9 m/s.

## The confident-wrong signature

Reported σ_h (= √(σ_E·σ_N)) versus actual error:

| | σ_h p50 | σ_h p95 | σ_h max |
|---|---|---|---|
| good epochs (≤5 m) | 0.097 | 0.353 | 0.806 |
| bad epochs (>5 m) | **0.235** | 0.893 | 1.461 |

Bad epochs carry a median reported uncertainty of 0.235 m while being 22 m wrong — the
filter is **~100× overconfident**. Small residuals with a large constant position offset
is the classic signature of a wrong double-difference ambiguity set: the ambiguities
absorb the offset exactly, so the innovations stay small, information stays high, and
the covariance stays small while the position is wrong.

This also explains the frozen standstill plateau: with the wrong ambiguities the
solution is *self-consistent*, so standing still changes nothing and no amount of
additional data improves it. Moving changes the geometry, the wrong set becomes
inconsistent, and the filter snaps back.

## Negative result: gating cannot fix this

Both available confidence signals were tested directly. Neither separates good from bad.

**(a) `quality` is binary and uninformative within float.** Epoch 759 is `quality==2`
at 0.50 m; epoch 767 is `quality==2` at 31.31 m — identical label, 60× the error.
1937 of 2096 epochs carry the same label.

**(b) σ does not separate.** Sweeping a σ_h threshold:

| thr (m) | kept | of which bad | p50 | p95 | max |
|---|---|---|---|---|---|
| 0.04 | 1762 | 312 | 1.689 | 23.037 | 39.290 |
| 0.10 | 1143 | 269 | 1.960 | 23.466 | 39.290 |
| 0.20 | 579 | 170 | 2.666 | 26.134 | 39.290 |
| 0.30 | 300 | 142 | 4.363 | 29.159 | 39.290 |
| 0.50 | 95 | 75 | 12.491 | 31.041 | 39.290 |

There is **no usable operating point**: at σ≥0.50 m, 75 of 95 retained epochs are bad
and p95 is still 31 m.

**This retro-explains an earlier negative result.** A previous round measured a
per-epoch fallback selector at 57.6% and concluded it was "strictly dominated". That
was correct but mis-attributed: the selector had no signal to select on. `quality` is
binary across a 130× error range and σ does not separate. **The tail cannot be fixed by
gating; it must be fixed at the source, by preventing the confidently-wrong state.**

## Round 2 — internal instrumentation: root cause found

The gating analysis above is correct but incomplete. Instrumenting `GnssRtkIekf::process_epoch`
directly (temporary, reverted) and joining the filter's own ECEF positions to `reference.csv`
truth identifies the mechanism.

**The excursion is a constant ECEF position offset, not divergence.**

During tow 282060–282075 (rover stationary), the ECEF error vector `(pos − truth)` is:

```
(-10.0, +37.4, +14.2) m,  |b| = 41.2 m
angular deviation from the mean direction: p50 = 0.1 deg, p95 = 0.3 deg
```

A fixed vector to within 0.3° for ~90 s. The prior window (282036–282050) shows the
transition: the bias grows 11.9 → 50.9 m over 8 epochs, then holds.

This single fact explains every previously puzzling observation:

| observation | explanation |
|---|---|
| error grows at ≈ vehicle speed | a constant bias accumulates as the vehicle drives away from it |
| error frozen while stationary | the bias does not move, and neither does the vehicle |
| σ_h p50 = 0.235 m while 22 m wrong | the wrong state is *self-consistent* — small innovations, high information |
| instant recovery on motion resume | geometry change breaks the wrong state |
| DD count steady at 27–31 | no satellite is lost; the solution is not starved |

This is the **DD float ambiguity–position degeneracy**: position and double-difference
ambiguities are only weakly separable under the geometry present here, and the filter
settles into a mutually consistent wrong pair.

### Hypotheses tested and eliminated this round

| hypothesis | result |
|---|---|
| Phase-innovation slip gate re-seeding ambiguities (`mod.rs:243-263`) | **0 events** in 800 epochs; `widelane_ar` is enabled so the path is live |
| Gross-error screen rejecting satellites (`screen.rs:41`) | **0 events** in 800 epochs |
| Satellite dropout at onset | rover BeiDou/GPS counts *rise* at tow 282035 (6→7→8 GPS) |
| Raw filter σ_p as a gate | **falsified** — as the threshold rises the fraction of bad epochs kept *increases*, 28% → 50%; error 8 epochs after a σ firing has p50 only 3.2 m |
| Filter lag / dead reckoning | **falsified** — filter moved 95 m over 12 epochs vs ~66 m of truth; it overshoots, and the error vector is direction-stable |

Note on the σ_p result: σ_p does rise 0.082 → 0.311 four epochs before the position
departs, which looks like early warning, but it is not usable — most σ firings are
benign, so it has no precision as a gate.

### Practical implication

Because the wrong state is *self-consistent*, no residual- or covariance-based detector
can see it: the filter's own evidence supports it. Detection has to come from a
**different** source — cross-epoch position consistency (a stationary rover's position
must not translate) or an independent solution (PPP, or a second reference satellite set).

The principled long-term cure is partial ambiguity resolution: fixing only the
well-constrained ambiguities (long arc, high elevation, low variance) pins the position
and destroys the degeneracy by construction. The repository already has the machinery
(`mw/`, `ar_subsets/`, `widelane.rs`) but does not apply it at these epochs — every
excursion epoch reports `quality == 2` (float).

## Related defects found while tracing this

- `SwfgSolution` (`swfg/engine/mod.rs:406-409`) hardcodes `clock_bias_m: 0.0`,
  `solver_iterations: 1`, `error: None`, and exposes no uncertainty at all. σ plumbing
  exists end-to-end elsewhere (`StreamingEpochSolution.std_east/north/up`,
  `SmoothedEpoch.cov_position`), so this is a reporting gap, not a missing computation.
- `screen.rs:84` uses `.max_by(|(_, a), (_, b)| a.abs().total_cmp(&b.abs()))`. NaN
  ordering in `total_cmp` is the reason a single non-finite observation is currently
  masked rather than rejected — consistent with defect C4.
- `eval_f9p_rover.rs` hardcoded `with_env_filter("warn")`, ignoring `RUST_LOG`. Now
  reads `GNEISS_LOG`, defaulting to `warn`.

## Next step

Break the degeneracy at its source. Candidate detectors that do **not** rely on the
filter's own (demonstrably unreliable) residuals:

1. **Stationarity consistency** — while the rover is stationary, the estimated position
   must not translate. Run 1 is stationary for 70 s with a 41 m static offset; this is
   directly observable without truth.
2. **Cross-epoch position-rate consistency** — an implied ground speed inconsistent
   with the Doppler-derived speed is a model-inconsistency signal.
3. **Partial ambiguity resolution** during kinematic float epochs.

(1) is cheapest and needs no truth, but must not misfire on genuine motion; (3) is the
principled fix.