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

## Related defects found while tracing this

- `SwfgSolution` (`swfg/engine/mod.rs:406-409`) hardcodes `clock_bias_m: 0.0`,
  `solver_iterations: 1`, `error: None`, and exposes no uncertainty at all. σ plumbing
  exists end-to-end elsewhere (`StreamingEpochSolution.std_east/north/up`,
  `SmoothedEpoch.cov_position`), so this is a reporting gap, not a missing computation.
- `screen.rs:84` uses `.max_by(|(_, a), (_, b)| a.abs().total_cmp(&b.abs()))`. NaN
  ordering in `total_cmp` is the reason a single non-finite observation is currently
  masked rather than rejected — consistent with defect C4.

## Next step

Attack the ambiguity degeneracy directly — detect the confidently-wrong float state
(small residuals, large position drift) and break it, rather than gating on a signal
that has been measured not to exist.