# Network RTK — Verified State & Next Architecture Steps

Status as of branch `network-rtk-long-baseline` (@ 696657d). All numbers
below are reproducible via `scripts/check_network_benchmark.py` (full-day
CORS set: rover P224 + bases P181/OHLN/CAPO/P225/P222/SLAC, 15–50 km).

## Verified results (smoothed, full day)

| base | km | horiz p50 | fixed-only p50 | fix rate | vert RMS |
|------|----|-----------|----------------|----------|----------|
| P181 | 15.0 | 24 mm | 23 mm | 96.9% | 106 mm |
| OHLN | 16.5 | 23 mm | 22 mm | 95.1% | 323 mm |
| CAPO | 16.6 | 33 mm | 32 mm | 96.9% | 60 mm |
| P225 | 21.9 | 38 mm | 37 mm | 97.8% | 130 mm |
| P222 | 38.0 | 72 mm | 68 mm | 89.5% | 749 mm |
| SLAC | 49.7 | 101 mm | 94 mm | 82.5% | 779 mm |
| **NETWORK** | — | **28 mm** | **28 mm** | **97.6%** | **79 mm** |

Regression guard: `scripts/check_network_benchmark.py` (nine budgets,
fails closed on parse errors; vertical budget would have caught the
round-7 ZWD divergence at product level).

## What produced the gains (all opt-in via `widelane_ar`)

1. Cadence-aware cycle-slip gating — 30 s data used to trip the gap test
   every epoch, re-seeding ambiguities forever.
2. Base-stream slip detector — base counter resets were invisible.
3. Reverse-safe process noise — negative-dt diagonals corrupted long
   backward passes.
4. Innovation-gated slip re-seeding — undetected slips can no longer
   destroy the filter state.
5. Multi-base consensus fusion — component-wise median over agreeing
   fixed candidates; empirical-scatter covariance; continuity gate.

## Measured dead ends (do not retry as-is)

- **Per-pair DD-iono random walk**: iono/N slide together along the
  per-band null direction (ΔN/ΔI ≈ disp/λ₁); unobservable without
  spatial structure.
- **SD-MW UPD correction**: raw MW arc means carry per-(station,sat)
  code multipath; differencing satellites leaves it in, and applying it
  injected bias (P222 RMS +59%).
- **Single-station scalar ZWD**: diverges unconstrained; step-saturated
  version still net-negative everywhere (SLAC vertical 0.78 → 2.19 m).
  The DD wet residual contains BOTH station residuals mapped by nearly
  identical elevations — a rover-only state absorbs non-tropo error into
  the wrong bucket.

## Next architecture steps (in order)

### 1. State-space wet troposphere inside the DD filter
Two ZWD states (rover + active base) with NMF wet mapping differences,
random-walk Q ≈ 3e-7 m²/s each, seeded ~15 cm. Elevations between
stations differ by <0.5°, so the two columns are weakly separated — the
design must (a) estimate their DIFFERENCE as primary and (b) tie the sum
toward zero with a loose pseudo-observation, or interpolate base-side
residuals from the network instead of estimating them. Target: P222/SLAC
vertical RMS < 0.3 m. Requires the matrix-state plumbing added in
`enable_zwd` (currently dormant behind `widelane_ar`).

### 2. Phase-only network UPD estimation (activates the MW cascade)
Raw MW carries per-(station,sat) code multipath — measured. Phase-only
alternative now feasible BECAUSE fusion exists: take the fused network
trajectory as a known rover trajectory, recompute each base's wide-lane
floats against it (geometry fully determined → no code term needed),
then decompose fractional residuals across bases × satellites
(Σu_sat = 0) exactly as PPP-AR networks do. Feed corrected integers to
the dormant cascade (`widelane::resolve_cascade`) and the FAR veto
(`far_matches_widelanes`). Expected: fix-rate headroom at P222/SLAC and
honest fixes through disturbed windows.

### 3. Combiner trust model
`combiner.rs` blesses co-wrong pass agreement with hard-coded gates
(0.5/0.2/10 m). End-of-day windows show both passes agreeing on wrong
fixes ~14 m apart from truth. Per-base continuity gating (eval-level)
mitigates the product; moving the gate into combiner semantics would
make every consumer honest. Must stay gated to avoid walkthrough drift.

### Measured blocker: iono-free engagement starvation at long baselines
Instrumented outcome counts (full day, RUST_LOG=debug):
- P181 (15 km): 1,480 Solutions / 3,285 NotEngaged+Rejected.
- SLAC (49.7 km): **7** Solutions / 292 NotEngaged+Rejected.
The >=6-both-band-pairs floor plus gate failures starve the stage exactly
where it matters most: at 49.7 km the fixed set rarely carries six
both-band pairs simultaneously (partial per-band fixing). Until this is
addressed (lower floor with stricter per-pair validation, or partial-set
IF solving using whichever both-band pairs exist), long-baseline products
ride per-band projected positions and inherit their atmospheric error —
the dominant term in SLAC/P222 vertical and horizontal tails.

### Post-calibration measurement & DD separability limit
After baseline-aware variance calibration, iono-free engagement improved
13x at SLAC (7 -> 92 Solutions full-day) - the calibrated gate works.
However SLAC smoothed vertical RMS did NOT improve correspondingly,
which surfaces a hard architectural limit: **DD observations cannot
separate rover-ZWD from base-ZWD.** Their mapping-function coefficients
are nearly collinear (stations 40 km apart share satellite elevations
within ~0.3 deg), so only a linear combination is observable, and
round-8 evidence shows even that combination does not improve vertical
at long baselines when absorbed by a single state. Breaking the
degeneracy requires non-DD information: absolute tropo products (VMF1
grids), PPP solutions at reference stations, or a network-side
interpolation model. Until one of those exists, vertical accuracy at
40+ km baselines is bounded by the Saastamoinen model (~0.7 m RMS
observed) rather than by the estimator.

### Known benign anomalies
- OHLN vertical: 93 episodic excursions (|v| > 30 cm, RMS 1.4 m over
  them) with unbiased p50 — wet-tropo/multipath activity specific to
  that baseline; diluted by fusion. Environmental, not a defect.
