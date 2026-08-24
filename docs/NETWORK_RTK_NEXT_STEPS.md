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

## Measured null result: satellite PCO in DD RTK (do not re-attempt casually)

Implemented end-to-end (orbital-fixed body frame from sun ephemeris,
igs14.atx L1 PCOs for all 32 GPS PRNs, phase-centre translation applied
at satellite-position extraction) and A/B measured: fix rates identical
to 0.1%, all p50s unchanged, per-base vertical biases unchanged
(CAPO -53 mm persists with corrections on). Physics: satellite PCO is a
satellite-side error entering both stations' single differences, so it
cancels in the double difference to |d| * baseline/range ~ sub-mm at
these baselines. Code removed; this note is the durable record.

What WOULD move systematics: receiver-side antenna PCO/PCV. The set
spans three antenna families (TRM59800 rover+P181/P222/SLAC, ASH701945B
OHLN, LEIAR20 CAPO, TRM29659 P225); cross-family differential PCO/PCV
does NOT cancel between stations and plausibly explains CAPO's -53 mm
vertical bias. Requires matching RINEX antenna types to ANTEX receiver
entries AND base coordinates consistent with ARP conventions.

## Measured: north-south tropospheric gradient signature (afternoon)

Per-base mean signed vertical error by local (PDT) hour, outliers >0.5 m
excluded (WL_DUMP per-epoch CSVs). At hours 14-18 local the bases split
by bearing from the rover:

| base | bearing | h16 signed v |
|------|---------|-------------|
| P181 | NW  | +82 mm |
| OHLN | N   | +41 mm |
| CAPO | S   | -148 mm |
| P225 | ESE | -125 mm |
| P222 | SSE | -139 mm |
| SLAC | SSW | -142 mm |

Antisymmetric about the rover's latitude, coherent across four southern
stations at -125..-148 mm, peaking exactly in the marine-surge window.
This is the classical signature of a north-south wet-delay gradient:
single-scalar rover ZWD cannot represent it (no azimuth dependence in
its mapping), so the residual leaks into vertical position with sign set
by baseline azimuth.

Engineering consequence: add tropo NORTH/EAST gradient states to the
long-baseline filter (mapping m_grad(el)*[cos az, sin az], random walk,
gated behind widelane_ar like the ZWD state). Predicted effect: removes
the afternoon-signed vertical systematics on the four southern bases and
improves their float quality / fix rates. This replaces vague
"two-station tropo" as the concrete Phase-5 design.

## Tropo gradient states: implemented, gated, measured (opt-in)

Two-state [north, east] horizontal wet-delay gradient (cot(el) mapping,
satellite-minus-reference at the rover, random walk 5e-11 m^2/s, seed
sigma 2 mm), wired through float filter H AND the iono-free re-estimation
(as known correction). Enabled with `GNEISS_TROPO_GRAD=1`; default off.

Measured A/B (p95 lens, outliers >2 m excluded from shape):

- Pooled v_p95: -6.2 mm all six (163.7 -> 157.5); -7.5 mm excluding OHLN.
  P225 v_p95 alone: 249 -> 228 mm.
- Horizontal p95: P181/CAPO/P225 improve; OHLN regresses badly
  (h_p95 112 -> 181 mm, h_p50 +9 mm) — its marine-layer environment
  interacts poorly with gradient-driven float shifts.
- Long baselines (P222/SLAC): only ~7 of 2875 epochs change (rare AR
  decision flips in the tails) — corrections are sub-rounding there.

Verdict: mechanism partially validated; kept opt-in until OHLN is either
absorbed by a two-station model or conditioned per-base. Candidate next
tuning: tighter GRAD_RW or az-binned robustness for OHLN-class data.

Tooling added: scripts/compare_dumps.py (A/B quantile comparator over
WL_DUMP CSVs); eval Vertical Error line now carries p95 (appended after
RMS to preserve guard's positional parsing).

## Receiver PCO (differential datum form): CAPO vertical bias fixed

Two fixes landed together:
1. q_accel regression repaired: a lock-Q sweep had silently rewritten the
   eval's convergence phase to 1e-9, disabling the two-phase design's
   loose phase since it landed. Restored to 1e-6 -> lock 1e-8 after
   15 min. OHLN h_p50 recovers 32 -> 23 mm; P181/CAPO pay ~1-2 mm.
2. `GNEISS_RECV_PCO=1` applies the DIFFERENTIAL receiver L1 PCO
   (base minus rover antenna, both ARP-referenced) to the base position.
   Same-family baselines shift <=5 mm; CAPO's Leica LEIAR20 is +39.7 mm
   Up vs the rover Trimble family.

Measured: CAPO v_p50 -54 -> -14 mm (74% of the daily bias removed);
all other bases unchanged within noise. CAPO v_p95 rises 117 -> 156 mm —
the elevation-dependent PCV residual that PCO cannot address; that plus
receiver PCV generally is the next antenna-modeling step. Multi-GNSS
scoping closed: no Galileo in any base file (GLONASS partial), so DD
multi-constellation is data-blocked, not code-blocked.
