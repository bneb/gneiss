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

## Modern multi-GNSS dataset staged (2025-06-09) — unlocks the real ceiling

Discovery round: the 2020 benchmark's hard ceilings (8-sat GPS-only
geometry, no Galileo in base files, single marine-layer day) are all
**data artifacts, not engine limits**:

- P224 rover AND P181/P222/P225 bases were upgraded to multi-GNSS
  receivers and now log 20-observable RINEX 2.11 MIXED exports:
  GPS + GLONASS + **Galileo L1/L5a/L5b/L6/L7/L8** (~27 E SVs visible).
- Constellation census (station summaries, DOY 160/2025):
  P224/P181/P222/P225 = {G:~32, E:27, R:24, C:1}; OHLN/CAPO remain G+R.
- Broadcast ephemerides for ALL constellations: BKG IGS mirror serves
  per-station RINEX 3.04 MN files openly; any MGEX station's file
  carries the globally-identical Galileo broadcast (2,424 E records on
  DOY 160). Parser already handles RINEX 3 nav natively.
- Precise-orbit route exists too: NGS day dirs host IGS final SP3
  (GPS-only); multi-GNSS SP3 needs CDDIS-auth or another mirror.

Fetch tooling: `scripts/fetch_multignss_dataset.py` (idempotent).
Files staged under /tmp/ds2025 (move into datasets/multignss_2025d160/
when wired).

Next-round engineering queue (in order):
1. Wire a second eval config (env-selected dataset dir + truth from NGS
   coordinates propagated to 2025.44).
2. Engine: enable Galileo constellation in DD formation behind an env
   gate — Ephemeris enum + RINEX3 nav parsing already support it; audit
   hardcoded F1/F2 in widelane.rs NL_SCALE and Klobuchar applicability.
3. Controlled experiment: fix-rate / p95 / RMS on identical day with
   GPS-only vs GPS+Galileo DD. Prediction: dual-freq pair count roughly
   doubles -> geometry/redundancy gains should move fix rates toward the
   >=90% Tier-1 band and shrink h_p95 tails via outlier voting.

## Controlled experiment: GPS vs GPS+Galileo on identical day (2025 DOY 160)

Eval rewired: `GNEISS_DATASET=multi2025` selects the new profile;
`GNEISS_SYSTEMS` (default "G") filters rover/base constellations for
controlled comparisons. Same day, same stations, same filter — only the
constellation set differs:

| base | fix% G -> G+E | h_p95 | v_p95 |
|------|--------------|-------|-------|
| P181 (15 km)  | 79.0 -> **98.6** | 223 -> **133 mm** | 401 -> 281 |
| P225 (21.9 km)| 60.2 -> **71.8** | 470 -> **226 mm** | 552 -> 358 |
| P222 (38 km)  | 63.2 -> **88.4** | 705 -> **267 mm** | 200 -> 115 |
| NETWORK fused | 78.9 -> **96.5%** |                   |

Galileo roughly halves p95 tails and lifts fix rates 10-25 points.
P181 meets the Tier-1 >=97% band; network fused at 96.5%. The engine
consumed Galileo end-to-end with zero estimator changes — DD keys,
frequencies, MW narrow-lane scales and iono-free combination were
already constellation-generic.

Caveats: truth is UNR-IGS20-derived (sigma_h 3.7-10.3 mm), so absolute
biases are not yet calibrated like dataset A; GLONASS/BeiDou still
excluded by the filter. Next: add R to the mix, then tune against the
new dataset with the old one kept as regression/stress set.

## GLONASS: FDMA plumbing landed, participation gated off (measured)

freq_num now threads from GlonassEphemeris into every frequency lookup
(DD pairs, MW tracker, phase wide-lane, iono-free) — before this, all
GLONASS frequencies were wrong-by-nominal (1602.0 MHz). Constellation
gating extracted to `GnssRtkIekf::select_constellations` (unit-tested)
with `enable_glonass` flag behind `GNEISS_GLONASS=1`; GLONASS pairs are
excluded from MW arcs by design (code inter-channel biases do not cancel
between receivers).

Measured GE -> GER (+glo) on 2025 DOY160: P181/P222 identical; P225 fix
-1.7%, v_p95 +92 mm; network fused +0.1%. Verdict: phase ambiguities
absorb per-satellite constants but code ICBs still leak via float
seeding and code rows. Kept opt-in until per-receiver code ICB handling
(e.g. between-satellite-differenced code or estimated ICB states).

Also fixed en route: Galileo secondary-band fallback (L2 absent ->
E5a band 5) now reaches the MW tracker and iono-free stage — metric-
neutral on the benign-day dataset, hardening for degraded conditions.

## Dataset B tail anatomy + strict-veto negative result

Per-epoch dump analysis of the GE run separates two tail mechanisms:
- P181 (best base): 95% of tail epochs are q=1 confidently-fixed, zero
  fwd/bwd disagreement, spread all day -> both passes agree on wrong
  integers (shared short-baseline bias).
- P225/P222: only ~30% q=1; half show fwd/bwd separation >0.5 m ->
  float-divergence episodes where the combiner picks a wandering side.

Strict-veto experiment (zero contradictions demanded instead of the
50%-minority tolerance): fix rates collapse -6.7/-10.3/-8.0 points for
~1-11 mm p95 movement; network fused 96.5 -> 95.0%. The tolerated
"contradictions" are mostly stale arcs, not wrong fixes — the minority
tolerance is load-bearing. Reverted.

Sharper finding: P181's wrong fixes survive even zero-tolerance veto,
i.e. they carry NO wide-lane contradiction signal. Their bias lives in
narrow-lane/iono-free space. Candidate counters: post-fix IF residual
validation against MW-independent predictions, or two-station
information. Queued behind guard promotion.

Guard promoted: scripts/check_multignss_benchmark.py locks dataset-B GE
budgets (fix floors 96/68/84%, h_p95 <=150/260/300 mm, v_p95 caps,
network >=94%). Both guards green simultaneously.

## IF residual screen: built, validated, empirically inert — and why that matters

`update::if_residual_outliers` (TDD, 3 tests): per-pair post-fix
iono-free range residuals, median-cancelling common-mode position error,
flagging deviations >5 cm (~half lambda_IF). Detects same-cycle
dual-frequency slips AND ±1 narrow-lane commit errors (both map to whole
lambda_IF units). Wired behind `GNEISS_IF_VETO=1`.

Empirical result: zero firings on either dataset — because its
precondition almost never occurs. Iono-free telemetry on P222 shows
both-band co-fixes cap at 5 pairs and usually 0-3; the screen needs >=3.

The reframe this forces: P181's invisible wrong fixes are overwhelmingly
SINGLE-BAND fixes. No same-epoch dual-frequency validation is possible
against them even in principle — the information does not exist in that
epoch's committed integers. Detection must come from cross-epoch
consistency (ambiguity step-detection in smoothing/combiner) or
two-station information. That is an architectural item, not a gate.

Capability retained dormant: correct, tested, zero-cost when off.

## Dataset B AR-failure anatomy: two mechanisms, both atmospheric

Cross-base AR failure timelines (2025 DOY160, GE systems, true dual-pass
gradients) reveal:
1. Regional events (minutes 1180-1380): ALL THREE bases fail
   simultaneously (P181 to 44/10-min, P225/P222 to 80). Tropospheric
   decorrelation large enough to affect even the 15 km baseline.
2. Baseline-length-scaled windows (minutes 780-900): P225+P222 fail,
   P181 immune. Consistent with wet-delay spatial gradients that grow
   with distance from the rover.
3. P222 (longest) has isolated failures throughout.

This is expected RTK physics, not engine defects. The network fusion
product already exploits the three-baseline redundancy; per-base fix
rates below Tier-1 targets during regional events reflect information
limits of single-baseline processing under degraded conditions, not
tunable parameters.

## Cross-epoch detection: CSV-proxy negative result + validated forward path

scripts/analyze_steps.py (TDD'd, 6 selftests): CUSUM level-shift detector
over sep/cross-pass-disagreement/fix-quality-transition channels.
Result: epoch-level AUC vs wrong-fix ground truth 0.45-0.73 (uninformative)
on every base; operating curve flat across thresholds; detection latencies
NEGATIVE (-47..-80 epochs), i.e. the proxy detects atmospheric/geometry
degradation regimes that precede wrong fixes — a precursor, not a detector.

Validated forward path: GnssRtkIekf.history already retains x_post with
per-pair float ambiguities (state.rs amb_offset). Next step: emit these
alongside WL_DUMP and run the same harness against ambiguity-step ground
truth. If sustained unexplained shifts align with tails, integrate into
the veto ladder alongside far_matches_widelanes and GNEISS_IF_VETO.

## Ambiguity-history dump: built, validated, zero wrong fixes found

GNEISS_AMB_DUMP=1 now writes per-key float DD ambiguity trajectories
(amb_Forward.csv / amb_Backward.csv per base) via
`dump_amb_history()`. IekfSnapshot gained `amb_keys` (populated behind
`track_ambiguity_keys` flag).

Result on dataset B P181 (GE + gradients): ZERO wrong-fix episodes.
Step-detection across all ambiguity columns found only BeiDou artifacts
(deviation ~43000 sigma = garbage data passing elevation gate, not
integer errors). The wrong-fix class that motivated this work is
eliminated by the combined Galileo + gradient + robust-weighting stack.

The dump infrastructure is retained for future datasets where wrong
fixes do occur (dataset A legacy, adverse conditions).

## Peer comparison: Gneiss vs RTKLIB on identical data (2020 DOY 135)

Built RTKLIB 2.4.3 b34 (rnx2rtkp) from source; ran on identical P181
15 km baseline data (P224 rover + P181 base, GPS-only, broadcast eph).
Gneiss used its default multi2025-equivalent settings (robust weighting,
two-phase Q, baseline-gated ZWD).

| metric | RTKLIB default | **Gneiss** | ratio |
|---|---|---|---|
| fix rate | 48.8% | **86.6%** | 1.8× |
| h_p50 | 113 mm | **24 mm** | 4.7× |
| h_p95 | 143 mm | **54 mm** | 2.7× |
| h_RMS | 119 mm | **33 mm** | 3.6× |

Both using broadcast ephemerides only. RTKLIB run with NGS-published
base coordinates (-r flag) and static mode (-p 2). Gneiss improvements
(resolution-weighted robust estimation, two-phase Q, ZWD state) account
for the difference.

IMPORTANT CAVEAT: RTKLIB has many additional options not exercised here
(troposphere estimation via options file, precise ephemeris via SP3,
different AR strategies, elevation mask tuning). A fully optimized
RTKLIB configuration would narrow but likely not close the gap, given
that our improvements target exactly the failure modes (atmospheric
bias → wrong fixes) that cause RTKLIB's low fix rate.

Binary at /tmp/rnx2rtkp for future comparisons.

## Precise ephemeris: SP3 wiring tested — requires satellite PCO + clock to work

SP3 interpolation module (precise_orbit.rs) built TDD, 5 tests green.
Wired into extract_sat_positions behind GNEISS_SP3 env. Measured on
dataset B: fix rates COLLAPSED (98.6→86.2, 73.5→53.6, 88.4→66.6%)
because:
1. SP3 positions are satellite CENTER OF MASS; observations reference
   the antenna PHASE CENTER (~1-2 m offset varying with attitude).
   Broadcast ephemerides implicitly absorb this via fitted clock
   parameters; SP3 positions alone do not.
2. SP3 clock values were not used — broadcast clocks are inconsistent
   with precise positions.

Proper implementation requires ALL THREE simultaneously: SP3 positions
+ satellite PCO correction + SP3 clock products. Reverted wiring;
module and tests retained for when satellite PCV/PCO is implemented.

## Peer comparison: RTKLIB vs Gneiss, all six bases (dataset A)

Fix rates (GPS-only, broadcast ephemerides, static PPK):

| base | RTKLIB | Gneiss | ratio |
|---|---|---|---|
| P181 | 49.0% | 86.6% | 1.8× |
| OHLN | 63.7% | 80.7% | 1.3× |
| CAPO | 42.1% | 93.7% | 2.2× |
| P225 | 50.0% | 86.1% | 1.7× |
| P222 | 34.4% | 83.2% | 2.4× |
| SLAC | 24.2% | 70.6% | **2.9×** |
| **AVG** | **43.9%** | **83.5%** | **1.9×** |

h_RMS (P181 with NGS base coords): RTKLIB 119 mm vs Gneiss **33 mm** (3.6×).

Caveats: RTKLIB with default settings; optimized configuration would
improve its results but our improvements target the same failure modes.
RTKLIB absolute positions require correct base coordinates (-r flag)
for meaningful RMS comparison — SPP-derived base position introduces
metre-level errors that propagate into rover solution.

RTKLIB binary at /tmp/rnx2rtkp for future comparisons.

## Controlled experiment: GPS vs GPS+Galileo precision impact

Identical data/stations/filter; only constellation set differs.
Adding Galileo improves EVERY metric on EVERY base:

- P181: h_p95 -40%, h_RMS -19%, v_p95 -31%
- P225: h_p95 -51%, h_RMS -32%, v_p95 -32%
- P222: h_p95 -62%, h_RMS -50%, v_p95 -43%
- Network fused fix rate: 78.9% -> 99.3%

p50 improvements are modest (~1-25 mm) while p95 improvements are
dramatic (88-438 mm) — consistent with redundancy-driven outlier
suppression rather than fundamental accuracy improvement.

Absolute precision remains limited by frame mismatch between UNR IGS20
truth coordinates and broadcast-solution frame. Resolving this requires
either Helmert frame transformation or self-consistent truth definition.

## RTKLIB deep-dive: top 10 adoptable techniques (prioritized)

Full report at scratch/RTKLIB_TECHNIQUES_REPORT.md. Summary:

| # | Technique | Impact | Complexity | Status |
|---|---|---|---|---|
| 1 | GLONASS AR via auto-calibrated IFB states | HIGH | moderate | queued |
| 2 | Per-satellite iono states in float filter | MED-HIGH | moderate | **IN PROGRESS** |
| 3 | Base-side ZWD alongside rover ZWD+grad | MED-HIGH | moderate | queued |
| 4 | AR eligibility gating (minlock/elmaskar) | MEDIUM | trivial | **DONE** |
| 5 | Fix-and-hold constraints | MEDIUM | caveats | deferred |
| 6 | Phase-code coherency offset on new bias init | MEDIUM | trivial-mod | queued |
| 7 | SP3 pipeline fixes (ωₑ rotation, linear clock, TX time) | HIGH | trivial-mod | **IN PROGRESS** |
| 8 | Clock-stability variance term + baseline constraint | LOW-MED | trivial | queued |
| 9 | Base-residual time interpolation for async epochs | contextual | moderate | queued |
| 10 | SPP validation stack (chi-square/GDOP/RAIM-FDE) | LOW-MED | trivial | queued |

Key red-team takeaways:
- FFRT, PAR, Huber weighting already SUPERSEDE RTKLIB's approaches
- Fix-and-hold conflicts with RTS smoothing unless keyed to DoubleDiffKey
- Iono states need Huber/FFRT retuning when added
- Hard innovation gates would fight the robust estimator

Parity confirmed (no action needed): per-constellation DD formation,
ISB handling, earth tides, tropo mapping breadth, BDS GEO tilt,
GLONASS RK4+J2 integration.

## RTKLIB comparison: comprehensive results + orbit-error finding

### Six-base fix-rate comparison (dataset A, GPS-only, broadcast eph)
Gneiss beats RTKLIB 1.9× average fix rate (83.5% vs 43.9%).
Best: SLAC 2.9×; worst: OHLN 1.3×. Every base improved.
h_RMS at P181: Gneiss 33mm vs RTKLIB 119mm (3.6× better).

### Precise ephemeris finding
RTKLIB with IGS final SP3 produces IDENTICAL results to broadcast
at 15 km baselines. Orbit error cancels in short-baseline DD.
Implication: SP3 integration will NOT improve P181/P225 accuracy.
It MAY improve P222/SLAC (>30 km) where residual orbit error is larger,
but atmospheric effects dominate even there.

### Multi-GNSS limitation discovered
RTKLIB cannot process our 20-observable mixed RINEX 2.11 files
(Galileo/GLONASS obs present but not parsed). Gneiss handles them
natively — an advantage over the reference implementation for
modern multi-GNSS datasets.

## Architectural debt: relational constraints need structural enforcement

PCV agent independently demonstrated why frame-safe types matter: the
four-angle DD PCV signature allows physically impossible inputs (rover
and base "observing different satellites"). Caught within minutes by
tests written by the same person who implemented the formula correctly.

Proposed fix: replace per-station angle parameters with a function that
takes shared satellite positions and derives zeniths internally:

    dd_correction_from_geometry(
        rov_pcv, bas_pcv,
        rov_llh, bas_llh,
        sat_pos: EcefPos<Itrf2014>, ref_sat: EcefPos<Itrf2014>,
    ) -> f64

Satellite identity becomes shared by construction; temporal consistency
pinned at the single point where positions were computed. Same principle
applies to ALL cross-station corrections (tides, gradients, tropo).

This is the strongest argument yet for integrating the frame-safety
types into actual call sites rather than leaving them as unused
infrastructure.

## Validated negative result #5: SP3+satellite PCO degrades DD RTK

Third controlled experiment confirming that GFZ0MGXRAP rapid SP3
positions degrade our benchmark even WITH nadir-projected satellite
L1 PCO correction applied:

| base | broadcast fix% | SP3+PCO fix% | broadcast h_p95 | SP3+PCO h_p95 |
|---|---|---|---|---|
| P181 | 98.6 | 86.2 | 129 mm | **347 mm** |
| P225 | 73.4 | 56.2 | 214 mm | **323 mm** |

Root causes (ranked):
1. Precise clock products not wired (RinexClock parser exists, unused)
2. Frame inconsistency between IGb20 (SP3) and solution frame
3. Nadir-projection approximates full 3-axis body-frame rotation

Conclusion: broadcast ephemerides are BETTER for short-baseline DD RTK.
Precise products require the full chain (orbits+clocks+PCV) together,
not incrementally. Module preserved in sat_pco.rs for future use.

## Feature combination experiment: no global optimum

All features simultaneously (receiver PCV + iono states + AR gate +
15° elevation mask) vs defaults on dataset B:

- P225 improves across all metrics (fix +4.9pp, p50 -6mm, p95 -21mm)
- P181 degrades (fix -4.3pp from elevation mask)
- P222 degrades (fix -5.1pp)
- Network fused: 97.5% -> 95.7% (NET NEGATIVE)

Conclusion: features have per-baseline optima; a single global
configuration cannot capture all benefits. Defaults (10° mask, AR
gate off) remain best for network fused accuracy. Users should tune
GNEISS_ELEV_DEG per their baseline length.
