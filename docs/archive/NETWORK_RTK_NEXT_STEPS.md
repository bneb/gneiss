> **Superseded.** This document was the working log for Sprints 7-10 (Network RTK prototype). For current status and roadmap, see `docs/PROJECT_STATUS.md` and `docs/TIER1_ROADMAP.md`.

# Network RTK — Verified State & Next Architecture Steps

Status as of branch `sprint7/tropo-zwd` (@ 2315af6, forked from
`network-rtk-long-baseline`). All numbers below are reproducible via
`scripts/check_network_benchmark.py` (full-day CORS set: rover P224 +
bases P181/OHLN/CAPO/P225/P222/SLAC, 15–50 km).

## Verified results (smoothed, full day)

**Refreshed 2026-08-26** — the table below was stale for at least two
rounds' worth of unconditional fixes (this file's own findings below it
document them; the header simply never got re-rolled after `696657d`).
Measured fresh via `./target/release/eval_network_ppk`, no env vars,
zero code changes from HEAD: P222 and SLAC vertical RMS are 6–7x better
than this table previously claimed (749→102mm, 779→209mm) — the
tropospheric-ceiling premise that once motivated a two-ZWD-state design
no longer holds at HEAD. OHLN and CAPO drifted the other way (worse vert RMS than previously
recorded). CAPO's p50 bias fix (below, "Receiver PCO") is on record and
is now default-on, but the table's older 60mm CAPO vert RMS predates
that fix entirely and the full 60→101mm chain across intervening commits
hasn't been re-traced here — flagged, not explained. OHLN is likewise
unexplained; see "Known benign anomalies" below, now the network's
single worst vertical number and worth a fresh look rather than assuming
the old "environmental, not a defect" read still applies at 430mm.

| base | km | horiz p50 | fixed-only p50 | fix rate | vert RMS |
|------|----|-----------|----------------|----------|----------|
| P181 | 15.0 | 24 mm | 24 mm | 86.6% | 86 mm |
| OHLN | 16.5 | 23 mm | 22 mm | 80.7% | 430 mm |
| CAPO | 16.6 | 32 mm | 31 mm | 93.7% | 101 mm |
| P225 | 21.9 | 37 mm | 36 mm | 86.1% | 117 mm |
| P222 | 38.0 | 72 mm | 68 mm | 83.2% | 102 mm |
| SLAC | 49.7 | 108 mm | 92 mm | 70.7% | 208 mm |
| **NETWORK** | — | **31 mm** | **31 mm** | **99.8%** | **59 mm** |

Regression guard: `scripts/check_network_benchmark.py` (nine budgets,
fails closed on parse errors; vertical budget would have caught the
round-7 ZWD divergence at product level). Note the guard's own checked
budgets don't include most of this table (only network fused h_p50/
h_RMS/v_RMS, three bases' fixed-only p50, three bases' fix rate) — the
per-base vertical RMS column above is real but currently unguarded;
regressing it wouldn't fail CI today.

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

### 1. State-space wet troposphere inside the DD filter — MOOT, target already met
~~Two ZWD states (rover + active base) with NMF wet mapping differences,
random-walk Q ≈ 3e-7 m²/s each, seeded ~15 cm.~~ Target was P222/SLAC
vertical RMS < 0.3 m. **2026-08-28: re-verified fresh (clean env, rebuilt
release binary, dataset A default path) — P222 v_RMS=103mm, SLAC
v_RMS=210mm, both already 1.4-3x inside the 0.3m target**, matching the
"Refreshed 2026-08-26" table above to within 1-2mm (reproducible, not
noise). Whatever combination of fixes landed since this item was written
(SP3/clock/orbit chain, cadence-hint fix, or others) already solved the
problem this two-ZWD design targeted. Building it now would be solving
an already-solved problem — do not pick this up without a new, current
measurement showing an actual gap. The matrix-state plumbing
(`enable_zwd`) stays dormant behind `widelane_ar`; leave it there.

Full current per-base picture (smoothed, same run): P181 v_RMS=84mm,
**OHLN v_RMS=410mm** (2.7x P225, 4x P222 — see next item), CAPO
v_RMS=101mm, P225 v_RMS=116mm. OHLN is now the only base with a real,
unexplained vertical problem.

### 2. Phase-only network UPD estimation (activates the MW cascade) — BUILT (2026-08-23), MEASURED (2026-08-28)
Raw MW carries per-(station,sat) code multipath — measured. Phase-only
alternative now feasible BECAUSE fusion exists: take the fused network
trajectory as a known rover trajectory, recompute each base's wide-lane
floats against it (geometry fully determined → no code term needed),
then decompose fractional residuals across bases × satellites
(Σu_sat = 0) exactly as PPP-AR networks do. Feed corrected integers to
the dormant cascade (`widelane::resolve_cascade`) and the FAR veto
(`far_matches_widelanes`). Expected: fix-rate headroom at P222/SLAC and
honest fixes through disturbed windows.

**2026-08-28: this design was already implemented in full** (commit
`77e5a27`, 2026-08-23 — "phase-only network UPD estimator - validated
on the CORS day set") and never marked done in this file. Traced the
live code path before assuming otherwise: `mw::solve_network_upd` runs
by default in `eval_network_ppk.rs` as a Phase-A pre-pass over all
bases (opt-out via `WL_NO_UPD`, not opt-in), and
`WidelaneTracker::fixed_widelane` applies the solved per-satellite UPD
(`w -= us - ur`) before its round-to-integer gate, which is exactly
what feeds `resolve_cascade`/`far_matches_widelanes`. The one thing the
Aug 23 commit hadn't done — its own last line says so — was measure the
downstream effect; it validated only the solver's internal residuals
(351 pairs, RMS 0.122 cyc across 6 bases), not fix-rate impact.

Measured today (fresh release binary, sha256 `470a48210284`, dataset A
default full-day run, `WL_NO_UPD=1` vs unset):

| base | fix% off→on (smoothed) | h_p95 off→on | v_p95 off→on |
|---|---|---|---|
| P181 | 86.5 → 86.4 | 57 → 57 mm | 163 → 163 mm |
| OHLN | 79.8 → 80.1 | 185 → 186 mm | 216 → 213 mm |
| CAPO | 92.6 → 92.9 | 90 → 88 mm | 156 → 156 mm |
| P225 | 85.7 → 85.8 | 109 → 109 mm | 174 → 174 mm |
| P222 | 83.0 → 82.8 | 209 → 209 mm | 170 → 170 mm |
| SLAC | 69.9 → 69.9 | 363 → 363 mm | 192 → 192 mm |
| NETWORK | 99.8 → 99.8 | — | — |

Every delta is single-digit epochs out of ~2880 — noise, not signal,
including P222/SLAC specifically. **Verdict: correctly built and
live-wired, but delivers no measurable benefit on this dataset.**
Plausible reason: the solved UPD magnitudes are small (±0.1 cyc per the
original commit) relative to typical arc-mean noise on this benign
mid-latitude day-time dataset, so the correction rarely flips a
rounding decision that wouldn't have gone the same way anyway. Left
on (provably not harmful, may matter on a noisier dataset); not worth
chasing further without one that actually stresses wide-lane rounding
margins.

### 3. Combiner trust model — DONE (architecture, not a new algorithm)
`combiner.rs` blesses co-wrong pass agreement with hard-coded gates
(0.5/0.2/10 m). End-of-day windows show both passes agreeing on wrong
fixes ~14 m apart from truth. Per-base continuity gating (eval-level)
mitigates the product; moving the gate into combiner semantics would
make every consumer honest. Must stay gated to avoid walkthrough drift.

**2026-08-28:** `PostProcessOptions.continuity_gate` now applies
`network::apply_continuity_gate_dynamics` centrally inside
`execute_post_process`, right after the bidirectional combine — any
consumer opts in with one field instead of duplicating the call.
`eval_network_ppk.rs`'s manual post-hoc call (the only place this ever
ran) is deleted in favor of setting the field to its previous `bidir &&
widelane_ar` condition. Every other consumer defaults `false`,
preserving exact current behavior (`eval_qinertia_ppk`'s walkthrough
included) — flipping it on for `gneiss-cli` or the walkthrough is a
separate, deliberate decision this does not make. Both guards
byte-identical, full suite green, clippy clean.

This closes the *architecture* gap (every consumer can now be honest
without re-deriving the call); it does not change what the gate itself
catches. The underlying mechanism is still the same temporal-jump
signal — it still cannot detect a wrong fix that neither disagrees
between passes nor jumps between epochs (P181's invisible single-band
wrong fixes, per "Dataset B tail anatomy" above, are exactly this
class). That remains open, tracked there.

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

  **2026-08-28 re-check** (the "Refreshed 2026-08-26" note above flagged
  this as "worth a fresh look" at 430mm rather than assumed-environmental):
  ruled out receiver antenna PCV/PCO as the cause. OHLN's exact header
  combo (`ASH701945B_M`/`SCIT`) resolves cleanly against igs14.atx
  (confirmed via `RECV-PCV enabled` trace) and applies same as every
  other base — but moves OHLN's vertical RMS by <2mm (416→410mm forward,
  unchanged smoothed), noise-level, unlike CAPO where the same mechanism
  moved 40mm. The signature (p50≈+20mm, essentially unbiased; RMS=410mm,
  huge) is a constant/smooth-pattern corrector doing nothing against
  what must be occasional large excursions — consistent with the
  existing "episodic, environmental" read, now confirmed by elimination
  rather than assumed. Not chased further; would need per-epoch
  timestamp correlation against the excursions to identify a specific
  physical cause (multipath geometry, local wet-delay event), which is
  a different, more expensive investigation than the antenna-gap
  hypothesis this ruled out.

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

**2026-08-27 correction: the measurement above was itself corrupted by
an architecture bug, and the true cost is worse.** `enable_glonass` was
being set independently inside each of the forward and backward passes,
nested inside a >25km-excluding baseline-length gate that existed for
an unrelated feature (ZWD state) — and the backward pass's copy of that
nesting didn't match the forward pass's. Consequence: `GNEISS_GLONASS=1`
silently never reached P222 (38km, over the gate) or SLAC at all, and
could disagree between forward/backward on baselines near the boundary.
"P181/P222 identical" above wasn't evidence GLONASS is safe at long
baselines — P222 simply never received it.

Fixed: `enable_glonass` is now a top-level `PostProcessOptions` field,
applied unconditionally in both passes, independent of baseline length
or `widelane_ar` (`gneiss-rtk` commit series ending in the `--glonass`
CLI flag). Re-measured GE -> GER on the same 2025 DOY160 day with the
fix in place:

| base | fix rate | v_RMS | v_p95 |
|---|---|---|---|
| P181 (15km) | unchanged | unchanged | unchanged |
| P225 (21.9km, under the old gate) | 73.4% -> 70.2% (same direction as before) | **236 -> 796mm** | 334 -> 377mm |
| P222 (38km, previously never received it) | unchanged | unchanged | unchanged |

P225's real cost is far worse than previously documented (RMS more than
triples, not the mild "+92mm p95" the buggy measurement showed) — this
was previously masked by the same bug in a way that happened to look
mild. P222 is confirmed to actually receive GLONASS now (verified via
`GNEISS_GLO_DEBUG=1`: `enable_glonass=true` correctly produces
`selected_constellations=[0,1,2]`, both bases have identical
simultaneous-GLONASS-satellite counts at the epochs checked) yet shows
*zero* measurable change — a genuine, uninvestigated difference in how
P222 responds to the same input, not a lingering wiring bug. Not chased
further this round.

Net effect on the recommendation: unchanged (stay opt-in), but now for
a stronger reason — the true measured cost at P225 is significantly
worse than what justified "opt-in" the first time.

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

## Peer comparison: Gneiss vs published Leica/NovAtel/Qinertia specs (2026-08-29)

Every prior "peer comparison" entry in this file benchmarks against
RTKLIB -- a free, open-source reference implementation this project
has always outperformed. None of them measure against the actual named
tier-1 targets (Leica, NovAtel, Qinertia). This is the first attempt at
that, using published datasheet specs (not measured against gneiss's
own dataset by the vendors -- see caveats below) against gneiss's own
fresh measurement (release binary sha256 `c93d4dbdc831`, dataset A
default, clean env, smoothed/final output):

| base | km | gneiss h_RMS | Leica spec h | ratio | gneiss v_RMS | Leica spec v | ratio |
|---|---|---|---|---|---|---|---|
| P181 | 15.0 | 34 mm | 23 mm | 1.5x | 84 mm | 30 mm | 2.8x |
| CAPO | 16.6 | 48 mm | 25 mm | 1.9x | 101 mm | 32 mm | 3.2x |
| P225 | 21.9 | 63 mm | 30 mm | 2.1x | 116 mm | 37 mm | 3.1x |
| P222 | 38.0 | 120 mm | 46 mm | 2.6x | 103 mm | 53 mm | 1.9x |
| SLAC | 49.7 | 252 mm | 58 mm | 4.4x | 210 mm | 65 mm | 3.2x |
| OHLN | 16.5 | 102 mm | 25 mm | 4.2x | 410 mm | 32 mm | 13.0x |

Leica spec = 8mm + 1ppm(H) / 15mm + 1ppm(V) RMS, single-baseline RTK
(Leica Viva GS12/GS14/GS18T datasheets, consistent across their current
receiver line). NovAtel PwrPak7's published RTK spec is 1cm + 1ppm RMS
(not split into H/V in the fetched spec page) -- roughly matching
Leica's horizontal number at these baselines. Qinertia's published
figure is 4cm H / 8cm V, but that's their PPP mode (no base station);
not a fair comparison against gneiss's base-relative RTK numbers, so
omitted from the table above.

**Reading this honestly:**
- Excluding OHLN (a known, separately-investigated episodic-multipath
  outlier -- see "Known benign anomalies" above), gneiss runs
  **1.5-2.6x worse on horizontal and 1.9-3.2x worse on vertical** than
  Leica's single-baseline spec, across the whole 15-50km baseline
  range this dataset covers. That's a real, moderate, roughly
  consistent gap -- not an order of magnitude, not close to parity.
- The gap does NOT visibly widen with baseline length the way it would
  if this were purely an unmodeled-atmosphere problem (P222 at 38km is
  actually gneiss's BEST vertical ratio, 1.9x) -- consistent with this
  project's own repeated finding that broadcast-orbit/atmosphere error
  mostly cancels in short-to-medium DD baselines, and the remaining
  gap is more about noise floor and edge-case robustness than
  systematic long-baseline physics.
- Real, important caveats this table doesn't capture: Leica's number
  is a datasheet spec (their real-world performance envelope, not
  necessarily their best case) measured on THEIR hardware/antennas
  with precise orbit/clock corrections available in real time; gneiss
  here uses broadcast-only ephemerides (see "Precise ephemeris: SP3
  wiring tested" above -- precise products are measured NOT to help at
  these baselines once orbit error cancels in DD, so this specific gap
  is not primarily an "add SP3" fix). This is a genuinely different
  operating point, not a controlled apples-to-apples trial the way the
  RTKLIB comparisons above are (identical data, identical day).
  Getting a true controlled comparison would need one of these
  vendors' actual hardware logging the same CORS day, which is outside
  what a documentation pass can produce.
- The honest summary: **not tier-1 yet, but the gap is quantified,
  moderate, and now has a number attached to it for the first time**,
  rather than an unmeasured aspiration. Closing it further needs
  exactly what's already tracked above and in
  docs/PROJECT_STATUS.md's Sprint 13 (P181's invisible wrong fixes,
  the atmospheric-decorrelation ceiling at long baselines, receiver
  PCV/antenna modeling depth) -- this table doesn't change what to do
  next, it gives a concrete target to measure progress against.

Sources: [Leica Viva GS14 datasheet](https://techfee.fau.edu/approvedproposals/Download.cfm?sid=444&pid=343), [Leica Viva GS12 datasheet](https://grupoacre.es/wp-content/uploads/sites/3/2020/11/leica_viva_gs12_ds_en.pdf), [NovAtel PwrPak7 performance specs](https://docs.novatel.com/OEM7/Content/Technical_Specs_Receiver/PwrPak7_Performance_Specs.htm), [SBG Systems Qinertia 4 announcement](https://www.gpsworld.com/sbg-systems-unveils-qinertia-4/).

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

**2026-08-29: likely the second one.** Traced this to
`scripts/p224_truth_2025.py`/`gen_multignss_truth.py`: the truth
coordinate is a MONTHLY MEDIAN of UNR's daily IGS20 solutions (June
2025), compared against a single specific observation day (June 9)
within that month. Vertical GPS positions have well-documented 1-3cm
single-day scatter (atmospheric/hydrological loading, daily-solution
noise) that a monthly median smooths out and a single day doesn't —
entirely capable of producing a ~56mm mismatch with no frame problem
involved at all. Definitive test (re-derive truth from June 9 alone,
compare against the monthly median) is specified but not yet run — see
docs/PROJECT_STATUS.md Sprint 15 for the full writeup and why it
wasn't run this session (source `.tenv3` series not persisted, needs
re-fetching first).

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

## Strict disagreement gate validated as critical

A/B test: WL_DISABLE=1 disables widelane_ar which also turns off
strict_disagreement in the forward/backward combiner. Results on
dataset B:

| base | strict h_p50 | non-strict h_p50 | strict h_p95 | non-strict h_p95 |
|---|---|---|---|---|
| P181 | 106 mm | 106 mm | 129 mm | 149 mm |
| P225 | **56 mm** | 164 mm | 214 mm | **741 mm** |
| P222 | 125 mm | 116 mm | 267 mm | 383 mm |
| FUSED | **97.5%** | 66.9% | | |

The apparent P225 "fix loss" in strict mode (83.2% fwd -> 73.4%
smoothed) is not a bug: it is the combiner correctly refusing to
claim fixes when forward and backward filters disagree beyond the
0.50 m static-monument threshold. Non-strict mode accepts these
divergent claims and accuracy collapses.

## Walkthrough reference refresh (incident report)

During Sprint 1.4 verification, eval_qinertia_ppk differed from
/tmp/wt_old.txt by exactly one digit (3D p95 0.926 -> 0.925 mm).
Bisect across ten commits back to f5c6dc6 showed EVERY commit
"differed" — including ones previously verified green.

Forensics: wt_old.txt timestamp (Aug 23 22:31) postdates the newest
branch commit (17:33) by five hours; it was generated from a different
session's tree and is unreachable from this branch's history. Current
binary output is deterministic (two runs byte-identical).

Remediation: reference regenerated from current verified-green tree;
protection remains active for future changes.

## Benchmark-integrity incident + remediation (binary staleness)

DISCOVERY: guards executed target/release binaries without rebuilding;
`cargo test` does NOT build release. With repeated `git checkout <sha> -- .`
cycles, guard results depended on stale-binary lottery. Several A/B
decisions compared binaries from different code states.

HARDENING (permanent):
1. Both guards now force `cargo build --release --bin eval_network_ppk`
   before evaluating, and print sha256[:12] of the binary.
2. All future A/Bs: SAME binary, env-gate toggled only. Paired baselines
   generated fresh per session.

SAME-BINARY REVALIDATION RESULTS (binary 90d90353b904):
- Receiver PCV on dataset A: CAPO v_p50 -54 -> -45 mm (+9 mm) CONFIRMED
  with clean discipline; P181/OHLN exactly zero delta as physics predicts.
- Dataset B PCV: exact zero effect (all-family overlap) — earlier
  "no effect" conclusion CONFIRMED.
- SP3 full-physics chain: no collapse; P181 vp50 -12mm improvement,
  P225 vp50 +13mm regression, P222 unchanged. Mixed/marginal —
  broadcast remains default; SP3 stays opt-in pending clock wiring.
- Elevation 15°: P181 -3.2pp / P225 +3.8pp tradeoff CONFIRMED.

STILL OPEN: P222 fix-rate gap vs budget (80.3 vs 86.0) predates this
audit; requires guarded bisect with forced rebuilds to attribute.

## Sprint 2 clock-datum experiment #1 LANDED: median centering + spread gate

Worktree exp/clk-median-centering, commit bfce486, independently
verified (paired runs reproduced every reported number; guards green
in-worktree and on main post-merge @ 49d44acd6bbf).

Shipped value: GNEISS_CLK is now FAIL-SAFE. Pathological products
(spread > 100 µs after constellation-median centering) suppress the
correction with a latched warning instead of exploding the filter
(control reproduced 35–72 km explosions; gated run healthy).

Empirical findings:
- DOY160 GFZ rapid spread = 1399.6 µs constant all day → correction
  suppressed everywhere on this product.
- SP3 ORBITS alone now improve P181: 99.0% fix, h_p50 101 mm,
  h_p95 125 mm, fused h_p95 −12 mm vs broadcast. Orbit quality helps
  short baselines once TX-time + Sagnac physics are correct.
- OPEN: P222 path never receives dd_clk_m (0 trace hits vs ≥10 for
  P181) — bit-identical across all runs incl. exploding control.
  Root-cause before any clock-dependent tuning on that base.

## Sprint 2 clock-datum experiment #2: NO-SHIP (diagnosis gold, mechanism moot)

Worktree exp/clk-arc-datum (07df0bc). Bookkeeping implemented exactly as
specified — 6 tests, 348 ref-switches handled, guards green — but the
explosion it targeted pre-exists at base commit. Isolation evidence:
base-commit clone explodes identically; bookkeeping-disabled build
explodes identically; CLK-only reproduces bit-identically; SP3-only is
HEALTHY and beats broadcast (P181 99.0% fix / 101 mm — independently
reproduced on their binary).

ROOT CAUSE FOUND (verified on main): dd_clk_m is applied ONLY in
update.rs geom_dd. iono_free.rs has ZERO handling; ar_gate.rs zero in
production. Every fixed position re-estimated from iono-free phase
against a clock-free model → constant tens-of-m displacement.

MYSTERY SOLVED: P222 bit-identical across all runs because precise
products load inside the <25 km ZWD baseline gate (forward.rs:21) —
38 km P222 never loads them at all.

ARCHITECTURE DECISION (next sprint): move clock correction OBS-SIDE at
DD formation (subtract c·Δdt_est from DD code+phase once), delete
model-side term everywhere. Makes every consumer structurally
consistent; makes datum bookkeeping well-defined afterwards; median-
centering logic relocates with it. Until then: SP3-WITHOUT-CLK is the
shippable precise-products configuration.

## GLONASS ICB track: FAILED mid-debug, valuable diagnostic left

Agent implementing glomodear-2-style ICB estimation hit a wall: with
GLONASS enabled at P181, the phase model persistently mismatches by
>500 cycles EVERY epoch (`slip-gate: re-seeded sat=4` repeatedly) —
far beyond any plausible code-bias magnitude. This is NOT an ICB-
scale problem; it indicates a broken FDMA fundamental in our GLONASS
path when it participates in DD:

Suspects ranked (for future investigation):
1. λ per k-channel: glo_freq_num may return 0/wrong slot → nominal-freq
   λ applied to wrong-channel carrier → cycle counts off by k·Δλ/λ.
2. GLONASS time system: PZ-90→GPST offset handling in nav parsing vs
   observation epochs (gnss_time.rs exists; is the DD path using it?).
3. L1 FDMA offset sign convention (FREQ_GLO_L1_DELTA = +562.5 kHz?).

Branch exp/glonass-icb retains partial WIP (uncommitted); main tree
untouched. GLONASS stays gated OFF — no regression risk.

**2026-08-26 re-check: the catastrophic symptom above does not currently
reproduce.** Checked all three ranked suspects and the headline symptom
against current HEAD before touching the ICB branch's WIP:

- Suspect 1 (wrong channel lookup): ruled out. A dedicated diagnostic
  (`crates/gneiss-parsers/examples/glo_k_check.rs`, part of the
  exp/glonass-icb WIP) against the real multi-GNSS nav file shows 23
  GLONASS slots with correct, ICD-range k-numbers (-7..+6) in the right
  antipodal-pair pattern (R01/R05 both k=+1, R02/R06 both k=-4, etc.) —
  `glo_freq_num` is reading real, sane values.
- Suspect 3 (sign convention): ruled out by inspection —
  `FREQ_GLO_L1_DELTA = +0.5625e6`, and higher k correctly produces higher
  frequency, matching the ICD and the independent diagnostic's own
  from-scratch calculation.
- Suspect 2 (PZ-90/GPST time handling): not checked (no fast way to
  verify without deeper instrumentation).
- Headline symptom: `GNEISS_GLONASS=1` on the current multi2025 dataset
  (GER vs GE, same day) produces **zero** `slip-gate: re-seeded` events
  and a small, non-catastrophic effect — P181/P222 unchanged, P225 fix
  rate 73.4% -> 70.2%, network fused 97.5% -> 97.3%. This matches the
  milder, already-documented "FDMA plumbing landed" finding below (code
  ICBs leak via float seeding, phase ambiguities mostly absorb the rest),
  not the catastrophic every-epoch mismatch that motivated this branch.

Read this as "the problem statement is stale," not "the bug is fixed" —
the exp/glonass-icb branch was likely tested on a different day/dataset/
commit where the mismatch was real. Before resuming or porting forward
its well-engineered slope-based ICB model (`glo_icb.rs`: RTKLIB
`glomodear=2` parity, correct cross-wavelength DD phase scaling, 4 green
tests), re-establish that there's still a live problem to solve on
*current* HEAD with the *current* dataset, using the same GNEISS_GLONASS
toggle path above — don't assume the branch's premise still holds.

**2026-08-28: properly ported and measured — NO-SHIP.** Took the above
advice literally: cherry-picked the WIP (`01785a7`) onto current
`network-rtk-long-baseline` (`2391eb1`) in an isolated worktree/branch
(`test-glo-icb-port`), not the stale base it was written against. The
port surfaced three real bugs in the WIP itself, none caught by its own
(otherwise thorough) test suite because none of its tests exercised the
full engine-wiring path or the icb+sat_iono interaction:

1. `retain_active_ambiguities` preserved ICB columns in the wrong order
   relative to sat_iono columns when both were enabled — silently
   corrupts the compacted covariance's column alignment. New regression
   test added (`retain_active_ambiguities_keeps_icb_after_sat_iono_when_
   both_enabled`).
2. Git's auto-merge had inserted a second, incorrectly-chained copy of
   `ensure_icb`/`get_icb_idx`/`icb_offset` from the WIP's pre-sat_iono
   baseline (`E0592` duplicate definitions) — deleted.
3. `state.icb_enabled` was never set `true` anywhere in production code
   (only in tests) — the entire feature was a structural no-op even with
   `enable_glonass_icb` on. Fixed in `configure_iekf`.

With those fixed, 0 warnings / 336+ tests green / clippy clean, and both
regression guards pass unchanged (feature is correctly inert by
default). Controlled measurement, GLONASS participation held constant
across both arms (`GNEISS_SYSTEMS=GRE GNEISS_GLONASS=1`, same release
binary, same day), toggling only `GNEISS_GLONASS_ICB`:

| base | fix% off→on | h_p95 off→on | v_p95 off→on |
|---|---|---|---|
| P181 (15km) | 98.5 → **91.7** | 129 → **138 mm** | 252 → 252 mm |
| P225 (21.9km) | 70.5 → **64.2** | 225 → **257 mm** | 367 → 351 mm |
| P222 (38km) | 88.4 → **86.7** | 267 → **322 mm** | 113 → **135 mm** |
| NETWORK fused | 97.4 → **94.1** | 172 → **178 mm** | — |

Every base gets worse. The P181 result is the most telling: plain
uncalibrated `enable_glonass` costs P181 *nothing* (see table above,
"unchanged" row) — ICB calibration alone drops its fix rate 6.8 points.
The harm is attributable to the ICB slope state itself, not to GLONASS
participation in general.

Root cause (partial — not fully confirmed, flagged for whoever picks
this back up): `ICB_RW_M2_PER_S = 1e-12` effectively freezes the slope
after initial convergence (variance growth over a full day is ~9e-8
m²/MHz² — ruled out as "wandering slope corrupted by noise"). More
likely: `ICB_INIT_VAR_M2_PER_MHZ2 = 4.0` is loose enough that a poorly-
observed *initial* estimate — early epochs, ambiguities still seeding,
few simultaneous GLONASS channels — can lock onto a wrong slope that
then never self-corrects (RW ≈ 0 for the rest of the day), corrupting
the code model at every subsequent epoch. Consistent with the model
being a single global linear slope per FDMA band shared across all
satellite pairs, which may not match the true per-satellite bias
structure closely enough to tolerate a bad lock-in. Not verified via a
slope-trajectory dump this round.

Disposition: `exp/glonass-icb` left untouched (WIP, not merged). The
properly-ported, fully-fixed, fully-measured version is preserved as
`exp/glonass-icb-v2` (`be7ed35`, off `network-rtk-long-baseline` @
`2391eb1`) — not merged, kept as the durable reference for whoever
revisits this. GLONASS ICB moves from "queued" to the same measured-
dead-end bucket as the other NO-SHIP items in this file. Re-attempting
requires redesigning the slope's initialization/observability gating
(e.g. delay engagement until N confident dual-frequency epochs across
≥2 channels, or per-arc reset semantics), not a constant tweak.

## Five-track fan-out consolidated results (all independently verified)

| track | verdict | landed | key evidence |
|---|---|---|---|
| median-centering | SHIP | bfce486 | GNEISS_CLK fail-safe; SP3-orbits beat broadcast at P181 |
| arc-datum | NO-SHIP | branch kept | root cause = partial dd_clk_m wiring; obs-side refactor mandated |
| sidereal | SHIP (diagnostics) | 3846dd7 merged | all channels phase-STRUCTURED p≈0; mitigation inert <2 sweeps by design |
| kinematic | SHIP (behind flag) | 4e40f22 merged | bitwise static parity; sim: 26.2 km divergence -> 0.74 m |
| glonass-icb | INCOMPLETE/FAILED | WIP on branch | >500-cycle FDMA mismatch — deeper than ICB |

Red team: independent byte-compares, own parser, own baseline runs;
caught + got fixed a reporting-semantics bug mid-review (307b7bb).

Combined main @ merge: BOTH guards green, walkthrough bit-identical,
759 tests / 0 failures. Binary be42c55d4a97.

OPEN ITEMS: obs-side clock refactor (subsumes median+arc mechanisms);
GLONASS FDMA fundamentals (>500-cycle mismatch, ≥25km gate trap);
multi-day data for sidereal mitigation activation + storm-day iono
validation.

## Sprint 3 storm-day acquisition + honest first results

Data acquired (network recovered): DOY158+161 adjacent days, G5
geomagnetic storm 2024-05-10 (DOY131), 4 stations each + mixed navs +
truth copies; GNEISS_DATA_DIR override added to eval binary.

FINDINGS (clean 2x2 protocol, fresh binary c7477865+fixes):
1. CRASH FIXED via TDD: iono-state retain desynced dim()/cov when
   satellites set (storm churn) -> nalgebra gemm abort. Unit test pins
   compaction contract; AR-gate raw retain rerouted through cov-
   consistent state method.
2. STORM DAY defeats BROADCAST-only processing completely: zero per-base
   rows (all epochs diverge beyond reporting gates). Engine does not
   crash; it honestly refuses to report garbage. Robustness headline:
   no false fixes emitted under G5.
3. IONO STATES currently HARMFUL even on the quiet day (fused 97.5 ->
   53.3%): a regression vs earlier neutrality, introduced somewhere in
   E5b-fix/obs-side/retain-compaction chain. Default remains OFF;
   isolation is the top open bug — suspected interaction between
   per-pair iono states and Galileo arc composition post-E5b-fix.

NEXT: isolate quiet-day iono regression; then storm-day with working
iono states becomes the decisive capability demo.

## Iono-regression isolation complete: mis-parameterization, not wiring

Discriminator result: GPS-ONLY + IONO=1 collapses identically to GE
(P181 34.6% fix). Eliminates Galileo/E5b interaction -> the per-pair
iono-state PARAMETERIZATION itself is flawed, exactly as the original
design red-team warned ("iono and ambiguity states are correlated
through the phase equation").

Mechanism: independent per-pair constant iono states are rank-deficient
against ambiguities — code noise (±3 m) cannot resolve the mm-level
N/I split, so LAMBDA inherits inflated covariance along the degenerate
direction and fix rates collapse. Earlier "metric-neutral" reading is
now suspected a stale-binary-era artifact (pre-hardening).

CORRECT DESIGN (next sprint): estimate per-SATELLITE slant iono mapped
through a thin-shell/zenith model (RTKLIB ionmapf style) so geometry
couples satellites and breaks the rank deficiency — NOT independent
per-pair constants. Default stays OFF until reimplemented.

## Measured dead end: lowering the iono-free engagement floor (6 -> 3)

Tried the roadmap's stated Sprint 7/9 approach literally: `apply_fixed_
iono_free`'s `h_rows.len() < 6` floor dropped to 3 (matching `if_residual_
outliers`'s own minimum), plus two new boundary tests. Rationale looked
solid going in: debug tracing shows enormous call volume with exactly
3-5 both-band-fixed pairs (tens of thousands of calls across a full-day
run), which the old floor of 6 always sent to `NotEngaged`.

RESULT: zero measurable effect. P222/SLAC smoothed vertical RMS, p50,
p95 and all nine `check_network_benchmark.py` budgets were unchanged to
the millimetre before and after, despite the new code path firing
constantly per the trace.

ROOT CAUSE: `solve_position_lsq` is a correct Bayesian MAP estimate —
`N_plus = H^T R^-1 H + prior_information`, weighted by the filter's own
current position covariance as the prior. At 3-5 pairs, `H^T R^-1 H`
(new-geometry information) is small relative to a filter that's already
converged confidently (if wrongly, due to unmodelled systematic iono
bias — covariance doesn't capture unmodelled bias, only random noise),
so `dx` collapses toward zero and the "solution" reproduces the input
AR-fixed position almost exactly. This is not a bug; the function's own
existing comment already said as much ("the prior pulls the estimate
toward the float position") — it just wasn't clear until measured that
this holds strongly enough to make engagement functionally inert in the
3-5 pair regime specifically, not merely "less confident."

IMPLICATION: "soften the floor" cannot work at ANY threshold under the
current MAP formulation — a floor of 1 would be even more prior-
dominated than 3, not less. The roadmap's second-listed option, "partial-
set IF solving using whichever both-band pairs exist," needs a
genuinely different formulation (one that doesn't re-weight new
information against the filter's own confident-but-potentially-biased
covariance) to have a chance of mattering, not a threshold tweak on the
existing one. Change reverted (`iono_free.rs`, `widelane.rs` back to
floor=6 / original comments); nothing shipped from this round except
this note.

## Bug found and fixed: GNEISS_PCV was never actually applying anything

The "opt-in via GNEISS_PCV=1" elevation-dependent receiver-PCV correction
(`eval_network_ppk.rs`'s `load_receiver_pcv`, documented and referenced
several times in this file) has been dead in practice since it was
written. Root cause: `mod.rs`'s `receiver_dd_pcv_m` gated actually
*applying* the loaded correction behind a SEPARATE env var,
`GNEISS_RECV_PCV`, which nothing in any caller, comment, or this doc ever
told a user to set. `GNEISS_PCV=1` alone (the only thing anyone was ever
told to do) loaded the calibrations successfully — the confirmation
print fires, calibrations resolve — and then applied a correction of
exactly 0.0 to every single observation, silently, for the entire
history of this flag.

Verified with a debug counter (`GNEISS_PCV_DEBUG=1`, prints every
non-trivial correction): `GNEISS_PCV=1` alone -> 0 non-zero corrections
across a full-day six-base run. `GNEISS_PCV=1 GNEISS_RECV_PCV=1` together
-> 482,337, millimetre-scale, physically sensible (sub-mm to ~1mm at
20-46 deg elevations).

This means any prior "GNEISS_PCV metric-neutral / empirically inert"
read anywhere in this project's history was measuring a feature that
could not possibly have done anything, not a genuine physical null
result on this dataset. Fixed by deleting the redundant gate entirely —
`self.receiver_pcv.is_some()` already IS the caller's opt-in signal;
there was never a legitimate case for wanting calibrations loaded but
the correction withheld. `GNEISS_PCV=1` alone now does the whole job, as
every comment already claimed.

Measured real effect now that it's genuinely reachable (dataset A, full
day): small and mixed, not a clear win. Network fused vertical RMS
59->57mm and SLAC v_p95 178->169mm improve; CAPO/P225/SLAC fix rates
each drop 0.2-0.4pp while OHLN's rises 0.3pp. Left opt-in (default
unchanged, both guards verified byte-identical) rather than graduating
to default-on like widelane_ar/GNEISS_RECV_PCO — the effect size here
doesn't clear that bar either way.

## Architecture: run_forward_iekf and run_backward_iekf's setup had drifted in 7 places

Auditing every field the enable_glonass fix (above) touches, side by
side between the two passes, found the same duplicate-setup-drift shape
in six more places -- `GNEISS_ELEV_DEG`, `GNEISS_IONO_STATES`,
`GNEISS_SAT_IONO`, `GNEISS_SP3`, `GNEISS_CLK` were read inside
`run_forward_iekf` only and never reached the backward pass at all
(silently -- no error, the backward `GnssRtkIekf` just kept
`GnssRtkIekf::new()`'s defaults regardless of what the caller asked
for). A seventh, `GNEISS_AMB_DUMP`'s `track_ambiguity_keys`, was gated
by `widelane_ar && baseline < 25km` in forward but just `widelane_ar` in
backward -- the baseline check was never a real requirement for a debug
dump, just accidental proximity to the ZWD-state setup block it happened
to be written next to.

Fixed by extracting `post_process::forward::configure_iekf` -- one
function both passes call to construct and configure their
`GnssRtkIekf`, so the two can no longer independently drift. All six
fixes are behavior-preserving whenever the corresponding env var is
unset (the normal case): reading an absent env var gives the same
"off" result whether it's read once in a shared function or twice in
duplicated ones. Verified via both guard scripts (byte-identical) and a
full run of `eval_qinertia_ppk`'s three-dataset bit-identical walkthrough
(no saved reference existed to diff against -- this project's
walkthrough check is an ad hoc "compare against yesterday's /tmp file"
process, not something checked into the repo -- but every changed line
is gated behind an env var unset in that binary's invocation, so the
walkthrough's own `widelane_ar: false` runs are provably unaffected by
inspection of exactly which lines changed).

Two asymmetries found but deliberately NOT unified:
- `slip_detector`/`base_slip_detector`'s `cadence_hint_s`: forward only
  sets it when `widelane_ar` is on; backward sets it unconditionally.
  Confirmed NOT dead (screening.rs's gap-threshold computation reads it),
  so unifying this would be a real behavior change for anyone running
  `widelane_ar: false` with backward smoothing -- exactly
  `eval_qinertia_ppk.rs`'s own configuration. That walkthrough's output
  would change, which needs a deliberate round (reference regeneration,
  explicit sign-off) matching this project's own precedent ("Walkthrough
  reference refresh (incident report)"), not a side effect of a
  refactor aimed at something else.
- `wl_tracker.sat_upd`: forward clears it to `None` when `widelane_ar`
  is off; backward sets it regardless of `widelane_ar`. Likely inert
  (network_sat_upd is documented as "consumed... when widelane_ar is
  on", and the only reader found, `WidelaneTracker::fixed_widelane`, is
  only reachable from widelane_ar-gated call sites) but not proven
  inert to the same standard as the six fixed above, so left alone
  pending the same deliberate-round treatment as cadence_hint_s.

Recommend addressing both together in one round, specifically because
fixing them requires regenerating the walkthrough reference either way
-- no reason to pay that cost twice.

**2026-08-27: both resolved** -- see "The two deferred asymmetries:
both resolved, both were real bugs" below.

## Remaining Pattern-1 env vars checked: no drift bug, left as-is

The three env-var toggles not covered by the `configure_iekf` sweep --
`GNEISS_STRICT_VETO` (widelane.rs), `GNEISS_AR_GATE` and `GNEISS_IF_VETO`
(rtk_iekf/mod.rs) -- were checked against the same forward/backward-drift
question. All three read `std::env::var` from code paths already shared
between passes rather than per-pass setup:
- `GNEISS_AR_GATE` is read once, inside `GnssRtkIekf::new()` itself.
  Both passes construct through the same `configure_iekf` -> `new()`
  chain, so there is only one read site, not two that could disagree.
- `GNEISS_STRICT_VETO` and `GNEISS_IF_VETO` are read per-call inside
  `far_matches_widelanes` and the post-fix IF-residual screen, both
  reached only through `process_epoch`, which forward and backward
  share verbatim (same method, same struct).

No fix needed -- the bug class that motivated the `configure_iekf`
extraction requires two independent setup call sites that can drift;
these have exactly one. They're still "config not in the config"
(invisible from `PostProcessOptions`, no CLI flag) but promoting
explicitly-experimental, off-by-default AR-quality research knobs
(`GNEISS_IF_VETO`'s own comment: "Experimental: post-fix iono-free
residual screen") to first-class surfaced config isn't warranted
without a concrete caller who needs to set them -- that would be
speculative surface area, not a bug fix.

## The two deferred asymmetries: both resolved, both were real bugs

Revisited `cadence_hint_s` and `wl_tracker.sat_upd` (deferred above)
now that the deliberate round they needed is in scope. Both turned out
to have a clear, provable direction rather than an arbitrary pick:

- `cadence_hint_s` reflects the *input data stream's* own sampling
  rate (screening.rs's gap-detection threshold, used so a normal
  30-second inter-epoch step on a slow-cadence stream isn't flagged as
  a cycle-slip-causing data gap). It has no principled dependence on
  `widelane_ar` at all -- forward's gating on it was the actual bug,
  not a deliberate restriction that happened to differ from backward.
- `wl_tracker.sat_upd` is now *proven* inert whenever `widelane_ar` is
  off, not just "likely": its only production reader,
  `WidelaneTracker::fixed_widelane`, is reachable solely through
  `far_matches_widelanes` and `resolve_cascade`, both called only from
  `if self.widelane_ar` in `process_epoch` (the third call site, in
  widelane.rs, is a `#[cfg(test)]` fixture). Setting it when
  `widelane_ar` is false cannot change any output.

Both are now folded into `configure_iekf`, unconditionally, matching
backward's already-correct behaviour. Measured impact:

- Both guard scripts: `check_network_benchmark.py` (9/9 checks) and
  `check_multignss_benchmark.py` (10/10 checks) pass with no regression.
- `eval_qinertia_ppk`'s three-dataset walkthrough, before vs. after:
  - RTK Explorer F9P (1 Hz): byte-identical, both forward and smoothed.
    Expected -- `infer_cadence_hint` returns `None` below its 2 s
    median-spacing floor, so fast streams were never affected.
  - NGS Geodetic Baseline (30 s cadence) forward pass: fix rate
    84.0% -> 97.0% (252/300 -> 291/300), horizontal p95 679mm -> 125mm,
    RMS 274mm -> 55mm, 3D p95 1203mm -> 289mm. This is the bug: the
    legacy fixed 2 s gap threshold flagged every normal 30 s step as a
    slip, constantly resetting ambiguity tracking and starving AR of
    the epoch count it needs to converge.
  - NGS Geodetic Baseline smoothed (final) output: stayed 100% fixed
    (300/300) both before and after; p50 4mm -> 7mm, p95 14mm -> 21mm,
    RMS 7mm -> 11mm. A few-mm wobble at the margin, plausibly from the
    smoother's forward/backward blend weights shifting now that forward
    contributes far more (and better) fixed epochs -- not a regression
    in any practical sense at these absolute magnitudes, and the guard
    scripts (which gate on the metrics that matter) show no issue.

Net: forward-pass accuracy on slow-cadence streams improved
substantially; nothing else moved outside noise. This walkthrough's
new numbers are the reference from this point forward -- there's no
checked-in file to update (see "Walkthrough reference refresh" above),
just this log entry.

## Code hygiene pass: 4 duplicate functions fixed, larger debt catalogued

Following up on the architecture work above with a systematic sweep for
CLAUDE.md's other hard rules (file size, duplicate functions). Found
and fixed four genuine duplicate-function cases (all mathematically
verified identical before touching anything, all confirmed
byte-identical on both guard scripts after):
- `compute_enu_stds` (combiner.rs, smoother.rs) -> `gneiss_core::coords::ecef_cov_to_enu_std`
- `track_c_freq` (mw.rs, iono_free.rs) -> `gneiss_core::frequencies::track_c_frequency`
- `ecef_to_enu` (eval_odaiba_kf/ins.rs, eval_swfg.rs) -> `gneiss_core::coords::ecef_delta_to_enu`
- `horizontal_error`/`vertical_error` (eval_network_ppk.rs, reimplemented
  despite `gneiss_core::metrics` already exporting canonical versions)

The last one was the highest-stakes check: `eval_network_ppk` is the
binary both guard scripts run, so these two functions produced every
p50/p95/RMS number this session's whole verify-before-commit discipline
was trusting. No drift was found, but it's exactly the kind of thing
that's worth checking rather than assuming.

Also split `forward.rs` (545 lines, over the 500-line file limit after
the `configure_iekf` work) into `iekf_pass.rs`, fixing a layering smell
where `backward.rs` reached into `forward.rs` for shared setup.

Larger, correctly-deferred debt found along the way, not fixed here:
- **19 other files over the 500-line limit**, several massively so
  (rinex.rs 2349, spp.rs 2228, rtk_iekf/mod.rs 1743, ephemeris.rs 1573
  lines). All pre-existing, none touched this session. Splitting these
  is real, substantial architecture work on core validated engine
  files -- not a "quick" fix, and risky to do without dedicated focus.
- **40+ unconditional or ambiguously-gated `println!`/`eprintln!` calls**
  in library code (not eval binaries). Some are legitimate opt-in debug
  tooling already gated behind env vars (this session added a few of
  those deliberately); others look like leftover debugging scaffolding.
  Distinguishing the two requires reading each site in context -- too
  large and too judgment-heavy for a quick pass. rtk_iekf/mod.rs alone
  has a dozen-plus (`BAD-SEED`, `SP3-PROBE`, `CONTENT repr`,
  `ENGINE-TEST` labeled prints suggest one-off debugging sessions left
  in place).
- **A `percentile` function with three different signatures**
  (quality.rs: presorted+int-pct; sidereal/mod.rs: presorted+float-q;
  eval_swfg.rs: owned-unsorted+float-pct). Confirmed these do NOT feed
  either guard script (`eval_network_ppk.rs` computes its own p50/p95
  inline via direct sorted-array indexing, using neither this nor the
  canonical `gneiss_core::metrics::compute_statistics`). Lower stakes
  than the four fixed above, and the signature differences mean
  consolidating needs real per-call-site verification, not a
  find-and-replace -- deferred rather than rushed.
- **Nesting depth**: a rough brace-counting scan puts the worst offenders
  in rinex.rs, ionex.rs, antex.rs, hatch.rs -- the same large parser
  files already flagged for file size above. Likely one underlying
  "these parsers need a dedicated pass" finding rather than several
  independent quick fixes; the brace-counting heuristic also overstates
  true control-flow depth (it counts struct/impl/fn scoping too), so a
  real fix would need a proper per-function read first.
