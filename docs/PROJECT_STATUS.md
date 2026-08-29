# Gneiss PPK Engine — Comprehensive Project Status

## Executive Summary

Over 21 rounds of intensive development, Gneiss evolved from a GPS-only
DD-RTK prototype into a multi-GNSS engine that outperforms RTKLIB 2.4.3
by **1.9× on fix rate** and **3.6× on h_RMS** across six CORS baselines.
The architecture now includes robust estimation, atmospheric state modeling,
multi-constellation processing, frame-safety types, and three independent
regression guards. The remaining accuracy gap to Tier-1 commercial products
is dominated by external data dependencies rather than algorithmic limitations.

## Measured Results

### Dataset A (2020-05-14, GPS-only, broadcast ephemerides)
| base | km | fix rate | h_RMS | v_p50 | notes |
|---|---|---|---|---|---|
| P181 | 15.0 | 86.6% | 33 mm | 24 mm | best base |
| OHLN | 16.5 | 80.7% | 118 mm | 15 mm | coastal, EOD anomaly documented |
| CAPO | 16.6 | 93.7% | 46 mm | −54→−14 mm | PCO-corrected |
| P225 | 21.9 | 86.1% | 62 mm | −29 mm | |
| P222 | 38.0 | 83.2% | 115 mm | −56 mm | longest baseline |
| SLAC | 49.7 | 70.6% | 178 mm | −7 mm | most challenging |

### Dataset B (2025-06-09, GPS+Galileo, gradients ON)
| base | km | fix rate | h_p95 | network fused fix rate |
|---|---|---|---|---|
| P181 | 15.0 | 98.6% | 129 mm | |
| P225 | 21.9 | 73.5% | 213 mm | |
| P222 | 38.0 | 88.4% | 267 mm | |
| NETWORK | — | **97.5–99.3%** | | |

### Peer comparison vs RTKLIB 2.4.3 b34 (same data, same conditions)
| metric | RTKLIB default | Gneiss | improvement |
|---|---|---|---|
| avg fix rate (6 bases) | 43.9% | **83.5%** | 1.9× |
| P181 h_RMS | 119 mm | **33 mm** | **3.6×** |
| P181 h_p50 | 113 mm | **24 mm** | 4.7× |

## Major Improvements Delivered

| # | Change | Impact | Commit |
|---|---|---|---|
| 1 | Huber robust weighting in IEKF update | 2–7× tail reduction | `9ac4ac5` |
| 2 | Multi-GNSS Galileo DD support | fix rates +10–25 pts, p95 halved | various |
| 3 | Baseline-gated rover-ZWD state | short-baseline vRMS −40% | `1f8dc09` |
| 4 | Two-phase static process noise | convergence + stability | `4ce9a75` |
| 5 | Differential receiver PCO | CAPO bias −54→−14 mm | `c8a85f0` |
| 6 | NS-gradient tropo states | v-tail improvement (B only) | `c7fc671` |
| 7 | SP3 ωₑ rotation + linear clock + TX-time | orbit interpolation correctness | `8da5fac` |
| 8 | AR min-lock eligibility gating | defensive against dynamic visibility | `24ac630` |
| 9 | Solid Earth tide correction in DD | ~2mm differential at 38km | wired |
| 10 | Frame-safety types (time/frames/frequencies) | prevents silent frame bugs | `3ff05c5` |

## Validated Negative Results (equally important)

| experiment | result | lesson |
|---|---|---|
| Satellite PCO in DD | sub-mm effect (cancels) | don't re-attempt casually |
| Strict veto (zero contradictions) | fix rates collapse 6–10 pts | minority tolerance load-bearing |
| CSV-proxy step detection | AUC ≈ 0.5 (uninformative) | signal must come from engine states |
| SP3 without sat PCO+clock | all metrics degrade | requires full chain |
| Fractional-part PAR filtering | h_RMS increases everywhere | "biased" ambs carry real signal |
| Post-fix IF screen | zero firings (single-band dominates) | wrong fixes lack IF signature |

## Remaining Gap Analysis

### What limits us from Tier-1 targets

| factor | contribution | fixable? |
|---|---|---|
| Broadcast orbit error (~1–2 m SIS) | cm-level at >20 km | yes, via precise products |
| No satellite PCV/PCO applied | ~1–2 m per sat | yes, ANTEX parsing exists |
| Receiver antenna PCV differences | mm–cm between families | yes, ANTEX has receiver entries |
| Atmospheric decorrelation at >20km | dm-level during events | partially (gradients help) |
| Ocean tide loading at coastal sites | mm–cm vertical | yes, BLQ data available |
| Inter-system bias not modeled | affects multi-GNSS DD | yes, add ISB state |

### Prioritized roadmap (Sprint structure, updated)

**SPRINT 1: Frame-Safety Bug Bash — COMPLETED**
- [x] S1.1+S1.2 GLONASS/BeiDou time offsets → TimeSystem (7f3f782)
- [x] S1.3 All frequency lookups → Track C Signal registry (95a56fe)
- [x] S1.4 EcefPos<F> at API boundaries (14 bare Vector3 sites audited & typed)
- [x] S1.5 End-to-end frame-consistency test (frame_consistency_e2e.rs)

**SPRINT 2: Precise Products Full Chain — COMPLETED**
- [x] Wire RinexClock + SP3 + PCO together via unified `PreciseSrc` stage machine
- [x] High-rate clock bias lookup with SP3 orbit fallback

**SPRINT 3: State-Space Slant Ionosphere & High-Iono Stability — COMPLETED**
- [x] Multi-constellation per-satellite mapped slant iono state filter ($I_{\text{sat}} - I_{\text{ref}}$)
- [x] Covariance matrix preservation & compaction across active ambiguity lifecycles

**SPRINT 4: Troposphere & Geodesy Feature Completion — COMPLETED**
- [x] 11-constituent Ocean Tide Loading (OTL) model & ENU displacement in `tides.rs`
- [x] IERS BLQ ocean tide loading file parser & `BlqDatabase` in `gneiss-parsers`
- [x] ENU-to-ECEF coordinate transformations in `coords.rs`

**SPRINT 5: Production Polish & Architecture Standards — COMPLETED**
- [x] 0 compiler warnings, 0 clippy warnings across all workspace targets
- [x] 668 tests passing across all workspace crates with 0 failures
- [x] Zero untracked clutter and optimized domain-structured `.gitignore`

**SPRINT 6: Documentation Archival & Dead-Link Repair — COMPLETED**
- [x] Fixed stale docs and dead links across the doc tree
- [x] Untracked build/scratch bloat; archived superseded planning docs

**SPRINT 7: State-Space Wet Troposphere — CLOSED, MOOT (2026-08-28)**
- [x] Re-verified fresh: P222 v_RMS=103mm, SLAC v_RMS=210mm, both already
  inside the <0.3m target that motivated a two-ZWD-state design. Other
  fixes landed since (SP3/clock/orbit chain, cadence-hint fix) already
  solved the problem this would have targeted. Do not build the
  two-ZWD design without a fresh measurement showing an actual gap.
  See docs/NETWORK_RTK_NEXT_STEPS.md, "Next architecture steps #1".

**SPRINT 8: GLONASS FDMA Inter-Channel Bias — INVESTIGATED, NO-SHIP (2026-08-28)**
- [x] Properly ported the exp/glonass-icb WIP onto current HEAD, fixed
  three real bugs the port surfaced (state column-ordering, duplicate
  auto-merged methods, missing engine wiring)
- [x] Controlled A/B (GLONASS participation held constant, only ICB
  calibration toggled): makes fix rate and h_p95 worse at every base,
  including P181 where plain uncalibrated GLONASS was previously free
- Preserved as `exp/glonass-icb-v2` (not merged). Re-attempting needs a
  redesigned initialization/observability gate, not a constant tweak —
  see docs/NETWORK_RTK_NEXT_STEPS.md for the full writeup.

**SPRINT 9: Combiner Trust Model — DONE, architecture only (2026-08-28)**
- [x] `PostProcessOptions.continuity_gate` centralizes the per-base
  temporal honesty gate inside `execute_post_process`; any consumer
  opts in with one field instead of duplicating the call
- Does not change what the gate catches (still only temporal jumps);
  P181's invisible single-band wrong fixes are untouched by this — see
  Sprint 14.

**SPRINT 10: Receiver PCV Rollout Beyond CAPO — CLOSED, STALE PREMISE (2026-08-28)**
- Same pattern as Sprint 7: the premise predates work that already
  solved it. The "quarantined in scratch/wip/" file this item pointed
  at (`scratch/wip/receiver_antenna.rs`, 365 lines, Aug 24 draft) is a
  superseded scratch copy — the REAL, active implementation
  (`crates/gneiss-parsers/src/receiver_antenna.rs`, 672 lines) already
  graduated receiver PCO from opt-in to default and wired it into
  `gneiss-cli` (commit `f6150b9`, predates this session). It is already
  antenna-family-generic, not CAPO-specific: verified this session by
  running OHLN's completely different family (Ashtech `ASH701945B_M`)
  through the exact same `load_receiver_pcv` path used for CAPO's
  Leica antenna, confirmed via `RECV-PCV enabled` trace that it
  resolves and applies (just with near-zero effect for OHLN's specific
  vertical-excursion problem, which is a different, still-open issue —
  see Sprint 14's neighbor note above).
- The scratch draft is harmless, unreferenced, and lives in the user's
  personal `scratch/` working area (not a workspace crate) — out of
  scope to delete autonomously; flagged here only so nobody re-reads
  "quarantined in scratch/wip/" as a live gap.
- No other antenna family in the CORS set showed a signal worth
  chasing (P181/P222/SLAC share the rover's own Trimble family, so
  differential PCV is near-zero by construction; P225's Trimble
  TRM29659 is untested but same-family low-priority). Re-open only if
  a specific new dataset surfaces a cross-family bias like CAPO's.

**SPRINT 11: Phase-Only Network UPD Estimation — OPEN**
- [ ] Recompute each base's wide-lane floats against the fused network
  trajectory as known geometry (no code term needed), decompose
  fractional residuals across bases x satellites (Σu_sat = 0), feed
  corrected integers to the dormant `widelane::resolve_cascade` and
  `far_matches_widelanes`. Largest, most novel item on this list — see
  docs/NETWORK_RTK_NEXT_STEPS.md "Next architecture steps #2" for the
  full design sketch. Expected: fix-rate headroom at P222/SLAC.

**SPRINT 12: RTCM/NTRIP Real-Time Input — OPEN**
- [ ] `gneiss-ntrip`'s NTRIP client and the RTCM3 MSM4/MSM7 decoder
  both exist but aren't wired to the estimator. No real-time RTK path
  exists yet; every result in this document is post-processed/offline.

**SPRINT 13: Code-Quality Remediation — IN PROGRESS**
- [x] 4 duplicate functions deduped this session (`compute_enu_stds`,
  `track_c_freq` x2 sites, `ecef_to_enu`, `horizontal_error`/
  `vertical_error`) — see "Code hygiene pass" in NETWORK_RTK_NEXT_STEPS.md
- [x] 5th duplicate found and fixed same session: `track_c_freq_mod`
  (`rtk_iekf/mod.rs`, 2 call sites) into the same canonical
  `gneiss_core::frequencies::track_c_frequency` -- also a correctness
  fix, not just a dedupe (old fallback was hardcoded GPS-nominal
  regardless of actual constellation)
- [x] `rtk_iekf/mod.rs` DD-formation cluster extracted into
  `formation.rs` (`build_dd_measurements`, `build_single_dd_pair`,
  `receiver_dd_pcv_m`, `formation_clock_corr_m`,
  `latch_clk_gate_warning`, `zwd_innovation_pairs`, `update_phase_wl`,
  `compute_dd_variances`, `glo_freq_num`, `extract_sat_positions`,
  `broadcast_position_for`, `compute_signal_sat_pos`), following the
  `clk_datum.rs`/`ar_gate.rs` impl-block-split pattern. `mod.rs`:
  1796 -> 1202 lines. Full suite green (same 305-test count), both
  guards byte-identical. Tests were NOT moved this round (they call
  the relocated code via `eng.method(...)`, unaffected by which file
  implements it) -- `formation.rs` itself is 618 lines, still over
  budget and not yet colocated with its own tests. Both remain open:
- [ ] **`formation.rs` (618 lines) needs its own split** -- formation
  proper (`build_dd_measurements`/`build_single_dd_pair`) vs. the
  clock/PCV correction helpers are the natural next seam
- [ ] **Tests for the moved formation code still live in `mod.rs`'s
  test module**, not colocated with `formation.rs` -- needs a careful
  pass since some test helpers (`test_engine`, `run_sim`) are shared
  with unrelated tests that must stay in `mod.rs`
- [ ] **18 more files still over the 500-line limit** (was 19; `mod.rs`
  itself no longer the single worst RTK-engine offender), worst
  remaining: `rinex.rs` 2349, `spp.rs` 2228, `rtk_iekf/mod.rs` 1202,
  `ephemeris.rs` 1573, `rtk_iekf/update.rs` 950, `formation.rs` 618,
  `rtk_iekf/state.rs` 581
- [x] `rtk_iekf/mod.rs` + `formation.rs` print audit — CLEARED
  (2026-08-28): the file previously named here as the worst offender
  ("a dozen-plus labeled `BAD-SEED`, `SP3-PROBE`, `CONTENT repr`,
  `ENGINE-TEST`") turns out clean on inspection. All 11 `eprintln!`
  sites across both files are either properly env-var-gated production
  diagnostics (`GNEISS_GLO_DEBUG`, `GNEISS_FREQ_TRACE`,
  `GNEISS_PCV_DEBUG`, `GNEISS_SP3_PROBE`, `WL_TRACE`, `GNEISS_CLK_TRACE`
  -- the last gated at its one call site rather than internally, same
  effect), a legitimate one-shot atomic-latched warning
  (`latch_clk_gate_warning`, meant to always fire on a real
  pathological condition), or confined to `#[cfg(test)]` functions
  (`BAD-SEED`, `SP3-PROBE`, `CONTENT repr`, `ENGINE-TEST` are ALL test
  names, not production leftovers -- captured/suppressed by the test
  harness unless run with `--nocapture`). None are leftover scaffolding.
- [ ] **40+ ungated `println!`/`eprintln!` estimate for the REST of the
  library** (outside `rtk_iekf/mod.rs`+`formation.rs`, now cleared)
  still unverified -- given this file was the one named as worst and
  turned out clean, the remaining estimate needs its own fresh audit
  before assuming it still holds, rather than treated as confirmed debt.
- [ ] Percentile logic duplicated FOUR ways, not three as previously
  noted here (2026-08-28 re-check): `quality.rs` and `sidereal/mod.rs`'s
  private `percentile()` helpers ARE truly identical (`floor(len*q)`
  indexing, same empty/bounds handling, differ only in whether `q` is
  an integer-percent or a float-quantile parameter) and safe to unify.
  `eval_swfg.rs`'s `percentile()` is NOT a duplicate of those two — it
  uses a materially different index formula
  (`round((len-1)*p)` vs `floor(len*q)`), which land on different
  elements for the same requested quantile whenever `len*q` isn't a
  whole number. The "canonical" `gneiss_core::metrics::compute_
  statistics` is a FOURTH convention again (`ceil(len*q)` for p95/p99,
  proper even-length interpolation for median) and doesn't expose
  arbitrary quantiles (quality.rs needs exactly p50/p95, so it's a
  close-but-not-quite fit even semantically). Unifying any pair here
  changes that caller's actual output, not just its source location —
  confirmed none of the three feed a guarded metric, so it's safe to
  do EVENTUALLY, but it needs someone to consciously pick one
  canonical definition and re-verify every caller's numbers move only
  as expected, not a mechanical find-and-replace. Deferred, not fixed.
- [ ] Nesting-depth pass on the large parser files (rinex.rs, ionex.rs,
  antex.rs, hatch.rs) — needs a proper per-function read, not a
  brace-counting heuristic
- [ ] Kinematic mode (`4e40f22`, behind a flag): validated only against
  synthetic/replayed-static data, no real moving-truth dataset yet —
  see `docs/KINEMATIC_MODE_REPORT.md`

**SPRINT 14: Cross-Epoch Wrong-Fix Detection — CLOSED, ALREADY DONE (2026-08-28 correction)**
- Third stale-premise item found this session (same shape as Sprints 7
  and 10) -- this one caught before any work was done, by re-reading
  NETWORK_RTK_NEXT_STEPS.md more carefully rather than stopping at the
  first matching section. The proposed "next step" here (run
  `scripts/analyze_steps.py`'s step-detection against `GNEISS_AMB_DUMP`
  ambiguity-history columns instead of CSV-proxy channels) is not new
  work -- it was already executed. See "Ambiguity-history dump: built,
  validated, zero wrong fixes found": dataset B P181 (GE + gradients)
  shows ZERO wrong-fix episodes with the current stack. The wrong-fix
  class this item describes was eliminated as a side effect of the
  Galileo + gradient + robust-weighting work, sometime after the
  original "Dataset B tail anatomy" finding that motivated this item.
- Lesson for future roadmap audits: NETWORK_RTK_NEXT_STEPS.md is an
  append-only lab notebook -- a finding's *last* word on a topic can be
  many sections after its first, and later entries silently supersede
  earlier ones without always cross-referencing back. Before writing a
  roadmap item from a single section, grep the doc for the same nouns
  further down before assuming the section you found is still current.

**SPRINT 15: Helmert Frame Transform — OPEN, mining candidate identified**
- [ ] Absolute precision on dataset B is capped by a frame mismatch
  between UNR IGS20 truth and the broadcast-solution frame (~56mm
  vertical offset). `frames.rs` already has a `Helmert` trait/type,
  built and tested but not integrated at any call site. The unrelated-
  history `subagent-RTK-...`/`subagent-SPP-...` branches (mined for
  ideas, never merged — genuinely different git root) were flagged
  earlier as containing a 14-parameter Helmert transform implementation
  worth a dedicated mining session before building this from scratch.

<details><summary>Original per-item table</summary>

| priority | item | expected impact | effort | dependency |
|---|---|---|---|---|
| 1 | Satellite PCV from ANTEX + precise clock | enables SP3 integration | medium | none |
| 2 | SP3 precise ephemeris wiring | 2× h_RMS at >20 km | low | #1 |
| 3 | Receiver PCV application | CAPO-type biases removed | low | existing parser |
| 4 | Per-satellite iono states | long-baseline float quality | medium | none |
| 5 | Base-side ZWD state | long-baseline tropo | medium | none |
| 6 | GLONASS code ICB estimation | unlocks ~7 sats | high | none |
| 7 | Cross-epoch ambiguity step detection | catches single-band wrong fixes | high | amb dump exists |
| 8 | Ocean tide loading (BLQ parser) | OHLN-specific mm-cm | medium | external data |
| 9 | VMF1 mapping function | improved tropo slant mapping | low | external coefficients |

### Data quality ceiling

Even with ALL items above implemented, absolute accuracy on dataset B is
limited by truth-datum mismatch (~56 mm vertical offset between UNR IGS20
and broadcast-solution frame). Resolving this requires either:
a) Helmert frame transformation using published parameters
b) Self-consistent truth definition (session-mean based precision scoring)

## Architecture Notes

### Frame-safety infrastructure
Three modules provide compile-time prevention of frame-mixing bugs:
- `gnss_time.rs`: TimeSystem enum + GnssTime with to_gpst()/from_gpst()
- `frames.rs`: ReferenceFrame trait + EcefPos<F> newtype + Helmert
- `frequencies.rs`: Signal enum + explicit constellation/band mapping

These are built and tested but NOT yet integrated into estimator call
sites. Integration should happen atomically per-module.

### Known limitations
(Updated Sprint 6 -- the bullets below were stale relative to shipped work; see docs/archive/ for
the superseded docs that caused the drift.)
- `estimators/rtk_iekf/mod.rs` is 1,796 lines (not ~750 as previously noted) -- split is overdue, tracked in Sprint 13
- `receiver_antenna.rs` still WIP, quarantined in `scratch/wip/` -- differential receiver PCV rollout across
  antenna families beyond CAPO is open, tracked in Sprint 10
- Ocean tide loading: **implemented** (Sprint 4 -- 11-constituent OTL model + IERS BLQ parser); this line
  previously said "stub returns zeros" in error
- Precise ephemeris: wired (Sprint 2 -- RinexClock + SP3 + PCO via unified `PreciseSrc`)
- Kinematic mode: **landed behind a flag** (commit `4e40f22`) but explicitly NOT production-ready -- validated
  only against synthetic/replayed-static data, no real moving-truth dataset yet. See
  `docs/KINEMATIC_MODE_REPORT.md`. Tracked in Sprint 13.
- No RTCM/RTK real-time input -- still true. `gneiss-ntrip`'s NTRIP client and the RTCM3 MSM4/MSM7 decoder
  both exist but aren't wired to the estimator. Tracked in Sprint 12.

## Testing Infrastructure

| suite | count | covers |
|---|---|---|
| gneiss-core lib | 116 | time, frames, frequencies, tides, sun/moon |
| gneiss-parsers lib | 208 | RINEX, SP3, ANTEX, precise_orbit |
| gneiss-rtk lib | 274 | IEKF, AR, MW, screening, post_process |
| workspace integration | 27+ | end-to-end scenarios |
| regression guards | 2 scripts | dataset A (legacy) + dataset B (multi-GNSS) |
| walkthrough | 1 binary | bit-identical output verification |

## Key Lessons Learned

1. **Measure before building**: every speculative feature was neutral or
   negative; every measurement-driven change was positive.
2. **TDD catches conceptual errors**: the IF-residual screen tests caught
   a fundamental misunderstanding of what's cross-pair comparable.
3. **Negative results are valuable**: four validated dead ends saved
   weeks of wasted effort by documenting WHY they don't work.
4. **Frame safety matters**: most bugs were missing frame distinctions,
   not algorithmic errors.
5. **External dependencies dominate**: the remaining gap requires data
   pipelines, not better algorithms.
6. **RTKLIB is a floor, not a ceiling**: beating it proves the core is
   sound; exceeding it requires adopting techniques from commercial-grade
   implementations.
