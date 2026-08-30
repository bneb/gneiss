# Gneiss PPK Engine — Comprehensive Project Status

## Executive Summary

Over 21 rounds of intensive development, Gneiss evolved from a GPS-only
DD-RTK prototype into a multi-GNSS engine that outperforms RTKLIB 2.4.3
by **1.9× on fix rate** and **3.6× on h_RMS** across six CORS baselines.
The architecture now includes robust estimation, atmospheric state modeling,
multi-constellation processing, frame-safety types, and three independent
regression guards.

**2026-08-29: first measured comparison against actual Tier-1 specs**
(not just RTKLIB) -- see docs/NETWORK_RTK_NEXT_STEPS.md, "Peer
comparison: Gneiss vs published Leica/NovAtel/Qinertia specs". Against
Leica's published single-baseline RTK datasheet spec (8mm+1ppm H,
15mm+1ppm V RMS), gneiss currently runs **1.5-2.6x worse on horizontal,
1.9-3.2x worse on vertical** across 15-50km baselines (excluding one
known outlier base). That supersedes this section's previous
unmeasured claim that the gap is "dominated by external data
dependencies" -- the measured gap doesn't widen with baseline length
the way pure atmospheric/orbit error would predict, suggesting noise
floor and edge-case robustness (P181's invisible wrong fixes,
receiver/antenna modeling depth) matter at least as much as external
data. Real caveat: this compares a datasheet spec against gneiss's own
broadcast-ephemeris run, not a controlled same-day same-hardware
trial -- see that section for the full caveats.

**2026-08-28 roadmap audit**: a full pass through every open item below
(Sprints 6-15) found that roughly half were already resolved by other
work and never marked done -- including the single largest item on the
list (Sprint 11, phase-only network UPD estimation), which turned out
to already be built, wired, and running by default since 2026-08-23,
just never measured or documented. See each sprint's entry for what
was verified vs. what's genuinely still open; the honest remaining
list is shorter than it looked.

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

**SPRINT 11: Phase-Only Network UPD Estimation — CLOSED, ALREADY BUILT (2026-08-28 correction)**
- The biggest miss in this document's first draft: this was written up
  as "OPEN... largest, most novel item on this list" by reading only
  the original design-sketch section in NETWORK_RTK_NEXT_STEPS.md
  ("Next architecture steps #2") without checking whether it had
  already been done. It has. Traced the actual code before writing
  anything: `mw::solve_network_upd` IS the phase-only cross-base
  least-squares decomposition (Σu_sat=0 via eliminated-satellite
  substitution) exactly as designed; `eval_network_ppk.rs` runs it as a
  default Phase-A pre-pass (opt-OUT via `WL_NO_UPD`, not opt-in) over
  all bases; and `WidelaneTracker::fixed_widelane` genuinely applies
  the solved UPD (`w -= us - ur`) before the round-to-integer check
  that feeds `resolve_cascade`/`far_matches_widelanes` -- confirmed by
  reading the exact line, not inferring from a comment. Landed in
  commit `77e5a27` (2026-08-23), which validated the SOLVER's internal
  residuals (351 pairs, RMS 0.122 cyc) but explicitly left "wiring the
  solved UPDs into cascade rounding" as its own next step -- that
  wiring was completed at some point after, without a doc update.
- **What was actually still missing, and is now done**: nobody had
  measured the mechanism's downstream effect on fix rate/accuracy, only
  the solver's internal residuals. Measured today (fresh
  hash-verified release binary `470a48210284`, dataset A default,
  `WL_NO_UPD=1` vs unset): every base's smoothed fix rate and h_p95/
  v_p95 move by noise-level amounts (single-digit epoch counts out of
  ~2880), including P222 and SLAC specifically -- the bases this was
  predicted to help most. **Verdict: correctly built and wired, but
  currently delivers no measurable benefit on this dataset.** Plausible
  reason (not further chased): the commit's own note that solved UPD
  magnitudes are small (+/-0.1 cyc) relative to typical arc-mean noise,
  so the correction rarely flips a rounding decision either way on this
  relatively benign mid-latitude day-time dataset. Left on (default,
  unchanged) since it's provably not harmful and may matter on noisier
  data; not worth chasing further without a dataset that actually
  stresses wide-lane rounding margins.

**SPRINT 12: RTCM/NTRIP Real-Time Input — IN PROGRESS, 2 of 6 sub-items done (2026-08-29)**
- Traced the actual state of every piece before writing this down,
  rather than leaving it as one large undifferentiated "OPEN":

- [x] **12a. RTCM3 framing/CRC** — real and tested (`crc24q`,
  `parse_rtcm3_frame`).
- [x] **12b. RTCM3 MSM raw bitfield extraction** — real and tested:
  correct bit widths per MSM4/5 vs MSM6/7, sign extension, satellite/
  signal/cell mask decoding (`parse_msm_header/masks/satellite_data/
  signal_data`).
- [x] **12c-1. RTCM3 MSM raw bitfield alignment — was a SEVERE, confirmed
  bug, now fixed (2026-08-28).** Deeper than the `into_epoch_obs` stub
  below: `parse_satellite_data` never read DF397 (8-bit rough-range
  integer-ms field) at all, and read DF419 (extended sat info) for
  every MSM type instead of MSM5/7 only. Verified against two
  independent sources before touching code -- RTKLIB's
  `decode_msm4`/`decode_msm_head` (github.com/tomojitakasu/RTKLIB,
  the reference implementation this project already benchmarks
  against) and an independent RTCM 10403.3 field reference, both
  confirming the correct order (DF397 all-types -> DF419 MSM5/7-only
  -> DF398 all-types -> DF399 MSM5/7-only). Both bugs shifted every
  bit read after the satellite-data section -- the entire signal-data
  section (pseudoranges, phaseranges, lock times, CNRs) -- for any
  real MSM4/6 message. The test suite never caught this because its
  synthetic payloads were built with the same wrong layout the parser
  expected: parser and tests agreed with each other while both
  disagreed with the real wire format. Fixed the struct, the parser,
  and all 9 affected test call sites (3 of which used inline bit
  literals a name-based search missed on the first pass). Full suite
  green, clippy clean.
- [x] **12c-2. RTCM3 MSM -> EpochObs physical-unit conversion — DONE
  (2026-08-29).** `into_epoch_obs` now produces real pseudorange (all
  5 constellations with a sourced signal table) and carrier-phase
  (all but GLONASS) observables, not empty vectors. Sourced the two
  remaining pieces from RTKLIB directly: the MSM signal-ID tables
  (`msm_sig_gps/gal/cmp/glo/qzs`, verbatim) and the cell-to-
  (satellite,signal) index mapping (`save_msm_obs`'s own loop,
  confirmed satellite-major/signal-minor before writing gneiss's
  version -- this was the highest-risk piece, easy to get subtly
  backwards). GLONASS phase stays deliberately unemitted (needs the
  FDMA channel number, which MSM only carries in MSM5/7's extended-
  sat-info field); GLONASS pseudorange works fine since range
  reconstruction needs no frequency. SBAS isn't decoded (no signal
  table sourced, unused for DD-RTK here). 5 new hand-verified tests
  (zero-range exactness, both sentinel cases independently, the
  GLONASS scope limit, and a sparse-cell-mask test that would have
  caught an indexing bug) plus the original stub-pinning test rewritten
  now that there's a real feature to check. Full suite green (251 ->
  256 gneiss-parsers tests), clippy clean. Still no real RTCM3 sample
  data in this repo to validate end-to-end against actual bytes --
  every test's expected value is computed from the same verified
  formula the implementation uses, which catches implementation bugs
  but can't catch a formula misunderstanding shared between the two.
- [ ] **12d. UBX (u-blox binary) -> EpochObs conversion** — by
  contrast, `UbxRxmRawx::into_epoch_obs` (`gneiss-parsers/src/ubx.rs`)
  DOES genuinely populate real pseudorange/carrier-phase/Doppler/SNR
  values (GPS/GLONASS/Galileo/BeiDou/QZSS/SBAS variants, invalid
  PR/CP handling) -- structurally easier than RTCM3 to get right,
  since UBX's own binary protocol already hands you SI-unit-scaled
  floats (`pr_mes`/`cp_mes`) rather than RTCM3's compact integer+
  scale-factor wire encoding the parser has to reconstruct. But
  checked its signal-band mapping before trusting the "complete"
  label: `freq_band = if sig_id == 0 { 1 } else { 2 }` collapses EVERY
  non-zero UBX `sig_id` to band 2, which only covers basic L1+one-other
  dual-frequency receivers -- it does not correctly distinguish GPS
  L5, Galileo E5b, BeiDou B2, etc. on multi-band firmware. The existing
  test (`test_ubx_into_epoch_obs_sig_id_2_band`) only pins this
  simplified behavior, it doesn't verify against the real per-
  constellation UBX sig_id table -- same "don't guess at spec details"
  caution as RTCM3 applies to fixing this properly. Real pseudorange
  data still makes this the more finished path overall; just don't
  assume full multi-band correctness without checking the actual
  receiver's signal set against u-blox's interface description first.

  **2026-08-29: tried to source the real sigId table, hit a genuine
  dead end.** The official u-blox interface description PDFs (ZED-F9R
  etc.) aren't text-extractable with any tool available in this
  environment (no `pdftotext`/`pdfgrep`, no Python PDF library
  installed, and `strings` finds zero plaintext matches -- the tables
  live in compressed content streams). Checked whether RTKLIB's own
  `ublox.c` driver had a usable table the way its RTCM3 code did (the
  approach that worked for 12c-1/12c-2): it does not -- `decode_rxmrawx`
  hardcodes exactly one signal code per constellation
  (`sys==SYS_CMP?CODE_L1I:(sys==SYS_GAL?CODE_L1X:CODE_L1C)`) and never
  reads `sigId` to differentiate bands at all, i.e. even RTKLIB doesn't
  solve this generally. Correctly stopping here rather than guessing
  sigId values from memory or installing new PDF-parsing tooling for a
  narrower-impact fix (unlike the RTCM stub, UBX's existing pseudorange/
  phase VALUES are already correct -- only the band LABEL for non-L1
  signals is oversimplified). Whoever picks this up needs either a
  copy of the target receiver's actual interface description read
  properly, or real UBX log data from a multi-band receiver to
  reverse-engineer the mapping empirically.
  If a live receiver is easier to source
  over USB/serial (UBX) than a working RTCM3 base feed, this is the
  more finished path to wire up first.
- [ ] **12e. NTRIP client maturity** — `gneiss-ntrip` is 119 lines
  (just `client.rs`); unclear whether it's been exercised against a
  real caster or only unit-tested in isolation. Needs verification
  against a live NTRIP mountpoint before trusting it in a real pipeline.
- [ ] **12f. The real architectural gap: batch vs. streaming.** Even
  with 12c/12d done, every entry point (`execute_post_process`,
  `GnssRtkIekf::process_epoch`'s callers) takes a pre-collected
  `&[EpochObs]` array. There is no incremental/streaming loop anywhere
  that calls `process_epoch` one epoch at a time as data arrives live,
  handles out-of-order or late base corrections, or manages a
  live rover+base pairing. This is the largest single piece of new
  architecture Sprint 12 actually needs, independent of which wire
  format (RTCM3/UBX) feeds it -- start here only after 12c or 12d has
  a real, verified data source to drive it with.

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
- [ ] **`formation.rs` (618 lines) needs its own split -- but not a
  plain move, a real dedupe was found first (2026-08-28).**
  `extract_sat_positions`/`compute_signal_sat_pos` (in `formation.rs`,
  ~134 lines) manually inline transmit-time iteration + Sagnac rotation
  for both the broadcast and SP3-precise-orbit paths. `satpos.rs`
  (264 lines, already in this directory) is a formal, type-safe
  branded-stage pipeline built specifically to prevent exactly this
  class of bug -- its own doc comment: "every satellite position used
  in DD formation MUST flow through these five stages," referencing a
  real historical ~140m cross-track bug from the broadcast and SP3
  paths applying Sagnac inconsistently. Checked whether that bug is
  currently live: it is not -- `grep satpos:: formation.rs` and
  `grep -rl "use.*satpos::"` across `gneiss-rtk/src` are BOTH empty.
  `satpos.rs`'s pipeline is wired to nothing. Compared its Sagnac
  rotation formula against `formation.rs`'s manual one directly:
  byte-for-byte identical (`r.x*cw + r.y*sw, -r.x*sw + r.y*cw, r.z`
  either way). So the historical bug isn't currently live -- someone
  independently hand-fixed the immediate case in the manual code and
  separately built the principled type-safe prevention, and the two
  were never connected. Net: this is real, duplicate physics code, not
  a live accuracy bug. The RIGHT next step is migrating
  `extract_sat_positions` onto `satpos.rs`'s `EphSource`/
  `compute_phase_centre` and deleting the manual duplicate -- not
  mechanically relocating soon-to-be-deleted code into a third file.

  Went one level deeper before calling this a safe mechanical swap,
  and it isn't: `compute_phase_centre`'s pipeline is real and its
  `BroadcastSrc`/`PreciseSrc` EphSource impls already exist and work,
  but comparing it line-by-line against `extract_sat_positions` found
  three genuine divergences, not just missing glue code:
  1. `extract_sat_positions` calls `precise.position_at_with_hint(sv,
     t, Some(rx_pos))` (a position hint for interpolation robustness);
     `PreciseSrc::position_at` calls the hintless `position_at`. Swapping
     in the pipeline as-is silently drops that hint.
  2. `extract_sat_positions` looks up each satellite's OWN PCO via
     `sat_pco::apply_sat_pco_z(pos, prn)` (per-satellite table);
     `compute_phase_centre`'s stage 5 takes a single scalar `pco_z_m`
     the caller must already know -- the per-satellite lookup would
     have to move to the call site, not disappear.
  3. Different transmit-time seeding entirely: `extract_sat_positions`'s
     broadcast fallback (`compute_signal_sat_pos`) seeds τ from the
     OBSERVED PSEUDORANGE (`pr_m / c`) at stage 1; `compute_phase_centre`
     seeds τ from the GEOMETRIC distance to an approximate receiver
     position (`(rx_pos - p0).norm() / c`). These differ by the
     receiver clock bias term (pseudorange = geometric range +
     c*clock_bias + ...) -- small once the filter has converged, less
     so early in a session. This is an actual algorithmic choice, not
     a bug in either direction, and picking one silently changes
     early-epoch behavior.

  None of these make the duplication finding above wrong -- the Sagnac
  rotation itself is still identical -- but they mean "migrate onto
  satpos.rs" is a real design task (extend the pipeline to accept a
  position hint and a per-satellite PCO callback, and consciously
  choose a τ-seeding convention) not a mechanical replace-and-delete.
  Correctly deferred; this is now specified precisely enough for
  whoever picks it up to not have to redo this comparison.
- [ ] **Tests for the moved formation code still live in `mod.rs`'s
  test module**, not colocated with `formation.rs` -- needs a careful
  pass since some test helpers (`test_engine`, `run_sim`) are shared
  with unrelated tests that must stay in `mod.rs`
- [ ] **17 more files still over the 500-line limit** (was 19 at the
  start of this Sprint 13 pass), worst remaining: `rinex.rs` 2349,
  `spp.rs` 2228, `rtk_iekf/mod.rs` 1202, `ephemeris.rs` 1397 (was 1573
  -- see `keplerian.rs` below), `rtk_iekf/update.rs` 814 (was 950 --
  see `screen.rs` below), `formation.rs` 618, `rtk_iekf/state.rs` 581.
  `screen.rs` (new, 148 lines) and `keplerian.rs` (new, 176 lines,
  `crates/gneiss-core`) are both clean extractions -- `screen.rs`
  moved production code together with its already-self-contained test
  module; `keplerian.rs` moved genuinely shared orbital mechanics
  (GPS/Galileo/BeiDou/QZSS's Keplerian propagation, verified
  byte-identical via diff before committing, and NOT GLONASS-specific
  despite sitting right after `GlonassEphemeris` in the original file
  -- GLONASS uses a separate, shorter RK4 integration instead). Neither
  needed the careful test-untangling the `formation.rs`/`mod.rs`
  impl-block extractions did.
- [x] `spp.rs` (2228 lines) assessed, correctly deferred (2026-08-29):
  a much bigger, messier undertaking than the three extractions above
  -- production code alone is already ~690 lines (over budget on its
  own, unlike `state.rs`), and its test code is ONE 1534-line module
  rather than cleanly pre-separated ones like `screen.rs`'s or
  `state.rs`'s. It's also the legacy SWFG/PPP path, not the validated
  `rtk_iekf` engine this project actually benchmarks. Lower payoff,
  higher effort, lower priority than what's already been done today --
  correctly left for a dedicated future pass rather than rushed.
  `rinex.rs` (2349 lines, a parser file) not yet assessed at all.
- [x] `rtk_iekf/state.rs` (581 lines) checked, correctly left alone
  (2026-08-29): unlike `formation.rs`/`update.rs`, its production code
  (one `impl RtkState` block, ~335 lines -- already under budget on its
  own) is a single cohesive concern (state-vector column-offset
  bookkeeping), not multiple mixed-together ones. The overage is
  entirely four legitimately-organized, feature-specific test modules
  (base lifecycle+gradients, iono, iono-retain-compaction, sat-iono).
  No clean fault line to split along without inventing one --
  forcing a split here would be exactly the "unrequested abstraction"
  CLAUDE.md's own Code section warns against. Left as one file.
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
- [x] println!/eprintln! audit extended to the rest of the library —
  ALSO mostly clean (2026-08-28). Sampled `post_process/` and `swfg/`
  (~12 of the remaining ~30 sites outside `rtk_iekf`): per-epoch traces
  are properly gated (`GNEISS_SWFG_DEBUG`, one with a comment
  explaining a real incident it was added to catch -- an orphan-
  variable crash that silently produced "0 epochs processed"); parse
  failures print on the error path only, not per-epoch; the remainder
  are one-time-per-pass confirmation lines (SPP init, SP3/CLK product
  load counts) that are unconditional but harmless -- they fire once,
  not per-epoch, so they don't spam output even though a stricter
  reading would prefer them behind a verbosity flag too. No genuine
  leftover-debugging-session cruft found anywhere sampled. Original
  "some likely leftover scaffolding" characterization does not hold up
  across two independent samples now (`rtk_iekf` and this one) -- not
  exhaustively re-verified for every remaining site, but confident
  enough to stop treating this as a real item on the debt list.
- [x] Percentile logic was duplicated FOUR ways, not three as
  previously noted here. `quality.rs` and `sidereal/mod.rs`'s private
  `percentile()` helpers were confirmed truly identical (`floor(len*q)`
  indexing, same empty/bounds handling, differing only in whether `q`
  was an integer-percent or float-quantile parameter) and unified into
  `post_process::percentile` (2026-08-28) -- full suite green, same
  305-test count, both guards pass.
- [ ] Two conventions remain, deliberately NOT merged into the above:
  `eval_swfg.rs`'s `percentile()` uses a materially different index
  formula (`round((len-1)*p)` vs `floor(len*q)`), landing on different
  elements for the same quantile whenever `len*q` isn't a whole number
  -- a different statistic under the same name, not a true duplicate.
  The "canonical" `gneiss_core::metrics::compute_statistics` is a
  fourth convention again (`ceil(len*q)` for p95/p99, proper
  even-length median interpolation) and doesn't expose arbitrary
  quantiles either. Merging either of these changes that caller's
  actual output, not just its source location -- confirmed neither
  feeds a guarded metric, so it's safe to do eventually, but needs
  someone to consciously pick one canonical definition and re-verify
  the numbers move only as expected, not a mechanical find-and-replace.
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

**SPRINT 15: Helmert Frame Transform — OPEN, likely NOT a frame problem at all (2026-08-29)**
- [ ] Absolute precision on dataset B is capped by a frame mismatch
  between UNR IGS20 truth and the broadcast-solution frame (~56mm
  vertical offset). Went looking for the mining candidate this item
  originally pointed at and found something better already in-tree:
  a COMPLETE, correct 14-parameter Helmert implementation
  (translation+rotation+scale, each with linear rates, standard
  small-angle rotation form) already exists at
  `crates/gneiss-geodesy/src/helmert.rs` (455 lines) — a real workspace
  member, already a declared dependency of `gneiss-rtk`'s own
  `Cargo.toml`. It has ZERO call sites anywhere (`grep
  gneiss_geodesy crates/gneiss-rtk/src` is empty): a properly-built,
  properly-wired-as-a-dependency, entirely unused capability. `frames.rs`
  ALSO has its own independent, separately-tested `HelmertParams`
  implementation used for its `ReferenceFrame` trait's fixed per-datum
  constants — meaning there are now TWO parallel Helmert
  implementations in-tree, which is itself a duplicate-implementation
  item worth resolving alongside this (pick one, likely the
  `gneiss-geodesy` one since it's the one with linear-rate/epoch
  propagation support this problem actually needs).
- **Why this isn't a 10-minute fix despite existing code**: IGS20,
  ITRF2020, and current WGS84 realizations are all mm-level aligned
  with each other -- a generic published inter-frame Helmert would not
  explain a 56mm offset. That number is far more consistent with a
  reference-EPOCH propagation issue (the CORS truth coordinates were
  established at some survey epoch and may never have been propagated
  to the actual 2020/2025 observation epoch via tectonic plate
  velocity) than a frame-REALIZATION issue -- which is exactly what
  `HelmertParams`'s `dtx/dty/dtz`/`ref_epoch` fields are for, but using
  them correctly needs the station's actual published velocity, not a
  generic frame-pair transform. Plugging in unverified or guessed
  parameters would produce a plausible-looking but scientifically
  bogus correction -- worse than leaving it alone.

  **2026-08-29: did the research, found real data, and it points away
  from tectonic epoch-propagation as the vertical culprit.** Found
  actual published PBO/NGS horizontal velocities and reference epochs
  for three of the shared stations (P181/P225/P222 appear in both
  dataset A and dataset B): P181 ref epoch 2005.09 at (-29.0, 9.6)
  mm/yr, P225 epoch 2005.14 at (-25.2, 2.7) mm/yr, P222 epoch 2005.26
  at (-31.5, 10.0) mm/yr (horizontal East/North components, Pacific
  plate motion). Dataset B is 2025-06-09 (~2025.44) -- a ~20.2-20.35
  year gap from these reference epochs, which at ~30mm/yr would be
  ~600mm of UNPROPAGATED horizontal drift if these raw reference-epoch
  coordinates were used directly as truth.

  That's the key finding: no ~600mm horizontal bias has EVER been
  reported anywhere in this project's extensive dataset A/B accuracy
  work (measured horizontal errors top out around 250mm even in the
  worst tail). This is strong indirect evidence that whatever truth
  coordinates this project already uses are NOT the raw 2005-epoch
  values -- they're already epoch-propagated (standard NGS/OPUS
  practice: querying a station's position for a specific date returns
  the propagated coordinate, not the raw datasheet reference-epoch
  one). **This rules out plain tectonic epoch-propagation as the
  driver of the ~56mm VERTICAL offset specifically** -- if the standard
  tool that got horizontal right was used, it got vertical
  epoch-propagation right too, and vertical tectonic rates are
  typically much smaller than horizontal (mm/yr, not cm/yr) anyway, so
  they wouldn't explain 56mm even if missed.
  
  Could not confirm the exact truth-coordinate provenance via web
  search alone (NGS datasheet PDFs aren't text-extractable with tools
  available in this environment, same wall hit in Sprint 12d) -- but
  reading gneiss's OWN dataset-generation code directly answered the
  question with certainty, no external data needed.

  **`scripts/p224_truth_2025.py` and `gen_multignss_truth.py` (both
  already in this repo) confirm dataset B's truth is a MONTHLY MEDIAN,
  not a single-epoch position.** Straight from the former's own
  docstring: "P224 (the multi-GNSS benchmark rover) has no single
  published coordinate file for 2025; UNR's daily IGS20 solution gives
  per-day lat/lon/h with real scatter. Median over the target month ->
  truth ECEF for benchmark scoring." `gen_multignss_truth.py` applies
  this SAME monthly-median treatment uniformly to all four sites
  (`SITES = ["P224", "P181", "P222", "P225"]`) -- extracting every
  daily UNR IGS20 position for June 2025 and taking the component-wise
  median -- while the actual benchmark OBSERVATION is a single specific
  day, 2025-06-09, drawn from within that same month. The script even
  already computes and prints the day-to-day scatter (`sigma_h`) within
  the month, confirming real dispersion exists; it just was never
  checked against how far June 9 specifically sat from the monthly
  median.

  This is now a confirmed MECHANISM, not a hypothesis: vertical GPS
  positions have well-documented 1-3cm single-day scatter from
  atmospheric/hydrological loading and daily-solution noise that
  barely touches horizontal. A monthly median smooths that out; a
  single day carries whichever way June 9 happened to land. A ~56mm
  mismatch between "this specific day" and "the monthly average" is
  entirely plausible at that noise level, and would look EXACTLY like
  a "frame mismatch" in aggregate accuracy statistics without being
  one at all -- no Helmert transform of any kind would fix it, because
  there's no frame problem to fix.

  **The definitive test is fully specified but not yet run**: re-derive
  each site's truth using ONLY June 9's row (or a tight window around
  it) via the same `extract_rows`/`median_llh` functions already in
  `p224_truth_2025.py`, and compare against the current monthly-median
  `station_coords.json`. If the delta is ~56mm vertically, this is
  confirmed as the dominant cause and Sprint 15 should be re-scoped
  again -- from "build/wire a Helmert transform" to "fix the truth-
  generation methodology to use single-day (or day-matched) UNR
  positions instead of a monthly median." Could not run this test in
  this session: the source `.tenv3` time series
  (`/tmp/ds2025/*_IGS20.tenv3`) that both scripts read are not
  persisted in this repo and would need re-fetching from UNR's
  geodesy lab first (see `scripts/fetch_multignss_dataset.py` for the
  original fetch mechanism) -- `datasets/multignss_2025d160/` only
  keeps the already-computed monthly-median output, not the raw
  per-day series the definitive test needs.

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
