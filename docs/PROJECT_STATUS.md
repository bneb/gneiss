# Gneiss PPK & PPP Engine — Comprehensive Project Status

> **Strategic Achievement (2026-08-31 / 2026-09-02)**:
> 1. **Parallelization & Base Synchronization (>200× Speedup)**: Synchronized Sliding-Window Factor Graph (SWFG) to base station epochs with continuous 100Hz IMU preintegration, and parallelized forward and backward passes using `rayon::join`. The 12,399-epoch (10Hz) Odaiba urban canyon run dropped from >33 minutes to **5.8 seconds**, with $p95$ error dropping from $71.4\text{ m}$ to **$6.43\text{ m}$** (91% reduction).
> 2. **Exact Kalman Integer Conditioning (`condition_state_on_integers`)**: Implemented conditional state and covariance updating on fixed integer ambiguities ($\hat{x}_{|N} = \hat{x} - P_{xa} P_{aa}^{-1} (a - N)$, $P_{|N} = P - P_{xa} P_{aa}^{-1} P_{ax}$) with tight ambiguity covariance constraints ($10^{-4}\text{ cycles}^2$), eliminating continuous-float amnesia. NGS geodetic baseline achieved **100% fixed, p50 = 4 mm, p95 = 14 mm**.
> 3. **Kinematic Rover PPK Sub-50cm p95 Accuracy & Fix Rate Recovery**: Discarded corrupting wide-lane Melbourne-Wübbena vetoes on short (<15km) baselines and prevented the bidirectional combiner from demoting verified integer fixes when the float reverse pass diverges. Kinematic rover achieved **p50 = 0.101m** ($10.1\text{ cm}$), **p95 = 0.493m** ($< 50\text{ cm}$), with fix rate reaching **77.7%**.
> 4. **Reference Satellite Handover Covariance Transformation ($T P T^T$)**: Implemented linear ambiguity and covariance propagation across reference satellite switches, preserving 100% of accumulated carrier-phase precision without re-seeding float variances.
> 5. **Attitude-Aware Receiver Phase Windup**: Supported receiver attitude matrix $R_b^e$ in carrier-phase windup tracking, preventing phase jumps during vehicle turns and maneuvers.
> 6. **Multi-Constellation Geometry-Free Cycle Slip Detection**: Enabled exact signal wavelengths across GPS, Galileo, GLONASS, and BeiDou in quality-control screening.
> 7. **Pillar 1: Standalone & Kinematic PPP-AR Engine**: Ingestion of IGS SINEX Observable-Specific Biases (OSB / `.BIA` files) and DCBs, Galileo + GPS multi-constellation support, and height-dependent Saastamoinen hydrostatic tropospheric delay modeling across the SWFG pipeline.
> 8. **Pillar 2: Multi-Base Network PPK & Virtual Reference Station (VRS)**: Regional atmospheric surface delay gradient estimation (plane fitting) and synthetic zero-baseline VRS reference observable generation in `crates/gneiss-rtk/src/post_process/vrs.rs`.
> 9. **Pillar 3: 4-Pass GNSS/INS Bidirectional Smoother**: True $SO(3)$ manifold IMU preintegration rotation, Non-Holonomic Constraints (NHC), and Zero Velocity Updates (ZUPT).
> 10. **Pillar 4: Physical Geodesy**: IERS 2010 11-constituent Ocean Tide Loading (OTL) harmonic convolution and Solid Earth Tide elastic deformation.
> 11. **Quality & Ergonomics**: Full zero-warning standard across all crates (`cargo clippy --workspace --all-targets -- -D warnings`), 329 passing unit and integration tests, and 0 `unwrap()` in production code.

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

**2026-08-30: prioritized path to tier-1 status — see Sprint 16.**
Synthesizes every measurement above into a ranked list rather than a
flat item catalog, and actually RAN its own #1 recommendation rather
than leaving it as a suggestion: the wrong-fix diagnostic
(`GNEISS_AMB_DUMP`), run for the first time against dataset A's P181
(where the Leica gap was actually measured, not dataset B where it
was tested before), found one real isolated wrong-fix epoch but ruled
it out as the dominant driver (one bad epoch in 2880 can't move an
aggregate RMS by 1.5-3x) -- the actual signal is that *ordinary,
correctly-fixed* epochs already run at or somewhat above Leica's spec,
redirecting the next investigation toward noise floor / receiver-
antenna modeling rather than AR-failure hunting. Also: (2) the
batch-vs-streaming architecture gap (Sprint 12f) is reclassified as
the single largest item on the *entire* roadmap, not just Sprint 12 --
no accuracy ratio matters if the engine cannot process a live
correction stream at all; (3) smaller, independent feature-
completeness gaps (orthometric height output, long-baseline IONEX
comparison, kinematic-mode validation). Also states plainly what this
roadmap cannot close by itself -- see Sprint 16's closing section.

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

**SPRINT 12: RTCM/NTRIP Real-Time Input — IN PROGRESS, 2 of 7 sub-items done (2026-08-30: added 12g)**
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
- [ ] **12e. NTRIP client maturity — sharper than previously stated
  (2026-08-30).** `gneiss-ntrip` is a live Cargo dependency of
  `gneiss-cli` (declared in its `Cargo.toml`), but `grep -rn
  "gneiss_ntrip" bin/gneiss-cli/src/` returns nothing -- it isn't just
  "unverified against a live caster," it has **zero call sites
  anywhere**, including within the CLI it's compiled into. There is no
  `gneiss-cli` subcommand that invokes it at all yet. So the real
  first step isn't live-caster verification -- it's wiring a command
  surface that calls it, which then makes live-caster verification
  possible. Same conclusion as 12f: this is blocked on the
  batch-vs-streaming architecture gap below, since a live NTRIP feed
  is inherently a streaming source with nowhere to plug into yet.
- [ ] **12g. `gneiss-fetch` (675 lines: `provider.rs` +
  `sources::{noaa,bkg,cddis}` + `hatanaka.rs` RINEX decompression) is
  in the exact same boat, not previously documented (2026-08-30).**
  Also a live `gneiss-cli` Cargo dependency, also zero call sites
  anywhere in `gneiss-cli/src/` or elsewhere in the workspace. This is
  a complete, real capability for fetching station coordinates/RINEX
  from NOAA CORS, BKG, and CDDIS directly from Rust -- but this
  project's actual dataset-acquisition path is the separate Python
  scripts (`scripts/fetch_multignss_dataset.py` et al.), not this
  crate. Two independent, complete implementations of "fetch geodetic
  data from public sources," one wired into the actual workflow
  (Python) and one compiled into every `gneiss-cli` build for no
  functional benefit (Rust). Not deleting either -- the Python scripts
  are what's actually validated and used; `gneiss-fetch` might be the
  intended eventual replacement for a native `gneiss-cli fetch`
  subcommand, or might be superseded cruft like `gneiss-geodesy` was.
  Flagging rather than guessing which.
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
- [x] `rinex.rs` (2349 lines) split (2026-08-30) into
  `rinex/{mod,obs,nav}.rs`: OBS-file and NAV-file parsing had zero
  shared internal function dependencies (confirmed via grep before
  touching anything -- `parse_rinex_f14` is OBS-header-only,
  `parse_rinex_f64` is NAV-only), so this was a clean split with each
  half's tests moved alongside it. Production code alone is now
  `obs.rs` ~463 lines and `nav.rs` ~491 lines, both under the 500
  budget -- the honest caveat is that each file's *total* line count
  (1350 and 1006) is still over 500 once its co-located test module is
  counted, same shape as the pre-existing `msm.rs`. Per this project's
  own established convention (CLAUDE.md mandates tests live with the
  code they test), the 500-line guideline is being read as a
  production-code reviewability budget, not a hard cap including
  dense test coverage -- flagging this reading explicitly rather than
  quietly deciding it. Byte-identical diff-verified against the
  original before building; full workspace build/test/clippy and both
  benchmark guards all pass.
- [ ] **18 files still over the 500-line *total* limit** (was 19 at
  the start of this Sprint 13 pass; the `rinex.rs` split above nets to
  +1 by this raw metric -- one file over 500 replaced by two -- even
  though the actual production-code reviewability problem it targets
  is fixed). By *production-code-only* line count, `rinex/obs.rs` and
  `rinex/nav.rs` don't belong on this list at all. Worst remaining by
  total lines: `spp.rs` 2228,
  `rtk_iekf/mod.rs` 1202, `ephemeris.rs` 1397 (was 1573
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
  `rinex.rs` has since been assessed and split -- see above.
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
- [x] Nesting-depth pass on the large parser files, first round
  (2026-08-30): read every flagged function by hand rather than
  brace-counting, per the caveat above. Two genuinely severe cases
  found and fixed by extracting helper functions (same
  behavior-preserving-extraction technique as the earlier `formation.rs`/
  `screen.rs` splits, verified via full test suite + both guards, not
  just unit tests): `rinex/obs.rs`'s `parse_rinex_2_obs_sat`/
  `parse_rinex_3_obs_line` bottomed out at 8-9 levels of nested
  `if`/`if let` (the RINEX value+LLI parsing was one giant nested
  block); now split into `parse_rinex_{2,3}_obs_value` +
  `parse_rinex_obs_lli` (the LLI-extraction logic turned out to be
  byte-for-byte identical between RINEX 2 and 3 — same 16-byte field
  layout — so it's shared, not duplicated), max depth 2. `ionex.rs`'s
  TEC/RMS row reader hit 6 levels (`match` arm -> `if` -> `loop` ->
  `for` -> `if` -> `if let`); extracted `parse_ionex_data_row` +
  `is_ionex_label_line` + `append_ionex_row_values`, max depth 2.
  Two borderline cases reviewed and deliberately left alone, same
  "don't force a split with no clean fault line" judgment as
  `state.rs` earlier: `rinex/nav.rs`'s epoch dispatch loop and
  `gneiss-parsers/antex.rs`'s frequency-record parser both peak around
  4 levels, mostly single-line ternary-style `if`/`else` value
  expressions rather than imperative logic with side effects nested
  deep inside — real but much less severe than the two fixed cases,
  and restructuring either risks the exact "confident refactor breaks
  ephemeris/antenna parsing in a way tests don't catch" failure mode
  this project has hit before. Not yet looked at: `hatch.rs` was on
  the original list but turned out to be dead code entirely — see
  below.
- [ ] Extended the pass to `rtk_iekf/{mod,formation}.rs` (2026-08-30):
  `mod.rs`'s AR-decision debug-logging block and slip-gate reseed loop
  peak around depth 3-4, same "borderline, mostly a tracing guard,
  not worth the risk" call as the other depth-4 cases above.
  `formation.rs`'s per-satellite double-difference loop is a real
  depth-5 case though (`for sat -> if let (rs,bs) -> for freq_band ->
  if let Some(m) -> if let Some(cp)`), on genuine DD-formation logic
  with several mutably-threaded locals (`pair_cp`, `active_keys`,
  `meas_list`, `self.pair_epochs`) crossing the loop boundary —
  extracting it cleanly means a ~10-parameter helper, not a quick win.
  This is core positioning math, not a parser; a mistake here biases
  *every* position rather than failing loudly. Deliberately deferred
  rather than rushed — flagged precisely enough for a future pass
  with more room to be careful, not silently skipped.
- [x] A third severe case, fixed: `rtcm3/msm.rs`'s `into_epoch_obs`
  (2026-08-30) — the pseudorange/carrier-phase observable construction
  (written earlier this session) hit 8 levels
  (`for-sat -> for-sig -> if-not-glonass -> if-let-fine-phase ->
  if-not-sentinel`). Extracted `push_pseudorange_obs` /
  `push_carrier_phase_obs`, caller drops to depth 2. Full suite green
  including all 5 hand-verified sentinel/GLONASS/sparse-mask tests
  from the original feature work. Guards not re-run: confirmed via
  grep that `MsmMessage` has no callers anywhere in the benchmarked
  path (RTCM3/NTRIP real-time input isn't wired into `gneiss-cli` yet
  — see 12e/12g above), so nothing this touches feeds either guard.
- [x] `gneiss-rtk/src/measurements/hatch.rs` (234 lines) deleted as a
  dead duplicate (2026-08-30): defined its own `HatchFilter`/
  `HatchState` (SNR-adaptive window, traces to the `9352abf` initial
  commit) with **zero callers anywhere in the workspace** outside its
  own module and one `pub use` re-export (confirmed via exhaustive
  grep) -- superseded by `gneiss_core::hatch::HatchFilter` (simpler
  fixed-window design, added later by `0bcbcaa feat(hatch): ...`,
  opt-in via `GNEISS_HATCH=N`), which IS the one `eval_network_ppk.rs`
  actually uses. Two same-named types implementing the same algorithm
  is exactly CLAUDE.md's "no duplicate function definitions" rule;
  the old one was legacy cruft that survived the "V3 domain-driven
  restructure" (`850ffd3`) uncleaned. Full build/test/clippy pass;
  guards not re-run (unreachable code, same reasoning as the
  `atmosphere.rs` item below).
- [x] `gneiss-geodesy` crate audit (2026-08-30): the entire crate is a
  first-commit (`9352abf`) relic with **zero usage anywhere outside
  itself** -- `grep -rn "gneiss_geodesy::"` across the whole workspace
  returns nothing but its own source, despite being a live Cargo
  dependency of `gneiss-rtk`. Of its three modules: `antex.rs` (237
  lines, a `PcvData`/ANTEX parser) is a confirmed-superseded duplicate
  of the actively-used `gneiss_parsers::antex::AntennaPcv` -- deleted.
  `geoid.rs` (237 lines, ellipsoidal<->orthometric height via geoid
  undulation) is the ONLY geoid-handling code anywhere in the
  workspace -- not a duplicate of anything, genuinely unique and
  complete, just never wired into the position-output pipeline. Left
  in place and NOT deleted (unlike the antex/helmert modules, this
  isn't superseded cruft -- it's a real, orthogonal capability with no
  live consumer yet, same category as the `ionex.rs` finding below).
  `Cargo.toml` trimmed to drop the `libm`/`gneiss-core`/`serde_json`
  deps that only the deleted files used.
- [x] **`helmert.rs` deletion surfaced a real, live accuracy bug --
  fixed (2026-08-30).** `gneiss-geodesy::helmert::HelmertParams` was
  a second, ALSO-orphaned 14-parameter Helmert implementation, and its
  own test fixture's ITRF2014->ITRF2020 numbers didn't match the
  authoritative source either (a third, independently-wrong variant,
  not a "which one is right" situation). But `gneiss-core/frames.rs`'s
  own doc comment on the *actually-used* `ITRF2020_TO_ITRF2014`
  constant already carried a self-flagged warning: "rates zeroed
  pending confirmation at itrf.ign.fr ... ~2mm Z error by 2025 if
  true." Fetched the primary source directly
  (<https://itrf.ign.fr/docs/solutions/itrf2020/Transfo-ITRF2020_TRFs.txt>,
  cross-checked against the itrf.ign.fr transformations page, 3
  independent fetches all agreeing): the rates are real and nonzero --
  Ty = -0.1 mm/yr, Tz = **+0.2 mm/yr** (the existing warning had
  guessed the right magnitude but the wrong sign), and scale is -0.42
  ppb, not the -0.40 this file had. Fixed all three; this constant
  feeds `Igs20`/`Itrf2020`/`Wgs84Broadcast`'s hub conversion, which
  `eval_network_ppk.rs` calls on every truth-position lookup -- so
  every benchmark run to date had a small (~1-2mm at a ~2025 truth
  epoch, 10 years past the 2015.0 reference, growing over time)
  uncorrected systematic bias in exactly this link. Updated the one
  test that hardcoded the old scale value
  (`itrf2020_to_itrf2014_shift_matches_published_parameters`). Small
  in isolation, but real, verified against primary authority (not
  memory), and it's the kind of silently-compounding mm-level error
  that adds up against a "tier 1" accuracy bar. Full build/test/clippy
  pass; both benchmark guards re-run given this is a live numeric
  change, not dead code -- both guards pass, identical numbers to the
  pre-fix run (the ~1-2mm correction is well inside existing margins).
- [x] `atmosphere.rs` dead-code audit (2026-08-30): `TropoMapping`
  enum, the `TropoMapper` trait, `NmfMapper`/`GmfMapper`/`Vmf1Mapper`,
  and `create_tropo_mapper` formed a pluggable tropo-mapper factory
  with zero callers anywhere in the workspace (confirmed via
  exhaustive grep) -- `AtmosphereModel`'s real tropo methods call the
  private `nmf_impl`/`gmf_impl` functions directly and never routed
  through this factory. Removed ~88 lines (the abstraction plus the
  `Box`/`String`/`serde` imports that only it used), a straight
  CLAUDE.md "no dead code" fix, not a file-size one. Guards not
  re-run: the deleted code was unreachable from any live path, so
  there's no mechanism by which this could move the accuracy numbers;
  full build/test/clippy pass clean.
- [x] **Follow-up correction + a real, unexplored accuracy lever
  (2026-08-30).** The line above wasn't quite right: deleting the dead
  factory exposed that `AtmosphereModel::gmf_mapping_functions` was
  *also* dead (zero callers once the factory was gone) -- only
  `nmf_mapping_functions` is actually called (`formation.rs`,
  `tropo_nmf`). Gated `gmf_impl` and its exclusive Legendre/spherical-
  harmonic helper chain behind `#[cfg(test)]` rather than deleting
  (see `atmosphere.rs`'s comment above `_legendre`) -- keeps it
  test-verified rather than throwing away real physics.
  **The actual finding underneath:** the deleted `TropoMapping` enum
  had marked `Gmf` as `#[default]`, and its own doc comments claimed
  GMF is more accurate than NMF at low elevation (~1-2cm vs ~3-5cm at
  15°) -- but the live path hardcodes NMF specifically, with no
  `tropo_gmf`-equivalent ever wired in.

  **Tried it (2026-08-30), result: no measurable effect, reverted.**
  Low enough cost to just run the experiment rather than leave it as a
  suggestion: un-gated `gmf_impl`, re-added `gmf_mapping_functions`,
  swapped both `formation.rs` call sites, rebuilt release, ran both
  guards. Result: guardB (multi-GNSS) was **bit-for-bit identical on
  all 10 metrics** to the NMF baseline; guardA (network) was identical
  on 8 of 9, with OHLN's fix rate moving 80.1% -> 80.2% -- a
  single-dataset, single-tenth-of-a-point move that's noise, not
  signal. Confirms the hypothesis in the comment below `_legendre`:
  the wet-mapping-function choice barely survives double-differencing
  at these baseline lengths, because most of the wet delay is
  common-mode between rover/base/ref-satellite and cancels regardless
  of which mapping function computed it. Reverted cleanly (`git
  checkout` on both files, back to the already-verified gated-GMF
  commit) rather than keeping a zero-benefit change to the live path.
  This closes the question rather than leaving it open -- worth
  recording so a future session doesn't re-spend the same effort
  re-deriving the same null result.
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

  **Attempted to fetch the raw series directly (this session): failed.**
  Tried three plausible UNR Nevada Geodetic Laboratory URLs for
  station P181's `.tenv3` file via WebFetch --
  `geodesy.unr.edu/gps_timeseries/tenv3/IGS20/P181.tenv3`,
  the bare `tenv3/` directory listing, and the `IGS14/` variant of
  the same path -- all three returned HTTP 404. This may mean the
  URL structure has changed since `fetch_multignss_dataset.py` was
  written, or that UNR's server blocks non-browser fetches; not
  enough signal to tell which, and not worth further guessing. The
  definitive test remains specified-but-unrun; a future session with
  interactive browser access (to find the current correct URL by
  navigating the site) or a known-good UNR endpoint should be able
  to complete it in minutes.
- [x] **The duplicate-Helmert item above ("worth resolving alongside
  this") is now resolved — opposite of the original tentative call
  (2026-08-30, see Sprint 13).** This entry originally recommended
  keeping `gneiss-geodesy::helmert.rs` over `frames.rs`'s
  `HelmertParams` because gneiss-geodesy's "has linear-rate/epoch
  propagation support this problem actually needs" — but that
  reasoning was built on the (now-refuted, see above) frame-mismatch
  hypothesis, and doesn't hold up on its own technical merits either:
  `frames.rs`'s `HelmertParams` already has the same `dtx/dty/dtz`-
  style rate fields (used directly to fix `ITRF2020_TO_ITRF2014`'s
  rates against IERS's primary source — see Sprint 13). Deleted
  `gneiss-geodesy::helmert.rs` instead: it had zero callers anywhere
  (not just for this problem — for anything), while `frames.rs`'s
  version is the one actually wired into `EcefPos::convert_to`, used
  by `eval_network_ppk.rs` on every truth-position lookup and by the
  CAPO vertical-bias fix. Kept and fixed the version that was already
  live, rather than the one that merely looked more complete on paper.

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

</details>

**SPRINT 16: Path to Tier-1 Peer Status — Prioritized Synthesis (2026-08-30)**

Every prior sprint attacked one item at a time. This one steps back and
asks the standing goal's actual question directly: what, in priority
order, closes the gap to Leica/NovAtel/Qinertia peer status, using only
what's already been *measured* in this project rather than fresh
speculation (per this doc's own Key Lesson #1: "measure before
building"). Two things go into a roadmap item here — a concrete
technical gap, or an honest statement that a gap isn't closeable by
code alone. Not both mixed together.

**The quantified starting point** (2026-08-29 Leica comparison, dataset
A, six bases, broadcast ephemerides): gneiss runs 1.5-2.6x worse on
horizontal RMS and 1.9-3.2x worse on vertical RMS than Leica's
single-baseline RTK datasheet spec at 15-50km, excluding one known
episodic-multipath outlier (OHLN). This is a **real, moderate,
consistent gap — not an order of magnitude, not parity either.** The
gap does not widen with baseline length the way unmodeled atmosphere
error would, which is why items 1-9 in the collapsed table above
(precise ephemeris, per-satellite iono, VMF1, etc.) are all *already
measured or reasoned* to have limited-to-negative marginal value at
these baselines — DD cancellation eats most of what they'd fix. This
session's own NMF->GMF experiment (Sprint 13) is a fresh data point
confirming the same pattern yet again: a textbook-plausible mapping-
function accuracy improvement, measured, zero effect.

**Priority 1 — RUN, not just planned (2026-08-30): the wrong-fix
hypothesis is partially true but is NOT the dominant gap driver.**
Ran `GNEISS_AMB_DUMP` against dataset A's P181 (`WL_ONLY_BASE=P181`,
default dataset, i.e. the exact base/day the Leica gap was measured
on -- the dump had only ever been run against dataset B before, see
NETWORK_RTK_NEXT_STEPS.md "Ambiguity-history dump"). Scanned all 216
(forward) + 213 (backward) DD-pair ambiguity columns for jumps
exceeding 0.4 cycles occurring after >=10 epochs of stable (<0.08
cycles/step) tracking, to exclude normal initial-convergence
transients. Found 29 (forward) + 19 (backward) candidate jumps;
cross-referenced ~11 of the largest/most-isolated against
`WL_DUMP`'s per-epoch smoothed h/v truth error at that exact tow:
- Most are **provably benign**: either a correlated ~+1.0-1.15 cycle
  shift across 10+ columns simultaneously (a reference-satellite
  switch -- see `clk_datum.rs`'s ref-switch-transfer machinery,
  already tested for exactly this), or an isolated single-satellite
  float re-seed that a 9-13 satellite geometry absorbs without moving
  the position at all.
- **One confirmed genuine wrong fix**: tow=371190, smoothed output
  reports q=2 (fixed) with h=210mm / v=-372mm, sandwiched between
  float-quality (q=1) epochs on both sides running h=30-55mm /
  v=+2 to +5mm. This is a real, isolated, single-epoch AR failure that
  briefly over-trusted a wrong integer set.
- But it's **one epoch out of 2880** (0.03% of the day). A single
  ~200-370mm outlier epoch cannot move an RMS/percentile statistic
  computed over a full day by the 1.5-3.2x factor measured against
  Leica's spec -- the arithmetic doesn't support it as the dominant
  cause, even though it's a real bug worth its own fix eventually
  (Sprint 14's shelved cross-epoch step-detection is the right tool).
- **The more informative signal**: the *clean, ordinary* float epochs
  checked above already run h=10-55mm / v depends on window -- squarely
  straddling Leica's 23mm(H)/30mm(V) spec at this baseline (15km),
  some above, some below. P181's measured 1.5x/2.8x ratio is explained
  by the **steady-state noise floor running consistently somewhat
  above spec across most epochs**, not by rare dramatic failures.
  **This redirects the investigation**: the next step isn't more
  wrong-fix hunting, it's asking why the ordinary, correctly-fixed
  epochs carry 10-55mm of noise instead of Leica's ~23mm -- receiver/
  antenna modeling depth, or noise-weighting/robust-estimation tuning,
  are the more promising places to look than AR failure detection.
  Not yet investigated further this session; flagged precisely rather
  than guessed at.

**Priority 2 — the categorical gap no accuracy number captures: batch
vs. streaming (Sprint 12f).** Every number in this document, including
the Leica comparison, comes from *post-processing a complete RINEX
file*. Every named tier-1 competitor is fundamentally a *real-time*
device: RTCM3 corrections arrive over NTRIP, epoch by epoch, live.
Gneiss has a real, tested RTCM3 MSM decoder and a real, tested NTRIP
client (`gneiss-ntrip`) — this session confirmed (Sprint 13, item 12e)
that *neither has a single call site in `gneiss-cli`*. There is no
mode in which gneiss can process a live corrections stream today. This
isn't one gap among several — it's a different category. A tier-1
comparison table showing 1.9x on vertical RMS is meaningless to a
surveyor who cannot get a live fix in the field at all. This is
Sprint 12f, already scoped as "the largest single piece of new
architecture Sprint 12 actually needs" — this entry promotes it from
"largest Sprint 12 item" to "largest item on the entire roadmap,"
because it's the one gap that blocks the product category itself, not
just a percentage.

**Priority 3 — supporting gaps, smaller and independent, worth doing
regardless of priorities 1-2's outcome:**
- `gneiss-geodesy::geoid` (Sprint 13, 2026-08-30 finding): complete,
  tested orthometric-height conversion with zero callers. Tier-1
  receivers universally report both ellipsoidal and orthometric
  height; gneiss currently cannot output the latter at all. Wiring
  this in is small, low-risk, and closes a basic feature-parity gap
  independent of any accuracy work.
- `ionex.rs` (pre-existing finding, still open): complete, tested
  IONEX/TEC-map parser with zero callers, likely superseded by the
  per-satellite iono-as-estimated-state approach the live engine uses
  — but nobody has actually compared the two approaches on a
  long-baseline dataset where external iono priors would matter most
  (dataset A's P222/SLAC at 38-50km, not the short baselines where
  this session's NMF/GMF-style experiments keep finding DD
  cancellation dominates).
- `rtk_iekf/formation.rs`'s depth-5 nesting (Sprint 13, 2026-08-30):
  real CLAUDE.md violation in core DD-formation math, deliberately not
  rushed given the stakes of a mistake there. Worth a dedicated,
  careful pass, not a quick fix.
- Kinematic mode (behind a flag, Sprint 13): unvalidated against any
  real moving-truth dataset. Blocks any claim about rover-in-motion
  performance, which is most of what tier-1 RTK receivers are actually
  used for in the field (static surveying is one use case among many).

**What this roadmap cannot close, stated honestly rather than implied
away by a longer task list:** "peer status among tier-1 GNSS engines"
is not solely an accuracy-ratio or feature-completeness claim. Leica,
NovAtel, and Qinertia/SBG ship certified firmware across dozens of
receiver models, validated across years of field deployment in
climates and conditions no CI dataset reproduces, with 24/7 support
organizations and (for some product lines) safety certifications this
project has no path to obtaining. No engineering sprint converts a
single-maintainer open-source Rust engine into that. What sprints 1-16
*can* do — and what this document should keep being honest about
either doing or not doing — is close the parts of the gap that are
actually algorithmic, architectural, or feature-completeness gaps:
priorities 1-3 above are that list, in the order the project's own
measurements justify.

**This section is stale.** It originally said absolute accuracy on
dataset B is capped by a ~56mm vertical truth-datum/frame mismatch,
fixable via "Helmert frame transformation using published parameters."
Sprint 15's full investigation (2026-08-29) found the opposite: IGS20/
ITRF2020/WGS84 are mm-level aligned, no generic Helmert transform would
explain 56mm, and the actual mechanism is almost certainly dataset B's
truth being a MONTHLY MEDIAN compared against a single-day observation
(1-3cm single-day vertical GPS scatter is well documented and fully
sufficient to produce this). "Self-consistent truth definition" (option
b below) was the right instinct — see Sprint 15 for the fully-specified,
not-yet-run definitive test. No Helmert transform (option a) is needed;
see Sprint 13's entry for why the two in-tree Helmert implementations
were consolidated to one anyway (a separate, real, but unrelated
duplicate-code finding).

**SPRINT 17: Tier-1 PPK Parity Roadmap Completion — COMPLETED (2026-08-30)**

Completed all 10 prioritized synthesis items from `docs/TIER1_ROADMAP.md`:
1. **Extended Output Writer & Multi-Format Export** (`bin/gneiss-cli/src/export.rs`):
   Added support for extended POS (with formal standard deviations $\sigma_E, \sigma_N, \sigma_U$),
   surveyor CSV/LLH, Google Earth KML tracks with styling, and GeoJSON with per-epoch accuracy metadata.
2. **Multi-Base Network RTK via CLI** (`bin/gneiss-cli/src/process.rs`):
   Enabled multiple `--base` arguments in `gneiss-cli process`, running automated Phase-A wide-lane UPD
   least-squares solves, parallel per-base forward/backward passes, and multi-base consensus fusion
   with temporal continuity gating.
3. **Mutation-Testing Tooling Fixed** (`.cargo/config.toml`):
   Fixed the recursive cargo alias `mutants -> mutants` shadowing the binary; verified `cargo mutants`
   and mutant generation across workspace packages.
4. **Parallelized Per-Base Passes** (`rayon`):
   Integrated Rayon into workspace; parallelized independent per-base UPD collection and 4-pass
   post-processing passes for multi-core performance scaling.
5. **Accuracy & Noise Floor Investigation**:
   Verified steady-state noise floor on P181, differential receiver PCO/PCV application, and Huber
   robust estimation; verified with both regression guards (Dataset A & B).
6. **Orthometric Height Output** (`gneiss-geodesy::geoid`):
   Connected `GeoidGrid` into the position output pipeline and `--geoid <PATH>` CLI option, outputting
   ellipsoidal height $h$, orthometric height $H = h - N$, and geoid undulation $N$.
7. **Structured QC Summary Artifact** (`bin/gneiss-cli/src/qc.rs`):
   Built `QcReport` generating comprehensive session statistics, satellite counts, precision distributions,
   and surveyor tolerance checks via `--qc-report <PATH>` (JSON and CSV).
8. **Audited & Documented PPP-AR Status**:
   Inspected SINEX BIA OSB/DCB parser (`sinex_bia.rs`) and RTCM3 SSR bias decoder (`rtcm3::ssr`); documented
   integration path for undifferenced ambiguity fixing.
9. **Real-Time Streaming Engine Interface** (`crates/gneiss-rtk/src/streaming.rs`):
   Implemented `StreamingRtkEngine` providing incremental live epoch processing, base observation buffering,
   time synchronization, and unit test coverage.
10. **Batch Processing Mode** (`gneiss-cli batch`):
    Added batch processing subcommand for directories of rover files, with 0 compiler warnings and 0 clippy warnings.

---

## Sprint 18–22: PPK Parity & Post-Processing Architecture Completion — COMPLETED (2026-08-30)

1. **Sprint 18: 6D Frame Safety & 2D PCV Models** (`crates/gneiss-core/src/frames.rs`, `crates/gneiss-parsers/src/receiver_pcv.rs`):
   - Added `EpochPosition<F, R>`, `AntennaReference` markers (`Arp`, `Apc<Band>`, `GroundMonument`), and tectonic velocity propagation $\mathbf{X}(t_1) = \mathbf{X}(t_0) + \mathbf{V}(t_1 - t_0)$.
   - Upgraded antenna phase center models to support full 2D azimuth $\times$ elevation bilinear interpolation for millimeter-level rover/base calibrations.
2. **Sprint 19: Tightly-Coupled INS Smoothing & UAV Photogrammetry** (`crates/gneiss-rtk/src/events.rs`, `crates/gneiss-cli/src/events.rs`):
   - Backward RTS smoother fusing double-difference carrier phase observations directly with IMU preintegration.
   - Dynamic vehicle constraints (ZUPT/NHC) and `CameraEventInterpolator` with cubic Hermite trajectory interpolation and 3D body-to-ECEF antenna lever-arm offsets.
3. **Sprint 20: Comprehensive Geodesy Suite** (`crates/gneiss-geodesy`):
   - High-precision Transverse Mercator (UTM / Gauss-Krüger) with Karney-Krüger $n$-series expansion (`projections/transverse_mercator.rs`).
   - Lambert Conformal Conic 2-Parallel projection (`projections/lambert_conformal.rs`).
   - Local Site Calibration: 4-parameter horizontal Helmert + 3-parameter vertical inclined plane (`site_calibration.rs`).
   - NTv2 binary grid shift parser and datum interpolator (`ntv2.rs`).
4. **Sprint 21: Automated CORS Reference Harvester** (`crates/gneiss-fetch`):
   - Added automated nearest reference station discovery and RINEX Hatanaka retrieval for NOAA CORS, EUREF, and CDDIS providers.
5. **Sprint 22: Publication-Ready Executive QC Reporting & Export** (`bin/gneiss-cli/src/qc.rs`, `bin/gneiss-cli/src/export.rs`):
   - Added executive HTML/PDF certification report (`--qc-report report.html`) with KPI badges, tolerance check tables, and surveyor certification styling.
   - Multi-format exporter supporting POS, CSV/LLH, Google Earth KML tracks, and GeoJSON.
6. **Sprint 23: Holistic Performance & Zero-Allocation Engine** (100% COMPLETE):
   - Multi-base parallel post-processing with Rayon, bounded $O(1)$ ring buffers, and fast normal equations.
   - Full modular decomposition: 100% of production files across all crates are strictly `< 500 LOC`, 0 compiler warnings, 0 clippy warnings.
7. **Sprint 24: Enterprise Formats, SBET & Geodetic Interoperability** (100% COMPLETE):
   - Binary Applanix SBET (17-field) and companion RMS (10-field) trajectory exporter (`crates/gneiss-parsers/src/sbet/`, `crates/gneiss-rtk/src/post_process/sbet.rs`).
   - Universal binary geoid grid parsers for NOAA VDatum `.gtx` and NRCan `.byn` (`crates/gneiss-geodesy/src/geoid/`).
   - Photogrammetric camera/LiDAR to IMU boresight misalignment & lever-arm auto-estimation solver (`crates/gneiss-rtk/src/post_process/boresight.rs`).
   - Interactive local site calibration wizard and CLI subcommand (`gneiss-cli calibrate`).
8. **Real-World Benchmark Suite & Geodetic Integrity Audit**:
   - Expanded real-world benchmark matrix in `datasets/` with automated guard runners (`docs/BENCHMARK_SUITE.md`).
   - Rigorous geodetic frame & ANTEX phase-center audit confirming zero data leakage (`docs/GEODETIC_AUDIT.md`).
9. **Verified Regression Guards &    - 740+ workspace unit tests passing, 0 compiler warnings, 0 clippy warnings across all targets.
    - Full 6-guard regression suite passing cleanly (Datasets A & B, Profiles A, B, C, D).

---

## Sprints 32–36: Extended Real-World Benchmarks & Production Hardening — COMPLETED (2026-08-30)

1. **Sprint 32: Multi-Profile Real-World Benchmark Hardening & Scintillation Resilience**:
   - Enforced automated regression guards evaluating actual estimator output across 6 real and simulated profiles:
     - `check_network_benchmark.py`: Full-day 24h NOAA CORS network (6 bases, 15–50km, network fused $p50 = 3.1\text{ cm}$, $\text{RMS} = 4.8\text{ cm}$).
     - `check_multignss_benchmark.py`: Multi-GNSS GPS+Galileo network ($97.6\%$ network fix rate).
     - `check_f9p_benchmark.py`: Real low-cost u-blox ZED-F9P kinematic tracking ($p50 = 20.3\text{ cm}$, $\text{RMS} = 26.7\text{ cm}$, $88.5\%$ fixed).
     - `check_kinematic_uav_benchmark.py`: High-dynamic circular vehicle ($\text{RMS} = 6.0\text{ mm}$) and open-sky kinematic PPK ($\text{RMS} = 5.9\text{ mm}$).
     - `check_storm_benchmark.py`: Severe cycle slip recovery ($\text{RMS} = 5.9\text{ mm}$) and satellite outage continuity ($\text{RMS} = 7.0\text{ mm}$).
     - `check_mgex_benchmark.py`: Real IGS tracking on Wettzell WTZR ($p50 = 38.8\text{ cm}$, final $dU = 1.5\text{ cm}$) and ALIC ($p50 = 1.05\text{ m}$).
2. **Sprint 33: Tightly-Coupled GNSS/INS Field Validation & Urban Dynamics**:
   - Tuned Non-Holonomic Constraints (NHC) and Zero Velocity Updates (ZUPT) in `crates/gneiss-rtk/src/swfg/pipeline/factors/dynamics.rs`.
   - Verified photogrammetric boresight auto-estimation on multi-pass flight lines.
3. **Sprint 34: Live Hardware-in-the-Loop (HITL), Serial RTCM3/UBX & Edge Execution**:
   - Supported direct serial port binary streams (`/dev/ttyUSB*`, COM) and RTCM3 MSM decoding.
   - High-concurrency async NTRIP v1/v2 client with automatic NMEA GGA feedback for VRS networks.
4. **Sprint 35: Visual GUI Workspace & Telemetry Diagnostics Polish**:
   - Embedded Web diagnostic GUI (`gneiss-cli gui`) with responsive polar skyplot, multi-channel carrier-phase residual charts, and live trajectory overlays.
   - One-click export to KML, GeoJSON, and Applanix SBET trajectory formats directly from UI.
5. **Sprint 36: Total CI Zero-Warning / Zero-Unwrap Cleanliness & Mutation Testing Gate**:
   - 0 compiler warnings, 0 clippy warnings across ALL workspace targets (`cargo clippy --workspace --all-targets -- -D warnings`).
   - All tests passing with 0 unwrap in production code (`unwrap_used = "deny"` enforced).

---

## Architecture Notes

### Frame-safety infrastructure
Three modules provide compile-time prevention of frame-mixing bugs:
- `gnss_time.rs`: TimeSystem enum + GnssTime with to_gpst()/from_gpst()
- `frames.rs`: ReferenceFrame trait + EcefPos<F> newtype + Helmert + EpochPosition<F, R> with 6D tectonic velocity propagation
- `frequencies.rs`: Signal enum + explicit constellation/band mapping

### Code Quality & Modular Decomposition Tracking
Strict enforcement of AGENTS.md rules (< 500 LOC/file, < 32 LOC/func, < 3 nesting, > 95% coverage, 0 mutation survivors with `cargo-mutants`).
100% of production source files in `crates/` and `bin/` are strictly `< 500 LOC`.

## Testing Infrastructure

| suite | count | covers |
|---|---|---|
| gneiss-core lib | 116 | time, frames, frequencies, tides, sun/moon |
| gneiss-parsers lib | 259 | RINEX, SP3, ANTEX, precise_orbit, RTCM3, SBF, UBX, BLQ, SINEX, SBET |
| gneiss-rtk lib | 316 | IEKF, AR, MW, screening, post_process, INS, events, streaming, SWFG, UDUC |
| gneiss-geodesy lib | 10 | projections, NTv2, site calibration, geoid (GTX/BYN) |
| gneiss-fetch lib | 25 | CORS discovery, Hatanaka uncompression |
| gneiss-cli | 8 | process, export, qc, batch, events, calibrate, gui, live |
| workspace integration | 36+ | end-to-end scenarios, benchmark matrices, real IGS CORS (WTZR, ALIC) |
| regression guards | 6 scripts | Datasets A & B + Profiles A (UAV), B (Storm), C (MGEX), D (F9P) |
| walkthrough | 1 binary | bit-identical output verification |

## Sprint 38: Geodetic Normalizations, UDUC Decomposition & Benchmark Integrity Audit (2026-08-31)
- **Geodetic Normalizations Added & Tested**:
  - Relativistic periodic orbit eccentricity range correction ($-2 \mathbf{r}\cdot\mathbf{v}/c$).
  - Gravitational Shapiro time delay range correction.
  - IERS 2010 Solid Earth Tide degree-2/3 elastic crustal deformation.
  - Continuous Wu et al. (1993) RHCP carrier phase windup tracking across $2\pi$ turns.
  - Antenna Reference Point (ARP) to Antenna Phase Center (APC) ANTEX receiver PCO projection.
- **Architectural & Numerical Improvements**:
  - Replaced $5\text{ cm}$ Huber threshold on carrier phase factors with $50\text{ cm}$ unclipped float pull, preventing premature downweighting of initial convergence.
  - Included marginal prior quadratic cost in `evaluate_total_error`, guaranteeing monotonic Levenberg-Marquardt step acceptance consistency.
  - Added persistent `StaticPose` formulation in SWFG for static sessions, accumulating information continuously across sliding-window marginalizations.
  - Modularized factor graph builder into `builder.rs` and `uduc_builder.rs` (< 500 LOC per file, 0 warnings).
- **Benchmark Integrity & Anti-Reward-Hacking Audit**:
  - Removed misleading synthetic unit mock (`truth_station + 2mm`) and renamed benchmarks to accurately state their physical regime.
  - **Real CORS DD-RTK / PPK** (`TMG2-TMGO`): **$p50 = 3.9\text{ mm}$ (39/40 fixed)** — genuine sub-centimeter performance on real raw observations.
  - **Real IGS Float PPP** (`WTZR`): **$p50 = 38.8\text{ cm}$, final $dU = +1.5\text{ cm}$** — standalone float PPP convergence on 600 epochs with CODE SP3/CLK.
  - **Real IGS Perturbed PPP** (`ALIC`): Filter recovers from 3m perturbed seed to **$< 1.50\text{ m}$**.
  - **High-Dynamic Simulations**: $\text{RMS} = 6.0\text{ mm}$ (circular motion), $5.9\text{ mm}$ (cycle slips), $7.0\text{ mm}$ (5s outage).
- **Roadmap to < 5cm Standalone PPP**:
  - Ingestion of IGS SINEX `.bia` satellite Observable-Specific Biases (OSB / FCB) to enable uncorrupted integer ambiguity resolution on standalone carrier phase arcs.
  - Application of $P_1-C_1$ Differential Code Biases (DCB).

## Key Lessons Learned

1. **Guard Against Reward Hacking**: Never name a test `sub_centimeter` if the assertion is `< 3.0m` or if the position was synthetically offset by 2 mm. Tests must test actual estimator output against ground truth.
2. **Measure before building**: every speculative feature was neutral or negative; every measurement-driven change was positive.
3. **TDD catches conceptual errors**: the IF-residual screen tests caught a fundamental misunderstanding of what's cross-pair comparable.
4. **Negative results are valuable**: documented dead ends saved weeks of wasted effort by recording WHY they don't work.
5. **Frame safety matters**: most bugs were missing frame distinctions, not algorithmic errors.
6. **External dependencies dominate**: the remaining gap in standalone PPP requires satellite phase bias (OSB) ingestion, not artificial tuning.
7. **RTKLIB is a floor, not a ceiling**: beating it proves the core is sound; exceeding it requires adopting techniques from commercial-grade implementations.
