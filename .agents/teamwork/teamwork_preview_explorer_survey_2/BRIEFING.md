# BRIEFING — 2026-09-25T21:12:00Z

## Mission
Survey Gneiss codebase for R2: Temporal Frame & Epoch Alignment Safety, leap seconds, epoch matching, and propose strictly typed epoch architecture.

## 🔒 My Identity
- Archetype: explorer
- Roles: investigator, surveyor
- Working directory: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_2
- Original parent: db66ae0c-b21b-4e14-ac97-93509c51c4b0
- Milestone: Survey 2 - Temporal Frame & Epoch Alignment Safety (R2)

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- Never read entire files into context — grep first, read line slices
- Target <500 LOC files, <=32 LOC functions, 0 unwraps in production, frame safety guarantees

## Current Parent
- Conversation ID: db66ae0c-b21b-4e14-ac97-93509c51c4b0
- Updated: 2026-09-25T21:12:00Z

## Investigation State
- **Explored paths**:
  - `crates/gneiss-core/src/time.rs`, `gnss_time.rs`, `ephemeris/keplerian.rs`, `keplerian.rs`, `glonass.rs`, `obs.rs`, `imu.rs`
  - `crates/gneiss-parsers/src/rinex/nav/`, `rinex/obs/`, `rtcm3/msm/`, `sp3.rs`, `antex.rs`, `csrs_pos.rs`
  - `crates/gneiss-fetch/src/sources/{bkg, noaa, cddis}.rs`
  - `crates/gneiss-rtk/src/post_process/{network, iekf_pass, combiner, backward}.rs`, `streaming.rs`, `estimators/rtk_iekf/predict.rs`, `swfg/engine/epoch.rs`
  - `crates/gneiss-rtk/src/bin/eval_odaiba_ins/` (`main.rs`, `odaiba_helpers.rs`), `eval_f9p_rover.rs`, `eval_qinertia_ppk.rs`, `eval_network_ppk.rs`, `eval_ppp.rs`, `eval_odaiba.rs`
- **Key findings**:
  - Discovered 14s satellite clock offset bug in `BeidouEphemeris` due to mixing BDT `toe` with GPST `toc`.
  - Discovered 16 sites performing direct `tow - tow` float subtractions that fail on week rollover (evaluating to `-604799s`).
  - Discovered leap-second omissions in `antex.rs` and `gneiss-fetch` converting `GpsTime` to UTC/Unix without subtracting 18s.
  - Cataloged 16 ad-hoc float rounding / truncation sites for epoch matching.
  - Analyzed `eval_odaiba_ins` IMU/GNSS sync and found fragile exact integer millisecond base station lookup.
  - Designed zero-cost typestate `Epoch<Scale: TimeScale>` and integer nanosecond `TimeDelta` architecture.
- **Unexplored areas**: None (R2 survey scope complete).

## Key Decisions Made
- Completed full R2 temporal survey and authored `survey_r2_temporal.md`.
- Authored self-contained `handoff.md` (Hard Handoff).

## Artifact Index
- DISPATCH.md — record of incoming dispatch instructions
- plan.md — plan summary
- progress.md — liveness and completion checklist
- survey_r2_temporal.md — comprehensive technical report for Requirement R2
- handoff.md — self-contained 5-component handoff report for the orchestrator
