# BRIEFING — 2026-06-21T20:06:04Z

## Mission
Investigate Bug 15: Incorrect Broadcast Clock TGD Correction in gneiss-core.

## 🔒 My Identity
- Archetype: Explorer
- Roles: Read-only investigator
- Working directory: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_bug15_3
- Original parent: 2edbcec8-b8bb-45cc-b32c-9a9af7206f15
- Milestone: milestone_1_tier_1_bugs

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- Limit edits to reports and analysis files inside working directory

## Current Parent
- Conversation ID: 2edbcec8-b8bb-45cc-b32c-9a9af7206f15
- Updated: 2026-06-21T20:09:35Z

## Investigation State
- **Explored paths**:
  - `crates/gneiss-core/src/ephemeris.rs` - Contains BeidouEphemeris and Keplerian calculation methods.
  - `crates/gneiss-parsers/src/rinex.rs` - Location of `build_beidou_ephemeris`.
  - `crates/gneiss-rtk/src/estimators/spp.rs` - Instantiations of BeidouEphemeris in SPP tests.
- **Key findings**:
  1. `BeidouEphemeris` is missing `tgd2` field.
  2. RINEX parser incorrectly parses TGD2 into `aodc`, and ignores actual `AODC`.
  3. `BeidouEphemeris::position_iono_free()` passes `0.0` as `tgd` to Keplerian calculation, missing correct B1I/B2I iono-free Timing Group Delay correction.
- **Unexplored areas**: None.

## Key Decisions Made
- Confirmed correct B1I/B2I iono-free TGD formula for Beidou.
- Planned updates to BeidouEphemeris struct, RINEX parser, position_iono_free calculation, and tests.

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_bug15_3/analysis.md — Bug 15 Analysis and Proposed Fix Strategy
- /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_bug15_3/handoff.md — Handoff Report
