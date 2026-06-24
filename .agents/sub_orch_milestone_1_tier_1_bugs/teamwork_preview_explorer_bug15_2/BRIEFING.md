# BRIEFING — 2026-06-21T20:06:04Z

## Mission
Investigate Bug 15 (Incorrect Broadcast Clock TGD Correction) in gneiss-core and propose a fix.

## 🔒 My Identity
- Archetype: Explorer
- Roles: Teamwork explorer
- Working directory: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_bug15_2
- Original parent: 2edbcec8-b8bb-45cc-b32c-9a9af7206f15
- Milestone: Bug 15: Incorrect Broadcast Clock TGD Correction

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- Code-only network mode (no external HTTP calls)
- Follow Handoff Protocol with 5-component report

## Current Parent
- Conversation ID: e2192921-f0a3-4e28-a7cb-96bd652471d1
- Updated: 2026-06-21T20:08:45Z

## Investigation State
- **Explored paths**:
  - `crates/gneiss-core/src/ephemeris.rs` - Location of BeidouEphemeris and TGD/BGD correction logic.
  - `crates/gneiss-parsers/src/rinex.rs` - Location of RINEX parser for Beidou and GPS ephemerides.
  - `crates/gneiss-rtk/src/estimators/spp.rs` - SPP measurement and dual-frequency/iono-free logic.
  - `crates/gneiss-rtk/src/engine/ppp.rs` - PPP measurement processing and iono-free clock logic.
- **Key findings**:
  1. `BeidouEphemeris` struct lacks a `tgd2` field.
  2. RINEX navigation parser (`rinex.rs`) incorrectly parses `TGD2` (at `vals[23]`) into `aodc` for Beidou, and completely ignores the actual `aodc` (at `vals[25]`).
  3. `BeidouEphemeris::position_iono_free` passes `0.0` as `tgd` to `calc_keplerian`, which is incorrect because BDS legacy clock parameters are referenced to the B3I frequency, not the ionosphere-free combination. For BDS B1I/B2I ionosphere-free dual frequency combination, the correct clock correction is $a_{f0} - T_{GD\_IF}$, where $T_{GD\_IF} = \frac{f_{B1I}^2 T_{GD1} - f_{B2I}^2 T_{GD2}}{f_{B1I}^2 - f_{B2I}^2}$.
  4. The Beidou test in `test_broadcast_clock_tgd_correct` expects `clk_if - clk_pos = tgd1` because it wrongly assumed `position_iono_free` should subtract `0.0` for BDS.
- **Unexplored areas**: None. The bug is fully located and characterized.

## Key Decisions Made
- Identified root cause in `BeidouEphemeris` fields, `rinex.rs` parsing offsets, and `position_iono_free` math.
- Formulated fix strategy including modifying the parser, adding `tgd2` field, correcting the iono-free TGD formula for Beidou, and updating the regression tests.

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_bug15_2/analysis.md — Final analysis report and proposed fix strategy
