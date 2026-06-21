# BRIEFING — 2026-06-20T20:17:30-07:00

## Mission
Investigate Bug 17: GLONASS Time Scale Discrepancy in gneiss-core.

## 🔒 My Identity
- Archetype: explorer
- Roles: Teamwork explorer
- Working directory: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_glonass
- Original parent: f16afb25-c177-42fe-985d-6840e173046f
- Milestone: milestone_1

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- CODE_ONLY network mode: MUST NOT access external websites/services, MUST NOT use run_command to run curl/wget/lynx.

## Current Parent
- Conversation ID: f16afb25-c177-42fe-985d-6840e173046f
- Updated: 2026-06-20T20:17:30-07:00

## Investigation State
- **Explored paths**:
  - `crates/gneiss-parsers/src/rinex.rs` (navigation parsing, tests)
  - `crates/gneiss-core/src/time.rs` (GPS time representation, calendar conversions)
  - `crates/gneiss-core/src/ephemeris.rs` (GLONASS orbit propagation)
  - `crates/gneiss-rtk/src/estimators/spp.rs` (Single Point Positioning measurement residuals)
- **Key findings**:
  - Found a 3-hour offset discrepancy in `rinex.rs` line 688: GLONASS navigation epochs are in GLONASST (Moscow Time = UTC+3h), but were parsed directly as GPST and only adjusted by adding leap seconds (+18.0s), missing the -3h offset.
  - The correct GPST conversion is `GPST = GLONASST - 3 hours + 18 seconds = GLONASST - 10782 seconds`.
  - Propagating the GLONASS orbit with a 3-hour timing error results in a satellite coordinate error of ~27,000 km, preventing convergence or usage of GLONASS in RTK/SPP.
  - The unit test `test_parse_rinex_3_nav_date` was written to assert the incorrect TOW of 422118.0 instead of 411318.0.
- **Unexplored areas**: None, the root cause has been precisely pinpointed and analyzed.

## Key Decisions Made
- Conformed to the read-only constraint: documented the findings and a complete, precise fix strategy in `handoff.md` without editing codebase source files.
- Recommended using `GpsTime + f64` addition for the fix to ensure robust week-boundary normalization.

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_glonass/handoff.md — Detailed findings, logic chain, and proposed fix strategy.
