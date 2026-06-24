# BRIEFING — 2026-06-22T01:03:42Z

## Mission
Investigate Galileo BGD correction mismatch where bgd_e1_e5b should be used for band 7 (E5b) observations instead of bgd_e1_e5a.

## 🔒 My Identity
- Archetype: teamwork_preview_explorer
- Roles: Teamwork explorer
- Working directory: /Users/kevin/projects/gneiss/.agents/explorer_bug_16
- Original parent: c1e1438e-1aa8-4425-a0b8-6dc0b1d37f21
- Milestone: explorer_bug_16

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- CODE_ONLY network mode: no external web access

## Current Parent
- Conversation ID: c1e1438e-1aa8-4425-a0b8-6dc0b1d37f21
- Updated: 2026-06-22T01:03:42Z

## Investigation State
- **Explored paths**:
  - `crates/gneiss-core/src/ephemeris.rs`: Definition of `GalileoEphemeris::position_e5b` and `Ephemeris::tgd`/`bgd_e5b`.
  - `crates/gneiss-core/src/signal.rs`: Mapping of frequency bands and constants.
  - `crates/gneiss-rtk/src/estimators/spp.rs`: Measurement construction and clock/position calculations.
  - `crates/gneiss-rtk/src/engine/spp_tight.rs`: Tight single-point positioning calculations.
- **Key findings**:
  - `GalileoEphemeris::position_e5b` is implemented but never called in the engine or estimators.
  - `Ephemeris::position` always delegates to constellation-specific `position`, which for Galileo uses `bgd_e1_e5a`.
  - `SppMeasurement` lacks a `freq_band` field, meaning the estimators cannot differentiate between Galileo E1 (band 1) and E5b (band 7) single-frequency observations.
  - `get_frequency` in `signal.rs` fails to map band 7 to Galileo E5b frequency, defaulting to GPS L1 instead.
- **Unexplored areas**: None, the bug boundary has been fully explored and defined.

## Key Decisions Made
- Scoped fix specifically to Galileo E5b (band 7) single-frequency observations.
- Drafted a diff patch `fix_bug_16.patch` containing all required updates.

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/explorer_bug_16/handoff.md — Handoff report containing observations, logic chain, caveats, conclusion, and verification method.
- /Users/kevin/projects/gneiss/.agents/explorer_bug_16/fix_bug_16.patch — Precise, machine-applicable diff patch file.
