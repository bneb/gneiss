# BRIEFING — 2026-09-25T21:13:15Z

## Mission
Implement Worker M2 deliverables (Features 6–10): typed epoch systems, compile-time cross-scale prevention, BeiDou broadcast ephemeris scale alignment, leap second fixes, EpochKey matching, and unit/compile-fail tests adhering strictly to AGENTS.md rules.

## 🔒 My Identity
- Archetype: Worker
- Roles: implementer, qa, specialist
- Working directory: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_worker_m2
- Original parent: db66ae0c-b21b-4e14-ac97-93509c51c4b0
- Milestone: M2 - Features 6–10 (Epoch Safety & Alignment)

## 🔒 Key Constraints
- File size strictly < 500 LOC
- Function size strictly <= 32 LOC
- Nesting depth strictly < 3 levels
- 0 compiler and clippy warnings under cargo clippy --workspace --all-targets -- -D warnings
- 0 unwrap() in production code (unwrap() in #[cfg(test)] is acceptable)
- No dead code, no duplicate functions
- Tests co-located in source or tests/
- Own files: crates/gneiss-core/src/time/, crates/gneiss-parsers/src/rinex/nav/builder.rs, crates/gneiss-parsers/src/antex.rs, crates/gneiss-fetch/src/sources/bkg.rs, crates/gneiss-core/tests/
- Genuine implementation only, no cheating or facades

## Current Parent
- Conversation ID: db66ae0c-b21b-4e14-ac97-93509c51c4b0
- Updated: not yet

## Task Summary
- **What to build**:
  1. Typed epoch systems (TimeScale trait, GpsScale, BdtScale, GstScale, GlonassScale, UtcScale, Epoch<Scale>, TimeDelta).
  2. Disallow cross-scale subtraction at compile time.
  3. Fix BeiDou broadcast ephemeris scale discrepancy in rinex/nav/builder.rs (and keplerian.rs if within scope or interface).
  4. Fix leap second handling in antex.rs and gneiss-fetch/sources/bkg.rs.
  5. EpochKey and tolerance-based matching is_within(tolerance).
  6. Unit tests and compile-fail tests (e.g. trybuild or similar compile-fail tests).
- **Success criteria**: All tests pass, 0 warnings, clean clippy, compliant LOC and nesting.
- **Interface contracts**: PROJECT.md, survey_r2_temporal.md
- **Code layout**: crates/gneiss-core/src/time/, crates/gneiss-parsers/, crates/gneiss-fetch/

## Key Decisions Made
- [TBD]

## Artifact Index
- DISPATCH.md — Assignment instructions
- BRIEFING.md — Working memory index
- progress.md — Heartbeat and progress tracking

## Change Tracker
- **Files modified**: None yet
- **Build status**: Untested
- **Pending issues**: None

## Quality Status
- **Build/test result**: Not yet run
- **Lint status**: Clean
- **Tests added/modified**: None yet

## Loaded Skills
- None
