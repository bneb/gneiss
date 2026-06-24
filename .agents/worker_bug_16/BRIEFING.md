# BRIEFING — 2026-06-22T01:04:30Z

## Mission
Implement the fix for Bug 16: Mismatched Galileo BGD Correction.

## 🔒 My Identity
- Archetype: worker
- Roles: implementer, qa, specialist
- Working directory: /Users/kevin/projects/gneiss/.agents/worker_bug_16
- Original parent: c1e1438e-1aa8-4425-a0b8-6dc0b1d37f21
- Milestone: Bug 16 Fix

## 🔒 Key Constraints
- CODE_ONLY network mode. No external HTTP requests.
- DO NOT CHEAT. All implementations must be genuine. Do not hardcode test results.
- Minimal change principle. Only modify what is necessary.

## Current Parent
- Conversation ID: c1e1438e-1aa8-4425-a0b8-6dc0b1d37f21
- Updated: not yet

## Task Summary
- **What to build**: Define `position_e5b` on `Ephemeris` to delegate to `GalileoEphemeris::position_e5b` when constellation is Galileo, else `position`. Handle band 7 in `get_frequency`. Add `freq_band` to `SppMeasurement`. Set `freq_band` in `build_single_measurement`. Update `compute_sat_state` and `process_measurement` to call `m.eph.position_e5b` if `m.freq_band == 7`. Add regression unit test.
- **Success criteria**: All code compiles, tests pass, regression test passes and covers the fix.
- **Interface contracts**: Rust crates `gneiss-core` and `gneiss-rtk`.
- **Code layout**: Source files and tests inside crates.

## Key Decisions Made
- [TBD]

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/worker_bug_16/handoff.md - Handoff report documenting the work.
