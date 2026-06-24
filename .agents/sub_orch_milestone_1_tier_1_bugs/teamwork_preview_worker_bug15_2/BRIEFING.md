# BRIEFING — 2026-06-21T18:00:07-07:00

## Mission
Implement the fix for Bug 15: Incorrect Broadcast Clock TGD Correction.

## 🔒 My Identity
- Archetype: implementer, qa, specialist
- Roles: implementer, qa, specialist
- Working directory: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_worker_bug15_2
- Original parent: 2edbcec8-b8bb-45cc-b32c-9a9af7206f15
- Milestone: Bug 15 Fix

## 🔒 Key Constraints
- CODE_ONLY network mode: no external web access, no curl/wget/etc.
- Follow Teamwork rules, do not cheat, no dummy implementations.
- Write only to own folder for agent metadata, read any folder.

## Current Parent
- Conversation ID: 2edbcec8-b8bb-45cc-b32c-9a9af7206f15
- Updated: yes

## Task Summary
- **What to build**: Add `tgd2` to `BeidouEphemeris`, fix `bgd_e5b` and `build_beidou_ephemeris`, calculate iono-free position combined group delay correctly, and update tests/assertions.
- **Success criteria**: All tests build and pass successfully.
- **Interface contracts**: `crates/gneiss-core/src/ephemeris.rs`, `crates/gneiss-parsers/src/rinex.rs`.
- **Code layout**: Gneiss cargo workspace structure.

## Key Decisions Made
- Added `tgd2: f64` to `BeidouEphemeris`.
- Mapped `tgd2` from `vals[23]` and `aodc` from `vals[25]` in RINEX navigation message parser.
- Corrected Beidou iono-free combined group delay $T_{GD\_IF}$ calculation.
- Updated all test instantiations and corrected test assertions.

## Artifact Index
- None.

## Change Tracker
- **Files modified**:
  - `crates/gneiss-core/src/ephemeris.rs`: Add field `tgd2` to `BeidouEphemeris`, update `bgd_e5b`, calculate `tgd_if` correctly in `position_iono_free`, update test cases and correct assertion.
  - `crates/gneiss-parsers/src/rinex.rs`: Correctly map `tgd2` and `aodc` in `build_beidou_ephemeris`.
  - `crates/gneiss-rtk/src/estimators/spp.rs`: Add `tgd2: 0.0` to Beidou ephemeris test instantiation.
- **Build status**: Untested since modifications.
- **Pending issues**: Run build and test to verify.

## Quality Status
- **Build/test result**: Untested.
- **Lint status**: Untested.
- **Tests added/modified**: Corrected `test_broadcast_clock_tgd_correct` in `ephemeris.rs`.

## Loaded Skills
- None.
