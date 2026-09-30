## 2026-09-25T21:13:05Z

You are Worker M2 for the Gneiss Frame Safety & Epoch Alignment Refactoring.
Your working directory is: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_worker_m2
Authoritative request file: /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md
Survey findings: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_2/survey_r2_temporal.md and /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_2/handoff.md
Project specification: /Users/kevin/projects/gneiss/PROJECT.md

MANDATORY: You MUST read /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md and /Users/kevin/projects/gneiss/AGENTS.md before starting work.

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Scope & File Boundaries (M2 - Features 6–10):
You exclusively own:
- crates/gneiss-core/src/time/ (e.g. epoch.rs, scales.rs, mod.rs, time.rs)
- crates/gneiss-parsers/src/rinex/nav/builder.rs (BeiDou scale alignment)
- crates/gneiss-parsers/src/antex.rs (leap second fix)
- crates/gneiss-fetch/src/sources/bkg.rs (leap second fix)
- Unit and compile-fail tests in crates/gneiss-core/tests/ (e.g. tests/test_epoch_safety.rs)

Deliverables:
1. Formalize typed epoch systems:
   - TimeScale trait and markers: GpsScale, BdtScale, GstScale, GlonassScale, UtcScale
   - Epoch<Scale: TimeScale> with integer nanoseconds (week: u32, tow_nanos: u64)
   - TimeDelta { nanos: i64 }
   - Week-rollover-safe epoch subtraction and comparison
2. Disallow cross-scale subtraction at compile time (Epoch<GpsScale> - Epoch<BdtScale> fails to compile).
3. Fix BeiDou broadcast ephemeris scale discrepancy in rinex/nav/builder.rs and keplerian.rs so toe and toc are in consistent time systems, eliminating the 14-second satellite clock offset error.
4. Fix leap-second omissions in antex.rs and gneiss-fetch to explicitly account for 18 leap seconds when converting to UTC/Unix.
5. Provide EpochKey and tolerance-based matching is_within(tolerance).
6. Comprehensive unit tests and compile-fail tests verifying:
   - Cross-system subtraction fails at compile time
   - Subtraction across week rollover computes exact elapsed seconds
   - BeiDou ephemeris clock bias is correct
7. Verify compliance with AGENTS.md:
   - File size strictly < 500 LOC
   - Function size strictly <= 32 LOC
   - Nesting depth strictly < 3
   - 0 compiler and clippy warnings under cargo clippy --workspace --all-targets -- -D warnings
   - 0 unwrap() in production code
   - All tests pass: cargo test --workspace
8. Write your handoff to handoff.md and send a completion message when done.
