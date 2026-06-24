## 2026-06-21T18:00:07Z
You are a Worker agent. Your task is to implement the fix for Bug 15: Incorrect Broadcast Clock TGD Correction.
Working directory: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_worker_bug15_2

Please read the Explorer analysis report at `/Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_bug15_2/analysis.md` for full context, files, and formulas.

Tasks:
1. Add `pub tgd2: f64` to `BeidouEphemeris` in `crates/gneiss-core/src/ephemeris.rs`.
2. Fix `bgd_e5b` in `Ephemeris` in `crates/gneiss-core/src/ephemeris.rs` to return `tgd2` for Beidou.
3. Update `build_beidou_ephemeris` in `crates/gneiss-parsers/src/rinex.rs` to map `tgd2` from `vals[23]` and `aodc` from `vals[25] as u32` (instead of mapping `vals[23]` to `aodc` and discarding the real `aodc`).
4. In `BeidouEphemeris::position_iono_free` in `crates/gneiss-core/src/ephemeris.rs`, calculate the combined group delay $T_{GD\_IF}$ using the formula:
   tgd_if = (f1_sq * self.tgd1 - f2_sq * self.tgd2) / (f1_sq - f2_sq)
   where f1 = gneiss_core::signal::FREQ_BDS_B1I and f2 = gneiss_core::signal::FREQ_GAL_E5B.
   Pass this `tgd_if` as the `tgd` argument to `calc_keplerian`.
5. Update all occurrences and instantiations of `BeidouEphemeris` in the project's codebase (such as in tests in `crates/gneiss-core/src/ephemeris.rs` and other crates) to include `tgd2` initialization.
6. Correct the assertion in `test_broadcast_clock_tgd_correct` in `crates/gneiss-core/src/ephemeris.rs` to verify that `clk_if - clk_pos = tgd1 - tgd_if`.
7. Verify your implementation by running build and test commands (e.g. `cargo build`, `cargo test --workspace`) and ensure all tests pass.

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A Forensic Auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Write your changes and verification outcomes to `changes.md` and a final handoff report to `handoff.md` in your working directory. Report back with the paths to these files when done.
