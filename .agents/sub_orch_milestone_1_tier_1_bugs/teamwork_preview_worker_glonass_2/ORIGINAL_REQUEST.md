## 2026-06-21T03:25:28Z
You are teamwork_preview_worker. Your working directory is /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_worker_glonass_2.
Your mission is to implement the fix for Bug 17: GLONASS Time Scale Discrepancy.
Please follow the proposed fix strategy:
1. In `crates/gneiss-parsers/src/rinex.rs`:
   Modify GLONASS epoch time conversion to GPST to subtract the 3-hour Moscow Time offset:
   ```rust
   let mut toc_gpst = parse_rinex_nav_epoch_time(&line, is_rinex_3);
   match current_constellation {
       Constellation::Glonass => toc_gpst = toc_gpst + (18.0 - 10800.0),
       Constellation::Beidou => toc_gpst = toc_gpst + 14.0,
       _ => {}
   }
   ```
2. Update the unit test `test_parse_rinex_3_nav_date` in `crates/gneiss-parsers/src/rinex.rs` to assert the corrected TOW (`411318.0`). This test acts as a regression test: it passes with the fix but would fail if the buggy code is restored.
3. Run `cargo test -p gneiss-parsers` to verify the fix and tests pass.
4. Run `cargo test --workspace` to ensure no other tests are broken.

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A Forensic Auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Write your changes and verification details in `changes.md` in your working directory and write a `handoff.md` with:
- Description of the fix
- Compilation and test results (with the exact cargo command and output)
- Verification that layout and regression tests behave correctly.

When complete, send a message to conversation ID f16afb25-c177-42fe-985d-6840e173046f.
