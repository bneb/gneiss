## 2026-06-21T04:39:48Z
You are teamwork_preview_worker. Your working directory is /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_worker_vel_att.
Your mission is to implement the fix for Bug 2: Velocity-Attitude Transition Sign Mismatch.
Please follow the proposed fix strategy:
1. In `crates/gneiss-rtk/src/engine/predictor.rs`:
   Negate the attitude-to-velocity coupling block computation:
   ```rust
   let vel_att = -f_e_skew * dt;
   ```
2. Add the unit test `test_transition_matrix_velocity_attitude_coupling` in `crates/gneiss-rtk/src/engine/tests_predictor.rs` as suggested. This unit test acts as a regression test: it passes with your fix but fails (the computed coupling has positive sign instead of negative) if the buggy code is restored.
3. Verify that the build succeeds and tests pass by running:
   - `cargo test -p gneiss-rtk`
   - `cargo test --workspace`

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A Forensic Auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Write your changes and verification details in `changes.md` in your working directory and write a `handoff.md` with:
- Description of the fix
- Compilation and test results (with the exact cargo command and output)
- Verification that layout and regression tests behave correctly.

When complete, send a message to conversation ID f16afb25-c177-42fe-985d-6840e173046f.
