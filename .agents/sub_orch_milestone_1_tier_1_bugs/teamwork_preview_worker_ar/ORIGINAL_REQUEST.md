## 2026-06-20T21:29:19-07:00
You are teamwork_preview_worker. Your working directory is /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_worker_ar.
Your mission is to implement the fix for Bug 9: Sequential AR Covariance Mismatch.
Please follow the proposed fix strategy:
1. In `crates/gneiss-rtk/src/engine/ppp_iekf.rs`:
   - Modify the signature of `resolve_widelane_ar` to accept `p: &DMatrix<f64>` as a parameter.
   - Replace all usages of `state.covariance` with `p` inside `resolve_widelane_ar`.
   - In `resolve_narrowlane_ar`, replace `state.covariance.nrows()` with `p_wl.nrows()`.
   - In `resolve_cascade_ar`:
     - Remove the call to `apply_state_vector(state, &x_current, p_current.clone());` inside the `for (constellation, group_cands) in &const_groups` loop.
     - Update the calls to `resolve_widelane_ar` inside the loop to pass `&p_current`.
     - Update the call to `resolve_widelane_ar` in the fallback block to pass `&p_current`.
     - Call `apply_state_vector` only at the end of the method after all validations (including global position validation) pass.
2. Implement a unit test in `crates/gneiss-rtk/src/engine/ppp_iekf.rs` (or `crates/gneiss-rtk/tests/`) that mocks sequential AR succeeding but final global validation failing (e.g. producing position jump > 20m). Assert that:
   - `resolve_cascade_ar` returns `Err`.
   - The input `state` covariance and state vector remain completely unchanged.
   This unit test acts as a regression test: it passes with your fix but fails (the state's covariance is mutated) if the buggy code is restored.
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
