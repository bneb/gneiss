# Plan - Phase Wind-up correction sign fix (Bug 18)

## Steps

1. **Verify original behaviour**:
   - Locate the functions/lines to be changed.
   - Run existing tests to ensure clean workspace.

2. **Implement sign corrections in `crates/gneiss-rtk/src/engine/measurement.rs`**:
   - Locate `apply_windup_to_obs` (around lines 42-49).
   - Change `*cp += windup;` and `*cp2 += windup;` to `-=`.
   - Add unit test `test_phase_windup_correction_sign_rtk` under `mod tests`.

3. **Implement sign corrections in `crates/gneiss-rtk/src/engine/ppp_iekf.rs`**:
   - Locate `predict_carrier_phase` or where carrier phase predictions are made (line 904). Change `cp1 + windup` to `cp1 - windup`.
   - Locate UDUC phase measurement residual construction (lines 1032, 1049). Change `sat.cp1.unwrap() + windup` and `sat.cp2.unwrap() + windup` to subtract `windup`.

4. **Verify correctness**:
   - Run `cargo test -p gneiss-rtk` and verify it passes.
   - Run `cargo test --workspace` and verify all tests pass.
   - Run `cargo fmt --check`.

5. **Regression Verification**:
   - Revert the changes temporarily and ensure that the added unit test fails.
   - Re-apply the changes and ensure it passes.

6. **Documentation and Handoff**:
   - Document changes in `changes.md`.
   - Write a detailed handoff report in `handoff.md`.
   - Send completion message to parent agent.
