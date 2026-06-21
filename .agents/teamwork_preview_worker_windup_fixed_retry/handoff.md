# Handoff Report - Phase Wind-up Correction Sign Fix (Bug 18)

## 1. Observation
The following observations were made on the codebase and git status:
- In `crates/gneiss-rtk/src/engine/measurement.rs`:
  - `apply_windup_to_obs` originally added the wind-up correction:
    ```rust
    fn apply_windup_to_obs(obs: &mut DdObservation, windup: f64) {
        if let Some(cp) = &mut obs.cp_l1 {
            *cp += windup;
        }
        if let Some(cp2) = &mut obs.cp_l2 {
            *cp2 += windup;
        }
    }
    ```
- In `crates/gneiss-rtk/src/engine/ppp_iekf.rs`:
  - `push_cp_measurement` originally added the windup:
    ```rust
    let l_meas = (cp1 + windup) * sat.lam1;
    ```
  - UDUC carrier phase residual calculations added the windup:
    ```rust
    let res_l1 = (sat.cp1.unwrap() + windup) * sat.lam1 - (expected_base - i1 + n1);
    ...
    let res_l2 = (sat.cp2.unwrap() + windup) * sat.lam2 - (expected_base - gamma * i1 + n2);
    ```
- In `crates/gneiss-rtk/src/engine/ppp.rs`:
  - Git history showed that commit `61bf5a63502483aea209ac5a3898ca67b9a4c941` ("fix(ppp, atmosphere): Bugs 18, 6, 12, 25 + regression tests") had already correctly updated `ppp.rs` call sites to subtract `wup` (`cp1 - wup`, `sat.cp2.unwrap() - wup`).
- Run `cargo test -p gneiss-rtk` compiled and ran all tests successfully after making the corrections.
- Reverting the sign change in `apply_windup_to_obs` to `+=` caused the newly added unit test `test_phase_windup_correction_sign_rtk` to fail with:
  ```
  assertion `left == right` failed
    left: 10.25
   right: 9.75
  ```

## 2. Logic Chain
1. **Mathematical Correction**: The phase wind-up correction `wup` represents the physical rotation of the effective phase in cycles. Correcting the carrier phase observation to obtain the physical distance requires subtracting the wind-up correction: `corrected_cp = raw_cp - windup`. Adding the wind-up correction doubles the sign error.
2. **Implementation Errors**: In `measurement.rs` and `ppp_iekf.rs`, the windup was incorrectly added. Changing the addition operators (`+=` or `+`) to subtraction (`-=` or `-`) corrects the mathematical definition.
3. **Regression Proof**: The added unit test `test_phase_windup_correction_sign_rtk` verifies that `apply_windup_to_obs` correctly subtracts the phase wind-up. When the fix is reverted, the test fails, confirming that the test is sensitive to the bug.

## 3. Caveats
No caveats. All proposed sites in `measurement.rs`, `ppp_iekf.rs`, and `ppp.rs` have been addressed and validated.

## 4. Conclusion
The phase wind-up correction sign fix (Bug 18) has been fully implemented across `measurement.rs` and `ppp_iekf.rs`, and verified in `ppp.rs`. The code builds successfully, formatting check passes, and all 258 unit tests in `gneiss-rtk` plus all workspace integration tests pass cleanly.

## 5. Verification Method
To independently verify the changes, run:
1. **Clean Workspace and Run Tests**:
   ```bash
   cargo test -p gneiss-rtk
   cargo test --workspace
   ```
   All tests should compile and pass cleanly.
2. **Inspect the Unit Test**:
   Open `crates/gneiss-rtk/src/engine/measurement.rs` and inspect `test_phase_windup_correction_sign_rtk` under `mod tests`.
3. **Format Check**:
   ```bash
   cargo fmt --check
   ```
   Should exit with status 0.
