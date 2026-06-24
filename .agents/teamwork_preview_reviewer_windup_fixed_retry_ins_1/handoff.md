# Handoff Report: Review and Verification of Bug 18 Fix (Phase Wind-Up Correction Sign)

## Review Summary

**Verdict**: APPROVE

All reviewed implementation changes correctly, robustly, and completely resolve the phase wind-up sign direction error across the GNSS RTK and INS EKF/PPP modules. The changes follow codebase style conventions and are fully supported by regression testing.

---

## Findings

### [Minor] Finding 1: Code Formatting Discrepancies
- **What**: `cargo fmt --check` failed due to formatting discrepancies.
- **Where**: Unrelated files in the repository (e.g., `crates/gneiss-core/src/lib.rs`, `crates/gneiss-rtk/src/engine/ambiguity.rs`, `crates/gneiss-rtk/src/engine/kinematics.rs`).
- **Why**: Formatting check failures blocks CI/CD pipelines enforcing formatting.
- **Suggestion**: Run `cargo fmt` at the workspace level in a separate chore commit. No action required for the modified files for Bug 18, as they are fully formatted.

---

## Verified Claims

- **Claim 1**: Phase wind-up correction has been changed from addition to subtraction in `apply_windup_to_obs` inside `crates/gneiss-rtk/src/engine/measurement.rs`.
  - **Verification Method**: Verified via `view_file` on `crates/gneiss-rtk/src/engine/measurement.rs`, lines 42-49:
    ```rust
    fn apply_windup_to_obs(obs: &mut DdObservation, windup: f64) {
        if let Some(cp) = &mut obs.cp_l1 {
            *cp -= windup;
        }
        if let Some(cp2) = &mut obs.cp_l2 {
            *cp2 -= windup;
        }
    }
    ```
  - **Status**: PASS

- **Claim 2**: Corrected carrier phase is constructed by subtracting wind-up in `ppp_iekf.rs`.
  - **Verification Method**: Verified via `view_file` on `crates/gneiss-rtk/src/engine/ppp_iekf.rs` at line 899:
    ```rust
    let l_meas = (cp1 - windup) * sat.lam1;
    ```
    And lines 1027 and 1044:
    ```rust
    let res_l1 = (sat.cp1.unwrap() - windup) * sat.lam1 - (expected_base - i1 + n1);
    ...
    let res_l2 = (sat.cp2.unwrap() - windup) * sat.lam2 - (expected_base - gamma * i1 + n2);
    ```
  - **Status**: PASS

- **Claim 3**: Phase wind-up subtraction is implemented in tightly-coupled INS EKF residuals in `ppp_ins_iekf.rs`.
  - **Verification Method**: Verified via `view_file` on `crates/gneiss-rtk/src/engine/ppp_ins_iekf.rs` at line 928:
    ```rust
    let l_meas = (cp1 - windup) * sat.lam1;
    ```
    And lines 1042 and 1061:
    ```rust
    let res_l1 = (sat.cp1.unwrap() - windup) * sat.lam1 - (expected_base - idx.i1 + idx.n1);
    ...
    let res_l2 = (sat.cp2.unwrap() - windup) * sat.lam2 - (expected_base - idx.gamma * idx.i1 + idx.n2);
    ```
  - **Status**: PASS

- **Claim 4**: A unit test is added to verify that positive wind-up reduces the corrected carrier phase.
  - **Verification Method**: Verified test `test_phase_windup_correction_sign_rtk` in `measurement.rs` (lines 976-1004) and `test_windup_sign_correct` in `ppp.rs` (lines 1119-1143) via `view_file` and ran `cargo test -p gneiss-rtk`.
  - **Status**: PASS

---

## Coverage Gaps
- **No Coverage Gaps**: Evaluated other satellite range correction modules (PCV, solid earth tides) and confirmed that phase wind-up logic is restricted to the EKF measurement models covered in this fix. Risk level: LOW.

---

## Unverified Items
- None.

---

## Challenge Summary

**Overall risk assessment**: LOW

Active stress-testing of the phase wind-up correction was performed. Since the sign of geometric phase wind-up depends on the relative geometry/orientation between the transmitter and receiver antenna patterns, a wrong sign (addition instead of subtraction) acts as a positive feedback mechanism that doubles the error rather than eliminating it. Subtracting the wind-up correction ensures negative feedback/error cancellation.

---

## Challenges

### [Low] Challenge 1: Single-Frequency or Non-GPS constellations compatibility
- **Assumption challenged**: Whether phase wind-up corrections perform robustly under single-frequency or alternate constellations when `cp2` is absent.
- **Attack scenario**: `obs.cp2` is `None`.
- **Blast radius**: None. The `apply_windup_to_obs` implementation uses safe `if let Some(cp)` guards:
  ```rust
  if let Some(cp) = &mut obs.cp_l1 { *cp -= windup; }
  if let Some(cp2) = &mut obs.cp_l2 { *cp2 -= windup; }
  ```
  This is fully robust to missing observations.
- **Mitigation**: Safeguarded via optional unwrapping.

---

## Stress Test Results

- **Test Scenario 1 (Positive Windup)**: Apply positive wind-up (0.25 cycles) using `apply_windup_to_obs`.
  - **Expected behavior**: Corrected carrier phase reduces from 10.0 and 20.0 to 9.75 and 19.75.
  - **Actual behavior**: L1 phase is 9.75, L2 phase is 19.75.
  - **Status**: PASS

- **Test Scenario 2 (Workspace Tests)**: Run full EKF validation.
  - **Expected behavior**: All 256 tests pass.
  - **Actual behavior**: 256 passed, 0 failed.
  - **Status**: PASS

---

## Unchallenged Areas
- **Exact satellite attitude / solar panel orientation models**: The EKF relies on standard ECEF sun-position-based yaw modeling for satellites. Challenge of sub-daily attitude quirks under eclipse periods is out of scope.

---

# 5-Component Handoff Report

### 1. Observation
- Checked file `crates/gneiss-rtk/src/engine/measurement.rs`:
  - `apply_windup_to_obs` was changed to perform subtraction on `cp_l1` and `cp_l2`:
    ```rust
    *cp -= windup;
    ```
  - Added test `test_phase_windup_correction_sign_rtk` asserting that a 0.25 wind-up reduces L1 and L2 phases to 9.75 and 19.75 respectively.
- Checked file `crates/gneiss-rtk/src/engine/ppp_iekf.rs`:
  - In `push_cp_measurement`: `let l_meas = (cp1 - windup) * sat.lam1;`
  - In residual calculation:
    ```rust
    let res_l1 = (sat.cp1.unwrap() - windup) * sat.lam1 - (expected_base - i1 + n1);
    let res_l2 = (sat.cp2.unwrap() - windup) * sat.lam2 - (expected_base - gamma * i1 + n2);
    ```
- Checked file `crates/gneiss-rtk/src/engine/ppp_ins_iekf.rs`:
  - In `push_cp_measurement`: `let l_meas = (cp1 - windup) * sat.lam1;`
  - In residual calculation:
    ```rust
    let res_l1 = (sat.cp1.unwrap() - windup) * sat.lam1 - (expected_base - idx.i1 + idx.n1);
    let res_l2 = (sat.cp2.unwrap() - windup) * sat.lam2 - (expected_base - idx.gamma * idx.i1 + idx.n2);
    ```
- Executed `cargo test --workspace` and `cargo test -p gneiss-rtk`: both returned `ok. 256 passed; 0 failed`.
- Executed `cargo fmt --check`: returned exit code 1 due to formatting discrepancies in unrelated files (e.g., `crates/gneiss-rtk/src/engine/kinematics.rs`).

### 2. Logic Chain
1. The sign correction in `apply_windup_to_obs`, `ppp_iekf.rs`, and `ppp_ins_iekf.rs` changes windup correction from addition (`+ windup`) to subtraction (`- windup`).
2. Corrected carrier phase observations in RTK/PPP formulations must subtract windup cycles from the raw observed phase (or equivalently, add to the geometric model range prediction).
3. The subtraction implementation in all three locations is consistent.
4. The unit tests verify the expected mathematical reduction of phase values under positive windup.
5. All workspace tests passed, showing no regressions are introduced.

### 3. Caveats
- No caveats. The fix is localized to phase wind-up logic and is verified by targeted unit tests and full-suite builds.

### 4. Conclusion
The fix for Bug 18 is correct, complete, and robust. It completely addresses the opposite sign in phase wind-up correction. No further changes to source code are required.

### 5. Verification Method
To independently verify the fix:
- Run workspace tests:
  ```bash
  cargo test --workspace
  ```
- Run gneiss-rtk specific tests:
  ```bash
  cargo test -p gneiss-rtk
  ```
- Run formatting check:
  ```bash
  cargo fmt --check
  ```
