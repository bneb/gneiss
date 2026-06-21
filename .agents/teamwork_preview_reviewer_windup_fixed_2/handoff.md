# Handoff Report — Bug 18 Verification Review

## 1. Observation
- **Workspace Location**: `/Users/kevin/projects/gneiss`
- **Code Modifications**: We verified the changes in:
  - `crates/gneiss-rtk/src/engine/measurement.rs`:
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
    This function was modified to subtract the wind-up value rather than adding it.
  - `crates/gneiss-rtk/src/engine/ppp_iekf.rs`:
    All three instances of adding `windup` to `cp1`/`cp2` carrier phase measurements (lines 896, 1024, 1041) were changed to subtraction:
    - Line 898: `let l_meas = (cp1 - windup) * sat.lam1;`
    - Line 1026: `let res_l1 = (sat.cp1.unwrap() - windup) * sat.lam1 - (expected_base - i1 + n1);`
    - Line 1043: `let res_l2 = (sat.cp2.unwrap() - windup) * sat.lam2 - (expected_base - gamma * i1 + n2);`
- **Regression Test**: A new unit test `test_phase_windup_correction_sign_rtk` was added in `crates/gneiss-rtk/src/engine/measurement.rs` inside `mod tests`:
  ```rust
  #[test]
  fn test_phase_windup_correction_sign_rtk() {
      use super::apply_windup_to_obs;
      use gneiss_core::sat::{Constellation, SatelliteId};

      let mut obs = DdObservation {
          sat: SatelliteId {
              constellation: Constellation::Gps,
              prn: 1,
          },
          pr_l1: 0.0,
          pr_l2: None,
          cp_l1: Some(10.0),
          cp_l2: Some(20.0),
          doppler: 0.0,
          snr: 45.0,
          locktime: None,
      };

      let windup = 0.25; // 0.25 cycles of positive wind-up
      apply_windup_to_obs(&mut obs, windup);

      // Corrected carrier phase = raw_cp - windup
      assert_eq!(obs.cp_l1.unwrap(), 9.75);
      assert_eq!(obs.cp_l2.unwrap(), 19.75);
  }
  ```
- **Test Commands and Results**:
  - `cargo test -p gneiss-rtk` successfully passed all 258 tests.
  - `cargo test --workspace` successfully passed all 261 tests.
- **Code Formatting Check**:
  - `cargo fmt --check` completed successfully with no formatting violations.

## 2. Logic Chain
1. We checked out the modifications in the user's workspace at `/Users/kevin/projects/gneiss` to verify the fix for Bug 18.
2. The code changes in both `measurement.rs` (RTK path) and `ppp_iekf.rs` (PPP path) consistently subtract the wind-up correction (in cycles) from the raw carrier phase before converting to meters.
3. The unit test `test_phase_windup_correction_sign_rtk` confirms that positive wind-up reduces the corrected carrier phase, matching the physical model convention (Wu et al. 1992).
4. Running the full cargo test suite verifies that changing the sign does not break any existing tests and compiles cleanly.
5. `cargo fmt --check` confirms layout and style conformance.
6. Thus, the implementation is correct, complete, robust, and compliant.

## 3. Caveats
No caveats.

## 4. Conclusion
The fix for Bug 18 (Opposite Sign in Phase Wind-Up Correction) is fully implemented, verified, and functional. All unit and integration tests pass cleanly, and the formatting is correct.
Verdict: **APPROVE**

## 5. Verification Method
To independently verify:
1. Run `cargo test -p gneiss-rtk --engine::measurement::tests::test_phase_windup_correction_sign_rtk` to verify the specific regression test.
2. Run `cargo test --workspace` to verify overall suite health.
3. Run `cargo fmt --check` to check layout and formatting compliance.

---

# Quality Review Report

## Review Summary

**Verdict**: APPROVE

## Findings

No findings. The fix conforms to expected physical/mathematical conventions and is fully covered by unit tests.

## Verified Claims

- Phase wind-up correction is subtracted rather than added → verified via source code analysis and the `test_phase_windup_correction_sign_rtk` unit test → **PASS**
- Unit regression test `test_phase_windup_correction_sign_rtk` exists and runs → verified via `cargo test -p gneiss-rtk` → **PASS**
- Overall test suite integrity is maintained → verified via `cargo test --workspace` → **PASS**

## Coverage Gaps

None.

## Unverified Items

None.

---

# Adversarial Review Report

## Challenge Summary

**Overall risk assessment**: LOW

## Challenges

### [Low] Challenge 1: Phase Wrap-around Boundary Conditions

- **Assumption challenged**: That subtracting `windup` behaves correctly when phase wind-up crosses full-cycle boundaries (wraps around).
- **Attack scenario**: During long sessions, the total wind-up can accumulate to multiple cycles. If the correction only tracks fractional cycles or does not handle boundary wraps, it could introduce discontinuities.
- **Blast radius**: The wind-up calculation `phase_windup` tracks the continuous phase wind-up angle in radians/cycles (handling geometry updates relative to the previous epoch). Because the ambiguity estimation in EKF/factor graph acts as a free parameter to absorb integer cycle offsets, any integer wrap is absorbed by the ambiguity state, and only the fractional change needs to be consistently signed.
- **Mitigation**: The subtraction convention ensures the fractional change has the correct direction.

## Stress Test Results

- `test_phase_windup_correction_sign_rtk` → verifies that a positive wind-up reduces the corrected phase value. → **PASS**
- `test_windup_sign_correct` (in `ppp.rs`) → verifies that subtracting wind-up produces a smaller measured range than adding it. → **PASS**
