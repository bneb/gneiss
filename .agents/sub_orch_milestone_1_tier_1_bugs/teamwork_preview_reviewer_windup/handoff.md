# Handoff Report — Bug 18: Opposite Sign in Phase Wind-Up Correction Review

## 1. Observation
- We inspected the codebase for the changes related to Bug 18: Opposite Sign in Phase Wind-Up Correction.
- We checked `crates/gneiss-rtk/src/engine/measurement.rs` and found that the phase wind-up correction is still being added at lines 44 and 47:
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
- We checked `crates/gneiss-rtk/src/engine/ppp.rs` and found that the wind-up correction `wup` is still added at lines 497, 498, 503, 545, 546, 588, and 589:
  ```rust
  (cp1 + wup) * sat.lam1
  (sat.cp2.unwrap() + wup) * sat.lam2
  ```
- We checked `crates/gneiss-rtk/src/engine/ppp_iekf.rs` and found that the wind-up correction `windup` is still added at lines 898, 1026, and 1043:
  ```rust
  let res_l1 = (sat.cp1.unwrap() + windup) * sat.lam1 - (expected_base - i1 + n1);
  ```
- We checked `crates/gneiss-rtk/src/engine/ppp_ins_iekf.rs` and found that the wind-up correction `windup` is still added at lines 928, 1042, and 1061:
  ```rust
  let res_l2 =
      (sat.cp2.unwrap() + windup) * sat.lam2 - (expected_base - idx.gamma * idx.i1 + idx.n2);
  ```
- We ran a git log search (`git log --all -S "test_phase_windup_correction_sign_rtk"`) and a file search for the unit regression test `test_phase_windup_correction_sign_rtk` across the repository and found **zero occurrences** in source files. The only references exist in the explorer's `analysis.md` and `handoff.md`.
- We ran `cargo test -p gneiss-rtk` and it failed because of an unrelated test failure:
  ```
  failures:
      engine::ppp_iekf::mutant_killer_tests::test_resolve_widelane_ar_insufficient
  test result: FAILED. 254 passed; 1 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.11s
  ```

## 2. Logic Chain
1. The orchestrator dispatched us to verify the changes and tests for Bug 18: Opposite Sign in Phase Wind-Up Correction.
2. Based on our observations, no phase wind-up correction changes have been implemented or checked in on `main` or any other branch in the git history. All relevant files still contain the incorrect addition of wind-up to the phase observations.
3. The specified regression test `test_phase_windup_correction_sign_rtk` does not exist in any codebase source files.
4. The test suite has an existing failure in `engine::ppp_iekf::mutant_killer_tests::test_resolve_widelane_ar_insufficient`.
5. Therefore, the work product cannot be verified or approved at this time.

## 3. Caveats
- We assume that no other local uncommitted stashes or unpushed branches contain the changes, which is supported by our comprehensive git history searches.
- The failure of `test_resolve_widelane_ar_insufficient` is likely caused by other uncommitted changes currently present in the working tree, which are unrelated to Bug 18.

## 4. Conclusion
The changes and tests for Bug 18: Opposite Sign in Phase Wind-Up Correction are **absent** from the codebase. The files under review still apply the wind-up correction as addition (`+`) rather than subtraction (`-`), and the unit regression test `test_phase_windup_correction_sign_rtk` is missing.
Verdict: **REQUEST_CHANGES** is issued because the implementation and tests do not exist in the workspace.

## 5. Verification Method
To verify the fix once implemented:
1. Inspect `measurement.rs`, `ppp.rs`, `ppp_iekf.rs`, and `ppp_ins_iekf.rs` to ensure all wind-up additions have been replaced with subtractions.
2. Check that the regression test `test_phase_windup_correction_sign_rtk` is present in `measurement.rs` and passes.
3. Run `cargo test -p gneiss-rtk` and verify all tests pass.

---

# Quality Review Report

## Review Summary

**Verdict**: REQUEST_CHANGES

## Findings

### [Critical] Finding 1: Implementation of Wind-Up Sign Correction is Missing

- **What**: The phase wind-up corrections are still additions (`+` or `+=`) rather than subtractions (`-` or `-=`).
- **Where**: 
  - `crates/gneiss-rtk/src/engine/measurement.rs` (lines 44, 47)
  - `crates/gneiss-rtk/src/engine/ppp.rs` (lines 497, 498, 503, 545, 546, 588, 589)
  - `crates/gneiss-rtk/src/engine/ppp_iekf.rs` (lines 898, 1026, 1043)
  - `crates/gneiss-rtk/src/engine/ppp_ins_iekf.rs` (lines 928, 1042, 1061)
- **Why**: Adding the wind-up correction results in the opposite sign correction, violating the physical models (Wu et al. 1992).
- **Suggestion**: Change all additions of phase wind-up correction to subtractions.

### [Critical] Finding 2: Unit Regression Test is Missing

- **What**: The unit regression test `test_phase_windup_correction_sign_rtk` is not present in `measurement.rs` or elsewhere.
- **Where**: `crates/gneiss-rtk/src/engine/measurement.rs`
- **Why**: There is no test verifying the correct sign of the wind-up correction.
- **Suggestion**: Implement the `test_phase_windup_correction_sign_rtk` unit test as designed in the explorer's `analysis.md`.

## Verified Claims
- None.

## Coverage Gaps
- **Wind-Up Application in Factor Graphs** — risk level: Low — recommendation: Accept risk for now. The factor graph solver does not yet implement carrier phase wind-up, but once it does, it must follow the subtraction convention.

## Unverified Items
- All items are unverified because the changes are not implemented.

---

# Adversarial Review Report

## Challenge Summary

**Overall risk assessment**: HIGH (due to missing implementation)

## Challenges

### [High] Challenge 1: Phase Wind-up Sign Assumption Failure

- **Assumption challenged**: That the wind-up correction is additive.
- **Attack scenario**: In high-precision positioning (PPP/RTK) over long sessions, satellite-receiver geometry changes significantly. An incorrect wind-up sign (+ instead of -) causes a systematic error of up to 0.5 cycles (~10cm on L1) that varies over time, leading to EKF divergence or incorrect integer ambiguity resolution.
- **Blast radius**: Degraded convergence time and accuracy in both RTK and PPP engines.
- **Mitigation**: Change the implementation to subtract the wind-up correction.

## Stress Test Results
- None performed since changes are absent.

## Unchallenged Areas
- Physical calculation in `gneiss_core::windup::phase_windup` — out of scope.
