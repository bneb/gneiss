# Handoff Report: Melbourne-Wübbena Slip Detection Verification

## 1. Observation
- In `crates/gneiss-rtk/src/engine/ppp_math.rs`, line 184 contains:
  ```rust
  let mw = (cp1 - cp2) - (p1 / lam1 + p2 / lam2) * (lam2 - lam1) / (lam1 + lam2);
  ```
- The unit test `test_mw_slip_detection` is defined at lines 646–668 in `crates/gneiss-rtk/src/engine/ppp_math.rs`:
  ```rust
      #[test]
      fn test_mw_slip_detection() {
          let p1 = 20000000.0;
          let p2 = p1;
          let cp1 = p1 / 0.19;
          let cp2 = p2 / 0.24;
          
          let (_, mw) = detect_mw_slip(cp1, 0.19, cp2, 0.24, p1, p2, 0.0, false, 2.0);
          
          // Move by 1000m (normal geometry change) -> should NOT trigger a slip
          let dist_change = 1000.0;
          let p1_new = p1 + dist_change;
          let p2_new = p2 + dist_change;
          let cp1_new = p1_new / 0.19;
          let cp2_new = p2_new / 0.24;
          
          let (slip, _) = detect_mw_slip(cp1_new, 0.19, cp2_new, 0.24, p1_new, p2_new, mw, true, 2.0);
          assert!(!slip, "MW should cancel geometry changes");
          
          // Now introduce a 5-cycle slip on L1 phase
          let cp1_slip = cp1_new + 5.0;
          let (slip2, _) = detect_mw_slip(cp1_slip, 0.19, cp2_new, 0.24, p1_new, p2_new, mw, true, 2.0);
          assert!(slip2, "MW should detect phase cycle slips");
      }
  ```
- Running the unit test with the correct code:
  `cargo test --package gneiss-rtk --lib -- engine::ppp_math::tests::test_mw_slip_detection`
  Output:
  ```
  running 1 test
  test engine::ppp_math::tests::test_mw_slip_detection ... ok

  test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 253 filtered out; finished in 0.00s
  ```
- Changing line 184 in `crates/gneiss-rtk/src/engine/ppp_math.rs` to:
  ```rust
  let mw = (cp1 - cp2) - (p1 / lam1 + p2 / lam2) * (lam1 * lam2) / (lam1 + lam2);
  ```
  And running `cargo test --package gneiss-rtk --lib -- engine::ppp_math::tests::test_mw_slip_detection`
  Output:
  ```
  running 1 test
  test engine::ppp_math::tests::test_mw_slip_detection ... FAILED

  failures:

  ---- engine::ppp_math::tests::test_mw_slip_detection stdout ----

  thread 'engine::ppp_math::tests::test_mw_slip_detection' (17511347) panicked at crates/gneiss-rtk/src/engine/ppp_math.rs:662:9:
  MW should cancel geometry changes
  ```
- Running the full workspace test suite with the correct restored code:
  `cargo test --workspace`
  Output:
  ```
  test result: ok. 253 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.09s
  ...
  test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
  ```

## 2. Logic Chain
1. **Mathematical Correction Verification**: Under geometry change $d$, the pseudorange term `(p1 / lam1 + p2 / lam2)` changes by $d / \lambda_1 + d / \lambda_2 = d (\lambda_1 + \lambda_2) / (\lambda_1 \lambda_2)$.
   - Multiplying this by the corrected scaling factor `(lam2 - lam1) / (lam1 + lam2)` yields $d (\lambda_2 - \lambda_1) / (\lambda_1 \lambda_2)$, which equals $d / \lambda_{WL}$ (where wide-lane wavelength $\lambda_{WL} = \lambda_1 \lambda_2 / (\lambda_2 - \lambda_1)$).
   - This cancels the geometry change term in the phase difference `cp1 - cp2` (which is in cycles, changing by $d / \lambda_{WL}$).
   - If the buggy factor `(lam1 * lam2) / (lam1 + lam2)` is used, the pseudorange term evaluates to $d$, which does not match the unit of cycles of `cp1 - cp2` (changing by $d / \lambda_{WL}$). This causes a large remaining geometry dependent term under movement, falsely triggering a cycle slip.
2. **Regression Test Verification**:
   - With the correct code, the test passes because the MW combination cancels geometry change and does not report a false cycle slip.
   - When the buggy code is restored, the geometry-dependent term is not cancelled, causing a large jump in the Melbourne-Wübbena value on geometry changes. This causes `slip` to be `true` and the test panics at `MW should cancel geometry changes`, proving the test acts as an effective regression test.
3. **Workspace Integrity**: Re-running all tests in the workspace shows everything passes successfully, confirming that the fix does not break any other module or feature in the workspace.

## 3. Caveats
- No caveats.

## 4. Conclusion
- The Melbourne-Wübbena combination scaling factor typo is completely and correctly fixed in the codebase using `(lam2 - lam1) / (lam1 + lam2)`.
- The regression test `test_mw_slip_detection` is functional, robust, and correctly fails when the bug is reintroduced and passes when the fix is present.
- All workspace tests pass.

## 5. Verification Method
1. Open `crates/gneiss-rtk/src/engine/ppp_math.rs` and verify the MW combination calculation:
   ```rust
   let mw = (cp1 - cp2) - (p1 / lam1 + p2 / lam2) * (lam2 - lam1) / (lam1 + lam2);
   ```
2. Run the specific unit test:
   ```bash
   cargo test --package gneiss-rtk --lib -- engine::ppp_math::tests::test_mw_slip_detection
   ```
3. Run all tests in the workspace:
   ```bash
   cargo test --workspace
   ```
