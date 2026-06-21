# Melbourne-Wübbena Dimensional Typo Verification Report

## 1. Observation

- **Exact File Path & Code**:
  In `crates/gneiss-rtk/src/engine/ppp_math.rs`, lines 174-190 contain the Melbourne-Wübbena slip detection logic:
  ```rust
  pub fn detect_mw_slip(
      cp1: f64, lam1: f64,
      cp2: f64, lam2: f64,
      p1: f64,
      p2: f64,
      prev_mw: f64,
      has_prev: bool,
      threshold_cycles: f64,
  ) -> (bool, f64) {
      let _wl = lam1 * lam2 / (lam2 - lam1); // widelane wavelength
      let mw = (cp1 - cp2) - (p1 / lam1 + p2 / lam2) * (lam2 - lam1) / (lam1 + lam2);
      if !has_prev {
          return (false, mw);
      }
      let jump = (mw - prev_mw).abs();
      (jump > threshold_cycles, mw)
  }
  ```

- **Git Commit Diff**:
  Commit `da013e27a4be9319e98e5389afa792140f4b49f4` resolved the typo:
  ```diff
  -    let mw = (cp1 - cp2) - (p1 / lam1 + p2 / lam2) * (lam1 * lam2) / (lam1 + lam2);
  +    let mw = (cp1 - cp2) - (p1 / lam1 + p2 / lam2) * (lam2 - lam1) / (lam1 + lam2);
  ```

- **Verbatim Unit Test**:
  Lines 646-669 in `crates/gneiss-rtk/src/engine/ppp_math.rs`:
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

- **Command Results**:
  - `cargo test --package gneiss-rtk --lib -- engine::ppp_math::tests::test_mw_slip_detection`
    ```
    running 1 test
    test engine::ppp_math::tests::test_mw_slip_detection ... ok
    test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 253 filtered out; finished in 0.00s
    ```
  - `cargo test --workspace`
    ```
    test result: ok. 253 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.26s
    ...
    test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.01s
    ```

## 2. Logic Chain

1. **Dimensional Consistency**:
   - The carrier phase difference term `cp1 - cp2` is expressed in **cycles**.
   - The pseudorange term `p1 / lam1 + p2 / lam2` is expressed in **cycles**.
   - In the buggy implementation, the scaling factor was `(lam1 * lam2) / (lam1 + lam2)`, which has units of **meters** ($m^2 / m = m$). Multiplying cycles by meters results in meter-cycles, which cannot be subtracted from cycles.
   - In the corrected implementation, the scaling factor is `(lam2 - lam1) / (lam1 + lam2)`, which is a ratio of meters to meters and thus **dimensionless**. Multiplying cycles by a dimensionless ratio preserves the unit of **cycles**, making the formula dimensionally consistent.

2. **Geometry Cancellation**:
   - The Melbourne-Wübbena combination must completely remove the geometry term ($\rho$), clock biases ($dt, dT$), and tropospheric delay ($T$).
   - When the satellite-receiver distance changes by $1000\text{ m}$ (a geometry change), the carrier phases and pseudoranges change accordingly.
   - Under the corrected formula, this geometric change is scaled by the dimensionless factor $\frac{\lambda_2 - \lambda_1}{\lambda_1 + \lambda_2} \approx 0.116$ (for GPS frequencies), resulting in a perfect cancellation (residual jump of $0.0\text{ cycles}$).
   - Under the buggy formula, the scaling factor was $\frac{\lambda_1 \lambda_2}{\lambda_1 + \lambda_2} \approx 0.106\text{ m}$. A geometry change of $1000\text{ m}$ would yield a residual jump of $\approx 96.5\text{ cycles}$, which vastly exceeds the $2.0\text{ cycles}$ threshold and triggers false positive slips.

3. **Regression Test Authenticity**:
   - The previous test used inconsistent pseudorange ($20,000,000\text{ m}$) and carrier phase ($1000\text{ cycles}$) values, meaning it failed to model physical relationships.
   - The updated test computes physically consistent phase values (`cp = p / lam`) and verifies both geometry cancellation (no false positives for $1000\text{ m}$ movement) and slip detection (correct detection of a $5\text{ cycle}$ phase slip).

4. **Integrity Validation**:
   - No dummy implementations, hardcoded test result short-circuits, or verification facades exist in `ppp_math.rs` or any of the tests.

## 3. Caveats

- **Satellite-Specific Wavelengths**: The correctness of the MW combination relies on passing the exact wavelengths (`lam1`, `lam2`) corresponding to the frequency band of the satellite. For GLONASS satellites, which use Frequency Division Multiple Access (FDMA), these wavelengths differ per satellite channel. The caller must ensure that satellite-specific wavelengths are passed correctly.
- **Other Fixes**: Commit `da013e2` also modified other mathematical routines (velocity correction in earth rotation, OSB corrections, IMU updater asserts). This review focused specifically on the Melbourne-Wübbena slip detection math (Bug 1).

## 4. Conclusion

**Verdict**: **APPROVE**

The correction to the Melbourne-Wübbena combination in `crates/gneiss-rtk/src/engine/ppp_math.rs` is mathematically and dimensionally correct. The updated unit test `test_mw_slip_detection` is authentic, physically consistent, and serves as a highly robust regression test. All workspace builds and tests pass successfully.

## 5. Verification Method

- **Command to Execute**:
  ```bash
  cargo test --package gneiss-rtk --lib -- engine::ppp_math::tests::test_mw_slip_detection
  ```
- **Files to Inspect**:
  - `crates/gneiss-rtk/src/engine/ppp_math.rs` (lines 174-190, 646-669)
- **Invalidation Conditions**:
  - Reversion of the scaling factor in `detect_mw_slip` to the buggy version.
  - Failure of the unit test `test_mw_slip_detection`.

---

## 6. Quality Review Report

**Verdict**: **APPROVE**

### Findings
- **Minor Formatting Issue**: Running `cargo fmt --check` flags minor styling differences in `crates/gneiss-rtk/src/engine/ppp_math.rs` (e.g., parameter wrapping formatting). However, this does not affect correctness or build stability.

### Verified Claims
- **Melbourne-Wübbena Geometry Cancellation** $\rightarrow$ verified via unit test `test_mw_slip_detection` $\rightarrow$ **PASS**
- **Workspace Build and Test Suite Success** $\rightarrow$ verified via `cargo test --workspace` $\rightarrow$ **PASS**

### Coverage Gaps
- None. The unit test covers both the cancellation of geometry changes and the detection of actual phase slips.

### Unverified Items
- None.

---

## 7. Challenge Report (Adversarial Review)

**Overall Risk Assessment**: **LOW** (post-fix)

### Challenges

#### [Medium] Challenge 1: Pseudorange Noise Susceptibility (False Positives)
- **Assumption Challenged**: The fixed threshold of $2.0\text{ cycles}$ is robust against pseudorange noise.
- **Attack Scenario**: Under low elevation angles or severe multipath environments, pseudorange noise $\sigma_P$ can easily exceed $2.5\text{ meters}$. Because pseudorange noise propagates into the MW combination scaled by $\approx 0.83 \sigma_P$, a $2.5\text{ m}$ noise spike translates to $\approx 2.08\text{ cycles}$ of noise in the MW calculation. This will exceed the fixed $2.0\text{ cycles}$ threshold, triggering a false cycle slip detection.
- **Blast Radius**: Temporary loss of carrier phase ambiguity lock, causing the EKF to reset ambiguity states and degrade positioning accuracy.
- **Mitigation**: Implement a time-smoothed MW combination (moving average) or dynamically scale the threshold based on satellite elevation and measured pseudorange variance.

#### [High] Challenge 2: Cycle Slip Detection Blind Spots (False Negatives)
- **Assumption Challenged**: The combination of Geometry-Free (GF) and Melbourne-Wübbena (MW) slip detectors will catch all cycle slips.
- **Attack Scenario**: Specific combinations of cycle slips on L1 and L2 will pass through both detectors unnoticed.
  - **Case A: $(9, 7)$ Cycle Slip**: A simultaneous slip of $+9\text{ cycles}$ on L1 and $+7\text{ cycles}$ on L2 yields:
    - $\Delta GF = 9 \times 0.19029 - 7 \times 0.24421 = 0.003\text{ m}$ (below the $0.05\text{ m}$ threshold).
    - $\Delta MW = 9 - 7 = 2\text{ cycles}$ (does not exceed the `> 2.0` threshold).
    - Result: The $(9,7)$ cycle slip is completely missed by both detectors.
  - **Case B: $(5, 4)$ Cycle Slip**: A simultaneous slip of $+5\text{ cycles}$ on L1 and $+4\text{ cycles}$ on L2 yields:
    - $\Delta GF = 5 \times 0.19029 - 4 \times 0.24421 = -0.025\text{ m}$ (below the $0.05\text{ m}$ threshold).
    - $\Delta MW = 5 - 4 = 1\text{ cycle}$ (below the $2.0\text{ cycles}$ threshold).
    - Result: The $(5,4)$ cycle slip is completely missed by both detectors.
- **Blast Radius**: Unchecked cycle slips leak into the EKF, corrupting state estimates and causing positioning drifts.
- **Mitigation**: Supplement the combined slip detector with EKF innovation-based checks (residual testing) in the measurement update step.

### Stress Test Results

- **Scenario 1**: Geometry change of $1000\text{ m}$ $\rightarrow$ expected behavior: $0.0\text{ cycles}$ jump $\rightarrow$ actual behavior: $0.0\text{ cycles}$ jump $\rightarrow$ **PASS**
- **Scenario 2**: L1 cycle slip of $+5\text{ cycles}$ $\rightarrow$ expected behavior: detected via MW ($5.0\text{ cycles}$ jump) $\rightarrow$ actual behavior: detected ($5.0\text{ cycles}$ jump) $\rightarrow$ **PASS**
- **Scenario 3**: L1/L2 slip of $(5, 4)\text{ cycles}$ $\rightarrow$ expected behavior: detected $\rightarrow$ predicted/actual behavior: **MISSED** (GF jump = $2.5\text{ cm}$, MW jump = $1\text{ cycle}$) $\rightarrow$ **FAIL** (Blind spot confirmed)

### Unchallenged Areas
- Triple-frequency cycle slip combinations (out of scope, as the current engine is dual-frequency).
