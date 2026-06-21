# Handoff Report — Bug 1: Melbourne-Wübbena Dimensional Typo

## 1. Observation

In `crates/gneiss-rtk/src/engine/ppp_math.rs`, the Melbourne-Wübbena (MW) combination function `detect_mw_slip` (lines 174–190) calculates `mw` in cycles.

In the current workspace version (at git commit `da013e27a4be9319e98e5389afa792140f4b49f4`), the formula is:
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

Prior to commit `da013e27a4be9319e98e5389afa792140f4b49f4` (introduced in commit `ebba760df51eb5b503d8a396c0de77a83ba79586`), the formula was:
```rust
    let mw = (cp1 - cp2) - (p1 / lam1 + p2 / lam2) * (lam1 * lam2) / (lam1 + lam2);
```

## 2. Logic Chain

To evaluate the dimensions (units) of the terms in the original formula:
1. `cp1 - cp2`:
   - `cp1` and `cp2` are carrier phase observations in **cycles**.
   - Thus, `cp1 - cp2` has units of **cycles**.
2. `p1 / lam1 + p2 / lam2`:
   - `p1` and `p2` are pseudoranges in **meters**.
   - `lam1` and `lam2` are wavelengths in **meters**.
   - Therefore, `p1 / lam1` and `p2 / lam2` are pseudoranges converted to **cycles**.
   - Their sum `(p1 / lam1 + p2 / lam2)` has units of **cycles**.
3. `(lam1 * lam2) / (lam1 + lam2)`:
   - This term represents the narrowlane wavelength $\lambda_{NL} = \frac{\lambda_1 \lambda_2}{\lambda_1 + \lambda_2}$.
   - Since $\lambda_1, \lambda_2$ are in meters, $\lambda_{NL}$ has units of **meters**.
4. Dimensional Mismatch:
   - Multiplying the cycles term `(p1 / lam1 + p2 / lam2)` by $\lambda_{NL}$ (meters) results in a term with units of **meters * cycles** (or **meters** if cycles are treated as dimensionless).
   - Subtracting this from `cp1 - cp2` (cycles) represents a dimensional mismatch (subtracting meters from cycles).

To calculate the Melbourne-Wübbena combination in cycles, we must convert the narrowlane pseudorange (which is in narrowlane cycles) to widelane cycles. We do this by scaling it by the ratio of narrowlane to widelane wavelengths:
$$\frac{\lambda_{NL}}{\lambda_{WL}} = \frac{\frac{\lambda_1 \lambda_2}{\lambda_1 + \lambda_2}}{\frac{\lambda_1 \lambda_2}{\lambda_2 - \lambda_1}} = \frac{\lambda_2 - \lambda_1}{\lambda_1 + \lambda_2}$$

This ratio is dimensionless.
Thus, replacing the narrowlane wavelength $\lambda_{NL}$ with the dimensionless ratio $\lambda_{NL} / \lambda_{WL}$ gives the dimensionally correct formula:
$$\phi_{MW} = (\phi_1 - \phi_2) - \frac{\lambda_2 - \lambda_1}{\lambda_1 + \lambda_2} \left( \frac{P_1}{\lambda_1} + \frac{P_2}{\lambda_2} \right)$$
In Rust code:
`let mw = (cp1 - cp2) - (p1 / lam1 + p2 / lam2) * (lam2 - lam1) / (lam1 + lam2);`

This matches the fixed formula currently in `ppp_math.rs`.

## 3. Caveats

- The bug has already been resolved in the workspace under commit `da013e27a4be9319e98e5389afa792140f4b49f4`.
- This report explains the mathematical and dimensional origin of the bug that was fixed in that commit.

## 4. Conclusion

The Melbourne-Wübbena cycle slip calculation previously had a dimensional typo where the narrowlane wavelength $\lambda_{NL}$ was used as the scaling factor instead of the dimensionless ratio $\lambda_{NL} / \lambda_{WL}$. This resulted in a cycles-vs-meters mismatch.

The correct fix strategy is to scale the narrowlane pseudorange cycle sum by the dimensionless ratio of wavelengths:
$$\frac{\lambda_2 - \lambda_1}{\lambda_1 + \lambda_2}$$
which translates in code to:
`let mw = (cp1 - cp2) - (p1 / lam1 + p2 / lam2) * (lam2 - lam1) / (lam1 + lam2);`

This is already implemented and verified in the current codebase state.

## 5. Verification Method

To verify the correct execution and logic of the slip detection math:
1. Run the targeted test suite:
   ```bash
   cargo test --package gneiss-rtk --lib -- engine::ppp_math::tests::test_mw_slip_detection
   ```
2. Verify that the test passes successfully:
   ```
   running 1 test
   test engine::ppp_math::tests::test_mw_slip_detection ... ok
   ```
