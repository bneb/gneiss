# Handoff Report — Bug 15 Fix Integrity Audit

## 1. Observation

- **Implementation files modified**:
  - `crates/gneiss-core/src/ephemeris.rs` (lines 47-55, 517-546, 579-608, 647-680, 713-742)
  - `crates/gneiss-rtk/src/engine/ppp.rs` (lines 427, 430-433, 450)
  - `crates/gneiss-rtk/src/estimators/spp.rs` (lines 185-189, 194-198)
- **Tests file modified**:
  - `crates/gneiss-core/src/ephemeris.rs` (lines 1383-1515, unit test `test_broadcast_clock_tgd_correct`)
- **Commands executed and output**:
  - Running `cargo test ephemeris::tests::test_broadcast_clock_tgd_correct` succeeded:
    ```
    running 1 test
    test ephemeris::tests::test_broadcast_clock_tgd_correct ... ok
    test result: ok. 1 passed; 0 failed; 0 ignored; 47 filtered out; finished in 0.00s
    ```
  - Running `cargo test` succeeded:
    ```
    test result: ok. 257 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.30s
    ```
- **Integrity Level**: `demo` (read from `.agents/ORIGINAL_REQUEST.md`)

---

## 2. Logic Chain

1. **Bug Identification**: Bug 15 pertains to incorrect broadcast clock Timing Group Delay (TGD) corrections. Broadcast clocks are referenced to the dual-frequency ionosphere-free (L1/L2) combination. For dual-frequency or ionosphere-free observations, TGD cancels out in the combination, and thus the clock error should not be corrected by TGD.
2. **Implementation Check**:
   - The fix implements `position_iono_free` for the `Ephemeris` enum and individual constellation ephemeris types (`GpsEphemeris`, `GalileoEphemeris`, `BeidouEphemeris`, `QzssEphemeris`).
   - These methods call `calc_keplerian` passing `0.0` as the group delay argument (instead of `tgd`, `bgd_e1_e5a`, etc.), which correctly prevents subtraction of the group delay.
   - For `GlonassEphemeris`, which has no group delay, it delegates directly to `position`.
   - The PPP engine (`ppp.rs`) computes satellite states using `position_iono_free` because PPP is iono-free, and removes the manual `tgd()` addition hacks.
   - The SPP estimator (`spp.rs`) conditionally calls `position_iono_free` when `m.is_iono_free` is true, and `position` when it is false.
   - Since these changes correctly bypass TGD subtraction when appropriate, the mathematical implementation is genuine and correct (Observation: git diff).
3. **Integrity Violations Check**:
   - **Hardcoded test results**: The unit test `test_broadcast_clock_tgd_correct` dynamically checks that `position_iono_free(t)` minus `position(t)` equals the ephemeris's group delay. There are no static hardcoded assertions of expected positions or clock errors (Observation: git diff).
   - **Facade implementations**: Every implemented method actually computes Keplerian orbit parameters and clock corrections using full formulas. None return constant dummy values (Observation: git diff).
   - **Pre-populated logs/results**: No pre-populated logs or cheat artifacts exist in the workspace (Observation: `git status --ignored`).
   - **Layout Compliance**: All modifications are made to files in designated directories under `crates/gneiss-core/src/` and `crates/gneiss-rtk/src/`. No source or test files are placed inside `.agents/`.
4. **Conclusion Support**: The logic chain shows that the fix is correct, genuine, builds, passes all tests, and respects the layout and integrity constraints.

---

## 3. Caveats

- No caveats.

---

## 4. Conclusion

The Bug 15 fix is authentic, correct, and fully compliant with all integrity requirements.

---

## Forensic Audit Report

**Work Product**: Bug 15 fix implementation in gneiss project
**Profile**: General Project (Integrity mode: demo)
**Verdict**: CLEAN

### Phase Results
- **Hardcoded output detection**: PASS — No hardcoded test results or expected values found in the test suite or source code.
- **Facade detection**: PASS — Full and genuine implementations of `position_iono_free` delegate correctly to `calc_keplerian`.
- **Pre-populated artifact detection**: PASS — No pre-populated log files, result files, or verification artifacts exist.
- **Build and run**: PASS — The workspace successfully compiles, and the test suite executes with 257 passing tests.
- **Output verification**: PASS — Mathematical verification confirms that passing `0.0` as the group delay in `calc_keplerian` correctly prevents subtracting TGD/BGD.
- **Dependency audit**: PASS — No third-party packages or wrappers were introduced.

### Evidence
#### Git Diff (Summary of core changes)
```rust
// crates/gneiss-core/src/ephemeris.rs
pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
    match self {
        Ephemeris::Gps(e) => e.position_iono_free(t),
        Ephemeris::Galileo(e) => e.position_iono_free(t),
        Ephemeris::Beidou(e) => e.position_iono_free(t),
        Ephemeris::Qzss(e) => e.position_iono_free(t),
        Ephemeris::Glonass(e) => e.position(t),
    }
}
```

---

## 5. Verification Method

To verify the audit findings:
1. Run the specific unit test added for Bug 15:
   ```bash
   cargo test ephemeris::tests::test_broadcast_clock_tgd_correct
   ```
2. Run the entire workspace test suite:
   ```bash
   cargo test
   ```
3. Inspect `crates/gneiss-core/src/ephemeris.rs` to verify that `position_iono_free` calls `calc_keplerian` with `0.0` as the `tgd` argument.
