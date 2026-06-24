# Handoff Report - Bug 16 Fix: Mismatched Galileo BGD Correction

## 1. Observation
In `crates/gneiss-core/src/ephemeris.rs`, the Galileo ephemeris broadcast clock error calculation subtracted `bgd_e1_e5a` (BGD for E1/E5a) unconditionally via `calc_keplerian`. However, Galileo E5b tracking (band 7) requires the use of `bgd_e1_e5b` group delay correction.

The following files were identified as needing modifications:
- `crates/gneiss-core/src/ephemeris.rs`: `position` on `Ephemeris` and `GalileoEphemeris` did not support E5b specific clock calculations.
- `crates/gneiss-core/src/signal.rs`: `get_frequency` did not handle band 7 frequency lookup.
- `crates/gneiss-rtk/src/estimators/spp.rs`: `SppMeasurement` had no `freq_band` field, meaning the SPP state estimation could not distinguish when band 7 was being tracked.
- `crates/gneiss-rtk/src/engine/spp_tight.rs`: `process_measurement` evaluated satellite clock corrections via `position` unconditionally, ignoring band 7.

## 2. Logic Chain
1. Added `position_e5b` method to `GalileoEphemeris` which calls `calc_keplerian` passing `self.bgd_e1_e5b` instead of `self.bgd_e1_e5a`.
2. Implemented `position_e5b` on `Ephemeris` which delegates to `GalileoEphemeris::position_e5b` when the constellation is Galileo, and falls back to `position` for other constellations.
3. Updated `get_frequency` in `crates/gneiss-core/src/signal.rs` to return `FREQ_GAL_E5B` for band 7 (for Galileo and Beidou).
4. Added `pub freq_band: u8` field to `SppMeasurement`.
5. Updated `build_single_measurement` in `crates/gneiss-rtk/src/estimators/spp.rs` to compute `freq_band` based on observations:
   - If `p2_opt` is present: returns `7` for Galileo/Beidou if band 7 is matched, `5` for Galileo if band 5 is matched, `6` for Beidou if band 6 is matched, else `2`.
   - If only `p1_opt` is present: returns `1`.
6. Updated `compute_sat_state` in `crates/gneiss-rtk/src/estimators/spp.rs` and `process_measurement` in `crates/gneiss-rtk/src/engine/spp_tight.rs` to invoke `position_e5b` if `m.freq_band == 7`, ensuring that the correct `bgd_e1_e5b` correction is subtracted from the broadcast clock error when tracking band 7.
7. Added a regression test `test_galileo_position_e5b_regression` verifying that the clock error difference between `position` and `position_e5b` matches the difference between `bgd_e1_e5a` and `bgd_e1_e5b`.

## 3. Caveats
- Beidou uses `tgd1` as defined on `BeidouEphemeris` for all bands since there is no separate BGD field for other Beidou bands currently defined in the ephemeris struct. Therefore, delegating to `position` for non-Galileo constellations under `position_e5b` is correct and safe.

## 4. Conclusion
The implementation of the mismatched Galileo BGD correction is complete. The correct group delay (`bgd_e1_e5b`) is now used whenever band 7 is tracked in both SPP state calculation and tight EKF measurement processing.

## 5. Verification Method
Verify that all workspace tests compile and pass successfully by running:
```bash
cargo test --workspace
```
And verify that the specific regression test is run and passes:
```bash
cargo test --package gneiss-core --lib -- ephemeris::tests::test_galileo_position_e5b_regression
```
If the buggy code is restored (by removing the delegate in `Ephemeris::position_e5b` and making it call `position`), the regression test will fail as the calculated clock difference will be `0` instead of `bgd_e1_e5b - bgd_e1_e5a`.
