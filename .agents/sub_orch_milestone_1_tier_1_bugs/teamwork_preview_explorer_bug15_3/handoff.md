# Handoff Report: Bug 15 — Incorrect Broadcast Clock TGD Correction

## 1. Observation
- In `crates/gneiss-core/src/ephemeris.rs`:
  - `BeidouEphemeris` struct definition (lines 155-180) does not contain `tgd2`:
    ```rust
    pub struct BeidouEphemeris {
        pub sat: SatelliteId,
        pub toe: GpsTime,
        pub toc: GpsTime,
        pub af0: f64,
        pub af1: f64,
        pub af2: f64,
        // ...
        pub tgd1: f64,
        pub aode: u32,
        pub aodc: u32,
    }
    ```
  - `BeidouEphemeris::position_iono_free` passes `0.0` as `tgd` (line 707):
    ```rust
    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        // ...
        calc_keplerian(
            t_bdt,
            self.toe,
            self.toc,
            self.af0,
            self.af1,
            self.af2,
            // ...
            0.0,
            MU_BDS,
            OMEGA_E_BDS,
            is_bds_geo,
        )
    }
    ```
  - The unit test `test_broadcast_clock_tgd_correct` (line 1524) asserts:
    ```rust
    assert!((clk_if - clk_pos - bds_eph.tgd1).abs() < 1e-15);
    ```

- In `crates/gneiss-parsers/src/rinex.rs`:
  - `build_beidou_ephemeris` maps fields from `vals` as:
    ```rust
    tgd1: vals[22],
    aodc: vals[23] as u32,
    ```

- In `crates/gneiss-rtk/src/estimators/spp.rs`:
  - `eph_bds` is instantiated in a unit test (line 1054) without `tgd2`.

---

## 2. Logic Chain
- **Step 1**: The BeiDou broadcast clock parameters ($a_{f0}, a_{f1}, a_{f2}$) are referenced to the B3I frequency. Users of other legacy frequencies (e.g. B1I, B2I) or combinations thereof must apply group delay corrections (TGD1 for B1I, TGD2 for B2I).
- **Step 2**: The dual-frequency ionosphere-free combination of B1I and B2I requires a combined Timing Group Delay correction:
  $$T_{GD\_IF} = \frac{f_1^2 \cdot TGD1 - f_2^2 \cdot TGD2}{f_1^2 - f_2^2}$$
- **Step 3**: Currently, `position_iono_free` passes `0.0` as the Timing Group Delay correction for BeiDou. This ignores the TGD correction entirely for BDS iono-free users, resulting in up to ~18 meters of satellite clock bias error.
- **Step 4**: The RINEX parser maps `vals[23]` to `aodc`. But according to the RINEX 3.x specification, `vals[23]` represents TGD2, and `vals[25]` represents AODC. Therefore, TGD2 is lost/ignored, and AODC is populated with the incorrect value.
- **Step 5**: To apply the correct TGD correction, `BeidouEphemeris` must have `tgd2`, the RINEX parser must correctly parse `tgd2` and `aodc`, and `position_iono_free()` must calculate and apply $T_{GD\_IF}$.

---

## 3. Caveats
- No caveats. The root cause is fully identified and conforms to both the BeiDou ICD and RINEX 3.x specifications.

---

## 4. Conclusion
- A read-only investigation has located three issues causing incorrect Broadcast Clock TGD correction for BeiDou:
  1. `tgd2` is missing from `BeidouEphemeris`.
  2. RINEX parser maps `vals[23]` (TGD2) to `aodc` instead of parsing `tgd2` and using `vals[25]` for `aodc`.
  3. `BeidouEphemeris::position_iono_free()` applies `0.0` Timing Group Delay, ignoring the B1I/B2I iono-free Timing Group Delay correction.

The proposed fix strategy is detailed in `analysis.md`.

---

## 5. Verification Method
1. Run `cargo test -p gneiss-core` and verify all tests pass.
2. In particular, check `test_broadcast_clock_tgd_correct` to ensure the corrected BeiDou assertion behaves according to the B1I/B2I iono-free group delay formula.
3. Check `cargo test -p gneiss-parsers` and `cargo test -p gneiss-rtk` to verify no compilation errors occur due to adding the new struct field.
