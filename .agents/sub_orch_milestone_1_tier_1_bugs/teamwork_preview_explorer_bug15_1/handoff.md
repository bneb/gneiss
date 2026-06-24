# Handoff Report — Bug 15: Incorrect Broadcast Clock TGD Correction

## 1. Observation
The following code structures were directly observed in the Gneiss codebase:
* **`crates/gneiss-parsers/src/rinex.rs` (lines 573-574)**:
  ```rust
  tgd1: vals[22],
  aodc: vals[23] as u32,
  ```
* **`crates/gneiss-core/src/ephemeris.rs` (lines 680-689)**:
  ```rust
  pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
      let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
      let is_bds_geo = self.sat.prn <= 5 || self.sat.prn >= 59;
      calc_keplerian(
          t_bdt,
          ...
          0.0,
          MU_BDS,
          OMEGA_E_BDS,
          is_bds_geo,
      )
  }
  ```
* **`crates/gneiss-core/src/ephemeris.rs` (lines 1517-1524)**:
  ```rust
          let bds_eph = BeidouEphemeris {
              ...
              tgd1: 4.0e-9,
              aode: 1,
              aodc: 1,
          };
          let bds = Ephemeris::Beidou(bds_eph.clone());
          let (_, _, clk_pos, _) = bds.position(t);
          let (_, _, clk_if, _) = bds.position_iono_free(t);
          assert!((clk_if - clk_pos - bds_eph.tgd1).abs() < 1e-15);
  ```

## 2. Logic Chain
1. **Clock reference**: The legacy D1/D2 Beidou broadcast clock coefficients (af0, af1, af2) are referenced to the B3I frequency.
2. **Single frequency corrections**: 
   - Single-frequency B1I: $\Delta t_{SV}(B1I) = \Delta t_{SV}(B3I) - T_{GD1}$
   - Single-frequency B2I: $\Delta t_{SV}(B2I) = \Delta t_{SV}(B3I) - T_{GD2}$
3. **Dual frequency combination**: For the B1I/B2I ionosphere-free dual-frequency combination, the clock correction is:
   $$\Delta t_{SV\_IF} = \frac{f_1^2 \Delta t_{SV}(B1I) - f_2^2 \Delta t_{SV}(B2I)}{f_1^2 - f_2^2} = \Delta t_{SV}(B3I) - T_{GD\_IF}$$
   where $T_{GD\_IF} = \frac{f_{B1I}^2 T_{GD1} - f_{B2I}^2 T_{GD2}}{f_{B1I}^2 - f_{B2I}^2}$.
4. **Current implementation bugs**: 
   - `BeidouEphemeris` is missing `tgd2` and incorrectly maps `TGD2` (`vals[23]`) to `aodc`, discarding the actual `aodc` (`vals[25]`).
   - `BeidouEphemeris::position_iono_free` applies a timing group delay correction of `0.0`, neglecting $T_{GD\_IF}$.
   - The test assertion in `test_broadcast_clock_tgd_correct` verifies this incorrect `0.0` delay correction instead of expecting `tgd1 - tgd_if`.

## 3. Caveats
* This applies to legacy D1/D2 Beidou broadcast messages. Modern Beidou CNAV navigation messages have different clock references, but the current Gneiss engine parses and works with legacy D1/D2 format broadcast ephemerides from RINEX navigation files.
* B1I and B2I/E5b are assumed to be the frequencies used for the BDS dual-frequency combination, as implemented in `spp.rs` and `ppp.rs`.

## 4. Conclusion
To correct Bug 15, the following read-only proposed modifications must be implemented:
1. Update the `BeidouEphemeris` struct to add `pub tgd2: f64` in `crates/gneiss-core/src/ephemeris.rs`.
2. Update the `bgd_e5b` method in `Ephemeris` to return `self.tgd2` for Beidou.
3. Update `build_beidou_ephemeris` in `crates/gneiss-parsers/src/rinex.rs` to map `tgd2: vals[23]` and `aodc: vals[25] as u32`.
4. Modify `BeidouEphemeris::position_iono_free` in `crates/gneiss-core/src/ephemeris.rs` to compute $T_{GD\_IF}$ and pass it to `calc_keplerian`.
5. Update all instantiations of `BeidouEphemeris` across the tests to initialize `tgd2`, and correct the assertion in `test_broadcast_clock_tgd_correct`.

## 5. Verification Method
* Build the project: `cargo build`
* Run Gneiss-core tests: `cargo test -p gneiss-core`
* Verify that `test_broadcast_clock_tgd_correct` passes successfully.
