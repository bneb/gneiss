# Handoff Report — Bug 15: Incorrect Broadcast Clock TGD Correction

## 1. Observation
In Gneiss, Beidou broadcast clock correction logic has the following code definitions:
* **`crates/gneiss-parsers/src/rinex.rs:573-574`**:
  ```rust
  tgd1: vals[22],
  aodc: vals[23] as u32,
  ```
  Here `vals[23]` (which is `TGD2` according to the RINEX 3 format specs for Beidou) is parsed into `aodc`, and the actual `aodc` value (at `vals[25]`) is discarded.
* **`crates/gneiss-core/src/ephemeris.rs:680-712`**:
  ```rust
  pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
      let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
      let is_bds_geo = self.sat.prn <= 5 || self.sat.prn >= 59;
      calc_keplerian(
          t_bdt,
          self.toe,
          self.toc,
          self.af0,
          self.af1,
          self.af2,
          ...
          0.0,
          MU_BDS,
          OMEGA_E_BDS,
          is_bds_geo,
      )
  }
  ```
  Here `0.0` is passed as the `tgd` argument to `calc_keplerian`, applying no timing group delay correction for the ionosphere-free combination.
* **`crates/gneiss-core/src/ephemeris.rs:1517-1524`**:
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
  This test asserts that `clk_if - clk_pos = tgd1`, validating the incorrect behavior where `clk_if` has no TGD correction.

## 2. Logic Chain
* BDS legacy D1/D2 navigation message broadcast clock parameters (af0, af1, af2) are referenced to the B3I frequency.
* For single-frequency B1I, the satellite clock correction is:
  $$\Delta t_{SV}(B1I) = \Delta t_{SV\_brdc} - T_{GD1}$$
* For single-frequency B2I (which shares frequency with Galileo E5b), the satellite clock correction is:
  $$\Delta t_{SV}(B2I) = \Delta t_{SV\_brdc} - T_{GD2}$$
* For the B1I/B2I dual-frequency ionosphere-free combination, the satellite clock correction must be the ionosphere-free combination of the two single-frequency clock corrections:
  $$\Delta t_{SV\_IF} = \frac{f_1^2 \Delta t_{SV}(B1I) - f_2^2 \Delta t_{SV}(B2I)}{f_1^2 - f_2^2} = \Delta t_{SV\_brdc} - T_{GD\_IF}$$
  where:
  $$T_{GD\_IF} = \frac{f_{B1I}^2 T_{GD1} - f_{B2I}^2 T_{GD2}}{f_{B1I}^2 - f_{B2I}^2}$$
* Therefore, `BeidouEphemeris` must store both `tgd1` and `tgd2`, which must be correctly parsed from the RINEX file (offsets `vals[22]` and `vals[23]`, respectively, with `aodc` at `vals[25]`).
* `BeidouEphemeris::position_iono_free` must calculate and apply $T_{GD\_IF}$ as the `tgd` argument.
* The test assertion in `test_broadcast_clock_tgd_correct` must be corrected to expect `clk_if - clk_pos = T_GD1 - T_GD_IF`.

## 3. Caveats
* This analysis applies to the legacy D1/D2 broadcast navigation message. Modern Beidou CNAV navigation messages have different clock references, but the current Gneiss engine parses and works with legacy D1/D2 format broadcast ephemerides from RINEX navigation files.
* Assumes B1I and B2I are the frequencies used for the BDS dual-frequency combination, as implemented in `spp.rs` and `ppp.rs`.

## 4. Conclusion
To fix Bug 15, the following changes are required:
1. Add `tgd2: f64` to `BeidouEphemeris` in `ephemeris.rs`.
2. Fix `bgd_e5b` in `Ephemeris` to return `tgd2` for Beidou.
3. Update `build_beidou_ephemeris` in `rinex.rs` to map `tgd2` from `vals[23]` and `aodc` from `vals[25]`.
4. In `BeidouEphemeris::position_iono_free`, calculate $T_{GD\_IF} = \frac{f_{B1I}^2 T_{GD1} - f_{B2I}^2 T_{GD2}}{f_{B1I}^2 - f_{B2I}^2}$ and pass it to `calc_keplerian`.
5. Update test fixtures to initialize `tgd2` and adjust the assertion in `test_broadcast_clock_tgd_correct`.

## 5. Verification Method
* Build the project: `cargo build`
* Run Gneiss-core unit tests: `cargo test -p gneiss-core`
* Verify that the updated `test_broadcast_clock_tgd_correct` passes.
