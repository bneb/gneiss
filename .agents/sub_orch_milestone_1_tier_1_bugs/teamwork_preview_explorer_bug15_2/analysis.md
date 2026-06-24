# Bug 15 Analysis: Incorrect Broadcast Clock TGD Correction

## Overview
In GNSS positioning engines, the timing group delay (TGD) or broadcast group delay (BGD) corrections must be correctly applied to the satellite broadcast clock calculations to align single-frequency and dual-frequency measurements.
For GPS, Galileo, and QZSS, the legacy broadcast clock coefficients (af0, af1, af2) are referenced to the dual-frequency ionosphere-free (IF) combination. Thus:
- Single-frequency users must subtract TGD or BGD.
- Dual-frequency ionosphere-free users do not apply any TGD correction (pass `0.0`).

However, for Beidou, the legacy D1/D2 navigation message broadcast clock parameters (af0, af1, af2) are referenced to the B3I frequency. Therefore:
- Single-frequency B1I users must subtract $T_{GD1}$.
- Single-frequency B2I/B2b users must subtract $T_{GD2}$.
- Dual-frequency B1I/B2I ionosphere-free users must subtract the combined ionosphere-free group delay:
  $$T_{GD\_IF} = \frac{f_{B1I}^2 T_{GD1} - f_{B2I}^2 T_{GD2}}{f_{B1I}^2 - f_{B2I}^2}$$

Currently, Gneiss-core incorrectly assumes Beidou works like GPS, passing `0.0` as `tgd` in `BeidouEphemeris::position_iono_free`. Furthermore, `BeidouEphemeris` is completely missing the `tgd2` field, and the RINEX parser parses `TGD2` into `aodc` while discarding the real `aodc` value.

---

## Code Observations

### 1. Missing `tgd2` and Incorrect Parser Mapping in `crates/gneiss-parsers/src/rinex.rs`
In `crates/gneiss-parsers/src/rinex.rs` (lines 541-577):
```rust
fn build_beidou_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<gneiss_core::ephemeris::Ephemeris> {
    Some(gneiss_core::ephemeris::Ephemeris::Beidou(
        gneiss_core::ephemeris::BeidouEphemeris {
            sat,
            toc,
            toe: GpsTime::new(toc.week, vals[8]),
            af0,
            af1,
            af2,
            aode: vals[0] as u32,
            crs: vals[1],
            delta_n: vals[2],
            m0: vals[3],
            cuc: vals[4],
            e: vals[5],
            cus: vals[6],
            sqrt_a: vals[7],
            cic: vals[9],
            omega0: vals[10],
            cis: vals[11],
            i0: vals[12],
            crc: vals[13],
            omega: vals[14],
            omega_dot: vals[15],
            idot: vals[16],
            tgd1: vals[22],
            aodc: vals[23] as u32,
        },
    ))
}
```
In RINEX 3 Navigation files, the Beidou broadcast message parameters are mapped to `vals` where:
- `vals[22]` is `TGD1` (Timing Group Delay 1, B1I relative to B3I)
- `vals[23]` is `TGD2` (Timing Group Delay 2, B2I relative to B3I)
- `vals[25]` is `AODC` (Age of Data Clock)

Currently, Gneiss parses `TGD2` (at `vals[23]`) into `aodc` of `BeidouEphemeris`, and completely discards `TGD2` and the actual `aodc` (at `vals[25]`).

### 2. Incorrect `position_iono_free` for Beidou in `crates/gneiss-core/src/ephemeris.rs`
In `crates/gneiss-core/src/ephemeris.rs` (lines 680-712):
```rust
impl BeidouEphemeris {
    ...
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
}
```
By passing `0.0` as `tgd`, it does not apply any group delay correction. But the B1I/B2I ionosphere-free clock correction should subtract the combined $T_{GD\_IF}$.

### 3. Invalid Assertion in `test_broadcast_clock_tgd_correct`
In `crates/gneiss-core/src/ephemeris.rs` (lines 1521-1524):
```rust
        let bds = Ephemeris::Beidou(bds_eph.clone());
        let (_, _, clk_pos, _) = bds.position(t);
        let (_, _, clk_if, _) = bds.position_iono_free(t);
        assert!((clk_if - clk_pos - bds_eph.tgd1).abs() < 1e-15);
```
This test asserts that `clk_if - clk_pos = tgd1`, validating the incorrect behavior where `clk_if` has no TGD correction.

---

## Proposed Fix Strategy

1. **Update `BeidouEphemeris` Structure**
   In `crates/gneiss-core/src/ephemeris.rs`, add `pub tgd2: f64` to `BeidouEphemeris`.
   Update all struct instantiations of `BeidouEphemeris` in the project (including tests in `ephemeris.rs`, `spp.rs`, etc.) to initialize `tgd2: 0.0` or appropriate dummy values.

2. **Fix `bgd_e5b` Method in `Ephemeris`**
   In `crates/gneiss-core/src/ephemeris.rs` (lines 89-95), modify `bgd_e5b` to return `tgd2` for Beidou:
   ```rust
   pub fn bgd_e5b(&self) -> f64 {
       match self {
           Ephemeris::Galileo(e) => e.bgd_e1_e5b,
           Ephemeris::Beidou(e) => e.tgd2,
           other => other.tgd(),
       }
   }
   ```

3. **Correct RINEX Parser Mapping**
   In `crates/gneiss-parsers/src/rinex.rs` (lines 541-577), fix `build_beidou_ephemeris`:
   ```rust
   fn build_beidou_ephemeris(...) {
       Some(gneiss_core::ephemeris::Ephemeris::Beidou(
           gneiss_core::ephemeris::BeidouEphemeris {
               ...
               tgd1: vals[22],
               tgd2: vals[23],
               aodc: vals[25] as u32,
           },
       ))
   }
   ```

4. **Apply Correct Ionosphere-Free Clock Correction for Beidou**
   In `crates/gneiss-core/src/ephemeris.rs`, modify `BeidouEphemeris::position_iono_free` to calculate and apply $T_{GD\_IF}$:
   ```rust
   pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
       let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
       let is_bds_geo = self.sat.prn <= 5 || self.sat.prn >= 59;

       let f1 = gneiss_core::signal::FREQ_BDS_B1I;
       let f2 = gneiss_core::signal::FREQ_GAL_E5B; // shares BDS B2I/E5b frequency
       let f1_sq = f1 * f1;
       let f2_sq = f2 * f2;
       let tgd_if = (f1_sq * self.tgd1 - f2_sq * self.tgd2) / (f1_sq - f2_sq);

       calc_keplerian(
           t_bdt,
           self.toe,
           self.toc,
           self.af0,
           self.af1,
           self.af2,
           ...
           tgd_if,
           MU_BDS,
           OMEGA_E_BDS,
           is_bds_geo,
       )
   }
   ```

5. **Update Regression Tests**
   Update `test_broadcast_clock_tgd_correct` in `crates/gneiss-core/src/ephemeris.rs` to initialize `tgd2` and assert the correct relationship between `clk_if` and `clk_pos`.
   Specifically, if `tgd1 = 4.0e-9` and `tgd2 = 2.0e-9`, then:
   $$T_{GD\_IF} \approx 2.487168 \cdot 4.0\text{e-}9 - 1.487168 \cdot 2.0\text{e-}9 = 6.974336\text{e-}9\text{ s}$$
   `clk_if` should be $\Delta t_{SV\_brdc} - T_{GD\_IF}$.
   `clk_pos` should be $\Delta t_{SV\_brdc} - T_{GD1}$.
   So:
   $$\text{clk\_if} - \text{clk\_pos} = T_{GD1} - T_{GD\_IF}$$
   We should assert this correct difference:
   `assert!((clk_if - clk_pos - (bds_eph.tgd1 - tgd_if)).abs() < 1e-15);`
