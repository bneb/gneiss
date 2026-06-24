# Bug 15 Analysis: Incorrect Broadcast Clock TGD Correction

## Summary
In the Gneiss navigation engine, the timing group delay (TGD) or broadcast group delay (BGD) corrections align single-frequency and dual-frequency broadcast clocks. For GPS, Galileo, and QZSS, legacy broadcast clock coefficients are referenced to the dual-frequency ionosphere-free (IF) combination, so single-frequency users subtract TGD/BGD, while dual-frequency users apply no TGD correction. However, Beidou legacy D1/D2 broadcast clock parameters are referenced to the B3I frequency. Consequently, single-frequency B1I users must subtract $T_{GD1}$, single-frequency B2I users must subtract $T_{GD2}$, and dual-frequency B1I/B2I ionosphere-free users must subtract the combined ionosphere-free timing group delay $T_{GD\_IF}$. Currently, Gneiss-core incorrectly ignores $T_{GD\_IF}$ for Beidou, missing the $T_{GD2}$ field entirely, and the RINEX parser incorrectly maps $T_{GD2}$ to `aodc`.

---

## Detailed Observations

### 1. Code Location of Beidou Ephemeris Definition
In `crates/gneiss-core/src/ephemeris.rs`, the `BeidouEphemeris` struct is defined as follows:
```rust
pub struct BeidouEphemeris {
    pub sat: SatelliteId,
    pub toe: GpsTime,
    pub toc: GpsTime,
    pub af0: f64,
    pub af1: f64,
    pub af2: f64,
    ...
    pub tgd1: f64,
    pub aode: u32,
    pub aodc: u32,
}
```
It only has `tgd1` and lacks `tgd2`.

### 2. Incorrect Parser Mapping in `crates/gneiss-parsers/src/rinex.rs`
In `crates/gneiss-parsers/src/rinex.rs` (lines 573-574), Beidou ephemeris is constructed as follows:
```rust
            tgd1: vals[22],
            aodc: vals[23] as u32,
```
According to RINEX 3/4 broadcast ephemeris specifications for Beidou:
- `vals[22]` contains `TGD1` (B1I relative to B3I)
- `vals[23]` contains `TGD2` (B2I relative to B3I)
- `vals[25]` contains `AODC` (Age of Data Clock)

Currently, the parser incorrectly maps `TGD2` to `aodc` and discards the actual `aodc` value.

### 3. Incorrect `position_iono_free` Correction for Beidou
In `crates/gneiss-core/src/ephemeris.rs` (lines 680-712):
```rust
    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        ...
        calc_keplerian(
            t_bdt,
            ...
            0.0, // Incorrectly applies no TGD correction
            MU_BDS,
            OMEGA_E_BDS,
            is_bds_geo,
        )
    }
```
Passing `0.0` as `tgd` is only correct if the broadcast clock coefficients are referenced to the dual-frequency ionosphere-free combination (like GPS or Galileo). Since Beidou's reference is B3I, Beidou dual-frequency B1I/B2I clock correction must subtract $T_{GD\_IF}$.

### 4. Invalid Regression Test Assertion
In `crates/gneiss-core/src/ephemeris.rs` (lines 1521-1524):
```rust
        let bds = Ephemeris::Beidou(bds_eph.clone());
        let (_, _, clk_pos, _) = bds.position(t);
        let (_, _, clk_if, _) = bds.position_iono_free(t);
        assert!((clk_if - clk_pos - bds_eph.tgd1).abs() < 1e-15);
```
This test asserts that the difference between the dual-frequency clock and the single-frequency clock is exactly `tgd1`, confirming that the current implementation incorrectly applies no TGD correction for `position_iono_free`.

---

## Mathematical Logic

1. **Beidou Clock Corrections**:
   - $\Delta t_{SV}(B3I) = a_{f0} + a_{f1}(t - t_{oc}) + a_{f2}(t - t_{oc})^2 + \Delta t_r$
   - $\Delta t_{SV}(B1I) = \Delta t_{SV}(B3I) - T_{GD1}$
   - $\Delta t_{SV}(B2I) = \Delta t_{SV}(B3I) - T_{GD2}$

2. **Ionosphere-Free Combination**:
   $$\text{PR}_{IF} = \frac{f_1^2 \text{PR}_1 - f_2^2 \text{PR}_2}{f_1^2 - f_2^2}$$
   Substituting the satellite clocks gives the ionosphere-free satellite clock correction:
   $$\Delta t_{SV\_IF} = \frac{f_1^2 \Delta t_{SV}(B1I) - f_2^2 \Delta t_{SV}(B2I)}{f_1^2 - f_2^2} = \Delta t_{SV}(B3I) - T_{GD\_IF}$$
   where:
   $$T_{GD\_IF} = \frac{f_1^2 T_{GD1} - f_2^2 T_{GD2}}{f_1^2 - f_2^2}$$

3. **Beidou Frequency Settings in Gneiss**:
   In `crates/gneiss-core/src/signal.rs`:
   - $f_1 = f_{B1I} = 1561.098\text{ MHz}$
   - $f_2 = f_{B2I} = f_{E5b} = 1207.140\text{ MHz}$

---

## Proposed Fix Strategy

### 1. Update `BeidouEphemeris` Struct
Add `pub tgd2: f64` to `BeidouEphemeris` in `crates/gneiss-core/src/ephemeris.rs`. Update all instances of `BeidouEphemeris` construction in `crates/gneiss-core/src/ephemeris.rs` and `crates/gneiss-rtk/src/estimators/spp.rs` to include `tgd2` (initialized to `0.0` or appropriate test values).

### 2. Fix `bgd_e5b` in `Ephemeris` Enum
Update the `bgd_e5b` helper in `crates/gneiss-core/src/ephemeris.rs` to return `self.tgd2` for Beidou:
```rust
    pub fn bgd_e5b(&self) -> f64 {
        match self {
            Ephemeris::Galileo(e) => e.bgd_e1_e5b,
            Ephemeris::Beidou(e) => e.tgd2,
            other => other.tgd(),
        }
    }
```

### 3. Correct RINEX Parser Mapping
In `crates/gneiss-parsers/src/rinex.rs`, fix `build_beidou_ephemeris`:
```rust
            tgd1: vals[22],
            tgd2: vals[23],
            aodc: vals[25] as u32,
```

### 4. Implement Combined $T_{GD\_IF}$ in Beidou `position_iono_free`
Modify `position_iono_free` for `BeidouEphemeris` in `crates/gneiss-core/src/ephemeris.rs`:
```rust
    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
        let is_bds_geo = self.sat.prn <= 5 || self.sat.prn >= 59;

        let f1 = crate::signal::FREQ_BDS_B1I;
        let f2 = crate::signal::FREQ_GAL_E5B; // shares BDS B2I frequency
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
            self.crs,
            self.crc,
            self.cuc,
            self.cus,
            self.cic,
            self.cis,
            self.m0,
            self.e,
            self.sqrt_a,
            self.delta_n,
            self.omega0,
            self.omega_dot,
            self.i0,
            self.idot,
            self.omega,
            tgd_if,
            MU_BDS,
            OMEGA_E_BDS,
            is_bds_geo,
        )
    }
```

### 5. Update Test Cases and Assertions
Update the Beidou test fixture in `test_broadcast_clock_tgd_correct`:
```rust
        // Test Beidou
        let bds_eph = BeidouEphemeris {
            ...
            tgd1: 4.0e-9,
            tgd2: 2.0e-9,
            aode: 1,
            aodc: 1,
        };
        let bds = Ephemeris::Beidou(bds_eph.clone());
        let (_, _, clk_pos, _) = bds.position(t);
        let (_, _, clk_if, _) = bds.position_iono_free(t);
        
        let f1 = crate::signal::FREQ_BDS_B1I;
        let f2 = crate::signal::FREQ_GAL_E5B;
        let tgd_if = (f1 * f1 * bds_eph.tgd1 - f2 * f2 * bds_eph.tgd2) / (f1 * f1 - f2 * f2);
        
        assert!((clk_if - clk_pos - (bds_eph.tgd1 - tgd_if)).abs() < 1e-15);
```
