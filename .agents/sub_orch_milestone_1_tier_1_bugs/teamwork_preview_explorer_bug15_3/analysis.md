# Bug 15 Analysis: Incorrect Broadcast Clock TGD Correction

## Overview
In GNSS legacy navigation messages (Lnav/Dnav), satellite broadcast clock parameters ($a_{f0}, a_{f1}, a_{f2}$) are referenced to specific signal frequencies or combinations:
- **GPS/QZSS**: Clock parameters are referenced to the L1/L2 ionosphere-free (IF) combination.
- **Galileo**: Clock parameters are referenced to the E1/E5a ionosphere-free combination.
- **BeiDou (BDS)**: Clock parameters are referenced to the B3I frequency (band 7 / 1268.52 MHz).

Because BeiDou clock parameters are referenced to B3I, a group delay correction must be applied for both single-frequency and dual-frequency ionosphere-free (IF) users using other frequencies.
Specifically:
- **Single-frequency B1I user**: Satellite clock correction is $\Delta t_{SV} - T_{GD1}$.
- **Single-frequency B2I user**: Satellite clock correction is $\Delta t_{SV} - T_{GD2}$.
- **Dual-frequency B1I/B2I IF user**: Satellite clock correction is $\Delta t_{SV} - T_{GD\_IF}$, where:
  $$T_{GD\_IF} = \frac{f_{B1I}^2 T_{GD1} - f_{B2I}^2 T_{GD2}}{f_{B1I}^2 - f_{B2I}^2}$$

Currently, the `gneiss` engine has three critical issues:
1. `BeidouEphemeris` completely lacks a `tgd2` field.
2. The RINEX navigation parser (`rinex.rs`) incorrectly parses `TGD2` into the `aodc` field, and ignores the actual `AODC` value.
3. `BeidouEphemeris::position_iono_free()` passes `0.0` as the Timing Group Delay (TGD) correction, which ignores BDS group delay corrections for B1I/B2I users and causes errors up to ~18 meters in satellite clock corrections.

---

## Detailed Findings

### 1. Missing `tgd2` Field in `BeidouEphemeris`
In `crates/gneiss-core/src/ephemeris.rs` (lines 155–180), the `BeidouEphemeris` struct is defined as follows:
```rust
pub struct BeidouEphemeris {
    pub sat: SatelliteId,
    pub toe: GpsTime,
    pub toc: GpsTime,
    pub af0: f64,
    pub af1: f64,
    pub af2: f64,
    // ... orbit parameters ...
    pub tgd1: f64,
    pub aode: u32,
    pub aodc: u32,
}
```
There is no `tgd2` field to store the B2I/B3I group delay.

### 2. Incorrect RINEX Parser Mapping for BeiDou
In `crates/gneiss-parsers/src/rinex.rs` (lines 541–577), `build_beidou_ephemeris` maps the parsed fields:
```rust
            tgd1: vals[22],
            aodc: vals[23] as u32,
```
According to the RINEX 3.x specification, for BeiDou:
- `vals[22]` corresponds to `TGD1` (B1I/B3I group delay).
- `vals[23]` corresponds to `TGD2` (B2I/B3I group delay).
- `vals[25]` corresponds to the Age of Data Clock (`AODC`).
Thus, `vals[23]` is parsed as `aodc` instead of `tgd2`, and the actual `aodc` value is ignored.

### 3. Incorrect `position_iono_free()` logic for BeiDou
In `crates/gneiss-core/src/ephemeris.rs` (lines 680–713), the ionosphere-free position method for BeiDou is implemented as:
```rust
    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
        let is_bds_geo = self.sat.prn <= 5 || self.sat.prn >= 59;
        calc_keplerian(
            t_bdt,
            // ...
            0.0, // <-- INCORRECT: passes 0.0 as tgd
            MU_BDS,
            OMEGA_E_BDS,
            is_bds_geo,
        )
    }
```
Passing `0.0` is correct for GPS/Galileo (whose clock parameters are already referenced to the IF combination), but incorrect for BeiDou, which is referenced to B3I. For BeiDou, it must pass the B1I/B2I ionosphere-free timing group delay correction $T_{GD\_IF}$.

---

## Proposed Fix Strategy

### 1. Update `BeidouEphemeris` Struct
Add the `tgd2` field to the struct in `crates/gneiss-core/src/ephemeris.rs`:
```rust
pub struct BeidouEphemeris {
    // ...
    pub tgd1: f64,
    pub tgd2: f64,
    pub aode: u32,
    pub aodc: u32,
}
```

### 2. Correct `position_iono_free` for BeiDou
Compute and apply the B1I/B2I iono-free Timing Group Delay correction in `BeidouEphemeris::position_iono_free`:
```rust
    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
        let is_bds_geo = self.sat.prn <= 5 || self.sat.prn >= 59;

        let f1 = crate::signal::FREQ_BDS_B1I;
        let f2 = crate::signal::FREQ_GAL_E5B; // B2I frequency shares Galileo E5b
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

### 3. Correct the RINEX Parser
Update `build_beidou_ephemeris` in `crates/gneiss-parsers/src/rinex.rs`:
```rust
        gneiss_core::ephemeris::BeidouEphemeris {
            // ...
            tgd1: vals[22],
            tgd2: vals[23],
            aodc: vals[25] as u32,
        },
```

### 4. Update Unit Tests and Instantiations
Because `tgd2` is added to `BeidouEphemeris`, all instantiations of the struct in tests must be updated to include `tgd2: 0.0` (or appropriate mock values):
1. **In `crates/gneiss-core/src/ephemeris.rs`**:
   - `test_bds_boundaries`
   - `test_broadcast_clock_tgd_correct`
   - Initializations for `bds_geo_eph` and `bds_igso_eph`.
2. **In `crates/gneiss-rtk/src/estimators/spp.rs`**:
   - `eph_bds` initialization in tests.

Also update `test_broadcast_clock_tgd_correct` for BeiDou:
```rust
        // Test Beidou
        let bds_eph = BeidouEphemeris {
            // ...
            tgd1: 4.0e-9,
            tgd2: 3.0e-9,
            aode: 1,
            aodc: 1,
        };
        let bds = Ephemeris::Beidou(bds_eph.clone());
        let (_, _, clk_pos, _) = bds.position(t);
        let (_, _, clk_if, _) = bds.position_iono_free(t);

        let f1 = FREQ_BDS_B1I;
        let f2 = FREQ_GAL_E5B;
        let f1_sq = f1 * f1;
        let f2_sq = f2 * f2;
        let tgd_if = (f1_sq * bds_eph.tgd1 - f2_sq * bds_eph.tgd2) / (f1_sq - f2_sq);

        assert!((clk_if - clk_pos - (bds_eph.tgd1 - tgd_if)).abs() < 1e-15);
```
