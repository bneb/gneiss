# Bug 15 Analysis: Incorrect Broadcast Clock TGD Correction

## 1. Location of TGD Correction Subtraction
The broadcast satellite clock correction calculation is implemented in `calc_keplerian` in `crates/gneiss-core/src/ephemeris.rs`. Specifically, the subtraction of the group delay (`tgd`) occurs on line 465:
```rust
let clk_err = af0 + af1 * tc + af2 * tc * tc + F * e * sqrt_a * libm::sin(ek) - tgd;
```
This function is called by the `position` method of each constellation's ephemeris struct (`GpsEphemeris`, `GalileoEphemeris`, `BeidouEphemeris`, `QzssEphemeris`), which passes the corresponding group delay value:
- `GpsEphemeris`: passes `self.tgd` (line 500)
- `GalileoEphemeris`: passes `self.bgd_e1_e5a` (line 532)
- `BeidouEphemeris`: passes `self.tgd1` (line 570)
- `QzssEphemeris`: passes `self.tgd` (line 602)

## 2. Handling of Single-Frequency vs. Dual-Frequency Clock Corrections
- **GPS/Galileo/BeiDou/QZSS Broadcast Reference**: The broadcast clock parameters ($a_{f0}, a_{f1}, a_{f2}$) are referenced to the dual-frequency (ionosphere-free) combination of signals (e.g., L1/L2 for GPS). Therefore, the raw broadcast clock correction represents the dual-frequency clock error directly.
- **Single-Frequency Mode**: Requires subtracting the group delay ($T_{GD}$ or $BGD$) from the raw clock correction to match the single-frequency observation group delay.
- **Dual-Frequency/Ionosphere-Free Mode**: The group delay cancels out in the ionosphere-free combination. Therefore, $T_{GD}$ must not be subtracted from the satellite clock correction.
- **Codebase Discrepancies**:
  - **In `crates/gneiss-rtk/src/engine/ppp.rs` (PPP, dual-frequency)**: Since `eph.position()` returns a clock error with $T_{GD}$ already subtracted, `ppp.rs` performs a manual addition to undo the correction (lines 431-434):
    ```rust
    dt_s = brdc_clk + eph.tgd();
    ```
  - **In `crates/gneiss-rtk/src/estimators/spp.rs` (SPP, single or dual-frequency)**: `SppMeasurement` defines `is_iono_free` to indicate whether a measurement is a dual-frequency combination. However, `compute_sat_state` directly applies the clock error from `eph.position()`, meaning $T_{GD}$ is incorrectly subtracted even for dual-frequency combined measurements.

## 3. Proposed Fix Strategy
A clean, correct fix strategy involves adding a new method `position_iono_free` to the ephemeris implementations to return the raw broadcast clock bias (i.e., with `tgd = 0.0` passed to `calc_keplerian`). This preserves the default `position` behavior for single-frequency L1 callers while providing a clean API for dual-frequency/ionosphere-free estimators.

### Proposed Changes

#### A. In `crates/gneiss-core/src/ephemeris.rs`:
Add `position_iono_free` to `Ephemeris` enum and individual ephemeris structs:
```rust
impl Ephemeris {
    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        match self {
            Ephemeris::Gps(e) => e.position_iono_free(t),
            Ephemeris::Galileo(e) => e.position_iono_free(t),
            Ephemeris::Beidou(e) => e.position_iono_free(t),
            Ephemeris::Qzss(e) => e.position_iono_free(t),
            Ephemeris::Glonass(e) => e.position(t), // GLONASS has no TGD
        }
    }
}

impl GpsEphemeris {
    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(
            t, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n, self.omega0, self.omega_dot,
            self.i0, self.idot, self.omega,
            0.0, // Pass 0.0 to avoid subtracting TGD
            MU_GPS, OMEGA_E_GPS, false,
        )
    }
}

impl GalileoEphemeris {
    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(
            t, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n, self.omega0, self.omega_dot,
            self.i0, self.idot, self.omega,
            0.0, // Pass 0.0 to avoid subtracting BGD
            MU_GAL, OMEGA_E_GAL, false,
        )
    }
}

impl BeidouEphemeris {
    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
        let is_bds_geo = self.sat.prn <= 5 || self.sat.prn >= 59;
        calc_keplerian(
            t_bdt, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n, self.omega0, self.omega_dot,
            self.i0, self.idot, self.omega,
            0.0, // Pass 0.0 to avoid subtracting TGD
            MU_BDS, OMEGA_E_BDS, is_bds_geo,
        )
    }
}

impl QzssEphemeris {
    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(
            t, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n, self.omega0, self.omega_dot,
            self.i0, self.idot, self.omega,
            0.0, // Pass 0.0 to avoid subtracting TGD
            MU_GPS, OMEGA_E_GPS, false,
        )
    }
}
```

#### B. In `crates/gneiss-rtk/src/engine/ppp.rs`:
Update `compute_sat_state` to use `position_iono_free` and remove the TGD addition logic:
```rust
    let brdc_clk = eph.position_iono_free(t_nom).2;

    let mut clk_found = dt_s != 0.0;
    if !precise {
        dt_s = brdc_clk;
        clk_found = true;
    }
```

#### C. In `crates/gneiss-rtk/src/estimators/spp.rs`:
Update `compute_sat_state` to query `position_iono_free` when `m.is_iono_free` is true:
```rust
    let sat_clk_err_rough = if m.is_iono_free {
        m.eph.position_iono_free(t_tx_sat_gps).2
    } else {
        m.eph.position(t_tx_sat_gps).2
    };

    let t_tx_true = t_tx_sat - sat_clk_err_rough;
    let t_tx_true_gps = GpsTime::new(m.time.week, t_tx_true);

    let (sat_pos, _, sat_clk_err, _) = if m.is_iono_free {
        m.eph.position_iono_free(t_tx_true_gps)
    } else {
        m.eph.position(t_tx_true_gps)
    };
```

## 4. Recommended Regression Test Design
A new test `test_broadcast_clock_tgd_correct` should be added in `crates/gneiss-core/src/ephemeris.rs` to verify that `position` and `position_iono_free` compute clock biases with and without TGD subtracted, respectively.

```rust
    #[test]
    fn test_broadcast_clock_tgd_correct() {
        let t = GpsTime::new(2000, 100000.0);
        let tgd_val = 1.5e-8_f64;

        // 1. GPS
        let gps_eph = GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 5 },
            toe: t, toc: t,
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.0, idot: 0.0, omega: 0.0,
            tgd: tgd_val, iode: 1, iodc: 1,
        };
        let eph_gps = Ephemeris::Gps(gps_eph);
        let (_, _, clk_gps, _) = eph_gps.position(t);
        let (_, _, clk_gps_if, _) = eph_gps.position_iono_free(t);
        assert!((clk_gps - (-tgd_val)).abs() < 1e-15);
        assert!(clk_gps_if.abs() < 1e-15);
        assert!((clk_gps_if - clk_gps - tgd_val).abs() < 1e-15);

        // 2. Galileo
        let gal_eph = GalileoEphemeris {
            sat: SatelliteId { constellation: Constellation::Galileo, prn: 3 },
            toe: t, toc: t,
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.0, sqrt_a: 5440.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.0, idot: 0.0, omega: 0.0,
            bgd_e1_e5a: tgd_val, bgd_e1_e5b: 2.0 * tgd_val, iod_nav: 1,
        };
        let eph_gal = Ephemeris::Galileo(gal_eph);
        let (_, _, clk_gal, _) = eph_gal.position(t);
        let (_, _, clk_gal_if, _) = eph_gal.position_iono_free(t);
        assert!((clk_gal - (-tgd_val)).abs() < 1e-15);
        assert!(clk_gal_if.abs() < 1e-15);
        assert!((clk_gal_if - clk_gal - tgd_val).abs() < 1e-15);

        // 3. BeiDou
        let bds_eph = BeidouEphemeris {
            sat: SatelliteId { constellation: Constellation::Beidou, prn: 6 },
            toe: t, toc: t,
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.0, idot: 0.0, omega: 0.0,
            tgd1: tgd_val, tgd2: 2.0 * tgd_val, aode: 1, aodc: 1,
        };
        let eph_bds = Ephemeris::Beidou(bds_eph);
        let (_, _, clk_bds, _) = eph_bds.position(t);
        let (_, _, clk_bds_if, _) = eph_bds.position_iono_free(t);
        assert!((clk_bds - (-tgd_val)).abs() < 1e-15);
        assert!(clk_bds_if.abs() < 1e-15);
        assert!((clk_bds_if - clk_bds - tgd_val).abs() < 1e-15);

        // 4. QZSS
        let qzss_eph = QzssEphemeris {
            sat: SatelliteId { constellation: Constellation::Qzss, prn: 193 },
            toe: t, toc: t,
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.0, idot: 0.0, omega: 0.0,
            tgd: tgd_val, iode: 1, iodc: 1,
        };
        let eph_qzss = Ephemeris::Qzss(qzss_eph);
        let (_, _, clk_qzss, _) = eph_qzss.position(t);
        let (_, _, clk_qzss_if, _) = eph_qzss.position_iono_free(t);
        assert!((clk_qzss - (-tgd_val)).abs() < 1e-15);
        assert!(clk_qzss_if.abs() < 1e-15);
        assert!((clk_qzss_if - clk_qzss - tgd_val).abs() < 1e-15);
    }
```
