use gneiss_core::obs::{ObsCode, SatObs};
use gneiss_core::sat::Constellation;
use gneiss_core::time::GpsTime;
use nalgebra::Vector3;

pub const LIGHT_SPEED: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

#[derive(Default)]
pub struct OsbCorrections {
    pub p1: Option<f64>,
    pub p2: Option<f64>,
    pub cp1: Option<f64>,
    pub cp2: Option<f64>,
    pub osb_p1: f64,
    pub osb_p2: f64,
    pub osb_cp1: f64,
    pub osb_cp2: f64,
}

pub fn apply_osb_corrections(
    sinex_opt: Option<&gneiss_parsers::sinex_bia::SinexBias>,
    sat_obs: &SatObs,
    time: GpsTime,
    f1: f64,
    f2: f64,
    f1_b: u8,
    f2_b: u8,
) -> OsbCorrections {
    let get = |t, fb| {
        sat_obs
            .observations
            .iter()
            .find(|o| o.code.obs_type == t && o.code.signal.freq_band == fb)
            .map(|o| (o.value, o.code))
    };
    let mut out = OsbCorrections::default();

    let p1_obs = get(gneiss_core::obs::ObsType::Pseudorange, f1_b);
    let p2_obs = get(gneiss_core::obs::ObsType::Pseudorange, f2_b);
    let cp1_obs = get(gneiss_core::obs::ObsType::CarrierPhase, f1_b);
    let cp2_obs = get(gneiss_core::obs::ObsType::CarrierPhase, f2_b);

    out.p1 = p1_obs.map(|o| o.0);
    out.p2 = p2_obs.map(|o| o.0);
    out.cp1 = cp1_obs.map(|o| o.0);
    out.cp2 = cp2_obs.map(|o| o.0);

    if let Some(sinex) = sinex_opt {
        let convert = |b: f64| b * 1e-9 * LIGHT_SPEED;

        let get_osb = |obs: Option<(f64, ObsCode)>| -> f64 {
            obs.and_then(|(_, c)| sinex.get_bias(sat_obs.sat, c, time))
                .map(convert)
                .unwrap_or(0.0)
        };

        out.osb_p1 = get_osb(p1_obs);
        out.osb_p2 = get_osb(p2_obs);
        out.osb_cp1 = get_osb(cp1_obs);
        out.osb_cp2 = get_osb(cp2_obs);

        if let Some(v) = out.p1.as_mut() {
            *v -= out.osb_p1;
        }
        if let Some(v) = out.p2.as_mut() {
            *v -= out.osb_p2;
        }
        if let Some(v) = out.cp1.as_mut() {
            *v -= out.osb_cp1 / (LIGHT_SPEED / f1);
        }
        if let Some(v) = out.cp2.as_mut() {
            *v -= out.osb_cp2 / (LIGHT_SPEED / f2);

            let is_cp2_fb = cp2_obs
                .and_then(|(_, c)| sinex.get_exact_bias(sat_obs.sat, c, time))
                .is_none()
                && cp2_obs
                    .and_then(|(_, c)| sinex.get_bias(sat_obs.sat, c, time))
                    .is_some();
            if is_cp2_fb && sat_obs.sat.constellation == Constellation::Gps {
                if let Some((_, c)) = cp2_obs {
                    let s = c.to_string();
                    if s == "L2L" || s == "L2S" || s == "L2X" {
                        *v -= 0.25;
                    }
                }
            }
        }
    }
    out
}

pub fn apply_earth_rotation(
    raw_pos: Vector3<f64>,
    raw_vel: Vector3<f64>,
    rcv_pos: Vector3<f64>,
) -> (Vector3<f64>, Vector3<f64>) {
    let rng1 = (raw_pos - rcv_pos).norm();
    let tau1 = rng1 / LIGHT_SPEED;
    let th1 = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * tau1;
    let c1 = libm::cos(th1);
    let s1 = libm::sin(th1);
    let p1 = Vector3::new(
        raw_pos.x * c1 + raw_pos.y * s1,
        -raw_pos.x * s1 + raw_pos.y * c1,
        raw_pos.z,
    );
    let _v1 = Vector3::new(
        raw_vel.x * c1 + raw_vel.y * s1,
        -raw_vel.x * s1 + raw_vel.y * c1,
        raw_vel.z,
    );

    let rng2 = (p1 - rcv_pos).norm();
    let tau2 = rng2 / LIGHT_SPEED;
    let th2 = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * tau2;
    let c2 = libm::cos(th2);
    let s2 = libm::sin(th2);
    let p2 = Vector3::new(
        raw_pos.x * c2 + raw_pos.y * s2,
        -raw_pos.x * s2 + raw_pos.y * c2,
        raw_pos.z,
    );
    let v2 = Vector3::new(
        raw_vel.x * c2 + raw_vel.y * s2,
        -raw_vel.x * s2 + raw_vel.y * c2,
        raw_vel.z,
    );

    (p2, v2)
}

pub fn detect_cycle_slip(sat_obs: &SatObs, prev: u32) -> (bool, u32) {
    let mut slip = false;
    let mut lk = prev.saturating_add(1);

    for obs in sat_obs
        .observations
        .iter()
        .filter(|o| o.code.obs_type == gneiss_core::obs::ObsType::CarrierPhase)
    {
        if let Some(lli) = obs.lli {
            let slip_cond = (lli & 1) != 0 || (lli & 2) != 0;
            if slip_cond {
                return (true, 0);
            }
        }

        if let Some(l) = obs.lock_time {
            let l32 = l as u32;
            if l32 == 0 || l32 < prev {
                return (true, l32);
            } else {
                lk = lk.min(l32);
            }
        }
    }

    if lk == 0 && prev > 0 {
        slip = true;
    }
    (slip, lk)
}

/// Detect cycle slip using geometry-free (L1-L2) combination.
/// Removes geometry, clock, and troposphere; a jump indicates a slip or ionospheric spike.
/// Returns true if the GF combination jumps by more than `threshold_m` meters.
pub fn detect_gf_slip(
    cp1: f64, lam1: f64,
    cp2: f64, lam2: f64,
    prev_gf: f64,
    has_prev: bool,
    threshold_m: f64,
) -> (bool, f64) {
    let gf = cp1 * lam1 - cp2 * lam2;
    if !has_prev {
        return (false, gf);
    }
    let jump = (gf - prev_gf).abs();
    (jump > threshold_m, gf)
}

/// Detect cycle slip using Melbourne-Wübbena (MW) widelane combination.
/// The MW combination removes geometry, ionosphere, troposphere, and clock,
/// isolating widelane ambiguity. A jump indicates a cycle slip.
/// Returns true if MW jumps by more than `threshold_cycles` cycles.
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
    let mw = (cp1 - cp2) - (p1 / lam1 + p2 / lam2) * (lam1 * lam2) / (lam1 + lam2);
    if !has_prev {
        return (false, mw);
    }
    let jump = (mw - prev_mw).abs();
    (jump > threshold_cycles, mw)
}

/// Combined cycle slip detection using LLI, lock-time, geometry-free, and MW.
/// Returns (is_slip, new_lock_time).
/// `gf_prev` and `mw_prev` are updated in-place with the current values if no slip.
pub fn detect_slip_combined(
    sat_obs: &SatObs,
    prev_lock: u32,
    cp1: Option<f64>, lam1: f64,
    cp2: Option<f64>, lam2: f64,
    p1: Option<f64>,
    p2: Option<f64>,
    gf_prev: &mut Option<f64>,
    mw_prev: &mut Option<f64>,
) -> (bool, u32) {
    // LLI and lock-time checks (hardware-reported)
    let (hw_slip, lk) = detect_cycle_slip(sat_obs, prev_lock);
    if hw_slip {
        *gf_prev = None;
        *mw_prev = None;
        return (true, lk);
    }

    // Geometry-free check (requires dual-frequency phase)
    if let (Some(c1), Some(c2)) = (cp1, cp2) {
        if let Some(prev) = *gf_prev {
            let (gf_slip, new_gf) = detect_gf_slip(c1, lam1, c2, lam2, prev, true, 0.05);
            if gf_slip {
                *gf_prev = None;
                *mw_prev = None;
                return (true, 0);
            }
            *gf_prev = Some(new_gf);
        } else {
            // Initialize GF
            let (_, new_gf) = detect_gf_slip(c1, lam1, c2, lam2, 0.0, false, 0.0);
            *gf_prev = Some(new_gf);
        }
    }

    // Melbourne-Wübbena check (requires dual-frequency phase + pseudorange)
    if let (Some(c1), Some(c2), Some(pr1), Some(pr2)) = (cp1, cp2, p1, p2) {
        if let Some(prev) = *mw_prev {
            let (mw_slip, new_mw) = detect_mw_slip(c1, lam1, c2, lam2, pr1, pr2, prev, true, 2.0);
            if mw_slip {
                *gf_prev = None;
                *mw_prev = None;
                return (true, 0);
            }
            *mw_prev = Some(new_mw);
        } else {
            // Initialize MW
            let (_, new_mw) = detect_mw_slip(c1, lam1, c2, lam2, pr1, pr2, 0.0, false, 0.0);
            *mw_prev = Some(new_mw);
        }
    }

    (false, lk)
}

pub fn compute_tropo_dry(
    rcv_pos_llh: Vector3<f64>,
    el: f64,
    t: GpsTime,
    mapper: &dyn gneiss_core::atmosphere::TropoMapper,
) -> (f64, f64) {
    if el < 0.0 {
        return (0.0, 0.0);
    }
    let tp = gneiss_core::atmosphere::TropoParams::default();
    let alt_m = rcv_pos_llh.z;
    if alt_m > 20000.0 {
        return (0.0, 0.0);
    }

    let press_scale = libm::pow(1.0 - 0.0000226 * alt_m, 5.225);
    let p_z = tp.press_hpa * press_scale;

    let lat_rad = rcv_pos_llh.x;
    let z_dry_scale = 1.0 - 0.00266 * libm::cos(2.0 * lat_rad) - 0.00028 * alt_m / 1000.0;
    let z_dry = (0.0022768 * p_z) / z_dry_scale;

    let (mh, mw) = mapper.mapping_functions(rcv_pos_llh, el, t);
    (z_dry * mh, mw)
}

pub fn compute_iono_free(f1: f64, f2: f64, v1: f64, v2: f64) -> f64 {
    if f1 == f2 || f1 == 0.0 || f2 == 0.0 {
        return v1;
    }
    let g = (f1 * f1) / (f2 * f2);
    (g * v1 - v2) / (g - 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::sat::SatelliteId;
    use std::str::FromStr;

    #[test]
    fn test_apply_osb_corrections_no_sinex() {
        let mut obs = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![],
        };
        obs.sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        obs.observations.push(gneiss_core::obs::Observation {
            code: ObsCode::from_str("C1C").unwrap(),
            value: 100.0,
            lli: None,
            lock_time: None,
        });

        let res = apply_osb_corrections(None, &obs, GpsTime::new(0, 0.0), 1.0e9, 1.0e9, 1, 2);
        assert_eq!(res.p1, Some(100.0));
        assert_eq!(res.osb_p1, 0.0);
    }

    #[test]
    fn test_apply_earth_rotation_zero() {
        let (p, _v) = apply_earth_rotation(
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
        );
        assert!((p.x - 1.0).abs() < 1e-9);
        assert!((p.y - 0.0).abs() < 1e-9);
        assert!((p.z - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_detect_cycle_slip_no_slip() {
        let mut obs = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![],
        };
        obs.observations.push(gneiss_core::obs::Observation {
            code: ObsCode::from_str("L1C").unwrap(),
            value: 100.0,
            lli: Some(0),
            lock_time: Some(10),
        });
        let (slip, lk) = detect_cycle_slip(&obs, 5);
        assert_eq!(slip, false);
        assert_eq!(lk, 6);
    }
    #[test]
    fn test_detect_cycle_slip_locktime_equal() {
        let mut obs = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![],
        };
        obs.observations.push(gneiss_core::obs::Observation {
            code: ObsCode::from_str("L1C").unwrap(),
            value: 100.0,
            lli: Some(0),
            lock_time: Some(5),
        });
        let (slip, lk) = detect_cycle_slip(&obs, 5);
        assert_eq!(slip, false);
        assert_eq!(lk, 5);
    }

    #[test]
    fn test_detect_cycle_slip_lli_bit1() {
        let mut obs = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![],
        };
        obs.observations.push(gneiss_core::obs::Observation {
            code: ObsCode::from_str("L1C").unwrap(),
            value: 100.0,
            lli: Some(1),
            lock_time: Some(10),
        });
        let (slip, lk) = detect_cycle_slip(&obs, 5);
        assert_eq!(slip, true);
        assert_eq!(lk, 0);
    }

    #[test]
    fn test_detect_cycle_slip_lli_bit2() {
        let mut obs = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![],
        };
        obs.observations.push(gneiss_core::obs::Observation {
            code: ObsCode::from_str("L1C").unwrap(),
            value: 100.0,
            lli: Some(2),
            lock_time: Some(10),
        });
        let (slip, lk) = detect_cycle_slip(&obs, 5);
        assert_eq!(slip, true);
        assert_eq!(lk, 0);
    }

    #[test]
    fn test_detect_cycle_slip_locktime_zero() {
        let mut obs = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![],
        };
        obs.observations.push(gneiss_core::obs::Observation {
            code: ObsCode::from_str("L1C").unwrap(),
            value: 100.0,
            lli: Some(0),
            lock_time: Some(0),
        });
        let (slip, lk) = detect_cycle_slip(&obs, 5);
        assert_eq!(slip, true);
        assert_eq!(lk, 0);
    }

    #[test]
    fn test_detect_cycle_slip_locktime_decrease() {
        let mut obs = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![],
        };
        obs.observations.push(gneiss_core::obs::Observation {
            code: ObsCode::from_str("L1C").unwrap(),
            value: 100.0,
            lli: Some(0),
            lock_time: Some(3),
        });
        let (slip, lk) = detect_cycle_slip(&obs, 5);
        assert_eq!(slip, true);
        assert_eq!(lk, 3);
    }

    #[test]
    fn test_compute_tropo_dry_negative_el() {
        let mapper = gneiss_core::atmosphere::NmfMapper;
        let (d, w) = compute_tropo_dry(Vector3::zeros(), -0.1, GpsTime::new(0, 0.0), &mapper);
        assert_eq!(d, 0.0);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_compute_tropo_dry_high_alt() {
        let mapper = gneiss_core::atmosphere::NmfMapper;
        let (d, w) = compute_tropo_dry(Vector3::new(0.0, 0.0, 20001.0), 1.0, GpsTime::new(0, 0.0), &mapper);
        assert_eq!(d, 0.0);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_compute_tropo_dry_valid() {
        let mapper = gneiss_core::atmosphere::NmfMapper;
        let (d, w) = compute_tropo_dry(Vector3::new(0.5, 0.0, 100.0), 0.5, GpsTime::new(0, 0.0), &mapper);
        assert_eq!(d, 4.742703419689178);
        assert_eq!(w, 2.081891197333091);
    }

    #[test]
    fn test_compute_iono_free() {
        let val = compute_iono_free(1.5e9, 1.2e9, 10.0, 15.0);
        // g = (1.5 / 1.2)^2 = 1.5625
        // (1.5625 * 10.0 - 15.0) / 0.5625 = (15.625 - 15.0) / 0.5625 = 0.625 / 0.5625 = 1.1111111111111112
        assert!((val - 1.1111111111111112).abs() < 1e-9);
    }

    #[test]
    fn test_compute_iono_free_zero_f1() {
        let val = compute_iono_free(0.0, 1.2e9, 10.0, 20.0);
        assert_eq!(val, 10.0);
    }

    #[test]
    fn test_compute_iono_free_zero_f2() {
        let val = compute_iono_free(1.5e9, 0.0, 10.0, 20.0);
        assert_eq!(val, 10.0);
    }

    #[test]
    fn test_compute_iono_free_equal_f() {
        let val = compute_iono_free(1.5e9, 1.5e9, 10.0, 20.0);
        assert_eq!(val, 10.0);
    }
    #[test]
    fn test_apply_earth_rotation_nonzero() {
        let raw_pos = Vector3::new(20000000.0, 10000000.0, 5000000.0);
        let raw_vel = Vector3::new(1000.0, 2000.0, 3000.0);
        let rcv_pos = Vector3::new(10000000.0, -10000000.0, 2000000.0);
        let (p, v) = apply_earth_rotation(raw_pos, raw_vel, rcv_pos);

        let rng1 = (raw_pos - rcv_pos).norm();
        let tau1 = rng1 / LIGHT_SPEED;
        let th1 = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * tau1;
        let c1 = libm::cos(th1);
        let s1 = libm::sin(th1);
        let p1 = Vector3::new(
            raw_pos.x * c1 + raw_pos.y * s1,
            -raw_pos.x * s1 + raw_pos.y * c1,
            raw_pos.z,
        );
        let _v1 = Vector3::new(
            raw_vel.x * c1 + raw_vel.y * s1,
            -raw_vel.x * s1 + raw_vel.y * c1,
            raw_vel.z,
        );

        let rng2 = (p1 - rcv_pos).norm();
        let tau2 = rng2 / LIGHT_SPEED;
        let th2 = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * tau2;
        let c2 = libm::cos(th2);
        let s2 = libm::sin(th2);
        let p2 = Vector3::new(
            raw_pos.x * c2 + raw_pos.y * s2,
            -raw_pos.x * s2 + raw_pos.y * c2,
            raw_pos.z,
        );
        let v2 = Vector3::new(
            raw_vel.x * c2 + raw_vel.y * s2,
            -raw_vel.x * s2 + raw_vel.y * c2,
            raw_vel.z,
        );

        assert_eq!(p.x, p2.x);
        assert_eq!(p.y, p2.y);
        assert_eq!(p.z, p2.z);
        assert_eq!(v.x, v2.x);
        assert_eq!(v.y, v2.y);
        assert_eq!(v.z, v2.z);

        // Sanity checks
        assert!((p.x - raw_pos.x).abs() > 10.0);
        assert!((v.x - raw_vel.x).abs() > 0.005);
    }

    #[test]
    fn test_apply_osb_corrections_cp2_fb() {
        use gneiss_core::obs::{ObsCode, Observation, SatObs};
        use gneiss_core::sat::Constellation;
        use gneiss_core::time::GpsTime;
        use gneiss_parsers::sinex_bia::{BiasRecord, BiasType, SinexBias};

        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let record = BiasRecord {
            bias_type: BiasType::Osb,
            sat: sat_id,
            station: None,
            obs1: "L2W".parse::<ObsCode>().unwrap(),
            obs2: None,
            start_time: GpsTime::new(0, 0.0),
            end_time: GpsTime::new(0, 100.0),
            unit: "ns".to_string(),
            value: 2.0,
            std_dev: 0.0,
        };
        let sinex = SinexBias {
            records: vec![record],
        };
        // Check all three codes that trigger -0.25
        for code_str in ["L2L", "L2S", "L2X"] {
            let mut obs = SatObs {
                sat: sat_id,
                observations: vec![],
            };
            let code = code_str.parse::<ObsCode>().unwrap();
            obs.observations.push(Observation {
                code,
                value: 300.0,
                lli: None,
                lock_time: None,
            });

            let out = apply_osb_corrections(
                Some(&sinex),
                &obs,
                GpsTime::new(0, 50.0),
                1.5e9,
                1.2e9,
                1,
                2,
            );
            if let Some(cp2) = out.cp2 {
                assert_eq!(cp2, 297.35, "failed for code {}", code_str);
            } else {
                panic!("cp2 is None for {}", code_str);
            }
        }

        // Check a code that DOES NOT trigger -0.25 but still triggers fallback
        let mut obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let code = "L2C".parse::<ObsCode>().unwrap();
        obs.observations.push(Observation {
            code,
            value: 300.0,
            lli: None,
            lock_time: None,
        });
        let out = apply_osb_corrections(
            Some(&sinex),
            &obs,
            GpsTime::new(0, 50.0),
            1.5e9,
            1.2e9,
            1,
            2,
        );
        if let Some(cp2) = out.cp2 {
            assert_eq!(cp2, 300.0, "failed for code L2C");
        } else {
            panic!("cp2 is None for L2C");
        }
    }

    #[test]
    fn test_gf_slip_detection() {
        // L1=1000 cycles, L2=800 cycles, λ1=0.19, λ2=0.24
        // GF = 1000*0.19 - 800*0.24 = 190 - 192 = -2.0m
        let (_, gf) = detect_gf_slip(1000.0, 0.19, 800.0, 0.24, 0.0, false, 0.05);
        assert!((gf - (-2.0)).abs() < 0.01);
        // No slip when values are the same as previous
        let (slip, _) = detect_gf_slip(1000.0, 0.19, 800.0, 0.24, gf, true, 0.05);
        assert!(!slip);
        // Slip: changing L2 by 1 cycle (0.24m) exceeds 0.05m threshold
        let (slip2, _) = detect_gf_slip(1000.0, 0.19, 801.0, 0.24, gf, true, 0.05);
        assert!(slip2);
    }

    #[test]
    fn test_mw_slip_detection() {
        // MW with identical values → no slip
        let (_, mw) = detect_mw_slip(1000.0, 0.19, 800.0, 0.24, 20000000.0, 20000000.0, 0.0, false, 2.0);
        // First epoch: no previous → no slip
        let (slip, _) = detect_mw_slip(1000.0, 0.19, 800.0, 0.24, 20000000.0, 20000000.0, mw, true, 2.0);
        assert!(!slip);
        // Large MW jump of 5 cycles → slip
        let (slip2, _) = detect_mw_slip(1000.0, 0.19, 805.0, 0.24, 20000000.0, 20000000.0, mw, true, 2.0);
        assert!(slip2);
    }
}
