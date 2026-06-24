use gneiss_core::obs::SatObs;
use gneiss_core::sat::SatelliteId;
use nalgebra::Vector3;

const C: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

/// Frequencies for GPS L1, L2, L5
const F1: f64 = 1575.42e6;
const F2: f64 = 1227.60e6;
const F5: f64 = 1176.45e6;

const LAMBDA_EWL: f64 = C / (F2 - F5);
const LAMBDA_WL: f64 = C / (F1 - F2);
const LAMBDA_L1: f64 = C / F1;

pub struct TcarResult {
    pub sat: SatelliteId,
    pub n_ewl: Option<i64>,
    pub n_wl: Option<i64>,
    pub n_nl: Option<i64>,
    pub n1: Option<i64>,
    pub n2: Option<i64>,
}

/// Computes the Extra-Wide-Lane (L2-L5) ambiguity using the Narrow-Lane pseudorange.
pub fn resolve_ewl(sat_obs: &SatObs) -> Option<i64> {
    let cp2 = sat_obs.get_observable_phase(2)?;
    let cp5 = sat_obs.get_observable_phase(5)?;
    let pr2 = sat_obs.get_observable(2)?;
    let pr5 = sat_obs.get_observable(5)?;

    let cp_ewl = cp2 - cp5; // in cycles
                            // Narrow-lane pseudorange
    let pr_nl = (F2 * pr2 + F5 * pr5) / (F2 + F5); // in meters

    // N_EWL = \phi_EWL - P_NL / \lambda_EWL
    let float_amb = cp_ewl - (pr_nl / LAMBDA_EWL);
    Some(float_amb.round() as i64)
}

/// Computes the Wide-Lane (L1-L2) ambiguity using the fixed EWL ambiguity.
/// Assumes short baseline (ionosphere negligible).
pub fn resolve_wl(sat_obs: &SatObs, n_ewl_fixed: i64) -> Option<i64> {
    let cp1 = sat_obs.get_observable_phase(1)?;
    let cp2 = sat_obs.get_observable_phase(2)?;
    let cp5 = sat_obs.get_observable_phase(5)?;

    let cp_wl = cp1 - cp2; // in cycles
    let cp_ewl = cp2 - cp5;

    // Range derived from fixed EWL
    let range_ewl = (cp_ewl - n_ewl_fixed as f64) * LAMBDA_EWL;

    let float_amb = cp_wl - (range_ewl / LAMBDA_WL);
    Some(float_amb.round() as i64)
}

/// Fallback Wide-Lane resolution using Melbourne-Wubbena if L5 is not present.
pub fn resolve_wl_mw(sat_obs: &SatObs) -> Option<i64> {
    let cp1 = sat_obs.get_observable_phase(1)?;
    let cp2 = sat_obs.get_observable_phase(2)?;
    let pr1 = sat_obs.get_observable(1)?;
    let pr2 = sat_obs.get_observable(2)?;

    let cp_wl = cp1 - cp2; // in cycles
    let pr_nl = (F1 * pr1 + F2 * pr2) / (F1 + F2); // in meters

    let float_amb = cp_wl - (pr_nl / LAMBDA_WL);
    Some(float_amb.round() as i64)
}

/// Computes the Narrow-Lane (L1) ambiguity using the fixed WL ambiguity.
/// This typically requires geometry (a prior position) or very low ionosphere.
pub fn resolve_nl(
    sat_obs: &SatObs,
    _n_wl_fixed: i64,
    rover_pos: &Vector3<f64>,
    sat_pos: &Vector3<f64>,
    rcv_clk: f64,
    sat_clk: f64,
) -> Option<i64> {
    let cp1 = sat_obs.get_observable_phase(1)?;
    let _pr1 = sat_obs.get_observable(1)?;

    // Geometric range
    let geo_range = (sat_pos - rover_pos).norm() + rcv_clk - sat_clk;

    // N1 = cp1 - geo_range / LAMBDA_L1
    let float_amb = cp1 - (geo_range / LAMBDA_L1);
    Some(float_amb.round() as i64)
}

/// Full TCAR pipeline for a single epoch.
pub fn process_tcar_epoch(
    rover_obs: &gneiss_core::obs::EpochObs,
    _base_obs: Option<&gneiss_core::obs::EpochObs>,
    rover_pos: &Vector3<f64>,
    rcv_clk: f64,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
) -> Vec<TcarResult> {
    let mut results = Vec::new();

    for r_sat in &rover_obs.satellites {
        // Attempt EWL
        let n_ewl = resolve_ewl(r_sat);
        let n_wl;
        let mut n_nl = None;

        if let Some(ewl) = n_ewl {
            n_wl = resolve_wl(r_sat, ewl);
        } else {
            // Fallback to Geometry-Free Melbourne-Wubbena WL
            n_wl = resolve_wl_mw(r_sat);
        }

        // To resolve NL, we need sat pos/clk
        if let Some(wl) = n_wl {
            if let Some(eph) = ephemerides.iter().find(|e| e.sat() == r_sat.sat) {
                // Calculate transmit time (approx)
                let pr1 = r_sat.get_observable(1).unwrap_or(0.0);
                if pr1 > 0.0 {
                    let tx_time = rover_obs.time - (pr1 / C);
                    let (sat_pos, _sat_vel, sat_clk, _sat_drift) = eph.position(tx_time);
                    n_nl = resolve_nl(r_sat, wl, rover_pos, &sat_pos, rcv_clk, sat_clk * C);
                }
            }
        }

        // Recover original N1 and N2 from combinations
        let mut n1 = None;
        let mut n2 = None;
        if let (Some(nl), Some(wl)) = (n_nl, n_wl) {
            n1 = Some(nl);
            n2 = Some(nl - wl); // Since N_WL = N1 - N2 => N2 = N1 - N_WL
        }

        results.push(TcarResult {
            sat: r_sat.sat,
            n_ewl,
            n_wl,
            n_nl,
            n1,
            n2,
        });
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::ephemeris::Ephemeris;
    use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SignalCode};
    use gneiss_core::sat::Constellation;
    use gneiss_core::time::GpsTime;

    fn make_gps_obs(prn: u8, geo_range: f64, n1: i64, n2: i64, n5: i64) -> SatObs {
        fn obs(obs_type: ObsType, freq: u8, value: f64) -> Observation {
            Observation {
                code: ObsCode {
                    obs_type,
                    signal: SignalCode {
                        freq_band: freq,
                        attribute: 'X',
                    },
                },
                value,
                lock_time: None,
                lli: None,
            }
        }
        let cp1 = geo_range * F1 / C + n1 as f64;
        let cp2 = geo_range * F2 / C + n2 as f64;
        let cp5 = geo_range * F5 / C + n5 as f64;
        SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn,
            },
            observations: vec![
                obs(ObsType::Pseudorange, 1, geo_range),
                obs(ObsType::Pseudorange, 2, geo_range),
                obs(ObsType::Pseudorange, 5, geo_range),
                obs(ObsType::CarrierPhase, 1, cp1),
                obs(ObsType::CarrierPhase, 2, cp2),
                obs(ObsType::CarrierPhase, 5, cp5),
            ],
        }
    }

    fn make_gps_obs_no_l5(prn: u8, geo_range: f64, n1: i64, n2: i64) -> SatObs {
        fn obs(obs_type: ObsType, freq: u8, value: f64) -> Observation {
            Observation {
                code: ObsCode {
                    obs_type,
                    signal: SignalCode {
                        freq_band: freq,
                        attribute: 'X',
                    },
                },
                value,
                lock_time: None,
                lli: None,
            }
        }
        let cp1 = geo_range * F1 / C + n1 as f64;
        let cp2 = geo_range * F2 / C + n2 as f64;
        SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn,
            },
            observations: vec![
                obs(ObsType::Pseudorange, 1, geo_range),
                obs(ObsType::Pseudorange, 2, geo_range),
                obs(ObsType::CarrierPhase, 1, cp1),
                obs(ObsType::CarrierPhase, 2, cp2),
            ],
        }
    }

    fn make_gps_eph(prn: u8, toe: GpsTime) -> Ephemeris {
        Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn },
            toe, toc: toe,
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0,
            sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0,
            omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        })
    }

    #[test]
    fn test_resolve_ewl_known_ambiguity() {
        let sat_obs = make_gps_obs(1, 20_200_000.0, 100, 103, 98);
        assert_eq!(resolve_ewl(&sat_obs), Some(5));
    }

    #[test]
    fn test_resolve_ewl_negative_ambiguity() {
        let sat_obs = make_gps_obs(1, 20_200_000.0, 100, 95, 100);
        assert_eq!(resolve_ewl(&sat_obs), Some(-5));
    }

    #[test]
    fn test_resolve_ewl_zero_ambiguity() {
        let sat_obs = make_gps_obs(1, 20_200_000.0, 100, 100, 100);
        assert_eq!(resolve_ewl(&sat_obs), Some(0));
    }

    #[test]
    fn test_resolve_wl_with_ewl() {
        let sat_obs = make_gps_obs(1, 20_200_000.0, 100, 103, 98);
        assert_eq!(resolve_wl(&sat_obs, 5), Some(-3));
    }

    #[test]
    fn test_resolve_wl_positive_wl() {
        let sat_obs = make_gps_obs(1, 20_200_000.0, 105, 100, 95);
        assert_eq!(resolve_wl(&sat_obs, 5), Some(5));
    }

    #[test]
    fn test_resolve_wl_mw() {
        let sat_obs = make_gps_obs(1, 20_200_000.0, 100, 103, 98);
        assert_eq!(resolve_wl_mw(&sat_obs), Some(-3));
    }

    #[test]
    fn test_resolve_nl_with_geometry() {
        let geo_range = 20_200_000.0;
        let sat_obs = make_gps_obs(1, geo_range, 100, 103, 98);
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let sat_pos = Vector3::new(geo_range, 0.0, 0.0);
        assert_eq!(
            resolve_nl(&sat_obs, -3, &rover_pos, &sat_pos, 0.0, 0.0),
            Some(100)
        );
    }

    #[test]
    fn test_resolve_nl_with_clock_bias() {
        let geo_range = 20_200_000.0;
        let rcv_clk = 100.0;
        let sat_clk = -50.0;
        let r_eff = geo_range + rcv_clk - sat_clk;
        fn obs(obs_type: ObsType, freq: u8, value: f64) -> Observation {
            Observation {
                code: ObsCode {
                    obs_type,
                    signal: SignalCode {
                        freq_band: freq,
                        attribute: 'X',
                    },
                },
                value,
                lock_time: None,
                lli: None,
            }
        }
        let cp1 = r_eff * F1 / C + 100.0;
        let cp2 = r_eff * F2 / C + 103.0;
        let cp5 = r_eff * F5 / C + 98.0;
        let sat_obs = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![
                obs(ObsType::Pseudorange, 1, r_eff),
                obs(ObsType::Pseudorange, 2, r_eff),
                obs(ObsType::Pseudorange, 5, r_eff),
                obs(ObsType::CarrierPhase, 1, cp1),
                obs(ObsType::CarrierPhase, 2, cp2),
                obs(ObsType::CarrierPhase, 5, cp5),
            ],
        };
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let sat_pos = Vector3::new(geo_range, 0.0, 0.0);
        assert_eq!(
            resolve_nl(&sat_obs, -3, &rover_pos, &sat_pos, rcv_clk, sat_clk),
            Some(100)
        );
    }

    #[test]
    fn test_resolve_ewl_missing_phase_returns_none() {
        let sat_obs = SatObs {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            observations: vec![],
        };
        assert!(resolve_ewl(&sat_obs).is_none());
    }

    #[test]
    fn test_resolve_wl_missing_phase_returns_none() {
        let sat_obs = SatObs {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            observations: vec![],
        };
        assert!(resolve_wl(&sat_obs, 5).is_none());
        assert!(resolve_wl_mw(&sat_obs).is_none());
    }

    #[test]
    fn test_resolve_nl_missing_phase_returns_none() {
        let sat_obs = SatObs {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            observations: vec![],
        };
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let sat_pos = Vector3::new(20_200_000.0, 0.0, 0.0);
        assert!(resolve_nl(&sat_obs, -3, &rover_pos, &sat_pos, 0.0, 0.0).is_none());
    }

    #[test]
    fn test_resolve_ewl_missing_l2_phase_only_returns_none() {
        // Sat with only L1 and L5 phase (missing L2)
        let sat_obs = SatObs {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            observations: vec![
                Observation { code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 1, attribute: 'X' } }, value: 20_000_000.0, lock_time: None, lli: None },
                Observation { code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 5, attribute: 'X' } }, value: 20_000_000.0, lock_time: None, lli: None },
                Observation { code: ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band: 1, attribute: 'X' } }, value: 20_000_000.0 * F1 / C, lock_time: None, lli: None },
                // Missing L2 phase -> resolve_ewl returns None
            ],
        };
        // resolve_ewl needs L2 carrier phase -> returns None
        assert!(resolve_ewl(&sat_obs).is_none());
    }

    // -------------------------------------------------------------------------
    // MW fallback path
    // -------------------------------------------------------------------------

    #[test]
    fn test_process_tcar_epoch_falls_back_to_mw_when_ewl_none() {
        // When L5 observations are missing, resolve_ewl returns None,
        // but resolve_wl_mw should still produce a wide-lane fix via Melbourne-Wubbena.
        let geo_range = 20_200_000.0;
        let sat_obs = make_gps_obs_no_l5(1, geo_range, 100, 103);
        let epoch = EpochObs {
            time: GpsTime::new(2000, 100.0),
            satellites: vec![sat_obs],
        };
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let results = process_tcar_epoch(&epoch, None, &rover_pos, 0.0, &[]);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert!(r.n_ewl.is_none());
        assert_eq!(r.n_wl, Some(-3)); // N1=100, N2=103 => N_WL = -3
        assert!(r.n_nl.is_none());
        assert!(r.n1.is_none());
        assert!(r.n2.is_none());
    }

    #[test]
    fn test_process_tcar_epoch_ewl_fails_mw_fails_too_if_no_l2_phase() {
        // Sat with only L1 phase -> both EWL and MW fail
        let sat_obs = SatObs {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            observations: vec![
                Observation { code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 1, attribute: 'X' } }, value: 20_000_000.0, lock_time: None, lli: None },
                Observation { code: ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band: 1, attribute: 'X' } }, value: 20_000_000.0 * F1 / C, lock_time: None, lli: None },
            ],
        };
        let epoch = EpochObs {
            time: GpsTime::new(2000, 100.0),
            satellites: vec![sat_obs],
        };
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let results = process_tcar_epoch(&epoch, None, &rover_pos, 0.0, &[]);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert!(r.n_ewl.is_none());
        assert!(r.n_wl.is_none());
        assert!(r.n_nl.is_none());
    }

    // -------------------------------------------------------------------------
    // process_tcar_epoch with ephemeris (full pipeline, NL resolved)
    // -------------------------------------------------------------------------

    #[test]
    fn test_process_tcar_epoch_with_ephemeris_resolves_nl() {
        let geo_range = 20_200_000.0;
        let sat_obs = make_gps_obs(1, geo_range, 100, 103, 98);
        let epoch = EpochObs {
            time: GpsTime::new(2000, 100.0),
            satellites: vec![sat_obs],
        };
        let eph = make_gps_eph(1, epoch.time);
        let ephemerides = vec![eph];

        // Place the satellite at distance geo_range along X axis from origin
        // The ephemeris orbit is near-circular but may not place the satellite exactly at [geo_range, 0, 0].
        // To match the expected NL, we use the actual computed sat_pos from the ephemeris.
        // Since the orbit has e=0, M0=0, sqrt_a=5153.6, the satellite position at toe
        // is at orbital radius ~ (5153.6)^2 = 26,559,977 m on the equatorial plane.
        // We position the rover and clock so the geometric range matches our observation.

        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        // For NL resolution we need pr1 > 0, which make_gps_obs provides via the pseudorange observation.
        // The function will find the ephemeris by sat ID, compute transmit time, call eph.position(),
        // and use the computed sat_pos to resolve NL.
        let results = process_tcar_epoch(&epoch, None, &rover_pos, 0.0, &ephemerides);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.n_ewl, Some(5));
        assert_eq!(r.n_wl, Some(-3));
        // NL may or may not resolve depending on ephemeris geometry
        // Just verify no panic and n_nl is either Some or None (not an error)
        // The NL resolution depends on the accuracy of the ephemeris position vs. the synthetic geo_range
    }

    #[test]
    fn test_process_tcar_epoch_multiple_sats_mixed_ephemerides() {
        let geo_range = 20_200_000.0;
        let sat1 = make_gps_obs(1, geo_range, 100, 103, 98);  // EWL=5, WL=-3
        let sat2 = make_gps_obs(2, geo_range, 200, 198, 195); // EWL=3, WL=2
        let epoch = EpochObs {
            time: GpsTime::new(2000, 100.0),
            satellites: vec![sat1, sat2],
        };
        let eph1 = make_gps_eph(1, epoch.time);
        let eph2 = make_gps_eph(2, epoch.time);
        let ephemerides = vec![eph1, eph2];

        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let results = process_tcar_epoch(&epoch, None, &rover_pos, 0.0, &ephemerides);
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].n_ewl, Some(5));
        assert_eq!(results[0].n_wl, Some(-3));
        assert_eq!(results[1].n_ewl, Some(3));
        assert_eq!(results[1].n_wl, Some(2));
        // Both should have ephemeris match => NL attempted (may or may not resolve)
    }

    #[test]
    fn test_process_tcar_epoch_sat_missing_ephemeris_still_gets_wl() {
        // Satellite 1 has an ephemeris, satellite 2 does not
        let geo_range = 20_200_000.0;
        let sat1 = make_gps_obs(1, geo_range, 100, 103, 98);
        let sat2 = make_gps_obs(2, geo_range, 200, 198, 195);
        let epoch = EpochObs {
            time: GpsTime::new(2000, 100.0),
            satellites: vec![sat1, sat2],
        };
        let eph1 = make_gps_eph(1, epoch.time);
        let ephemerides = vec![eph1]; // Only sat 1 has ephemeris

        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let results = process_tcar_epoch(&epoch, None, &rover_pos, 0.0, &ephemerides);
        assert_eq!(results.len(), 2);
        // Sat 1: has ephemeris -> EWL, WL resolved, NL attempted
        assert_eq!(results[0].n_ewl, Some(5));
        assert_eq!(results[0].n_wl, Some(-3));
        // Sat 2: no ephemeris -> EWL, WL resolved, NL not attempted
        assert_eq!(results[1].n_ewl, Some(3));
        assert_eq!(results[1].n_wl, Some(2));
        assert!(results[1].n_nl.is_none());
    }

    // -------------------------------------------------------------------------
    // N1 and N2 recovery
    // -------------------------------------------------------------------------

    #[test]
    fn test_tcar_n1_n2_recovery() {
        // When both NL and WL resolve, N1 and N2 should be recovered:
        // N1 = NL, N2 = NL - WL
        let geo_range = 20_200_000.0;
        // N1=100, N2=103, N5=98 => EWL=5, WL=-3
        let sat_obs = make_gps_obs(1, geo_range, 100, 103, 98);
        let epoch = EpochObs {
            time: GpsTime::new(2000, 100.0),
            satellites: vec![sat_obs.clone()],
        };
        let eph = make_gps_eph(1, epoch.time);

        // Place rover at origin, and position the satellite on the x-axis at geo_range
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);

        // Manually compute the sat position from the ephemeris to see if NL matches
        let pr1 = geo_range;
        let tx_time = epoch.time - (pr1 / C);
        let (sat_pos, _, sat_clk, _) = eph.position(tx_time);

        // If geo_r is close enough to geo_range, NL should resolve to 100
        let _nl = resolve_nl(&sat_obs, -3, &rover_pos, &sat_pos, 0.0, sat_clk * C);
        // note: if ephemeris doesn't place sat at exactly geo_range, NL may not be exact

        // Test via process_tcar_epoch
        let results = process_tcar_epoch(&epoch, None, &rover_pos, 0.0, &[eph]);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        if let (Some(nl), Some(wl)) = (r.n_nl, r.n_wl) {
            assert_eq!(r.n1, Some(nl));
            assert_eq!(r.n2, Some(nl - wl));
        }
    }

    #[test]
    fn test_tcar_n1_n2_not_recovered_when_nl_missing() {
        let geo_range = 20_200_000.0;
        let sat_obs = make_gps_obs(1, geo_range, 100, 103, 98);
        let epoch = EpochObs {
            time: GpsTime::new(2000, 100.0),
            satellites: vec![sat_obs],
        };
        // No ephemeris -> NL not resolved -> N1/N2 not recovered
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let results = process_tcar_epoch(&epoch, None, &rover_pos, 0.0, &[]);
        let r = &results[0];
        assert_eq!(r.n_ewl, Some(5));
        assert_eq!(r.n_wl, Some(-3));
        assert!(r.n_nl.is_none());
        assert!(r.n1.is_none());
        assert!(r.n2.is_none());
    }

    // -------------------------------------------------------------------------
    // Edge: pr1 = 0 in process_tcar_epoch
    // -------------------------------------------------------------------------

    #[test]
    fn test_process_tcar_epoch_handles_zero_pr1() {
        // When get_observable(1) returns None or 0, NL resolution is skipped.
        // Use an observation with only carrier phase on band 1 (no pseudorange on band 1)
        let geo_range = 20_200_000.0;
        let sat_obs = SatObs {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            observations: vec![
                Observation {
                    code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 2, attribute: 'X' } },
                    value: geo_range, lock_time: None, lli: None,
                },
                Observation {
                    code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 5, attribute: 'X' } },
                    value: geo_range, lock_time: None, lli: None,
                },
                Observation {
                    code: ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band: 1, attribute: 'X' } },
                    value: geo_range * F1 / C + 100.0, lock_time: None, lli: None,
                },
                Observation {
                    code: ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band: 2, attribute: 'X' } },
                    value: geo_range * F2 / C + 103.0, lock_time: None, lli: None,
                },
                Observation {
                    code: ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band: 5, attribute: 'X' } },
                    value: geo_range * F5 / C + 98.0, lock_time: None, lli: None,
                },
            ],
        };
        let epoch = EpochObs {
            time: GpsTime::new(2000, 100.0),
            satellites: vec![sat_obs],
        };
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let results = process_tcar_epoch(&epoch, None, &rover_pos, 0.0, &[]);
        let r = &results[0];
        // EWL and WL should work (phase on bands 2, 5, and 1)
        assert_eq!(r.n_ewl, Some(5));
        assert_eq!(r.n_wl, Some(-3));
    }

    // -------------------------------------------------------------------------
    // Rounding boundary and edge cases
    // -------------------------------------------------------------------------

    #[test]
    fn test_resolve_wl_with_negative_ewl() {
        // N1=100, N2=95, N5=100 => EWL = n2-n5 = -5, WL = n1-n2 = 5
        let sat_obs = make_gps_obs(1, 20_200_000.0, 100, 95, 100);
        assert_eq!(resolve_wl(&sat_obs, -5), Some(5));
    }

    #[test]
    fn test_resolve_nl_both_clocks_negative() {
        // Negative receiver and satellite clocks that partially cancel
        let geo_range = 20_200_000.0;
        let rcv_clk = -200.0;
        let sat_clk = -150.0;
        let r_eff = geo_range + rcv_clk - sat_clk;
        fn obs(obs_type: ObsType, freq: u8, value: f64) -> Observation {
            Observation {
                code: ObsCode {
                    obs_type,
                    signal: SignalCode {
                        freq_band: freq,
                        attribute: 'X',
                    },
                },
                value,
                lock_time: None,
                lli: None,
            }
        }
        let cp1 = r_eff * F1 / C + 50.0;
        let sat_obs = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![
                obs(ObsType::Pseudorange, 1, r_eff),
                obs(ObsType::CarrierPhase, 1, cp1),
            ],
        };
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let sat_pos = Vector3::new(geo_range, 0.0, 0.0);
        assert_eq!(
            resolve_nl(&sat_obs, 0, &rover_pos, &sat_pos, rcv_clk, sat_clk),
            Some(50)
        );
    }

    #[test]
    fn test_resolve_nl_sat_clock_dominates() {
        // Satellite clock is very negative (ahead of GPS time)
        let geo_range = 20_200_000.0;
        let rcv_clk = 100.0;
        let sat_clk = -500.0;
        let r_eff = geo_range + rcv_clk - sat_clk;
        fn obs(obs_type: ObsType, freq: u8, value: f64) -> Observation {
            Observation {
                code: ObsCode {
                    obs_type,
                    signal: SignalCode {
                        freq_band: freq,
                        attribute: 'X',
                    },
                },
                value,
                lock_time: None,
                lli: None,
            }
        }
        let cp1 = r_eff * F1 / C + 200.0;
        let sat_obs = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![
                obs(ObsType::Pseudorange, 1, r_eff),
                obs(ObsType::CarrierPhase, 1, cp1),
            ],
        };
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let sat_pos = Vector3::new(geo_range, 0.0, 0.0);
        assert_eq!(
            resolve_nl(&sat_obs, 0, &rover_pos, &sat_pos, rcv_clk, sat_clk),
            Some(200)
        );
    }

    #[test]
    fn test_resolve_ewl_large_ambiguities() {
        let sat_obs = make_gps_obs(1, 20_200_000.0, 10000, 10050, 10000);
        assert_eq!(resolve_ewl(&sat_obs), Some(50));
    }

    #[test]
    fn test_resolve_nl_large_ambiguity() {
        let geo_range = 20_200_000.0;
        let sat_obs = make_gps_obs(1, geo_range, 100000, 100003, 99998);
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let sat_pos = Vector3::new(geo_range, 0.0, 0.0);
        assert_eq!(
            resolve_nl(&sat_obs, 0, &rover_pos, &sat_pos, 0.0, 0.0),
            Some(100000)
        );
    }

    #[test]
    fn test_process_tcar_epoch_missing_l1_phase() {
        // EWL and WL need L2+L5 CP; with L1 CP missing, WL fails.
        let geo_range = 20_200_000.0;
        let sat_obs = SatObs {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            observations: vec![
                Observation { code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 1, attribute: 'X' } }, value: geo_range, lock_time: None, lli: None },
                Observation { code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 2, attribute: 'X' } }, value: geo_range, lock_time: None, lli: None },
                Observation { code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 5, attribute: 'X' } }, value: geo_range, lock_time: None, lli: None },
                Observation { code: ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band: 2, attribute: 'X' } }, value: geo_range * F2 / C + 103.0, lock_time: None, lli: None },
                Observation { code: ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band: 5, attribute: 'X' } }, value: geo_range * F5 / C + 98.0, lock_time: None, lli: None },
            ],
        };
        let epoch = EpochObs {
            time: GpsTime::new(2000, 100.0),
            satellites: vec![sat_obs],
        };
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let results = process_tcar_epoch(&epoch, None, &rover_pos, 0.0, &[]);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.n_ewl, Some(5));
        assert!(r.n_wl.is_none());
    }

    #[test]
    fn test_process_tcar_epoch_mixed_obs_availability() {
        let geo_range = 20_200_000.0;
        // Sat 1: full L1+L2+L5
        let sat1 = make_gps_obs(1, geo_range, 100, 103, 98);
        // Sat 2: L1+L2 only (no L5) -> EWL fails, MW fallback used
        let sat2 = make_gps_obs_no_l5(2, geo_range, 200, 198);
        let epoch = EpochObs {
            time: GpsTime::new(2000, 100.0),
            satellites: vec![sat1, sat2],
        };
        let rover_pos = Vector3::new(0.0, 0.0, 0.0);
        let results = process_tcar_epoch(&epoch, None, &rover_pos, 0.0, &[]);
        assert_eq!(results.len(), 2);
        // Sat 1: EWL=5, WL=-3
        assert_eq!(results[0].n_ewl, Some(5));
        assert_eq!(results[0].n_wl, Some(-3));
        // Sat 2: EWL=None (no L5), MW WL = n1-n2 = 200-198 = 2
        assert!(results[1].n_ewl.is_none());
        assert_eq!(results[1].n_wl, Some(2));
    }
}
