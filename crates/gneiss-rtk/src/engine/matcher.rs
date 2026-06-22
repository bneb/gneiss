use crate::filter::DdObservation;
use gneiss_core::obs::{EpochObs, ObsType};

pub fn match_observations(
    rover_obs: &EpochObs,
    base_obs: &EpochObs,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
) -> Vec<(DdObservation, DdObservation)> {
    let mut matched_obs = Vec::new();
    for r_sat in &rover_obs.satellites {
        if !ephemerides.iter().any(|e| e.sat() == r_sat.sat) {
            continue;
        }
        if let Some(b_sat) = base_obs.satellites.iter().find(|s| s.sat == r_sat.sat) {
            let r_pr_l1 = r_sat
                .observations
                .iter()
                .find(|o| o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == 1);
            let r_pr_l2 = r_sat
                .observations
                .iter()
                .find(|o| o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == 2);
            let r_cp_l1 = r_sat
                .observations
                .iter()
                .find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 1);
            let r_cp_l2 = r_sat
                .observations
                .iter()
                .find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 2);
            let r_dop = r_sat
                .observations
                .iter()
                .find(|o| o.code.obs_type == ObsType::Doppler && o.code.signal.freq_band == 1)
                .map(|o| o.value)
                .unwrap_or(0.0);
            let b_pr_l1 = b_sat
                .observations
                .iter()
                .find(|o| o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == 1);
            let b_pr_l2 = b_sat
                .observations
                .iter()
                .find(|o| o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == 2);
            let b_cp_l1 = b_sat
                .observations
                .iter()
                .find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 1);
            let b_cp_l2 = b_sat
                .observations
                .iter()
                .find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 2);
            let b_dop = b_sat
                .observations
                .iter()
                .find(|o| o.code.obs_type == ObsType::Doppler && o.code.signal.freq_band == 1)
                .map(|o| o.value)
                .unwrap_or(0.0);
            let r_snr = r_sat
                .observations
                .iter()
                .find(|o| o.code.obs_type == ObsType::Snr && o.code.signal.freq_band == 1)
                .map(|o| o.value)
                .unwrap_or(25.0);
            let r_lock = r_sat
                .observations
                .iter()
                .find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 1)
                .and_then(|o| o.lock_time);

            if r_pr_l1.is_none() {
                tracing::debug!(
                    "Sat {:?} missing rover PR1. Rover Obs: {:?}",
                    r_sat.sat,
                    r_sat
                        .observations
                        .iter()
                        .map(|o| o.code)
                        .collect::<Vec<_>>()
                );
            }
            if b_pr_l1.is_none() {
                tracing::debug!(
                    "Sat {:?} missing base PR1. Base Obs: {:?}",
                    b_sat.sat,
                    b_sat
                        .observations
                        .iter()
                        .map(|o| o.code)
                        .collect::<Vec<_>>()
                );
            }
            if r_cp_l1.is_none() {
                tracing::debug!("Sat {:?} missing rover CP1", r_sat.sat);
            }
            if b_cp_l1.is_none() {
                tracing::debug!(
                    "Sat {:?} missing base CP1. Base Obs: {:?}",
                    b_sat.sat,
                    b_sat
                        .observations
                        .iter()
                        .map(|o| o.code)
                        .collect::<Vec<_>>()
                );
            }

            if let (Some(r_pr1), Some(b_pr1)) = (r_pr_l1, b_pr_l1) {
                tracing::debug!(
                    "Sat {:?} L1 PR: Rover {} (val: {}), Base {} (val: {})",
                    r_sat.sat,
                    r_pr1.code,
                    r_pr1.value,
                    b_pr1.code,
                    b_pr1.value
                );
                if let (Some(r_pr2), Some(b_pr2)) = (r_pr_l2, b_pr_l2) {
                    tracing::debug!(
                        "Sat {:?} L2 PR: Rover {} (val: {}), Base {} (val: {})",
                        r_sat.sat,
                        r_pr2.code,
                        r_pr2.value,
                        b_pr2.code,
                        b_pr2.value
                    );
                }
                matched_obs.push((
                    DdObservation {
                        sat: r_sat.sat,
                        pr_l1: r_pr1.value,
                        pr_l2: r_pr_l2.map(|o| o.value),
                        cp_l1: r_cp_l1.map(|o| o.value),
                        cp_l2: r_cp_l2.map(|o| o.value),
                        doppler: r_dop,
                        snr: r_snr,
                        locktime: r_lock,
                    },
                    DdObservation {
                        sat: b_sat.sat,
                        pr_l1: b_pr1.value,
                        pr_l2: b_pr_l2.map(|o| o.value),
                        cp_l1: b_cp_l1.map(|o| o.value),
                        cp_l2: b_cp_l2.map(|o| o.value),
                        doppler: b_dop,
                        snr: 25.0,
                        locktime: Some(1000),
                    },
                ));
            }
        }
    }
    matched_obs
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::ephemeris::{Ephemeris, GalileoEphemeris, GpsEphemeris};
    use gneiss_core::obs::{Observation, ObsCode, ObsType, SatObs, SignalCode};
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    fn make_obs_code(obs_type: ObsType, freq_band: u8, attribute: char) -> ObsCode {
        ObsCode {
            obs_type,
            signal: SignalCode {
                freq_band,
                attribute,
            },
        }
    }

    fn make_obs(val: f64, obs_type: ObsType, freq_band: u8, lock_time: Option<u16>) -> Observation {
        Observation {
            code: make_obs_code(obs_type, freq_band, 'C'),
            value: val,
            lock_time,
            lli: None,
        }
    }

    fn make_full_obs(sat: SatelliteId, pr_l1: f64, pr_l2: f64, cp_l1: f64, cp_l2: f64) -> SatObs {
        SatObs {
            sat,
            observations: vec![
                make_obs(pr_l1, ObsType::Pseudorange, 1, None),
                make_obs(pr_l2, ObsType::Pseudorange, 2, None),
                make_obs(cp_l1, ObsType::CarrierPhase, 1, Some(100)),
                make_obs(cp_l2, ObsType::CarrierPhase, 2, None),
                make_obs(100.0, ObsType::Doppler, 1, None),
                make_obs(45.0, ObsType::Snr, 1, None),
            ],
        }
    }

    fn make_ephemeris(sat: SatelliteId) -> Ephemeris {
        Ephemeris::Gps(GpsEphemeris {
            sat,
            toe: GpsTime::new(2000, 0.0),
            toc: GpsTime::new(2000, 0.0),
            af0: 0.0,
            af1: 0.0,
            af2: 0.0,
            crs: 0.0,
            crc: 0.0,
            cuc: 0.0,
            cus: 0.0,
            cic: 0.0,
            cis: 0.0,
            m0: 0.0,
            e: 0.0,
            sqrt_a: 0.0,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.0,
            idot: 0.0,
            omega: 0.0,
            tgd: 0.0,
            iode: 0,
            iodc: 0,
        })
    }

    fn make_epoch(satellites: Vec<SatObs>, week: u32, tow: f64) -> EpochObs {
        EpochObs {
            time: GpsTime::new(week, tow),
            satellites,
        }
    }

    fn sat(constellation: Constellation, prn: u8) -> SatelliteId {
        SatelliteId { constellation, prn }
    }

    // -----------------------------------------------------------------------
    // Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_basic_match() {
        let rover = make_epoch(
            vec![
                make_full_obs(sat(Constellation::Gps, 1), 100.0, 200.0, 300.0, 400.0),
                make_full_obs(sat(Constellation::Gps, 3), 500.0, 600.0, 700.0, 800.0),
            ],
            2000,
            1.0,
        );
        let base = make_epoch(
            vec![
                make_full_obs(sat(Constellation::Gps, 1), 101.0, 201.0, 301.0, 401.0),
                make_full_obs(sat(Constellation::Gps, 3), 501.0, 601.0, 701.0, 801.0),
            ],
            2000,
            1.0,
        );
        let ephemerides = vec![
            make_ephemeris(sat(Constellation::Gps, 1)),
            make_ephemeris(sat(Constellation::Gps, 3)),
        ];

        let result = match_observations(&rover, &base, &ephemerides);
        assert_eq!(result.len(), 2);

        // First match: G01
        assert_eq!(result[0].0.sat, sat(Constellation::Gps, 1));
        assert_eq!(result[0].0.pr_l1, 100.0);
        assert_eq!(result[0].0.pr_l2, Some(200.0));
        assert_eq!(result[0].0.cp_l1, Some(300.0));
        assert_eq!(result[0].0.cp_l2, Some(400.0));
        assert_eq!(result[0].0.doppler, 100.0);
        assert_eq!(result[0].0.snr, 45.0);
        assert_eq!(result[0].0.locktime, Some(100));

        assert_eq!(result[0].1.sat, sat(Constellation::Gps, 1));
        assert_eq!(result[0].1.pr_l1, 101.0);
        assert_eq!(result[0].1.pr_l2, Some(201.0));

        // Second match: G03
        assert_eq!(result[1].0.sat, sat(Constellation::Gps, 3));
        assert_eq!(result[1].0.pr_l1, 500.0);
        assert_eq!(result[1].1.pr_l1, 501.0);
    }

    #[test]
    fn test_empty_rover() {
        let rover = make_epoch(vec![], 2000, 1.0);
        let base = make_epoch(
            vec![make_full_obs(sat(Constellation::Gps, 1), 101.0, 201.0, 301.0, 401.0)],
            2000,
            1.0,
        );
        let ephemerides = vec![make_ephemeris(sat(Constellation::Gps, 1))];

        let result = match_observations(&rover, &base, &ephemerides);
        assert!(result.is_empty());
    }

    #[test]
    fn test_empty_base() {
        let rover = make_epoch(
            vec![make_full_obs(sat(Constellation::Gps, 1), 100.0, 200.0, 300.0, 400.0)],
            2000,
            1.0,
        );
        let base = make_epoch(vec![], 2000, 1.0);
        let ephemerides = vec![make_ephemeris(sat(Constellation::Gps, 1))];

        let result = match_observations(&rover, &base, &ephemerides);
        assert!(result.is_empty());
    }

    #[test]
    fn test_no_ephemeris() {
        let rover = make_epoch(
            vec![make_full_obs(sat(Constellation::Gps, 1), 100.0, 200.0, 300.0, 400.0)],
            2000,
            1.0,
        );
        let base = make_epoch(
            vec![make_full_obs(sat(Constellation::Gps, 1), 101.0, 201.0, 301.0, 401.0)],
            2000,
            1.0,
        );
        // No ephemerides at all
        let ephemerides: Vec<Ephemeris> = vec![];

        let result = match_observations(&rover, &base, &ephemerides);
        assert!(result.is_empty());
    }

    #[test]
    fn test_no_base_match() {
        let rover = make_epoch(
            vec![make_full_obs(sat(Constellation::Gps, 1), 100.0, 200.0, 300.0, 400.0)],
            2000,
            1.0,
        );
        let base = make_epoch(
            // G03 in base, not G01
            vec![make_full_obs(sat(Constellation::Gps, 3), 501.0, 601.0, 701.0, 801.0)],
            2000,
            1.0,
        );
        let ephemerides = vec![
            make_ephemeris(sat(Constellation::Gps, 1)),
            make_ephemeris(sat(Constellation::Gps, 3)),
        ];

        let result = match_observations(&rover, &base, &ephemerides);
        assert!(result.is_empty());
    }

    #[test]
    fn test_rover_missing_pr_l1() {
        let rover_sat = SatObs {
            sat: sat(Constellation::Gps, 1),
            observations: vec![
                // Only L2 pseudorange, no L1
                make_obs(200.0, ObsType::Pseudorange, 2, None),
                make_obs(300.0, ObsType::CarrierPhase, 1, Some(100)),
            ],
        };
        let rover = make_epoch(vec![rover_sat], 2000, 1.0);
        let base = make_epoch(
            vec![make_full_obs(sat(Constellation::Gps, 1), 101.0, 201.0, 301.0, 401.0)],
            2000,
            1.0,
        );
        let ephemerides = vec![make_ephemeris(sat(Constellation::Gps, 1))];

        let result = match_observations(&rover, &base, &ephemerides);
        assert!(result.is_empty());
    }

    #[test]
    fn test_base_missing_pr_l1() {
        let rover = make_epoch(
            vec![make_full_obs(sat(Constellation::Gps, 1), 100.0, 200.0, 300.0, 400.0)],
            2000,
            1.0,
        );
        let base_sat = SatObs {
            sat: sat(Constellation::Gps, 1),
            observations: vec![
                make_obs(201.0, ObsType::Pseudorange, 2, None),
                make_obs(301.0, ObsType::CarrierPhase, 1, Some(100)),
            ],
        };
        let base = make_epoch(vec![base_sat], 2000, 1.0);
        let ephemerides = vec![make_ephemeris(sat(Constellation::Gps, 1))];

        let result = match_observations(&rover, &base, &ephemerides);
        assert!(result.is_empty());
    }

    #[test]
    fn test_duplicate_rover_satellite() {
        // Two entries for G01 in rover; only the first pair should match
        // (base has only one G01, find will match the first)
        let rover = make_epoch(
            vec![
                make_full_obs(sat(Constellation::Gps, 1), 100.0, 200.0, 300.0, 400.0),
                make_full_obs(sat(Constellation::Gps, 1), 999.0, 888.0, 777.0, 666.0),
            ],
            2000,
            1.0,
        );
        let base = make_epoch(
            vec![make_full_obs(sat(Constellation::Gps, 1), 101.0, 201.0, 301.0, 401.0)],
            2000,
            1.0,
        );
        let ephemerides = vec![make_ephemeris(sat(Constellation::Gps, 1))];

        let result = match_observations(&rover, &base, &ephemerides);
        // Two rover entries for G01, both match base G01 → two output pairs
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].0.pr_l1, 100.0);
        assert_eq!(result[1].0.pr_l1, 999.0);
    }

    #[test]
    fn test_mixed_constellations() {
        let rover = make_epoch(
            vec![
                make_full_obs(sat(Constellation::Gps, 1), 100.0, 200.0, 300.0, 400.0),
                make_full_obs(sat(Constellation::Galileo, 2), 500.0, 600.0, 700.0, 800.0),
            ],
            2000,
            1.0,
        );
        let base = make_epoch(
            vec![
                make_full_obs(sat(Constellation::Gps, 1), 101.0, 201.0, 301.0, 401.0),
                make_full_obs(sat(Constellation::Galileo, 2), 501.0, 601.0, 701.0, 801.0),
            ],
            2000,
            1.0,
        );
        let ephemerides = vec![
            make_ephemeris(sat(Constellation::Gps, 1)),
            Ephemeris::Galileo(GalileoEphemeris {
                sat: sat(Constellation::Galileo, 2),
                toe: GpsTime::new(2000, 0.0),
                toc: GpsTime::new(2000, 0.0),
                af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.0, sqrt_a: 0.0, delta_n: 0.0,
                omega0: 0.0, omega_dot: 0.0, i0: 0.0, idot: 0.0, omega: 0.0,
                bgd_e1_e5a: 0.0, bgd_e1_e5b: 0.0, iod_nav: 0,
            }),
        ];

        let result = match_observations(&rover, &base, &ephemerides);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].0.sat, sat(Constellation::Gps, 1));
        assert_eq!(result[1].0.sat, sat(Constellation::Galileo, 2));
    }

    #[test]
    fn test_partial_observations() {
        // Only PR L1 and CP L1 — no L2
        let rover_sat = SatObs {
            sat: sat(Constellation::Gps, 1),
            observations: vec![
                make_obs(100.0, ObsType::Pseudorange, 1, None),
                make_obs(300.0, ObsType::CarrierPhase, 1, Some(100)),
                make_obs(50.0, ObsType::Doppler, 1, None),
            ],
        };
        let base_sat = SatObs {
            sat: sat(Constellation::Gps, 1),
            observations: vec![
                make_obs(101.0, ObsType::Pseudorange, 1, None),
                make_obs(301.0, ObsType::CarrierPhase, 1, None),
                make_obs(51.0, ObsType::Doppler, 1, None),
            ],
        };
        let rover = make_epoch(vec![rover_sat], 2000, 1.0);
        let base = make_epoch(vec![base_sat], 2000, 1.0);
        let ephemerides = vec![make_ephemeris(sat(Constellation::Gps, 1))];

        let result = match_observations(&rover, &base, &ephemerides);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0.pr_l1, 100.0);
        assert_eq!(result[0].0.pr_l2, None);
        assert_eq!(result[0].0.cp_l1, Some(300.0));
        assert_eq!(result[0].0.cp_l2, None);
        // Doppler present
        assert_eq!(result[0].0.doppler, 50.0);
        // SNR missing → default 25.0
        assert_eq!(result[0].0.snr, 25.0);
        assert_eq!(result[0].0.locktime, Some(100));

        // Base side
        assert_eq!(result[0].1.pr_l1, 101.0);
        assert_eq!(result[0].1.pr_l2, None);
        assert_eq!(result[0].1.cp_l1, Some(301.0));
        assert_eq!(result[0].1.cp_l2, None);
        assert_eq!(result[0].1.doppler, 51.0);
        // Base SNR always 25.0
        assert_eq!(result[0].1.snr, 25.0);
        // Base locktime always Some(1000)
        assert_eq!(result[0].1.locktime, Some(1000));
    }

    #[test]
    fn test_doppler_default_zero() {
        let rover_sat = SatObs {
            sat: sat(Constellation::Gps, 1),
            observations: vec![
                make_obs(100.0, ObsType::Pseudorange, 1, None),
                make_obs(200.0, ObsType::Pseudorange, 2, None),
                make_obs(300.0, ObsType::CarrierPhase, 1, Some(100)),
            ],
        };
        let base_sat = SatObs {
            sat: sat(Constellation::Gps, 1),
            observations: vec![
                make_obs(101.0, ObsType::Pseudorange, 1, None),
                make_obs(201.0, ObsType::Pseudorange, 2, None),
                make_obs(301.0, ObsType::CarrierPhase, 1, None),
            ],
        };
        let rover = make_epoch(vec![rover_sat], 2000, 1.0);
        let base = make_epoch(vec![base_sat], 2000, 1.0);
        let ephemerides = vec![make_ephemeris(sat(Constellation::Gps, 1))];

        let result = match_observations(&rover, &base, &ephemerides);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0.doppler, 0.0);
    }

    #[test]
    fn test_epoch_time_fields() {
        // Different epochs (different week/tow) but same satellites → still match
        let rover = make_epoch(
            vec![make_full_obs(sat(Constellation::Gps, 1), 100.0, 200.0, 300.0, 400.0)],
            2001,
            0.0,
        );
        let base = make_epoch(
            vec![make_full_obs(sat(Constellation::Gps, 1), 101.0, 201.0, 301.0, 401.0)],
            2000,
            5.0,
        );
        let ephemerides = vec![make_ephemeris(sat(Constellation::Gps, 1))];

        let result = match_observations(&rover, &base, &ephemerides);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0.pr_l1, 100.0);
        assert_eq!(result[0].1.pr_l1, 101.0);
    }
}
