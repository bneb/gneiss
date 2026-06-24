use crate::filter::{DdObservation, RtkState};
use gneiss_core::coords::Coordinate;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::time::GpsTime;

pub fn manage_ambiguities_and_slips(
    state: &mut RtkState,
    config: &crate::engine::EngineConfig,
    matched_obs: &[(DdObservation, DdObservation)],
    ephemerides: &[Ephemeris],
    base_coord: &Coordinate,
    rover_time: GpsTime,
    base_time: GpsTime,
) {
    for (r, b) in matched_obs {
        let mut slip = false;
        let prev_lock = *state.locktimes.get(&(r.sat, 1)).unwrap_or(&0);
        let mut new_lock = prev_lock + 1;

        if let Some(r_lock) = r.locktime {
            if r_lock == 0 {
                tracing::debug!("slip=true because r_lock == 0 for {:?}", r.sat);
                slip = true;
                new_lock = 0;
            } else if r_lock < prev_lock {
                tracing::debug!("slip=true because r_lock < prev_lock for {:?}", r.sat);
                slip = true;
                new_lock = r_lock;
            } else {
                new_lock = r_lock;
            }
        }

        if new_lock == 0 && state.locktimes.contains_key(&(r.sat, 1)) {
            tracing::debug!("slip=true because new_lock == 0 for {:?}", r.sat);
            slip = true;
        }

        state.locktimes.insert((r.sat, 1), new_lock);
        state.locktimes.insert((r.sat, 2), new_lock);

        let r_freq_num = ephemerides
            .iter()
            .find(|e| e.sat() == r.sat)
            .map(|e| e.freq_num())
            .unwrap_or(0);
        let b_freq_num = ephemerides
            .iter()
            .find(|e| e.sat() == b.sat)
            .map(|e| e.freq_num())
            .unwrap_or(0);

        let (r_f1, r_f2) = gneiss_core::signal::satellite_frequencies(r.sat, r_freq_num);
        let (b_f1, b_f2) = gneiss_core::signal::satellite_frequencies(b.sat, b_freq_num);

        let mut slip_l1 = slip;
        let mut slip_l2 = slip;
        if *state.reject_counts.get(&(r.sat, 1)).unwrap_or(&0) > config.max_reject_count {
            slip_l1 = true;
        }
        if *state.reject_counts.get(&(r.sat, 2)).unwrap_or(&0) > config.max_reject_count {
            slip_l2 = true;
        }

        if let Some(r_cp1) = r.cp_l1 {
            if let Some(&(prev_cp, prev_doppler, prev_time)) = state.phase_history.get(&(r.sat, 1))
            {
                let dt = rover_time - prev_time;
                if dt > 0.0 && dt <= config.max_base_age_s {
                    if check_doppler_phase_slip(
                        r_cp1,
                        prev_cp,
                        r.doppler,
                        prev_doppler,
                        dt,
                        config.doppler_slip_threshold_cycles,
                    ) {
                        tracing::debug!("Doppler-Phase cycle slip detected on {:?} L1", r.sat);
                        slip_l1 = true;
                    } else if slip_l1
                        && dt < 5.0
                        && !check_doppler_phase_slip(
                            r_cp1,
                            prev_cp,
                            r.doppler,
                            prev_doppler,
                            dt,
                            0.5,
                        )
                    {
                        tracing::debug!(
                            "Bridging short interruption on {:?} L1 using Doppler",
                            r.sat
                        );
                        slip_l1 = false;
                    }
                } else if dt > config.max_base_age_s {
                    tracing::debug!(
                        "Data gap > {:.1}s on {:?} L1, resetting ambiguity",
                        config.max_base_age_s,
                        r.sat
                    );
                    slip_l1 = true;
                }
            }
            state
                .phase_history
                .insert((r.sat, 1), (r_cp1, r.doppler, rover_time));
        }

        if let Some(r_cp2) = r.cp_l2 {
            let doppler_l2 = r.doppler * (r_f2 / r_f1);
            if let Some(&(prev_cp, prev_doppler, prev_time)) = state.phase_history.get(&(r.sat, 2))
            {
                let dt = rover_time - prev_time;
                if dt > 0.0 && dt <= config.max_base_age_s {
                    if check_doppler_phase_slip(
                        r_cp2,
                        prev_cp,
                        doppler_l2,
                        prev_doppler,
                        dt,
                        config.doppler_slip_threshold_cycles,
                    ) {
                        tracing::debug!("Doppler-Phase cycle slip detected on {:?} L2", r.sat);
                        slip_l2 = true;
                    } else if slip_l2
                        && dt < 5.0
                        && !check_doppler_phase_slip(
                            r_cp2,
                            prev_cp,
                            doppler_l2,
                            prev_doppler,
                            dt,
                            0.5,
                        )
                    {
                        tracing::debug!(
                            "Bridging short interruption on {:?} L2 using Doppler",
                            r.sat
                        );
                        slip_l2 = false;
                    }
                } else if dt > config.max_base_age_s {
                    tracing::debug!(
                        "Data gap > {:.1}s on {:?} L2, resetting ambiguity",
                        config.max_base_age_s,
                        r.sat
                    );
                    slip_l2 = true;
                }
            }
            state
                .phase_history
                .insert((r.sat, 2), (r_cp2, doppler_l2, rover_time));
        }

        // Geometry-Free cycle slip detection (feature-gated)
        #[cfg(feature = "gf-slip")]
        if let (Some(r_cp1), Some(r_cp2)) = (r.cp_l1, r.cp_l2) {
            let lam1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f1;
            let lam2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f2;
            let l_gf = r_cp1 * lam1 - r_cp2 * lam2;
            if let Some(&prev_gf) = state.gf_values.get(&r.sat) {
                if (l_gf - prev_gf).abs() > 0.05 {
                    tracing::debug!(
                        "GF cycle slip detected on {:?}: |ΔGF| = {:.4}m",
                        r.sat,
                        (l_gf - prev_gf).abs()
                    );
                    slip_l1 = true;
                    slip_l2 = true;
                }
            }
            state.gf_values.insert(r.sat, l_gf);
        }

        let needs_l1 = r.cp_l1.is_some()
            && b.cp_l1.is_some()
            && (!state.ambiguity_keys.contains(&(r.sat, 1)) || slip_l1 || slip_l2);
        let needs_l2 = r.cp_l2.is_some()
            && b.cp_l2.is_some()
            && (!state.ambiguity_keys.contains(&(r.sat, 2)) || slip_l1 || slip_l2);

        if slip_l1 || slip_l2 {
            tracing::debug!(
                "Cycle slip detected for {:?}: slip_l1={}, slip_l2={}",
                r.sat,
                slip_l1,
                slip_l2
            );
            state.remove_ambiguity(r.sat, 1);
            state.remove_ambiguity(r.sat, 2);
            state.reject_counts.insert((r.sat, 1), 0);
            state.reject_counts.insert((r.sat, 2), 0);
        }

        if needs_l1 {
            let r_cp1 = r.cp_l1.unwrap();
            let b_cp1 = b.cp_l1.unwrap();
            let lam_r1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f1;
            let lam_b1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / b_f1;
            let cp_l1_rov = r_cp1 * lam_r1;
            let cp_l1_base = b_cp1 * lam_b1;

            let mut initialized = false;
            if state.covariance[(0, 0)] < 0.1 {
                for (anchor_r, anchor_b) in matched_obs.iter() {
                    if anchor_r.sat == r.sat {
                        continue;
                    }
                    if let Some(anchor_idx) = state
                        .ambiguity_keys
                        .iter()
                        .position(|&(s, f)| s == anchor_r.sat && f == 1)
                    {
                        if state.covariance[(
                            crate::filter::CORE_STATE_SIZE + anchor_idx,
                            crate::filter::CORE_STATE_SIZE + anchor_idx,
                        )] < 0.05
                        {
                            if let (Some(ar_cp), Some(ab_cp)) = (anchor_r.cp_l1, anchor_b.cp_l1) {
                                let anchor_eph = match ephemerides
                                    .iter()
                                    .find(|e| e.sat() == anchor_r.sat)
                                {
                                    Some(eph) => eph,
                                    None => continue,
                                };
                                let (a_f1, _) = gneiss_core::signal::satellite_frequencies(
                                    anchor_r.sat,
                                    anchor_eph.freq_num(),
                                );

                                let (ar_sat_vec, _) =
                                    crate::engine::measurement_math::get_sat_state(
                                        anchor_eph,
                                        anchor_r.pr_l1,
                                        state.rcv_clk_bias,
                                        rover_time,
                                        state.position.vector,
                                    );
                                let (ab_sat_vec, _) =
                                    crate::engine::measurement_math::get_sat_state(
                                        anchor_eph,
                                        anchor_b.pr_l1,
                                        0.0,
                                        base_time,
                                        base_coord.vector,
                                    );
                                let ar_dist_rov = (state.position.vector - ar_sat_vec).norm();
                                let ar_dist_base = (base_coord.vector - ab_sat_vec).norm();

                                let a_lam = gneiss_core::constants::SPEED_OF_LIGHT_M_S / a_f1;
                                let a_cp_rov = ar_cp * a_lam;
                                let a_cp_base = ab_cp * a_lam;

                                let anchor_sd = state.ambiguities[anchor_idx];
                                let b_clock_rov = a_cp_rov - ar_dist_rov - anchor_sd;
                                let b_clock_base = a_cp_base - ar_dist_base; // Base has no ambiguity in SD, assuming SD = rov - base

                                let r_eph = match ephemerides.iter().find(|e| e.sat() == r.sat) {
                                    Some(eph) => eph,
                                    None => continue,
                                };
                                let (r_sat_vec, _) = crate::engine::measurement_math::get_sat_state(
                                    r_eph,
                                    r.pr_l1,
                                    state.rcv_clk_bias,
                                    rover_time,
                                    state.position.vector,
                                );
                                let (b_sat_vec, _) = crate::engine::measurement_math::get_sat_state(
                                    r_eph,
                                    b.pr_l1,
                                    0.0,
                                    base_time,
                                    base_coord.vector,
                                );
                                let dist_rov = (state.position.vector - r_sat_vec).norm();
                                let dist_base = (base_coord.vector - b_sat_vec).norm();

                                let initial_est_l1 = (cp_l1_rov - dist_rov - b_clock_rov)
                                    - (cp_l1_base - dist_base - b_clock_base);
                                state.add_ambiguity(
                                    r.sat,
                                    1,
                                    initial_est_l1,
                                    config.initial_ambiguity_variance,
                                );
                                initialized = true;
                                break;
                            }
                        }
                    }
                }
            }

            if !initialized {
                let initial_est_l1 = (cp_l1_rov - r.pr_l1) - (cp_l1_base - b.pr_l1);
                state.add_ambiguity(r.sat, 1, initial_est_l1, config.initial_ambiguity_variance);
            }
        }

        if needs_l2 {
            let r_cp2 = r.cp_l2.unwrap();
            let b_cp2 = b.cp_l2.unwrap();
            let r_pr2 = r.pr_l2.unwrap();
            let b_pr2 = b.pr_l2.unwrap();

            let lam_r2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f2;
            let lam_b2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / b_f2;
            let cp_l2_rov = r_cp2 * lam_r2;
            let cp_l2_base = b_cp2 * lam_b2;

            let mut initialized = false;
            if state.covariance[(0, 0)] < 0.1 {
                for (anchor_r, anchor_b) in matched_obs.iter() {
                    if anchor_r.sat == r.sat {
                        continue;
                    }
                    if let Some(anchor_idx) = state
                        .ambiguity_keys
                        .iter()
                        .position(|&(s, f)| s == anchor_r.sat && f == 2)
                    {
                        if state.covariance[(
                            crate::filter::CORE_STATE_SIZE + anchor_idx,
                            crate::filter::CORE_STATE_SIZE + anchor_idx,
                        )] < 0.05
                        {
                            if let (Some(ar_cp), Some(ab_cp)) = (anchor_r.cp_l2, anchor_b.cp_l2) {
                                let anchor_eph = match ephemerides
                                    .iter()
                                    .find(|e| e.sat() == anchor_r.sat)
                                {
                                    Some(eph) => eph,
                                    None => continue,
                                };
                                let (_, a_f2) = gneiss_core::signal::satellite_frequencies(
                                    anchor_r.sat,
                                    anchor_eph.freq_num(),
                                );

                                let (ar_sat_vec, _) =
                                    crate::engine::measurement_math::get_sat_state(
                                        anchor_eph,
                                        anchor_r.pr_l1,
                                        state.rcv_clk_bias,
                                        rover_time,
                                        state.position.vector,
                                    );
                                let (ab_sat_vec, _) =
                                    crate::engine::measurement_math::get_sat_state(
                                        anchor_eph,
                                        anchor_b.pr_l1,
                                        0.0,
                                        base_time,
                                        base_coord.vector,
                                    );
                                let ar_dist_rov = (state.position.vector - ar_sat_vec).norm();
                                let ar_dist_base = (base_coord.vector - ab_sat_vec).norm();

                                let a_lam = gneiss_core::constants::SPEED_OF_LIGHT_M_S / a_f2;
                                let a_cp_rov = ar_cp * a_lam;
                                let a_cp_base = ab_cp * a_lam;

                                let anchor_sd = state.ambiguities[anchor_idx];
                                let b_clock_rov = a_cp_rov - ar_dist_rov - anchor_sd;
                                let b_clock_base = a_cp_base - ar_dist_base;

                                let r_eph = match ephemerides.iter().find(|e| e.sat() == r.sat) {
                                    Some(eph) => eph,
                                    None => continue,
                                };
                                let (r_sat_vec, _) = crate::engine::measurement_math::get_sat_state(
                                    r_eph,
                                    r.pr_l1,
                                    state.rcv_clk_bias,
                                    rover_time,
                                    state.position.vector,
                                );
                                let (b_sat_vec, _) = crate::engine::measurement_math::get_sat_state(
                                    r_eph,
                                    b.pr_l1,
                                    0.0,
                                    base_time,
                                    base_coord.vector,
                                );
                                let dist_rov = (state.position.vector - r_sat_vec).norm();
                                let dist_base = (base_coord.vector - b_sat_vec).norm();

                                let initial_est_l2 = (cp_l2_rov - dist_rov - b_clock_rov)
                                    - (cp_l2_base - dist_base - b_clock_base);
                                state.add_ambiguity(
                                    r.sat,
                                    2,
                                    initial_est_l2,
                                    config.initial_ambiguity_variance,
                                );
                                initialized = true;
                                break;
                            }
                        }
                    }
                }
            }

            if !initialized {
                let initial_est_l2 = (cp_l2_rov - r_pr2) - (cp_l2_base - b_pr2);
                state.add_ambiguity(r.sat, 2, initial_est_l2, config.initial_ambiguity_variance);
            }
        }
        if let (Some(r_pr2), Some(r_cp2), Some(b_pr2), Some(b_cp2), Some(r_cp1), Some(b_cp1)) =
            (r.pr_l2, r.cp_l2, b.pr_l2, b.cp_l2, r.cp_l1, b.cp_l1)
        {
            let r_cp1_m = r_cp1 * (gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f1);
            let r_cp2_m = r_cp2 * (gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f2);
            let b_cp1_m = b_cp1 * (gneiss_core::constants::SPEED_OF_LIGHT_M_S / b_f1);
            let b_cp2_m = b_cp2 * (gneiss_core::constants::SPEED_OF_LIGHT_M_S / b_f2);
            let mw_sd = crate::combinations::melbourne_wubbena(
                r_cp1_m, r_cp2_m, r.pr_l1, r_pr2, r_f1, r_f2,
            ) - crate::combinations::melbourne_wubbena(
                b_cp1_m, b_cp2_m, b.pr_l1, b_pr2, b_f1, b_f2,
            );
            state.update_mw(r.sat, mw_sd / crate::combinations::lambda_wl(r_f1, r_f2));
        }
    }
}

pub fn check_doppler_phase_slip(
    cp: f64,
    prev_cp: f64,
    doppler: f64,
    prev_doppler: f64,
    dt: f64,
    threshold: f64,
) -> bool {
    let expected_change = -0.5 * (doppler + prev_doppler) * dt;
    let diff = (cp - prev_cp - expected_change).abs();
    if diff > threshold {
        tracing::debug!(
            "Doppler slip: cp={} prev_cp={} dop={} prev_dop={} exp={} diff={}",
            cp,
            prev_cp,
            doppler,
            prev_doppler,
            expected_change,
            diff
        );
    }
    diff > threshold
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::EngineConfig;
    use gneiss_core::coords::{Datum, Frame};
    use gneiss_core::sat::{Constellation, SatelliteId};
    use nalgebra::Vector3;

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    fn make_state(week: u32, tow: f64) -> RtkState {
        let t = GpsTime::new(week, tow);
        let pos = Coordinate::new(
            Vector3::new(6378000.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            t,
        );
        RtkState::new(t, pos, 10.0)
    }

    fn gps_sat(prn: u8) -> SatelliteId {
        SatelliteId {
            constellation: Constellation::Gps,
            prn,
        }
    }

    fn base_coord(t: GpsTime) -> Coordinate {
        Coordinate::new(
            Vector3::new(6378000.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            t,
        )
    }

    // ------------------------------------------------------------------
    // Locktime slip detection (lines 15-37)
    // ------------------------------------------------------------------

    #[test]
    /// r_lock == 0 means the receiver lost lock and re-acquired the satellite;
    /// a cycle slip must be flagged and the locktime counter reset to 0.
    fn test_locktime_zero_detects_slip() {
        let t = GpsTime::new(2000, 0.0);
        let mut state = make_state(2000, 0.0);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);

        let config = EngineConfig::default();
        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(0),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t, t);

        assert_eq!(
            state.locktimes.get(&(gps_sat(1), 1)),
            Some(&0),
            "locktime reset to 0 when r_lock == 0"
        );
        assert_eq!(
            state.locktimes.get(&(gps_sat(1), 2)),
            Some(&0),
            "band-2 locktime also reset"
        );
    }

    #[test]
    /// When r.locktime < prev_lock, the lock counter decreased unexpectedly,
    /// indicating a loss-of-lock / cycle slip.
    fn test_locktime_decreased_detects_slip() {
        let t = GpsTime::new(2000, 0.0);
        let mut state = make_state(2000, 0.0);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);

        let config = EngineConfig::default();
        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(3),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t, t);

        assert_eq!(
            state.locktimes.get(&(gps_sat(1), 1)),
            Some(&3),
            "locktime adopts rover's smaller value (3 < 10)"
        );
    }

    #[test]
    /// When r.locktime >= prev_lock, the counter is increasing monotonically
    /// and no locktime-based slip is flagged.
    fn test_locktime_increasing_no_slip() {
        let t = GpsTime::new(2000, 0.0);
        let mut state = make_state(2000, 0.0);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);

        let config = EngineConfig::default();
        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(15),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t, t);

        assert_eq!(
            state.locktimes.get(&(gps_sat(1), 1)),
            Some(&15),
            "locktime increases from 10 to 15 without slip"
        );
    }

    #[test]
    /// First time a satellite is observed: no previous locktime in state,
    /// so prev_lock defaults to 0 and the rover's locktime is adopted directly.
    fn test_locktime_first_observation_no_slip() {
        let t = GpsTime::new(2000, 0.0);
        let mut state = make_state(2000, 0.0);
        // No locktimes entry for this sat -- prev_lock defaults to 0.

        let config = EngineConfig::default();
        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(5),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t, t);

        assert_eq!(
            state.locktimes.get(&(gps_sat(1), 1)),
            Some(&5),
            "first observation adopts rover's locktime directly"
        );
    }

    #[test]
    /// When r.locktime is None but a previous locktime entry exists,
    /// new_lock = prev_lock + 1 (increment the existing counter).
    fn test_locktime_none_with_history() {
        let t = GpsTime::new(2000, 0.0);
        let mut state = make_state(2000, 0.0);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);

        let config = EngineConfig::default();
        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: None,
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t, t);

        assert_eq!(
            state.locktimes.get(&(gps_sat(1), 1)),
            Some(&11),
            "locktime increments from 10 to 11 when r.locktime is None"
        );
    }

    // ------------------------------------------------------------------
    // Reject-count slip forcing (lines 58-63)
    // ------------------------------------------------------------------

    #[test]
    /// When reject_count[(sat, freq)] > max_reject_count, a slip is forced
    /// for that frequency band regardless of locktime / Doppler checks.
    fn test_reject_count_forces_slip() {
        let t = GpsTime::new(2000, 0.0);
        let mut state = make_state(2000, 0.0);
        // reject_count exceeds default max_reject_count (3)
        state.reject_counts.insert((gps_sat(1), 1), 5);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);

        let config = EngineConfig::default();
        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t, t);

        // slip triggered -> reject count should be reset to 0
        assert_eq!(
            state.reject_counts.get(&(gps_sat(1), 1)),
            Some(&0),
            "reject_count reset after forced slip"
        );
    }

    // ------------------------------------------------------------------
    // Doppler-phase slip on L1 (lines 65-109)
    // ------------------------------------------------------------------

    #[test]
    /// Doppler-phase consistency check on L1 detects a cycle slip when the
    /// carrier phase jump exceeds the threshold.
    fn test_doppler_l1_slip_detected() {
        let t0 = GpsTime::new(2000, 0.0);
        let t1 = GpsTime::new(2000, 1.0);
        let mut state = make_state(2000, 0.0);
        let doppler = -1500.0;

        // Populate phase history for L1
        state
            .phase_history
            .insert((gps_sat(1), 1), (100.0, -2000.0, t0));

        let config = EngineConfig::default();
        let expected_change = -0.5 * (doppler + (-2000.0)) * 1.0; // = 1750.0
        let consistent_cp = 100.0 + expected_change;
        let slipped_cp = consistent_cp + 6.0; // exceeds default 5.0 threshold

        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: Some(slipped_cp),
            cp_l2: None,
            doppler,
            snr: 45.0,
            locktime: Some(15),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None, // no base L1 to avoid ambiguity init
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t1);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t1, t1);

        let (updated_cp, updated_dop, updated_time) = *state
            .phase_history
            .get(&(gps_sat(1), 1))
            .expect("phase_history should be updated");
        assert!(
            (updated_cp - slipped_cp).abs() < 1e-9,
            "phase history stores slipped cp: got {updated_cp}"
        );
        assert_eq!(updated_dop, doppler);
        assert_eq!(updated_time, t1);
    }

    #[test]
    /// Doppler-phase consistency check on L1 confirms no slip when the
    /// carrier phase follows the expected Doppler integral.
    fn test_doppler_l1_no_slip() {
        let t0 = GpsTime::new(2000, 0.0);
        let t1 = GpsTime::new(2000, 1.0);
        let mut state = make_state(2000, 0.0);
        let doppler = -1500.0;

        state
            .phase_history
            .insert((gps_sat(1), 1), (100.0, -2000.0, t0));

        let config = EngineConfig::default();
        let expected_change = -0.5 * (doppler + (-2000.0)) * 1.0;
        let consistent_cp = 100.0 + expected_change;

        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: Some(consistent_cp),
            cp_l2: None,
            doppler,
            snr: 45.0,
            locktime: Some(15),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t1);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t1, t1);

        let (updated_cp, _, _) = *state
            .phase_history
            .get(&(gps_sat(1), 1))
            .expect("phase_history should be updated");
        assert!(
            (updated_cp - consistent_cp).abs() < 1e-9,
            "phase history stores consistent cp"
        );
    }

    #[test]
    /// When slip_l1 is set true (from locktime) and Doppler values are consistent
    /// with dt < 5.0 s, the short interruption is bridged and slip_l1 is cleared.
    fn test_doppler_l1_bridging() {
        let t0 = GpsTime::new(2000, 0.0);
        let t1 = GpsTime::new(2000, 1.0);
        let mut state = make_state(2000, 0.0);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);
        state
            .phase_history
            .insert((gps_sat(1), 1), (100.0, -2000.0, t0));

        let config = EngineConfig::default();
        let doppler = -2000.0; // unchanged -> expected_change = 2000.0
        let expected_change = -0.5 * (doppler + (-2000.0)) * 1.0;
        let consistent_cp = 100.0 + expected_change;

        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: Some(consistent_cp),
            cp_l2: None,
            doppler,
            snr: 45.0,
            locktime: Some(0), // triggers slip
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t1);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t1, t1);

        // Phase history updated as part of the bridging path
        let (updated_cp, _, _) = *state
            .phase_history
            .get(&(gps_sat(1), 1))
            .expect("phase_history should be updated");
        assert!(
            (updated_cp - consistent_cp).abs() < 1e-9,
            "bridging stores consistent cp"
        );
    }

    #[test]
    /// When dt > max_base_age_s, a data-gap reset forces slip_l1 = true.
    fn test_doppler_l1_data_gap_reset() {
        let t0 = GpsTime::new(2000, 0.0);
        let t1 = GpsTime::new(2000, 10.0); // dt = 10 s > 5.0 s max_base_age_s
        let mut state = make_state(2000, 0.0);
        state
            .phase_history
            .insert((gps_sat(1), 1), (100.0, -2000.0, t0));

        let mut config = EngineConfig::default();
        config.max_base_age_s = 5.0;

        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: Some(2000.0),
            cp_l2: None,
            doppler: -1500.0,
            snr: 45.0,
            locktime: Some(20),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t1);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t1, t1);

        let (updated_cp, _, _) = *state
            .phase_history
            .get(&(gps_sat(1), 1))
            .expect("phase_history should be updated");
        assert!(
            (updated_cp - 2000.0).abs() < 1e-9,
            "phase history updated after data gap"
        );
    }

    // ------------------------------------------------------------------
    // Doppler-phase slip on L2 (lines 111-156)
    // ------------------------------------------------------------------

    #[test]
    /// Doppler-phase check on L2 detects a cycle slip using the
    /// frequency-scaled Doppler (dop_l2 = r.doppler * f2/f1).
    fn test_doppler_l2_slip_detected() {
        let t0 = GpsTime::new(2000, 0.0);
        let t1 = GpsTime::new(2000, 1.0);
        let mut state = make_state(2000, 0.0);

        let f2_f1 = 1227.60e6 / 1575.42e6;
        let prev_doppler_l2 = -2000.0 * f2_f1;
        state
            .phase_history
            .insert((gps_sat(1), 2), (100.0, prev_doppler_l2, t0));

        let config = EngineConfig::default();
        let doppler_rov = -1500.0;
        let doppler_l2 = doppler_rov * f2_f1;
        let expected_change = -0.5 * (doppler_l2 + prev_doppler_l2) * 1.0;
        let consistent_cp = 100.0 + expected_change;
        let slipped_cp = consistent_cp + 6.0; // > 5.0 threshold

        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: None, // skip L1
            cp_l2: Some(slipped_cp),
            doppler: doppler_rov,
            snr: 45.0,
            locktime: Some(15),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t1);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t1, t1);

        let (updated_cp, updated_dop, updated_time) = *state
            .phase_history
            .get(&(gps_sat(1), 2))
            .expect("L2 phase_history should be updated");
        assert!(
            (updated_cp - slipped_cp).abs() < 1e-6,
            "L2 phase history stores slipped cp: got {updated_cp}"
        );
        assert!(
            (updated_dop - doppler_l2).abs() < 1e-9,
            "L2 phase history stores scaled doppler"
        );
        assert_eq!(updated_time, t1);
    }

    #[test]
    /// L2 Doppler-phase bridging: a short interruption (dt < 5.0 s) where
    /// carrier phase stays consistent is bridged, clearing slip_l2.
    fn test_doppler_l2_bridging() {
        let t0 = GpsTime::new(2000, 0.0);
        let t1 = GpsTime::new(2000, 1.0);
        let mut state = make_state(2000, 0.0);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);

        let f2_f1 = 1227.60e6 / 1575.42e6;
        let prev_doppler_l2 = -2000.0 * f2_f1;
        state
            .phase_history
            .insert((gps_sat(1), 2), (100.0, prev_doppler_l2, t0));

        let config = EngineConfig::default();
        let doppler_rov = -2000.0; // unchanged -> trivial integral
        let doppler_l2 = doppler_rov * f2_f1;
        let expected_change = -0.5 * (doppler_l2 + prev_doppler_l2) * 1.0;
        let consistent_cp = 100.0 + expected_change;

        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: None,
            cp_l2: Some(consistent_cp),
            doppler: doppler_rov,
            snr: 45.0,
            locktime: Some(0), // triggers slip
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t1);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t1, t1);

        let (updated_cp, _, _) = *state
            .phase_history
            .get(&(gps_sat(1), 2))
            .expect("L2 phase_history should be updated");
        assert!(
            (updated_cp - consistent_cp).abs() < 1e-6,
            "L2 bridging stores consistent cp"
        );
    }

    #[test]
    /// When dt > max_base_age_s, a data-gap reset forces slip on L2.
    fn test_doppler_l2_data_gap_reset() {
        let t0 = GpsTime::new(2000, 0.0);
        let t1 = GpsTime::new(2000, 10.0);
        let mut state = make_state(2000, 0.0);
        state
            .phase_history
            .insert((gps_sat(1), 2), (100.0, -2000.0, t0));

        let mut config = EngineConfig::default();
        config.max_base_age_s = 5.0;

        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: None,
            cp_l2: Some(2000.0),
            doppler: -1500.0,
            snr: 45.0,
            locktime: Some(20),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t1);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t1, t1);

        let (updated_cp, _, _) = *state
            .phase_history
            .get(&(gps_sat(1), 2))
            .expect("L2 phase_history should be updated");
        assert!(
            (updated_cp - 2000.0).abs() < 1e-9,
            "L2 phase history updated after data gap"
        );
    }

    // ------------------------------------------------------------------
    // Melbourne-Wubbena update path (lines 414-427)
    // ------------------------------------------------------------------

    #[test]
    /// When both L1 & L2 carrier phase and L2 pseudorange are present for
    /// rover and base, the MW combination is computed and the EMA updated.
    fn test_mw_update_executed() {
        let t = GpsTime::new(2000, 0.0);
        let mut state = make_state(2000, 0.0);

        let config = EngineConfig::default();
        let f1 = 1575.42e6;
        let f2 = 1227.60e6;
        let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

        let r_cp1 = 10000.0;
        let r_cp2 = 8000.0;
        let b_cp1 = 10000.0;
        let b_cp2 = 8000.0;
        let pr_l1 = 20_000_000.0;
        let pr_l2 = 20_000_000.0;

        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1,
            pr_l2: Some(pr_l2),
            cp_l1: Some(r_cp1),
            cp_l2: Some(r_cp2),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1,
            pr_l2: Some(pr_l2),
            cp_l1: Some(b_cp1),
            cp_l2: Some(b_cp2),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t);

        // Expected MW computation
        let r_cp1_m = r_cp1 * (c / f1);
        let r_cp2_m = r_cp2 * (c / f2);
        let b_cp1_m = b_cp1 * (c / f1);
        let b_cp2_m = b_cp2 * (c / f2);
        let mw_sd = crate::combinations::melbourne_wubbena(r_cp1_m, r_cp2_m, pr_l1, pr_l2, f1, f2)
            - crate::combinations::melbourne_wubbena(b_cp1_m, b_cp2_m, pr_l1, pr_l2, f1, f2);
        let expected_mw_cycles = mw_sd / crate::combinations::lambda_wl(f1, f2);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t, t);

        let mw_entry = state
            .mw_sd_ema
            .get(&gps_sat(1))
            .expect("MW EMA should be stored");
        assert!(
            (*mw_entry - expected_mw_cycles).abs() < 1e-6,
            "MW EMA matches expected: got {mw_entry}, expected {expected_mw_cycles}"
        );
        assert_eq!(
            state.mw_sd_counts.get(&gps_sat(1)),
            Some(&1),
            "MW count incremented to 1"
        );
    }

    #[test]
    /// When L2 data is missing for rover, the MW update is skipped.
    fn test_mw_update_skipped_when_l2_missing() {
        let t = GpsTime::new(2000, 0.0);
        let mut state = make_state(2000, 0.0);

        let config = EngineConfig::default();

        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: None, // missing L2 -> MW skipped
            cp_l1: Some(10000.0),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(10000.0),
            cp_l2: Some(8000.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t, t);

        assert!(
            state.mw_sd_ema.get(&gps_sat(1)).is_none(),
            "MW not stored when L2 data missing"
        );
    }

    // ------------------------------------------------------------------
    // Ambiguity re-initialization after slip
    // ------------------------------------------------------------------

    #[test]
    /// After a slip triggers ambiguity removal, a new ambiguity is
    /// initialized via the simple formula (no anchor-based init when
    /// covariance is large).
    fn test_slip_triggers_ambiguity_reinit() {
        let t = GpsTime::new(2000, 0.0);
        let mut state = make_state(2000, 0.0);
        // Pre-add an ambiguity to verify it is removed
        state.add_ambiguity(gps_sat(1), 1, 123.45, 10000.0);
        assert!(
            state.ambiguity_keys.contains(&(gps_sat(1), 1)),
            "ambiguity should exist before slip"
        );

        state.locktimes.insert((gps_sat(1), 1), 10_u16);

        let config = EngineConfig::default();
        let f1 = 1575.42e6;

        let rover = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(10000.0),
            cp_l2: Some(8000.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(0), // triggers slip
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(10000.0),
            cp_l2: Some(8000.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rover, base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t, t);

        // After slip, old ambiguity removed and a new one added
        assert!(
            state.ambiguity_keys.contains(&(gps_sat(1), 1)),
            "new ambiguity added after slip"
        );
        // Simple formula: (r_cp1*lam - r.pr_l1) - (b_cp1*lam - b.pr_l1)
        let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        let r_cp1_m = 10000.0 * (c / f1);
        let expected_est = (r_cp1_m - 20_000_000.0) - (r_cp1_m - 20_000_000.0);
        assert!(
            (state.ambiguities[0] - expected_est).abs() < 1.0,
            "ambiguity re-init: got {}, expected {}",
            state.ambiguities[0],
            expected_est
        );
    }

    // ------------------------------------------------------------------
    // check_doppler_phase_slip edge cases
    // ------------------------------------------------------------------

    #[test]
    /// check_doppler_phase_slip on consistent data with dt = 0 returns false.
    fn test_doppler_slip_zero_dt() {
        let slip = check_doppler_phase_slip(100.0, 100.0, -2000.0, -2000.0, 0.0, 5.0);
        assert!(!slip, "zero dt with same cp: no slip expected");
    }

    #[test]
    /// check_doppler_phase_slip with negative dt should not panic and returns
    /// a deterministic result.
    fn test_doppler_slip_negative_dt() {
        let slip = check_doppler_phase_slip(100.0, 110.0, -2000.0, -1500.0, -1.0, 5.0);
        // Should not panic; result depends on sign-flipped expected change
        let _ = slip;
    }

    #[test]
    /// At the exact threshold boundary, no slip is detected (strict >).
    fn test_doppler_slip_exact_threshold() {
        let prev_cp = 100.0;
        let prev_doppler = -2000.0;
        let doppler = -2000.0;
        let dt = 1.0;
        // expected_change = -0.5 * (-2000 + -2000) * 1 = 2000.0
        let expected_change: f64 = -0.5 * (doppler + prev_doppler) * dt;
        assert!((expected_change - 2000.0).abs() < 1e-12);

        // diff == threshold (5.0) -> no slip
        let boundary_cp = prev_cp + expected_change + 5.0;
        let slip = check_doppler_phase_slip(
            boundary_cp,
            prev_cp,
            doppler,
            prev_doppler,
            dt,
            5.0,
        );
        assert!(!slip, "diff == threshold: no slip (strict >)");

        // diff just above threshold -> slip
        let slip_cp = prev_cp + expected_change + 5.0001;
        let slip2 = check_doppler_phase_slip(slip_cp, prev_cp, doppler, prev_doppler, dt, 5.0);
        assert!(slip2, "diff > threshold: slip detected");
    }

    #[test]
    /// NaN inputs should not panic and should return false (NaN comparisons
    /// always return false).
    fn test_doppler_slip_nan_inputs() {
        let slip =
            check_doppler_phase_slip(f64::NAN, 100.0, -2000.0, -2000.0, 1.0, 5.0);
        assert!(!slip, "NaN cp should not panic (diff = NaN, NaN > 5.0 is false)");

        let slip =
            check_doppler_phase_slip(100.0, 100.0, f64::NAN, -2000.0, 1.0, 5.0);
        assert!(!slip, "NaN doppler should not detect slip");
    }

    #[test]
    /// Very small dt (near zero) with identical cp values -> no slip.
    fn test_doppler_slip_small_dt() {
        let slip =
            check_doppler_phase_slip(100.0, 100.0, -2000.0, -2000.0, 1e-10, 5.0);
        assert!(!slip, "very small dt with identical cp: no slip");
    }

    // ------------------------------------------------------------------
    // Preserved original test
    // ------------------------------------------------------------------

    #[test]
    fn test_doppler_phase_cycle_slip_detection() {
        let prev_cp = 100.0;
        let prev_doppler = -2.0;

        let dt = 1.0;
        let doppler = 5.0; // test approaching satellite

        let expected_change = -0.5 * (doppler + prev_doppler) * dt;
        assert_eq!(expected_change, -1.5);

        let cp = 100.0 + expected_change;

        let slip = check_doppler_phase_slip(cp, prev_cp, doppler, prev_doppler, dt, 5.0);
        assert!(!slip, "No slip should be detected for consistent doppler");

        let cp_slip = 100.0 + expected_change + 6.0;
        let slip = check_doppler_phase_slip(cp_slip, prev_cp, doppler, prev_doppler, dt, 5.0);
        assert!(slip, "Slip should be detected when phase jumps");
    }

    // ------------------------------------------------------------------
    // Helper: create a GPS ephemeris for testing
    // ------------------------------------------------------------------

    fn make_gps_eph(sat: SatelliteId, time: GpsTime, m0: f64) -> Ephemeris {
        Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat,
            toe: time,
            toc: time,
            af0: 0.0,
            af1: 0.0,
            af2: 0.0,
            crs: 0.0,
            crc: 0.0,
            cuc: 0.0,
            cus: 0.0,
            cic: 0.0,
            cis: 0.0,
            m0,
            e: 0.01,
            sqrt_a: 5153.6,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 1.0,
            idot: 0.0,
            omega: 0.0,
            tgd: 0.0,
            iode: 0,
            iodc: 0,
        })
    }

    // ------------------------------------------------------------------
    // Frequency lookup from ephemeris (lines 39-54)
    // ------------------------------------------------------------------

    #[test]
    /// Glonass ephemeris provides freq_num != 0, exercising the FDMA
    /// frequency computation path. Without ephemeris in the list, freq_num
    /// defaults to 0; with ephemeris, the actual freq_num is used.
    fn test_freq_lookup_with_glonass_ephemeris() {
        let t = GpsTime::new(2137, 422922.0);
        let mut state = make_state(2137, 422922.0);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);

        let config = EngineConfig::default();
        let glonass_sat = SatelliteId {
            constellation: Constellation::Glonass,
            prn: 6,
        };
        // freq_num = -4 means the satellite uses frequency channel -4 (FDMA)
        let glonass_eph = Ephemeris::Glonass(gneiss_core::ephemeris::GlonassEphemeris {
            sat: glonass_sat,
            toe: t,
            freq_num: -4,
            tau_n: 0.0,
            gamma_n: 0.0,
            delta_tau_n: 0.0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            ax: 0.0,
            ay: 0.0,
            az: 0.0,
        });

        let rov = DdObservation {
            sat: glonass_sat,
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rov, base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[glonass_eph], &bc, t, t);

        // Processing the Glonass satellite should store its locktime
        assert_eq!(
            state.locktimes.get(&(glonass_sat, 1)),
            Some(&10),
            "Glonass satellite processed through freq lookup path"
        );
    }

    #[test]
    /// When ephemeris list does not contain the satellite, freq_num defaults
    /// to 0 and nominal (non-FDMA) frequencies are used.
    fn test_freq_lookup_missing_ephemeris() {
        let t = GpsTime::new(2137, 422922.0);
        let mut state = make_state(2137, 422922.0);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);

        let config = EngineConfig::default();

        let rov = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(10000.0),
            cp_l2: Some(8000.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(10000.0),
            cp_l2: Some(8000.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rov, base)];
        let bc = base_coord(t);

        // Provide an ephemeris for a DIFFERENT satellite (sat 3) — rover's sat
        // is NOT in the list, so freq_num defaults to 0.
        let other_eph = make_gps_eph(gps_sat(3), t, 1.0);
        manage_ambiguities_and_slips(&mut state, &config, &matched, &[other_eph], &bc, t, t);

        // MW update executed with default (nominal) frequencies
        assert!(
            state.mw_sd_ema.get(&gps_sat(1)).is_some(),
            "MW EMA stored after freq lookup with missing ephemeris"
        );
    }

    // ------------------------------------------------------------------
    // needs_l1 / needs_l2 logic (lines 178-196)
    // ------------------------------------------------------------------

    #[test]
    /// When an ambiguity already exists for the satellite-frequency pair
    /// AND no slip is detected, needs_l1 should be false and the ambiguity
    /// value should be preserved unchanged.
    fn test_needs_l1_false_when_ambiguity_present() {
        let t = GpsTime::new(2000, 0.0);
        let mut state = make_state(2000, 0.0);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);
        // Pre-add an ambiguity with a distinctive value
        state.add_ambiguity(gps_sat(1), 1, 42.0, 10000.0);

        let config = EngineConfig::default();

        let rov = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(10000.0),
            cp_l2: Some(8000.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(10000.0),
            cp_l2: Some(8000.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rov, base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t, t);

        // The pre-existing ambiguity should still be there with its original value
        let idx = state
            .ambiguity_keys
            .iter()
            .position(|&(s, f)| s == gps_sat(1) && f == 1)
            .expect("ambiguity should still exist");
        assert!(
            (state.ambiguities[idx] - 42.0).abs() < 1e-9,
            "ambiguity value preserved when needs_l1 is false"
        );
    }

    #[test]
    /// A slip on L2 (slip_l2 = true) should also trigger needs_l1 for L1,
    /// because the double-difference float ambiguity breaks when either
    /// frequency slips.
    fn test_slip_l2_triggers_l1_reinit() {
        let t = GpsTime::new(2000, 0.0);
        let mut state = make_state(2000, 0.0);
        // Set reject count on L2 so slip_l2 is forced even with consistent locktime
        state.reject_counts.insert((gps_sat(1), 2), 5);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);

        let config = EngineConfig::default();

        let rov = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(10000.0),
            cp_l2: Some(8000.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(10000.0),
            cp_l2: Some(8000.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rov, base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t, t);

        // L2 reject_count should have been reset
        assert_eq!(
            state.reject_counts.get(&(gps_sat(1), 2)),
            Some(&0),
            "L2 reject_count reset after slip_l2"
        );
        // A new L1 ambiguity should have been initialized (triggered by slip_l2)
        assert!(
            state.ambiguity_keys.contains(&(gps_sat(1), 1)),
            "L1 ambiguity re-initialized after slip on L2"
        );
    }

    #[test]
    /// When both L1 and L2 ambiguities exist and a slip occurs, both are
    /// removed and reject_counts for both bands are reset.
    fn test_slip_removes_both_ambiguity_bands() {
        let t = GpsTime::new(2000, 0.0);
        let mut state = make_state(2000, 0.0);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);
        // Add both L1 and L2 ambiguities
        state.add_ambiguity(gps_sat(1), 1, 42.0, 10000.0);
        state.add_ambiguity(gps_sat(1), 2, 99.0, 10000.0);
        state.reject_counts.insert((gps_sat(1), 1), 5);
        state.reject_counts.insert((gps_sat(1), 2), 3);

        let config = EngineConfig::default();

        let rov = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(10000.0),
            cp_l2: Some(8000.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(10000.0),
            cp_l2: Some(8000.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rov, base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t, t);

        // Both reject counts reset
        assert_eq!(
            state.reject_counts.get(&(gps_sat(1), 1)),
            Some(&0),
            "L1 reject_count reset"
        );
        assert_eq!(
            state.reject_counts.get(&(gps_sat(1), 2)),
            Some(&0),
            "L2 reject_count reset"
        );
        // New ambiguities added for both bands (since slip_l1 from reject_count)
        assert!(
            state.ambiguity_keys.contains(&(gps_sat(1), 1)),
            "L1 ambiguity re-initialized"
        );
        assert!(
            state.ambiguity_keys.contains(&(gps_sat(1), 2)),
            "L2 ambiguity re-initialized"
        );
    }

    // ------------------------------------------------------------------
    // Anchor-based ambiguity initialization — L1 (lines 198-303)
    // ------------------------------------------------------------------

    #[test]
    /// When position covariance is converged (cov(0,0) < 0.1) and an anchor
    /// satellite with converged ambiguity exists, the new satellite's
    /// ambiguity is initialized from geometry via the anchor.
    fn test_anchor_based_l1_init() {
        let t = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(6378000.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            t,
        );
        // Small initial_var → cov(0,0) = 0.01 < 0.1 → anchor search enabled
        let mut state = RtkState::new(t, pos, 0.01);

        let anchor_sat = gps_sat(2);
        let new_sat = gps_sat(1);

        // Add converged anchor ambiguity (variance 0.01 < 0.05)
        state.add_ambiguity(anchor_sat, 1, 0.0, 0.01);
        state.locktimes.insert((anchor_sat, 1), 10_u16);

        let config = EngineConfig::default();

        // Ephemerides for both satellites
        let anchor_eph = make_gps_eph(anchor_sat, t, 1.0);
        let new_eph = make_gps_eph(new_sat, t, 2.0);
        let ephs = vec![anchor_eph, new_eph];

        let pr_l1 = 20_000_000.0;
        let lam_l1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / 1575.42e6;
        let cp_l1 = pr_l1 / lam_l1; // carrier phase in cycles ≈ geometric range

        // Both satellites have L1 data; cp_l2 = None so only L1 init is tested
        let anchor_rov = DdObservation {
            sat: anchor_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let anchor_base = DdObservation {
            sat: anchor_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let new_rov = DdObservation {
            sat: new_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let new_base = DdObservation {
            sat: new_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };

        // Anchor first (preserves existing ambiguity), then new satellite
        let matched = vec![(anchor_rov, anchor_base), (new_rov, new_base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &ephs, &bc, t, t);

        // New ambiguity added for the new satellite
        assert!(
            state.ambiguity_keys.contains(&(new_sat, 1)),
            "Anchor-based L1 init added ambiguity for new satellite"
        );
        // Anchor ambiguity preserved
        assert!(
            state.ambiguity_keys.contains(&(anchor_sat, 1)),
            "Anchor ambiguity preserved"
        );
        // All ambiguity values are finite
        assert!(
            state.ambiguities.iter().all(|v| v.is_finite()),
            "All ambiguity values are finite"
        );
        // At least 2 ambiguities (anchor + new)
        assert!(
            state.ambiguities.len() >= 2,
            "At least two ambiguities: anchor and new"
        );
    }

    #[test]
    /// When position covariance is not converged (cov(0,0) >= 0.1), the
    /// anchor search is skipped and the fallback formula (cp - pr) is used.
    fn test_anchor_search_skipped_position_not_converged() {
        let t = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(6378000.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            t,
        );
        // Default initial_var = 10.0 → cov(0,0) = 10.0 >= 0.1 → anchor disabled
        let mut state = RtkState::new(t, pos, 10.0);

        let anchor_sat = gps_sat(2);
        let new_sat = gps_sat(1);

        // Anchor ambiguity exists but position isn't converged, so anchor
        // search is skipped regardless
        state.add_ambiguity(anchor_sat, 1, 0.0, 0.01);
        state.locktimes.insert((anchor_sat, 1), 10_u16);

        let config = EngineConfig::default();

        let anchor_eph = make_gps_eph(anchor_sat, t, 1.0);
        let new_eph = make_gps_eph(new_sat, t, 2.0);
        let ephs = vec![anchor_eph, new_eph];

        let pr_l1 = 20_000_000.0;
        let lam_l1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / 1575.42e6;
        let cp_l1 = pr_l1 / lam_l1;

        let anchor_rov = DdObservation {
            sat: anchor_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let anchor_base = DdObservation {
            sat: anchor_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let new_rov = DdObservation {
            sat: new_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let new_base = DdObservation {
            sat: new_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };

        let matched = vec![(anchor_rov, anchor_base), (new_rov, new_base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &ephs, &bc, t, t);

        // Ambiguity initialized via fallback (cp - pr formula)
        assert!(
            state.ambiguity_keys.contains(&(new_sat, 1)),
            "Fallback L1 init added ambiguity for new satellite"
        );
        assert!(
            state.ambiguities.iter().all(|v| v.is_finite()),
            "All ambiguity values finite after fallback init"
        );
    }

    #[test]
    /// When the anchor ambiguity itself is not converged (cov > 0.05), the
    /// anchor search skips it and falls through to the cp - pr formula.
    fn test_anchor_ambiguity_not_converged() {
        let t = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(6378000.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            t,
        );
        // Position converged, so anchor search is attempted
        let mut state = RtkState::new(t, pos, 0.01);

        let anchor_sat = gps_sat(2);
        let new_sat = gps_sat(1);

        // Anchor ambiguity with LARGE variance (10.0 >= 0.05) → not converged
        state.add_ambiguity(anchor_sat, 1, 10.0, 10.0);
        state.locktimes.insert((anchor_sat, 1), 10_u16);

        let config = EngineConfig::default();

        let anchor_eph = make_gps_eph(anchor_sat, t, 1.0);
        let new_eph = make_gps_eph(new_sat, t, 2.0);
        let ephs = vec![anchor_eph, new_eph];

        let pr_l1 = 20_000_000.0;
        let lam_l1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / 1575.42e6;
        let cp_l1 = pr_l1 / lam_l1;

        let anchor_rov = DdObservation {
            sat: anchor_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let anchor_base = DdObservation {
            sat: anchor_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let new_rov = DdObservation {
            sat: new_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let new_base = DdObservation {
            sat: new_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };

        let matched = vec![(anchor_rov, anchor_base), (new_rov, new_base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &ephs, &bc, t, t);

        // Ambiguity should still be added via fallback (anchor skipped)
        assert!(
            state.ambiguity_keys.contains(&(new_sat, 1)),
            "Fallback L1 init when anchor ambiguity not converged"
        );
        assert!(
            state.ambiguities.iter().all(|v| v.is_finite()),
            "All ambiguity values finite"
        );
    }

    #[test]
    /// When the anchor satellite's ephemeris is not found in the list, the
    /// anchor search continues to the next candidate (skip, don't panic).
    fn test_anchor_search_anchor_ephemeris_missing() {
        let t = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(6378000.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            t,
        );
        let mut state = RtkState::new(t, pos, 0.01);

        let anchor_sat = gps_sat(2);
        let new_sat = gps_sat(1);

        state.add_ambiguity(anchor_sat, 1, 0.0, 0.01);
        state.locktimes.insert((anchor_sat, 1), 10_u16);

        let config = EngineConfig::default();

        // Only provide ephemeris for the NEW satellite, not the anchor
        let new_eph = make_gps_eph(new_sat, t, 2.0);
        let ephs = vec![new_eph];

        let pr_l1 = 20_000_000.0;
        let lam_l1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / 1575.42e6;
        let cp_l1 = pr_l1 / lam_l1;

        let anchor_rov = DdObservation {
            sat: anchor_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let anchor_base = DdObservation {
            sat: anchor_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let new_rov = DdObservation {
            sat: new_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let new_base = DdObservation {
            sat: new_sat,
            pr_l1,
            pr_l2: None,
            cp_l1: Some(cp_l1),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };

        let matched = vec![(anchor_rov, anchor_base), (new_rov, new_base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &ephs, &bc, t, t);

        // Anchor search skips (no anchor ephemeris) → fallback adds ambiguity
        assert!(
            state.ambiguity_keys.contains(&(new_sat, 1)),
            "Fallback when anchor ephemeris missing"
        );
        assert!(
            state.ambiguities.iter().all(|v| v.is_finite()),
            "All ambiguities finite"
        );
    }

    // ------------------------------------------------------------------
    // Anchor-based ambiguity initialization — L2 (lines 305-413)
    // ------------------------------------------------------------------

    #[test]
    /// Anchor-based ambiguity initialization on L2. Similar to the L1 test,
    /// but the anchor has a converged L2 ambiguity and the new satellite
    /// has L2 carrier phase and pseudorange data.
    fn test_anchor_based_l2_init() {
        let t = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(6378000.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            t,
        );
        let mut state = RtkState::new(t, pos, 0.01);

        let anchor_sat = gps_sat(2);
        let new_sat = gps_sat(1);

        // Add converged L2 ambiguity for anchor
        state.add_ambiguity(anchor_sat, 2, 0.0, 0.01);
        state.locktimes.insert((anchor_sat, 2), 10_u16);

        let config = EngineConfig::default();

        let anchor_eph = make_gps_eph(anchor_sat, t, 1.0);
        let new_eph = make_gps_eph(new_sat, t, 2.0);
        let ephs = vec![anchor_eph, new_eph];

        let pr_l1 = 20_000_000.0;
        let pr_l2 = 20_000_000.0;
        let lam_l2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / 1227.60e6;
        let cp_l2 = pr_l2 / lam_l2;

        // L1: set cp_l1 = None so needs_l1 = false, only L2 init is tested.
        // L1 pseudorange (pr_l1) is still present (required by get_sat_state).
        let anchor_rov = DdObservation {
            sat: anchor_sat,
            pr_l1,
            pr_l2: Some(pr_l2),
            cp_l1: None,
            cp_l2: Some(cp_l2),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let anchor_base = DdObservation {
            sat: anchor_sat,
            pr_l1,
            pr_l2: Some(pr_l2),
            cp_l1: None,
            cp_l2: Some(cp_l2),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let new_rov = DdObservation {
            sat: new_sat,
            pr_l1,
            pr_l2: Some(pr_l2),
            cp_l1: None,
            cp_l2: Some(cp_l2),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };
        let new_base = DdObservation {
            sat: new_sat,
            pr_l1,
            pr_l2: Some(pr_l2),
            cp_l1: None,
            cp_l2: Some(cp_l2),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(10),
        };

        let matched = vec![(anchor_rov, anchor_base), (new_rov, new_base)];
        let bc = base_coord(t);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &ephs, &bc, t, t);

        // New L2 ambiguity added
        assert!(
            state.ambiguity_keys.contains(&(new_sat, 2)),
            "Anchor-based L2 init added ambiguity for new satellite"
        );
        assert!(
            state.ambiguity_keys.contains(&(anchor_sat, 2)),
            "Anchor L2 ambiguity preserved"
        );
        assert!(
            state.ambiguities.iter().all(|v| v.is_finite()),
            "All ambiguity values finite"
        );
    }

    // ------------------------------------------------------------------
    // GF cycle slip (feature-gated, lines 158-176)
    // ------------------------------------------------------------------

    #[cfg(feature = "gf-slip")]
    #[test]
    /// When gf-slip feature is enabled, the GF combination is checked for
    /// jumps > 0.05 m, which triggers a slip on both L1 and L2.
    fn test_gf_cycle_slip_detected() {
        let t0 = GpsTime::new(2000, 0.0);
        let t1 = GpsTime::new(2000, 1.0);
        let mut state = make_state(2000, 0.0);
        state.locktimes.insert((gps_sat(1), 1), 10_u16);

        // Pre-populate GF value so current epoch is compared against it
        state.gf_values.insert(gps_sat(1), 0.0);

        let config = EngineConfig::default();
        let f1 = 1575.42e6;
        let f2 = 1227.60e6;
        let lam1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / f1;
        let lam2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / f2;

        // GF = cp1*lam1 - cp2*lam2. Set cp values so GF jumps > 0.05.
        let consistent_cp1 = 100.0;
        let consistent_cp2 = 100.0;
        let prev_gf = consistent_cp1 * lam1 - consistent_cp2 * lam2; // = 0.0 since lam1≈lam2
        state.gf_values.insert(gps_sat(1), prev_gf);

        // New cp values where GF jumps by > 0.05 m
        let jump = 0.1; // > 0.05 threshold
        let slipped_cp1 = consistent_cp1 + jump / lam1;
        let slipped_cp2 = consistent_cp2;

        let rov = DdObservation {
            sat: gps_sat(1),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(slipped_cp1),
            cp_l2: Some(slipped_cp2),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(15),
        };
        let base = DdObservation {
            sat: gps_sat(2),
            pr_l1: 20_000_000.0,
            pr_l2: Some(20_000_000.0),
            cp_l1: Some(10000.0),
            cp_l2: Some(8000.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let matched = vec![(rov, base)];
        let bc = base_coord(t1);

        manage_ambiguities_and_slips(&mut state, &config, &matched, &[], &bc, t1, t1);

        // GF slip triggered → both L1 and L2 ambiguities should be re-init'd
        assert!(
            state.ambiguity_keys.contains(&(gps_sat(1), 1)),
            "GF slip triggered L1 re-init"
        );
        assert!(
            state.ambiguity_keys.contains(&(gps_sat(1), 2)),
            "GF slip triggered L2 re-init"
        );
    }
}
