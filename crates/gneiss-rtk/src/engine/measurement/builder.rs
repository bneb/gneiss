use crate::engine::measurement::types::{EkfUpdates, MeasurementEnvironment};
use crate::engine::measurement::compute_innovations;
use crate::engine::measurement_math;
use crate::filter::{DdObservation, RtkState};
use nalgebra::{DMatrix, DVector};

pub struct EkfMeasurementMatrices {
    pub z: DVector<f64>,
    pub h: DMatrix<f64>,
    pub r: DMatrix<f64>,
    pub mt: Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
}

fn group_measurements_by_constellation(
    matched_obs: &[(DdObservation, DdObservation)],
) -> std::collections::HashMap<gneiss_core::sat::Constellation, Vec<(DdObservation, DdObservation)>>
{
    let mut const_groups = std::collections::HashMap::<
        gneiss_core::sat::Constellation,
        Vec<(DdObservation, DdObservation)>,
    >::new();
    for obs in matched_obs {
        const_groups
            .entry(obs.0.sat.constellation)
            .or_default()
            .push(obs.clone());
    }
    const_groups
}

fn select_reference_satellite(
    group: &[(DdObservation, DdObservation)],
    state: &RtkState,
    env: &MeasurementEnvironment,
) -> usize {
    let mut best_score = -1.0;
    let mut ref_idx = 0;
    let rov_llh = gneiss_core::coords::ecef_to_llh(state.position.vector);

    for (i, (r, _)) in group.iter().enumerate() {
        if let Some(eph) = env
            .ephemerides
            .iter()
            .filter(|e| e.sat() == r.sat)
            .min_by(|a, b| {
                let da = (a.toe().tow - state.time.tow).abs();
                let db = (b.toe().tow - state.time.tow).abs();
                da.partial_cmp(&db).unwrap()
            })
        {
            let tau = r.pr_l1 / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
            let t_tx = gneiss_core::time::GpsTime::new(state.time.week, state.time.tow - tau);
            let (sat_pos, _, _, _): (nalgebra::Vector3<f64>, _, _, _) = eph.position(t_tx);
            let (_, el) = gneiss_core::coords::az_el(rov_llh, state.position.vector, sat_pos);

            let score = if r.cp_l1.is_some() { el + 100.0 } else { el };
            if score > best_score {
                best_score = score;
                ref_idx = i;
            }
        }
    }
    ref_idx
}

fn update_reject_counts(
    state: &mut RtkState,
    h_row: &DMatrix<f64>,
    state_size: usize,
    passed: bool,
) {
    for c in crate::filter::CORE_STATE_SIZE..state_size {
        if h_row[(0, c)] > 0.5 {
            let key = state.ambiguity_keys[c - crate::filter::CORE_STATE_SIZE];
            if passed {
                state.reject_counts.insert(key, 0);
            } else {
                let count = *state.reject_counts.get(&key).unwrap_or(&0) + 1;
                state.reject_counts.insert(key, count);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn filter_innovations_chi_squared(
    state: &mut RtkState,
    state_size: usize,
    chi_pr: f64,
    chi_cp: f64,
    z_all: &[f64],
    h_all: &[Vec<f64>],
    r_all: &[f64],
    type_all: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> Vec<usize> {
    let mut safe_indices = Vec::new();
    for i in 0..z_all.len() {
        let mut h_row = DMatrix::zeros(1, state_size);
        for c in 0..state_size {
            h_row[(0, c)] = h_all[i][c];
        }
        let s_ii = (&h_row * &state.covariance * h_row.transpose())[(0, 0)] + r_all[i];
        let chi2 = z_all[i] * z_all[i] / s_ii;

        let threshold = match type_all[i].1 {
            0 => chi_pr * chi_pr,
            1 | 2 => chi_cp * chi_cp,
            3 => chi_pr * 1000.0,
            _ => chi_pr * chi_pr,
        };

        let passed = chi2 <= threshold;
        if passed {
            safe_indices.push(i);
        } else {
            tracing::debug!(
                "Rejected meas type {} with inn: {:.3}, chi2: {:.1}",
                type_all[i].1,
                z_all[i],
                chi2
            );
        }

        if type_all[i].1 == 1 || type_all[i].1 == 2 {
            update_reject_counts(state, &h_row, state_size, passed);
        }
    }
    safe_indices
}

fn build_final_measurement_matrices(
    state_size: usize,
    safe_indices: Vec<usize>,
    z_all: &[f64],
    h_all: &[Vec<f64>],
    r_all: &[f64],
    type_all: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> Option<EkfMeasurementMatrices> {
    if safe_indices.len() >= 4 {
        let mut z_vec = DVector::zeros(safe_indices.len());
        let mut h_mat = DMatrix::zeros(safe_indices.len(), state_size);
        let mut r_diagonals = Vec::new();
        let mut t_vec = Vec::new();

        for (new_i, &old_i) in safe_indices.iter().enumerate() {
            z_vec[new_i] = z_all[old_i];
            for c in 0..state_size {
                h_mat[(new_i, c)] = h_all[old_i][c];
            }
            r_diagonals.push(r_all[old_i]);
            t_vec.push(type_all[old_i]);
        }
        let r_mat = measurement_math::build_dense_covariance_matrix(&r_diagonals, &t_vec);
        Some(EkfMeasurementMatrices {
            z: z_vec,
            h: h_mat,
            r: r_mat,
            mt: t_vec,
        })
    } else {
        tracing::warn!(
            "measurement model empty! all_z={}, safe_indices={}",
            z_all.len(),
            safe_indices.len()
        );
        None
    }
}

pub fn build_measurement_model(
    state: &mut RtkState,
    matched_obs: &[(DdObservation, DdObservation)],
    env: &MeasurementEnvironment,
    chi_square_pr_threshold: f64,
    chi_square_cp_threshold: f64,
) -> Option<EkfMeasurementMatrices> {
    let mut all = EkfUpdates::new();
    let const_groups = group_measurements_by_constellation(matched_obs);

    for (_, group) in const_groups {
        if group.len() < 2 {
            continue;
        }
        let ref_idx = select_reference_satellite(&group, state, env);
        let mut group_clone = group.clone();
        let (ref_rover, ref_base) = group_clone.remove(ref_idx);

        if let Some(updates) = compute_innovations(state, &group_clone, &ref_rover, &ref_base, env)
        {
            all.extend(updates);
        }
    }

    let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();
    tracing::trace!("Pre-filter z_all len: {}", all.z.len());

    let safe_indices = filter_innovations_chi_squared(
        state,
        state_size,
        chi_square_pr_threshold,
        chi_square_cp_threshold,
        &all.z,
        &all.h,
        &all.r,
        &all.mt,
    );

    build_final_measurement_matrices(state_size, safe_indices, &all.z, &all.h, &all.r, &all.mt)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::config::EkfTuningConfig;
    use crate::filter::CORE_STATE_SIZE;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::ephemeris::{Ephemeris, GpsEphemeris};
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;
    use nalgebra::Vector3;
    use std::collections::HashMap;

    // ------------------------------------------------------------------
    // Helper factories
    // ------------------------------------------------------------------

    fn gps(prn: u8) -> SatelliteId {
        SatelliteId { constellation: Constellation::Gps, prn }
    }

    fn gal(prn: u8) -> SatelliteId {
        SatelliteId { constellation: Constellation::Galileo, prn }
    }

    fn make_obs(sat: SatelliteId, pr: f64, cp: Option<f64>) -> DdObservation {
        DdObservation {
            sat,
            pr_l1: pr,
            pr_l2: None,
            cp_l1: cp,
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        }
    }

    fn pair(sat: SatelliteId, pr: f64) -> (DdObservation, DdObservation) {
        (make_obs(sat, pr, None), make_obs(sat, pr, None))
    }

    fn make_state(time: GpsTime) -> RtkState {
        let pos = Coordinate::new(
            Vector3::new(6378137.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        RtkState::new(time, pos, 10.0)
    }

    fn make_eph(sat: SatelliteId, time: GpsTime, omega0: f64) -> Ephemeris {
        Ephemeris::Gps(GpsEphemeris {
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
            m0: 0.0,
            e: 0.01,
            sqrt_a: 5153.6,
            delta_n: 0.0,
            omega0,
            omega_dot: 0.0,
            i0: 1.0,
            idot: 0.0,
            omega: 0.0,
            tgd: 0.0,
            iode: 0,
            iodc: 0,
        })
    }

    fn make_env<'a>(
        ephs: &'a [Ephemeris],
        base: &'a Coordinate,
        tuning: &'a EkfTuningConfig,
        time: GpsTime,
    ) -> MeasurementEnvironment<'a> {
        MeasurementEnvironment {
            ephemerides: ephs,
            base_coord: base,
            base_time: time,
            lever_arm: Vector3::zeros(),
            omega_b: Vector3::zeros(),
            tuning,
            gnn_variances: HashMap::new(),
        }
    }

    // ===================================================================
    // group_measurements_by_constellation
    // ===================================================================

    #[test]
    fn test_group_empty_input() {
        let groups = group_measurements_by_constellation(&[]);
        assert!(groups.is_empty());
    }

    #[test]
    fn test_group_single_constellation() {
        let input = vec![pair(gps(1), 100.0), pair(gps(2), 200.0), pair(gps(3), 300.0)];
        let groups = group_measurements_by_constellation(&input);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups.get(&Constellation::Gps).unwrap().len(), 3);
    }

    #[test]
    fn test_group_multi_constellation() {
        let input = vec![pair(gps(1), 100.0), pair(gal(10), 200.0)];
        let groups = group_measurements_by_constellation(&input);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[&Constellation::Gps].len(), 1);
        assert_eq!(groups[&Constellation::Galileo].len(), 1);
    }

    #[test]
    fn test_group_data_integrity() {
        // Verify that the observation data survives grouping unchanged
        let input = vec![pair(gps(1), 42.5)];
        let groups = group_measurements_by_constellation(&input);
        let group = &groups[&Constellation::Gps];
        assert_eq!(group[0].0.pr_l1, 42.5);
        assert_eq!(group[0].1.pr_l1, 42.5);
    }

    // ===================================================================
    // select_reference_satellite
    // ===================================================================

    #[test]
    fn test_select_ref_prefers_carrier_phase() {
        let time = GpsTime::new(2137, 422922.0);
        let st = make_state(time);
        let sats = [gps(1), gps(2), gps(3)];

        // Only sat[0] has carrier phase -> gets +100 elevation bonus
        let group = vec![
            (make_obs(sats[0], 26560000.0, Some(0.0)), make_obs(sats[0], 26560000.0, None)),
            (make_obs(sats[1], 26570000.0, None), make_obs(sats[1], 26570000.0, None)),
            (make_obs(sats[2], 26580000.0, None), make_obs(sats[2], 26580000.0, None)),
        ];

        let ephs: Vec<_> = sats
            .iter()
            .enumerate()
            .map(|(i, &s)| make_eph(s, time, i as f64 * 0.1))
            .collect();
        let base = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let tuning = EkfTuningConfig::default();
        let env = make_env(&ephs, &base, &tuning, time);

        let idx = select_reference_satellite(&group, &st, &env);
        assert_eq!(idx, 0, "satellite with cp_l1 should be selected over those without");
    }

    #[test]
    fn test_select_ref_fallback_on_missing_ephemeris() {
        let time = GpsTime::new(2137, 422922.0);
        let st = make_state(time);
        let group = vec![pair(gps(1), 26560000.0), pair(gps(2), 26570000.0)];

        let base = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let tuning = EkfTuningConfig::default();
        // No ephemerides provided
        let env = make_env(&[], &base, &tuning, time);

        // Without ephemerides the function falls through the inner loop
        // and returns the default index 0
        let idx = select_reference_satellite(&group, &st, &env);
        assert_eq!(idx, 0, "should return first index when no ephemeris found");
    }

    // ===================================================================
    // update_reject_counts
    // ===================================================================

    #[test]
    fn test_reject_passed_resets_count() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat = gps(1);
        state.add_ambiguity(sat, 1, 0.0, 1.0);
        state.reject_counts.insert((sat, 1), 5);

        let state_size = CORE_STATE_SIZE + 1;
        let mut h_row = DMatrix::zeros(1, state_size);
        h_row[(0, CORE_STATE_SIZE)] = 1.0;

        update_reject_counts(&mut state, &h_row, state_size, true);
        assert_eq!(state.reject_counts[&(sat, 1)], 0);
    }

    #[test]
    fn test_reject_failed_increments() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat = gps(1);
        state.add_ambiguity(sat, 1, 0.0, 1.0);

        let state_size = CORE_STATE_SIZE + 1;
        let mut h_row = DMatrix::zeros(1, state_size);
        h_row[(0, CORE_STATE_SIZE)] = 1.0;

        update_reject_counts(&mut state, &h_row, state_size, false);
        assert_eq!(state.reject_counts[&(sat, 1)], 1);

        update_reject_counts(&mut state, &h_row, state_size, false);
        assert_eq!(state.reject_counts[&(sat, 1)], 2);
    }

    #[test]
    fn test_reject_no_ambiguity_contribution() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat = gps(1);
        state.add_ambiguity(sat, 1, 0.0, 1.0);
        state.reject_counts.insert((sat, 1), 3);

        let state_size = CORE_STATE_SIZE + 1;
        // h_row is all zeros -> no ambiguity column > 0.5
        let h_row = DMatrix::zeros(1, state_size);

        update_reject_counts(&mut state, &h_row, state_size, false);
        // reject count unchanged because h_row doesn't connect to ambiguity
        assert_eq!(state.reject_counts[&(sat, 1)], 3);
    }

    #[test]
    fn test_reject_no_ambiguities_defined() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        // No ambiguities -> loop from CORE_STATE_SIZE..CORE_STATE_SIZE is empty
        let h_row = DMatrix::zeros(1, CORE_STATE_SIZE);

        update_reject_counts(&mut state, &h_row, CORE_STATE_SIZE, false);
        assert!(state.reject_counts.is_empty());
    }

    // ===================================================================
    // filter_innovations_chi_squared
    // ===================================================================

    #[test]
    fn test_chi2_all_pass() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sz = CORE_STATE_SIZE;

        // Zero innovations always pass (chi2 = 0)
        let z_all = vec![0.0, 0.0];
        let h_all = vec![vec![0.0; sz], vec![0.0; sz]];
        let r_all = vec![1.0, 1.0];
        let mt_all = vec![(gps(1), 0, 0.0), (gps(2), 0, 0.0)];

        let safe = filter_innovations_chi_squared(&mut state, sz, 3.0, 3.0, &z_all, &h_all, &r_all, &mt_all);
        assert_eq!(safe, vec![0, 1]);
    }

    #[test]
    fn test_chi2_some_rejected() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sz = CORE_STATE_SIZE;

        // Second measurement has large innovation -> fails
        let z_all = vec![0.0, 100.0];
        let h_all = vec![vec![0.0; sz], vec![0.0; sz]];
        let r_all = vec![1.0, 1.0];
        let mt_all = vec![(gps(1), 0, 0.0), (gps(2), 0, 0.0)];

        let safe = filter_innovations_chi_squared(&mut state, sz, 3.0, 3.0, &z_all, &h_all, &r_all, &mt_all);
        assert_eq!(safe, vec![0]);
    }

    #[test]
    fn test_chi2_all_rejected() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sz = CORE_STATE_SIZE;

        // Both measurements have large innovations -> both fail
        let z_all = vec![50.0, 100.0];
        let h_all = vec![vec![0.0; sz], vec![0.0; sz]];
        let r_all = vec![1.0, 1.0];
        let mt_all = vec![(gps(1), 0, 0.0), (gps(2), 0, 0.0)];

        let safe = filter_innovations_chi_squared(&mut state, sz, 3.0, 3.0, &z_all, &h_all, &r_all, &mt_all);
        assert!(safe.is_empty());
    }

    #[test]
    fn test_chi2_doppler_uses_larger_threshold() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sz = CORE_STATE_SIZE;

        // Type 3 (Doppler) threshold = chi_pr * 1000 = 3.0 * 1000.0 = 3000.0
        // z=50 -> chi2 = 2500 <= 3000 -> pass
        // z=100 -> chi2 = 10000 > 3000 -> fail
        let z_all = vec![50.0, 100.0];
        let h_all = vec![vec![0.0; sz], vec![0.0; sz]];
        let r_all = vec![1.0, 1.0];
        let mt_all = vec![(gps(1), 3, 0.0), (gps(2), 3, 0.0)];

        let safe = filter_innovations_chi_squared(&mut state, sz, 3.0, 3.0, &z_all, &h_all, &r_all, &mt_all);
        assert_eq!(safe, vec![0], "doppler with z=50 should pass wider threshold");
    }

    #[test]
    fn test_chi2_unknown_type_falls_back_to_pr_threshold() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sz = CORE_STATE_SIZE;

        // Type 99 (unknown) -> falls to _ => chi_pr * chi_pr = 9.0
        // z=3 -> chi2 = 9 <= 9 -> pass
        // z=4 -> chi2 = 16 > 9 -> fail
        let z_all = vec![3.0, 4.0];
        let h_all = vec![vec![0.0; sz], vec![0.0; sz]];
        let r_all = vec![1.0, 1.0];
        let mt_all = vec![(gps(1), 99, 0.0), (gps(2), 99, 0.0)];

        let safe = filter_innovations_chi_squared(&mut state, sz, 3.0, 3.0, &z_all, &h_all, &r_all, &mt_all);
        assert_eq!(safe, vec![0]);
    }

    #[test]
    fn test_chi2_cp_rejected_increments_reject_count() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat = gps(1);
        state.add_ambiguity(sat, 1, 0.0, 1.0);

        let sz = CORE_STATE_SIZE + 1;
        let mut h_row = vec![0.0; sz];
        h_row[CORE_STATE_SIZE] = 1.0;

        // Type 1 (CP L1) with large innovation -> rejected, increments reject count
        let z_all = vec![100.0];
        let h_all = vec![h_row];
        let r_all = vec![1.0];
        let mt_all = vec![(sat, 1, 0.0)];

        let safe = filter_innovations_chi_squared(&mut state, sz, 3.0, 3.0, &z_all, &h_all, &r_all, &mt_all);
        assert!(safe.is_empty());
        assert_eq!(state.reject_counts[&(sat, 1)], 1);
    }

    #[test]
    fn test_chi2_cp_accepted_resets_reject_count() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat = gps(1);
        state.add_ambiguity(sat, 1, 0.0, 1.0);
        state.reject_counts.insert((sat, 1), 5);

        let sz = CORE_STATE_SIZE + 1;
        let mut h_row = vec![0.0; sz];
        h_row[CORE_STATE_SIZE] = 1.0;

        // Zero innovation -> passes, resets reject count
        let z_all = vec![0.0];
        let h_all = vec![h_row];
        let r_all = vec![1.0];
        let mt_all = vec![(sat, 1, 0.0)];

        let safe = filter_innovations_chi_squared(&mut state, sz, 3.0, 3.0, &z_all, &h_all, &r_all, &mt_all);
        assert_eq!(safe, vec![0]);
        assert_eq!(state.reject_counts[&(sat, 1)], 0);
    }

    // ===================================================================
    // build_final_measurement_matrices
    // ===================================================================

    #[test]
    fn test_build_final_sufficient_measurements() {
        let sz = CORE_STATE_SIZE;
        let safe = vec![0, 1, 2, 3];
        let z_all = vec![10.0, 20.0, 30.0, 40.0];
        let h_all: Vec<Vec<f64>> = (0..4)
            .map(|i| {
                let mut row = vec![0.0; sz];
                row[0] = (i + 1) as f64;
                row
            })
            .collect();
        let r_all = vec![0.1, 0.2, 0.3, 0.4];
        let mt_all = vec![
            (gps(1), 0, 0.0),
            (gps(2), 0, 0.0),
            (gps(3), 0, 0.0),
            (gps(4), 0, 0.0),
        ];

        let result = build_final_measurement_matrices(sz, safe, &z_all, &h_all, &r_all, &mt_all);
        assert!(result.is_some());
        let m = result.unwrap();

        assert_eq!(m.z.len(), 4);
        assert_eq!(m.h.nrows(), 4);
        assert_eq!(m.h.ncols(), sz);
        assert_eq!(m.r.nrows(), 4);
        assert_eq!(m.r.ncols(), 4);
        assert_eq!(m.mt.len(), 4);

        assert_eq!(m.z[0], 10.0);
        assert_eq!(m.z[3], 40.0);
        assert_eq!(m.h[(1, 0)], 2.0);
        assert_eq!(m.h[(3, 0)], 4.0);
    }

    #[test]
    fn test_build_final_insufficient_measurements() {
        let sz = CORE_STATE_SIZE;
        let safe = vec![0, 1, 2]; // only 3, minimum is 4

        let z_all = vec![1.0; 5];
        let h_all = vec![vec![0.0; sz]; 5];
        let r_all = vec![0.1; 5];
        let mt_all = vec![(gps(1), 0, 0.0); 5];

        let result = build_final_measurement_matrices(sz, safe, &z_all, &h_all, &r_all, &mt_all);
        assert!(result.is_none());
    }

    #[test]
    fn test_build_final_empty_safe_indices() {
        let sz = CORE_STATE_SIZE;
        let safe = vec![];

        let z_all = vec![1.0; 5];
        let h_all = vec![vec![0.0; sz]; 5];
        let r_all = vec![0.1; 5];
        let mt_all = vec![(gps(1), 0, 0.0); 5];

        let result = build_final_measurement_matrices(sz, safe, &z_all, &h_all, &r_all, &mt_all);
        assert!(result.is_none());
    }

    #[test]
    fn test_build_final_noncontiguous_safe_indices() {
        let sz = CORE_STATE_SIZE;
        // Non-contiguous safe indices
        let safe = vec![0, 2, 4, 6];
        let z_all = vec![10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0];
        let h_all: Vec<Vec<f64>> = (0..7)
            .map(|i| {
                let mut row = vec![0.0; sz];
                row[0] = i as f64;
                row
            })
            .collect();
        let r_all = vec![0.1; 7];
        let mt_all: Vec<_> = (0..7).map(|i| (gps(i as u8 + 1), 0, 0.0)).collect();

        let result = build_final_measurement_matrices(sz, safe, &z_all, &h_all, &r_all, &mt_all);
        assert!(result.is_some());
        let m = result.unwrap();

        assert_eq!(m.z.len(), 4);
        assert_eq!(m.z[0], 10.0); // original index 0
        assert_eq!(m.z[1], 30.0); // original index 2
        assert_eq!(m.z[2], 50.0); // original index 4
        assert_eq!(m.z[3], 70.0); // original index 6
        assert_eq!(m.mt[2].0, gps(5)); // original index 4 -> sat 5
    }

    // ===================================================================
    // build_measurement_model  (public entry point)
    // ===================================================================

    #[test]
    fn test_build_model_empty_matched_obs() {
        let time = GpsTime::new(2137, 422922.0);
        let mut st = make_state(time);
        let base = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let tuning = EkfTuningConfig::default();
        let env = make_env(&[], &base, &tuning, time);

        let result = build_measurement_model(&mut st, &[], &env, 3.0, 3.0);
        assert!(result.is_none(), "empty observations should yield None");
    }

    #[test]
    fn test_build_model_single_constellation_pair_skipped() {
        // Two satellites in the same constellation form a group of 2,
        // which passes the len check.  After removing the reference one
        // satellite remains, and compute_innovations is called.  If
        // valid ephemerides are provided the function should not panic.
        let time = GpsTime::new(2137, 422922.0);
        let mut st = make_state(time);
        let sats = [gps(1), gps(2)];

        let obs = vec![
            (make_obs(sats[0], 26560000.0, Some(0.0)), make_obs(sats[0], 26560000.0, None)),
            (make_obs(sats[1], 26570000.0, None), make_obs(sats[1], 26570000.0, None)),
        ];

        let ephs: Vec<_> = sats.iter().map(|&s| make_eph(s, time, 0.0)).collect();
        let base = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let tuning = EkfTuningConfig::default();
        let env = make_env(&ephs, &base, &tuning, time);

        // Should not panic; the result may be None since synthetic
        // observations won't produce realistic updates.
        let _result = build_measurement_model(&mut st, &obs, &env, 3.0, 3.0);
    }
}
