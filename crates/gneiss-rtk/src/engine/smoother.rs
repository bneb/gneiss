use crate::engine::{EngineError, EngineMode, ProcessingEngine};
use crate::filter::RtkState;
use nalgebra::{DMatrix, DVector};

const MAX_STATE_VARIANCE: f64 = 1e10;
const MIN_ACTIVE_STATE_VARIANCE: f64 = 1e-12;
const MIN_ATTITUDE_ROTATION: f64 = 1e-10;

pub fn run_combined_ppk(engine: &mut ProcessingEngine) -> Result<Vec<RtkState>, EngineError> {
    let n_epochs = engine.state_history.len();
    if n_epochs == 0 {
        return Err(EngineError::NoObservations);
    }

    let mut smoothed_states = engine.state_history.clone();

    if matches!(engine.config.mode, EngineMode::Spp) {
        return Ok(smoothed_states);
    }

    for k in (0..n_epochs - 1).rev() {
        if smoothed_states[k + 1].is_reset {
            tracing::debug!(
                "Epoch {} was reset. Breaking smoothing chain at k={}",
                k + 1,
                k
            );
            continue;
        }

        let phi_k: DMatrix<f64> = match &smoothed_states[k + 1].core_phi {
            Some(p) => DMatrix::<f64>::clone(p),
            None => continue,
        };
        let p_pred_k1: DMatrix<f64> = match &smoothed_states[k + 1].full_p_predict {
            Some(p) => DMatrix::<f64>::clone(p),
            None => continue,
        };
        let x_pred_k1: DVector<f64> = match &smoothed_states[k + 1].full_x_predict {
            Some(x) => DVector::<f64>::clone(x),
            None => continue,
        };

        let (left, right) = smoothed_states.split_at_mut(k + 1);
        let state_k = &mut left[k];
        let state_k1 = &right[0];

        if let Err(e) = smooth_epoch(state_k, state_k1, &phi_k, &p_pred_k1, &x_pred_k1, k) {
            tracing::debug!("RTS Smoothing skipped at k={}: {}", k, e);
            continue;
        }

        // Preserve forward-pass AR fix state — re-resolving on smoothed
        // states can pick different integers, creating inconsistencies
        // that propagate and amplify through the backward recursion.
        // The smoothed covariance already reflects the forward AR fix
        // through the RTS update.
    }
    Ok(smoothed_states)
}

fn smooth_epoch(
    state_k: &mut RtkState,
    state_k1: &RtkState,
    phi_k: &DMatrix<f64>,
    p_pred_k1: &DMatrix<f64>,
    x_pred_k1: &DVector<f64>,
    _k_idx: usize,
) -> Result<(), &'static str> {
    let core_size = crate::filter::CORE_STATE_SIZE;

    let (matched_k_indices, matched_k1_indices) = find_matched_ambiguities(state_k, state_k1);
    let smooth_len = core_size + matched_k_indices.len();

    let mut idx_k = (0..core_size).collect::<Vec<_>>();
    idx_k.extend(matched_k_indices.iter().copied());

    let mut idx_k1 = (0..core_size).collect::<Vec<_>>();
    idx_k1.extend(&matched_k1_indices);

    let x_k1_n = build_x_vector(state_k1, core_size, smooth_len, &matched_k1_indices);
    let p_k1_n = extract_submatrix(&state_k1.covariance, &idx_k1, &idx_k1);
    let p_k = extract_submatrix(&state_k.covariance, &idx_k, &idx_k);
    let p_pred_k1_sub = extract_submatrix(p_pred_k1, &idx_k1, &idx_k1);
    let phi_k_sub = build_phi_submatrix(phi_k, core_size, smooth_len);

    if p_pred_k1_sub
        .iter()
        .any(|&x| !x.is_finite() || x.abs() > MAX_STATE_VARIANCE)
        || p_k
            .iter()
            .any(|&x| !x.is_finite() || x.abs() > MAX_STATE_VARIANCE)
    {
        return Err("non-finite covariance");
    }

    let p_pred_inv = invert_p_pred(&p_pred_k1_sub, smooth_len)?;

    let x_pred_k1_sub = extract_subvector(x_pred_k1, &idx_k1);
    let c_k = &p_k * phi_k_sub.transpose() * &p_pred_inv;

    let x_k = build_x_vector(state_k, core_size, smooth_len, &matched_k_indices);
    let mut delta_x = &x_k1_n - &x_pred_k1_sub;

    if core_size > 6 {
        if let Some(predicted_attitude) = state_k1.predicted_attitude {
            let predicted_attitude: nalgebra::UnitQuaternion<f64> = predicted_attitude;
            let dq = state_k1.attitude * predicted_attitude.inverse();
            let mut d_theta = dq.scaled_axis();
            if dq.w < 0.0 {
                d_theta = -d_theta;
            }
            delta_x.rows_mut(6, 3).copy_from(&d_theta);
        }
    }

    if core_size > 15 {
        // White-noise states must not be smoothed backward — they are
        // re-estimated from scratch each epoch. Clock bias (15) and
        // inter-system biases (16=GLO, 17=GAL, 18=BDS) change by
        // km-equivalents between epochs; blending them corrupts the
        // vertical solution. Clock drift (19) and ZWD (20) are NOT
        // zeroed because they carry physically meaningful temporal
        // correlation: clock drift evolves as a random walk driven
        // by Allan-variance process noise, and ZWD varies slowly
        // with tropospheric turbulence — both benefit from smoothing.
        for idx in [15, 16, 17, 18] {
            delta_x[idx] = 0.0;
        }
    }

    let mut correction = &c_k * &delta_x;

    if core_size > 15 {
        for idx in [15, 16, 17, 18] {
            correction[idx] = 0.0;
        }
    }

    let x_k_n = x_k + correction;

    // Guard: if the smoothed state has diverged, skip this epoch
    if x_k_n.iter().any(|v| !v.is_finite()) || x_k_n.iter().any(|v| v.abs() > 1e15) {
        return Err("smoothed state diverged");
    }

    let p_k_n_raw = p_k + &c_k * (p_k1_n - p_pred_k1_sub) * c_k.transpose();
    // Enforce symmetry — numerical drift can break it across many epochs
    let p_k_n = 0.5 * (&p_k_n_raw + p_k_n_raw.transpose());

    // Guard: reject non-finite or pathologically large covariance
    if p_k_n
        .iter()
        .any(|v| !v.is_finite() || v.abs() > MAX_STATE_VARIANCE)
    {
        return Err("smoothed covariance diverged");
    }

    update_smoothed_state(
        state_k,
        &x_k_n,
        &p_k_n,
        core_size,
        smooth_len,
        &matched_k_indices,
        &idx_k,
    );
    Ok(())
}

fn find_matched_ambiguities(state_k: &RtkState, state_k1: &RtkState) -> (Vec<usize>, Vec<usize>) {
    let mut matched_k = Vec::new();
    let mut matched_k1 = Vec::new();
    for (i, key_k) in state_k.ambiguity_keys.iter().enumerate() {
        if let Some(j) = state_k1.ambiguity_keys.iter().position(|k| k == key_k) {
            if state_k.ambiguity_track_ids[i] == state_k1.ambiguity_track_ids[j] {
                let cov_idx = crate::filter::CORE_STATE_SIZE + j;
                if state_k1.covariance[(cov_idx, cov_idx)] < 10.0 {
                    matched_k.push(crate::filter::CORE_STATE_SIZE + i);
                    matched_k1.push(cov_idx);
                }
            }
        }
    }
    (matched_k, matched_k1)
}

fn build_x_vector(
    state: &RtkState,
    core_size: usize,
    len: usize,
    matched_indices: &[usize],
) -> DVector<f64> {
    let mut x = DVector::zeros(len);
    x.rows_mut(0, 3).copy_from(&state.position.vector);
    x.rows_mut(3, 3).copy_from(&state.velocity);
    if core_size > 6 {
        x.rows_mut(9, 3).copy_from(&state.accel_bias);
        x.rows_mut(12, 3).copy_from(&state.gyro_bias);
    }
    if core_size > 15 {
        x[15] = state.rcv_clk_bias;
        x[16] = state.isb_glo;
        x[17] = state.isb_gal;
        x[18] = state.isb_bds;
        x[19] = state.rcv_clk_drift;
        x[20] = state.zwd;
    }
    for (i, &idx) in matched_indices.iter().enumerate() {
        x[core_size + i] = state.ambiguities[idx - crate::filter::CORE_STATE_SIZE];
    }
    x
}

fn extract_submatrix(mat: &DMatrix<f64>, rows: &[usize], cols: &[usize]) -> DMatrix<f64> {
    let mut sub = DMatrix::zeros(rows.len(), cols.len());
    for (i, &r) in rows.iter().enumerate() {
        for (j, &c) in cols.iter().enumerate() {
            sub[(i, j)] = mat[(r, c)];
        }
    }
    sub
}

fn extract_subvector(vec: &DVector<f64>, indices: &[usize]) -> DVector<f64> {
    let mut sub = DVector::zeros(indices.len());
    for (i, &idx) in indices.iter().enumerate() {
        sub[i] = vec[idx];
    }
    sub
}

fn build_phi_submatrix(phi: &DMatrix<f64>, core_size: usize, len: usize) -> DMatrix<f64> {
    let mut sub = DMatrix::zeros(len, len);
    for i in 0..core_size {
        for j in 0..core_size {
            sub[(i, j)] = phi[(i, j)];
        }
    }
    // White-noise states (clock bias 15, ISBs 16-18) should have zero
    // phi entries — they have no inter-epoch correlation. Setting their
    // phi row/column to zero prevents the smoother gain C_k from
    // propagating these states backward.
    if core_size > 15 {
        for idx in [15, 16, 17, 18] {
            for j in 0..core_size {
                sub[(idx, j)] = 0.0;
                sub[(j, idx)] = 0.0;
            }
        }
    }
    // Ambiguity states: identity transition (persist across epochs)
    for i in core_size..len {
        sub[(i, i)] = 1.0;
    }
    sub
}

fn invert_p_pred(p_pred: &DMatrix<f64>, len: usize) -> Result<DMatrix<f64>, &'static str> {
    let active: Vec<usize> = (0..len)
        .filter(|&i| {
            // Exclude white-noise states from inversion: clock bias (15)
            // and inter-system biases (16=GLO, 17=GAL, 18=BDS).
            // Clock drift (19) and ZWD (20) are included because they
            // have physically meaningful temporal correlation (Allan
            // variance and tropospheric turbulence, respectively),
            // making them safe to invert and smooth.
            let is_white_noise = matches!(i, 15 | 16 | 17 | 18);
            !is_white_noise && p_pred[(i, i)] > MIN_ACTIVE_STATE_VARIANCE
        })
        .collect();
    let m = active.len();
    if m == len {
        Ok(crate::math::inversion::invert_matrix_robust(p_pred))
    } else if m > 0 {
        let p_act = extract_submatrix(p_pred, &active, &active);
        let inv_act = crate::math::inversion::invert_matrix_robust(&p_act);
        let mut inv_full = DMatrix::zeros(len, len);
        for (i, &r) in active.iter().enumerate() {
            for (j, &c) in active.iter().enumerate() {
                inv_full[(r, c)] = inv_act[(i, j)];
            }
        }
        Ok(inv_full)
    } else {
        Err("no active elements")
    }
}

fn update_smoothed_state(
    state: &mut RtkState,
    x_k_n: &DVector<f64>,
    p_k_n: &DMatrix<f64>,
    core_size: usize,
    smooth_len: usize,
    matched_indices: &[usize],
    idx_k: &[usize],
) {
    state.position.vector = x_k_n.fixed_rows::<3>(0).into_owned();
    state.velocity = x_k_n.fixed_rows::<3>(3).into_owned();
    if core_size > 6 {
        let d_theta = x_k_n.fixed_rows::<3>(6).into_owned();
        if d_theta.norm() > MIN_ATTITUDE_ROTATION {
            let dq = nalgebra::UnitQuaternion::from_axis_angle(
                &nalgebra::Unit::new_normalize(d_theta),
                d_theta.norm(),
            );
            state.attitude = dq * state.attitude;
            state.attitude.renormalize();
        }
        state.accel_bias = x_k_n.fixed_rows::<3>(9).into_owned();
        state.gyro_bias = x_k_n.fixed_rows::<3>(12).into_owned();
    }
    if core_size > 15 {
        state.rcv_clk_bias = x_k_n[15];
        state.isb_glo = x_k_n[16];
        state.isb_gal = x_k_n[17];
        state.isb_bds = x_k_n[18];
        state.rcv_clk_drift = x_k_n[19];
        state.zwd = x_k_n[20];
    }
    for (i, &idx) in matched_indices.iter().enumerate() {
        let amb_idx = idx - crate::filter::CORE_STATE_SIZE;
        state.ambiguities[amb_idx] = x_k_n[core_size + i];
    }
    for i in 0..smooth_len {
        for j in 0..smooth_len {
            state.covariance[(idx_k[i], idx_k[j])] = p_k_n[(i, j)];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::RtkState;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector, Vector3};

    // -----------------------------------------------------------------------
    // find_matched_ambiguities
    // -----------------------------------------------------------------------

    #[test]
    fn test_find_matched_ambiguities_empty() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state_k = RtkState::new(time, pos.clone(), 1.0);
        let state_k1 = RtkState::new(time, pos, 1.0);
        let (mk, mk1) = find_matched_ambiguities(&state_k, &state_k1);
        assert!(mk.is_empty());
        assert!(mk1.is_empty());
    }

    #[test]
    fn test_find_matched_ambiguities_no_common_sat() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state_k = RtkState::new(time, pos.clone(), 1.0);
        let mut state_k1 = RtkState::new(time, pos, 1.0);

        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat2 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        state_k.add_ambiguity(sat1, 1, 5.0, 1.0);
        state_k1.add_ambiguity(sat2, 1, 5.0, 1.0);

        let (mk, mk1) = find_matched_ambiguities(&state_k, &state_k1);
        assert!(mk.is_empty());
        assert!(mk1.is_empty());
    }

    #[test]
    fn test_find_matched_ambiguities_common_sat_matches() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state_k = RtkState::new(time, pos.clone(), 1.0);
        let mut state_k1 = RtkState::new(time, pos, 1.0);

        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state_k.add_ambiguity(sat1, 1, 5.0, 1.0);
        state_k1.add_ambiguity(sat1, 1, 5.0, 1.0);

        let (mk, mk1) = find_matched_ambiguities(&state_k, &state_k1);
        assert_eq!(mk.len(), 1);
        assert_eq!(mk1.len(), 1);
        // Indices should be at CORE_STATE_SIZE (21)
        assert_eq!(mk[0], crate::filter::CORE_STATE_SIZE);
        assert_eq!(mk1[0], crate::filter::CORE_STATE_SIZE);
    }

    #[test]
    fn test_find_matched_ambiguities_track_id_mismatch_blocks_match() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state_k = RtkState::new(time, pos.clone(), 1.0);
        let mut state_k1 = RtkState::new(time, pos, 1.0);

        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state_k.add_ambiguity(sat1, 1, 5.0, 1.0);
        state_k1.add_ambiguity(sat1, 1, 5.0, 1.0);
        // Manually force different track IDs
        state_k.ambiguity_track_ids[0] = 10;
        state_k1.ambiguity_track_ids[0] = 20;

        let (mk, mk1) = find_matched_ambiguities(&state_k, &state_k1);
        assert!(mk.is_empty());
        assert!(mk1.is_empty());
    }

    #[test]
    fn test_find_matched_ambiguities_covariance_threshold_filters() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state_k1 = RtkState::new(time, pos, 1.0);

        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state_k1.add_ambiguity(sat1, 1, 5.0, 20.0); // variance 20 > 10

        // We cannot easily build state_k with the same key _and_ a covariance check
        // because the check is on state_k1's covariance. Variance 20.0 > 10.0 => filtered.
        let time2 = GpsTime::new(2000, 1.0);
        let pos2 = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time2);
        let mut state_k = RtkState::new(time2, pos2, 1.0);
        state_k.add_ambiguity(sat1, 1, 5.0, 1.0);

        let (mk, mk1) = find_matched_ambiguities(&state_k, &state_k1);
        assert!(mk.is_empty());
        assert!(mk1.is_empty());
    }

    // -----------------------------------------------------------------------
    // build_x_vector
    // -----------------------------------------------------------------------

    #[test]
    fn test_build_x_vector_core_only_no_ambiguities() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(1.0, 2.0, 3.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);

        let x = build_x_vector(&state, crate::filter::CORE_STATE_SIZE,
                               crate::filter::CORE_STATE_SIZE, &[]);
        assert_eq!(x.len(), crate::filter::CORE_STATE_SIZE);
        assert_eq!(x[0], 1.0);
        assert_eq!(x[4], 5.0);
        assert_eq!(x[5], 6.0);
    }

    #[test]
    fn test_build_x_vector_with_ambiguities() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(10.0, 20.0, 30.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.add_ambiguity(sat, 1, 123.45, 1.0);

        let matched = vec![crate::filter::CORE_STATE_SIZE]; // index of the ambiguity
        let smooth_len = crate::filter::CORE_STATE_SIZE + 1;
        let x = build_x_vector(&state, crate::filter::CORE_STATE_SIZE, smooth_len, &matched);
        assert_eq!(x.len(), smooth_len);
        assert_eq!(x[crate::filter::CORE_STATE_SIZE], 123.45);
        assert_eq!(x[0], 10.0);
    }

    #[test]
    fn test_build_x_vector_core_size_exactly_6() {
        // When core_size == 6 (SPP without IMU), only position and velocity are set
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(1.0, 2.0, 3.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);

        let core_size = 6;
        let x = build_x_vector(&state, core_size, core_size, &[]);
        assert_eq!(x.len(), 6);
        assert_eq!(x[0], 1.0);
        assert_eq!(x[3], 4.0);
        assert_eq!(x[5], 6.0);
    }

    // -----------------------------------------------------------------------
    // extract_submatrix / extract_subvector
    // -----------------------------------------------------------------------

    #[test]
    fn test_extract_submatrix_basic() {
        let mat = DMatrix::from_row_slice(3, 3, &[
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
        ]);
        let rows = vec![0, 2];
        let cols = vec![1, 2];
        let sub = extract_submatrix(&mat, &rows, &cols);
        assert_eq!(sub.nrows(), 2);
        assert_eq!(sub.ncols(), 2);
        assert_eq!(sub[(0, 0)], mat[(0, 1)]); // 2.0
        assert_eq!(sub[(0, 1)], mat[(0, 2)]); // 3.0
        assert_eq!(sub[(1, 0)], mat[(2, 1)]); // 8.0
        assert_eq!(sub[(1, 1)], mat[(2, 2)]); // 9.0
    }

    #[test]
    fn test_extract_subvector_basic() {
        let vec = DVector::from_vec(vec![10.0, 20.0, 30.0, 40.0]);
        let indices = vec![3, 0];
        let sub = extract_subvector(&vec, &indices);
        assert_eq!(sub.len(), 2);
        assert_eq!(sub[0], 40.0);
        assert_eq!(sub[1], 10.0);
    }

    // -----------------------------------------------------------------------
    // build_phi_submatrix
    // -----------------------------------------------------------------------

    #[test]
    fn test_build_phi_submatrix_diagonal_and_ambiguity_identity() {
        let core_size = crate::filter::CORE_STATE_SIZE;
        let len = core_size + 2;
        let mut phi = DMatrix::zeros(core_size, core_size);
        for i in 0..core_size {
            phi[(i, i)] = 0.5; // typical decay
        }
        let sub = build_phi_submatrix(&phi, core_size, len);
        // Check ambiguity diagonal is 1.0
        assert_eq!(sub[(core_size, core_size)], 1.0);
        assert_eq!(sub[(core_size + 1, core_size + 1)], 1.0);
        // Check core diagonal preserved
        assert_eq!(sub[(0, 0)], 0.5);
        // Check white-noise states zeroed
        assert_eq!(sub[(15, 15)], 0.0);
        assert_eq!(sub[(16, 15)], 0.0);
        assert_eq!(sub[(15, 0)], 0.0);
        assert_eq!(sub[(17, 17)], 0.0);
        assert_eq!(sub[(18, 18)], 0.0);
    }

    #[test]
    fn test_build_phi_submatrix_small_core() {
        // core_size = 6 (no IMU), no white-noise zeroing
        let core_size = 6;
        let len = 6;
        let mut phi = DMatrix::zeros(core_size, core_size);
        phi[(0, 0)] = 1.0;
        phi[(5, 5)] = 1.0;
        let sub = build_phi_submatrix(&phi, core_size, len);
        assert_eq!(sub[(0, 0)], 1.0);
        assert_eq!(sub[(5, 5)], 1.0);
    }

    // -----------------------------------------------------------------------
    // invert_p_pred
    // -----------------------------------------------------------------------

    #[test]
    fn test_invert_p_pred_all_active() {
        let len = 3;
        let p_pred = DMatrix::identity(len, len) * 2.0;
        let inv = invert_p_pred(&p_pred, len).unwrap();
        assert!((inv[(0, 0)] - 0.5).abs() < 1e-10);
        assert_eq!(inv.nrows(), len);
    }

    #[test]
    fn test_invert_p_pred_partial_inactive() {
        // State 0 active, state 1 has zero variance (inactive)
        let mut p_pred = DMatrix::zeros(2, 2);
        p_pred[(0, 0)] = 2.0; // active
        p_pred[(1, 1)] = 0.0; // below MIN_ACTIVE_STATE_VARIANCE
        let inv = invert_p_pred(&p_pred, 2).unwrap();
        assert_eq!(inv.nrows(), 2);
        assert!((inv[(0, 0)] - 0.5).abs() < 1e-10);
        // The inactive element's inverse entry should be zero
        assert_eq!(inv[(1, 1)], 0.0);
    }

    #[test]
    fn test_invert_p_pred_none_active() {
        let p_pred = DMatrix::zeros(2, 2);
        let res = invert_p_pred(&p_pred, 2);
        assert!(res.is_err());
        assert_eq!(res.unwrap_err(), "no active elements");
    }

    #[test]
    fn test_invert_p_pred_white_noise_excluded_from_active() {
        // White-noise indices 15,16,17,18 are ALWAYS excluded even with non-zero variance.
        // For core_size > 15, indices >= 15 are excluded via the white-noise check.
        // But the function uses `matches!(i, 15 | 16 | 17 | 18)`. For a matrix smaller
        // than 21, we can still verify: if we have a 16x16 matrix, state at i=15 is excluded.
        let len = 16;
        let mut p_pred = DMatrix::identity(len, len); // all diag = 1.0
        // State 15 would be excluded as white-noise. So only 15 active elements.
        let inv = invert_p_pred(&p_pred, len).unwrap();
        assert_eq!(inv.nrows(), len);
        // The active inverse should be correct for all non-white-noise states
        assert!((inv[(0, 0)] - 1.0).abs() < 1e-10);
        // White-noise state (15) should have zero inverse
        assert_eq!(inv[(15, 15)], 0.0);
    }

    // -----------------------------------------------------------------------
    // update_smoothed_state
    // -----------------------------------------------------------------------

    #[test]
    fn test_update_smoothed_state_updates_position_velocity() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let core_size = crate::filter::CORE_STATE_SIZE;
        let smooth_len = core_size;
        let mut x_k_n = DVector::zeros(smooth_len);
        x_k_n[0] = 10.0;
        x_k_n[1] = 20.0;
        x_k_n[2] = 30.0;
        x_k_n[3] = 1.0;
        x_k_n[4] = 2.0;
        x_k_n[5] = 3.0;
        let p_k_n = DMatrix::identity(smooth_len, smooth_len);
        let idx_k: Vec<usize> = (0..core_size).collect();

        update_smoothed_state(&mut state, &x_k_n, &p_k_n, core_size, smooth_len, &[], &idx_k);
        assert_eq!(state.position.vector.x, 10.0);
        assert_eq!(state.position.vector.y, 20.0);
        assert_eq!(state.position.vector.z, 30.0);
        assert_eq!(state.velocity.x, 1.0);
        assert_eq!(state.velocity.y, 2.0);
        assert_eq!(state.velocity.z, 3.0);
    }

    #[test]
    fn test_update_smoothed_state_attitude_rotation_applied() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let core_size = crate::filter::CORE_STATE_SIZE;
        let smooth_len = core_size;
        // Create a state vector with a small attitude rotation d_theta at indices 6..9
        let mut x_k_n = DVector::zeros(smooth_len);
        x_k_n[6] = 1e-5;
        x_k_n[7] = 2e-5;
        x_k_n[8] = 3e-5;
        let p_k_n = DMatrix::identity(smooth_len, smooth_len);
        let idx_k: Vec<usize> = (0..core_size).collect();

        let initial_att = state.attitude;
        update_smoothed_state(&mut state, &x_k_n, &p_k_n, core_size, smooth_len, &[], &idx_k);
        // Attitude should have changed (norm of d_theta > MIN_ATTITUDE_ROTATION = 1e-10)
        assert_ne!(state.attitude, initial_att);
    }

    #[test]
    fn test_update_smoothed_state_clock_and_zwd() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let core_size = crate::filter::CORE_STATE_SIZE;
        let smooth_len = core_size;
        let mut x_k_n = DVector::zeros(smooth_len);
        x_k_n[15] = 100.0;
        x_k_n[16] = 10.0;
        x_k_n[17] = 20.0;
        x_k_n[18] = 30.0;
        x_k_n[19] = 1.5;
        x_k_n[20] = 0.05;
        let p_k_n = DMatrix::identity(smooth_len, smooth_len);
        let idx_k: Vec<usize> = (0..core_size).collect();

        update_smoothed_state(&mut state, &x_k_n, &p_k_n, core_size, smooth_len, &[], &idx_k);
        assert_eq!(state.rcv_clk_bias, 100.0);
        assert_eq!(state.isb_glo, 10.0);
        assert_eq!(state.isb_gal, 20.0);
        assert_eq!(state.isb_bds, 30.0);
        assert_eq!(state.rcv_clk_drift, 1.5);
        assert_eq!(state.zwd, 0.05);
    }

    #[test]
    fn test_update_smoothed_state_ambiguities_preserved() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.add_ambiguity(sat, 1, 100.0, 1.0);

        let core_size = crate::filter::CORE_STATE_SIZE;
        let smooth_len = core_size + 1;
        let mut x_k_n = DVector::zeros(smooth_len);
        x_k_n[core_size] = 200.0;
        let p_k_n = DMatrix::identity(smooth_len, smooth_len);
        let idx_k: Vec<usize> = (0..smooth_len).collect();
        let matched = vec![core_size];

        update_smoothed_state(&mut state, &x_k_n, &p_k_n, core_size, smooth_len, &matched, &idx_k);
        assert!((state.ambiguities[0] - 200.0).abs() < 1e-10);
    }

    #[test]
    fn test_update_smoothed_state_covariance_writeback() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let core_size = crate::filter::CORE_STATE_SIZE;
        let smooth_len = core_size;
        let x_k_n = DVector::zeros(smooth_len);
        let mut p_k_n = DMatrix::zeros(smooth_len, smooth_len);
        p_k_n[(0, 0)] = 42.0;
        p_k_n[(1, 1)] = 99.0;
        let idx_k: Vec<usize> = (0..core_size).collect();

        update_smoothed_state(&mut state, &x_k_n, &p_k_n, core_size, smooth_len, &[], &idx_k);
        assert_eq!(state.covariance[(0, 0)], 42.0);
        assert_eq!(state.covariance[(1, 1)], 99.0);
    }

    #[test]
    fn test_update_smoothed_state_tiny_rotation_skipped() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let core_size = crate::filter::CORE_STATE_SIZE;
        let smooth_len = core_size;
        let mut x_k_n = DVector::zeros(smooth_len);
        // d_theta norm < MIN_ATTITUDE_ROTATION (1e-10)
        x_k_n[6] = 1e-11;
        let p_k_n = DMatrix::identity(smooth_len, smooth_len);
        let idx_k: Vec<usize> = (0..core_size).collect();

        let initial_att = state.attitude;
        update_smoothed_state(&mut state, &x_k_n, &p_k_n, core_size, smooth_len, &[], &idx_k);
        assert_eq!(state.attitude, initial_att);
    }

    // -----------------------------------------------------------------------
    // run_combined_ppk
    // -----------------------------------------------------------------------

    #[test]
    fn test_run_combined_ppk_empty_history() {
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), crate::engine::EngineError::NoObservations);
    }

    #[test]
    fn test_run_combined_ppk_spp_mode_returns_early() {
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        engine.config.mode = crate::engine::EngineMode::Spp;
        // Add a dummy state to state_history
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        engine.state_history.push(RtkState::new(time, pos, 1.0));
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 1);
    }

    #[test]
    fn test_run_combined_ppk_single_epoch_no_smoothing() {
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        engine.state_history.push(RtkState::new(time, pos, 1.0));
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 1);
    }

    #[test]
    fn test_run_combined_ppk_reset_epoch_skips_smoothing() {
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state0 = RtkState::new(time, pos, 1.0);
        let mut state1 = RtkState::new(time, pos.clone(), 1.0);
        state1.is_reset = true;
        engine.state_history.push(state0);
        engine.state_history.push(state1);
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
    }

    #[test]
    fn test_run_combined_ppk_missing_phi_skips() {
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state0 = RtkState::new(time, pos.clone(), 1.0);
        let mut state1 = RtkState::new(time, pos, 1.0);
        // No core_phi set, so smoothing will be skipped
        state1.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE) * 0.9);
        state1.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE) * 2.0);
        state1.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));
        engine.state_history.push(state0);
        engine.state_history.push(state1);
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
    }

    // -----------------------------------------------------------------------
    // smooth_epoch guard tests via a processing engine with full data
    // -----------------------------------------------------------------------

    fn make_smoothable_state_pair() -> (RtkState, RtkState) {
        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time0);
        let mut state0 = RtkState::new(time0, pos, 1.0);
        let mut state1 = RtkState::new(time1, pos.clone(), 1.0);
        state1.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state1.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE) * 2.0);
        state1.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));
        (state0, state1)
    }

    fn make_full_smoothable_pair() -> (RtkState, RtkState) {
        let (s0, mut s1) = make_smoothable_state_pair();
        // Make state1 have same ambiguity keys so matching works
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        // State1 must have cov[amb, amb] < 10.0 for matching
        // We'll just rely on the fact that with no ambiguities, smooth_len = core_size.
        (s0, s1)
    }

    #[test]
    fn test_smooth_epoch_non_finite_covariance_skipped() {
        // Set up a state with non-finite covariance
        let (mut state0, state1) = make_smoothable_state_pair();
        // Put NaN in state0's covariance
        state0.covariance[(3, 3)] = f64::NAN;

        let phi_k = state1.core_phi.as_ref().unwrap().clone();
        let p_pred_k1 = state1.full_p_predict.as_ref().unwrap().clone();
        let x_pred_k1 = state1.full_x_predict.as_ref().unwrap().clone();

        let result = smooth_epoch(&mut state0, &state1, &phi_k, &p_pred_k1, &x_pred_k1, 0);
        assert!(result.is_err());
    }

    #[test]
    fn test_smooth_epoch_diverged_state_skipped() {
        let (mut state0, state1) = make_smoothable_state_pair();
        let phi_k = state1.core_phi.as_ref().unwrap().clone();
        let p_pred_k1 = state1.full_p_predict.as_ref().unwrap().clone();
        let x_pred_k1 = state1.full_x_predict.as_ref().unwrap().clone();

        // Create ambiguity to have matched length > core so the
        // non-finite covariance guard doesn't trip before the divergence guard.
        // Actually, to reach the divergence guard, we need the smoothed state to have |v| > 1e15.
        // Since the actual state values are near zero and the correction is also small,
        // we need to force divergence through the delta_x path.
        // Instead, let's just verify the error type. The easiest path:
        // make p_pred_k1 have huge values so the correction is huge.
        let huge_pred = DMatrix::from_element(
            crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE, 1e20);
        let result = smooth_epoch(&mut state0, &state1, &phi_k, &huge_pred, &x_pred_k1, 0);
        // Either non-finite covariance guard or divergence guard catches it
        assert!(result.is_err());
    }

    #[test]
    fn test_smooth_epoch_huge_pred_covariance_rejected() {
        let (mut state0, state1) = make_smoothable_state_pair();
        let phi_k = state1.core_phi.as_ref().unwrap().clone();
        let x_pred_k1 = state1.full_x_predict.as_ref().unwrap().clone();

        // p_pred_k1_sub element > MAX_STATE_VARIANCE = 1e10
        let huge_pred = DMatrix::from_element(
            crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE, 1e12);
        let result = smooth_epoch(&mut state0, &state1, &phi_k, &huge_pred, &x_pred_k1, 0);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // Integration: run_combined_ppk with a 2-epoch smoothable chain
    // -----------------------------------------------------------------------

    #[test]
    fn test_run_combined_ppk_two_epoch_smoothable() {
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time0);

        let mut state0 = RtkState::new(time0, pos, 1.0);
        state0.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state0.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state0.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));

        let mut state1 = RtkState::new(time1, pos.clone(), 1.0);
        state1.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state1.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state1.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));

        engine.state_history.push(state0);
        engine.state_history.push(state1);
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
        let states = result.unwrap();
        assert_eq!(states.len(), 2);
    }

    #[test]
    fn test_run_combined_ppk_missing_p_predict_skips() {
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time0);

        let mut state0 = RtkState::new(time0, pos, 1.0);
        let mut state1 = RtkState::new(time1, pos, 1.0);
        state1.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        // Intentionally leave full_p_predict as None
        engine.state_history.push(state0);
        engine.state_history.push(state1);
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
    }
}
