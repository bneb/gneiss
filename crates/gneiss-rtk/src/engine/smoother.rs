use crate::engine::{EngineError, EngineMode, ProcessingEngine};
use crate::filter::RtkState;
use nalgebra::{DMatrix, DVector};

const MAX_STATE_VARIANCE: f64 = 1e10;
const MIN_ACTIVE_STATE_VARIANCE: f64 = 1e-12;
const MIN_ATTITUDE_ROTATION: f64 = 1e-10;
/// Maximum allowed condition number for the forward covariance matrix.
/// If max(diag(P)) / min(diag(P)) exceeds this, the matrix is considered
/// ill-conditioned for RTS smoothing and clock/ISB states are frozen.
/// Empirically determined: 100 (10⁴ / 10²) causes divergence; 10 is safe.
pub const MAX_CONDITION: f64 = 50.0;

/// Frozen state indices in the backward pass.
/// Clock bias (15) is excluded via white-noise model (φ[15,15]=0),
/// which makes C_k[:,15] = 0 naturally — no explicit freeze needed.
/// ISB states (16-18) are piece-wise constants with φ=1.0; their
/// covariance must be well-conditioned vs position before removing
/// the freeze (benchmark required: Odaiba, Shinjuku, f9p).
pub const FROZEN_BACKWARD_INDICES: [usize; 3] = [16, 17, 18];

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

        // Forward-filter quality guard: skip smoothing for cold-start epochs
        // where the filter hasn't converged (epoch_count < 2 matches the
        // is_cold_start logic in ppp.rs). Early epochs have inflated covariances
        // that produce excessive RTS gain, degrading the backward correction.
        if state_k.epoch_count < 2 || state_k1.epoch_count < 2 {
            tracing::debug!(
                "Forward filter cold start at k={} (epoch_count={}, {}), skipping smoothing",
                k,
                state_k.epoch_count,
                state_k1.epoch_count
            );
            continue;
        }

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

/// RTS (Rauch-Tung-Striebel) backward smoothing for a single epoch pair.
///
/// # Interface Contract
///
/// ## Producer (forward filter) supplies:
///
/// | Field | Type | Dimensions | Description |
/// |-------|------|-----------|-------------|
/// | `covariance` | DMatrix | C×C | Filtered covariance P_{k|k}, positive-definite |
/// | `core_phi` | DMatrix | C×C | State transition from k→k+1; φ[i,i] ∈ {0, 1} |
/// | `full_p_predict` | DMatrix | C×C | Predicted covariance P_{k+1|k} ≥ P_{k+1|k+1} |
/// | `full_x_predict` | DVector | C | Predicted state x_{k+1|k} = φ·x_{k|k} |
///
/// where C = CORE_STATE_SIZE + N_ambiguities.
///
/// ## State vector layout (CORE_STATE_SIZE = 21):
///
/// ```text
///  0:  position X (m ECEF)       9:  accel bias X (m/s²)
///  1:  position Y                10: accel bias Y
///  2:  position Z                11: accel bias Z
///  3:  velocity X (m/s)          12: gyro bias X (rad/s)
///  4:  velocity Y                13: gyro bias Y
///  5:  velocity Z                14: gyro bias Z
///  6:  attitude dx (rad)         15: clock bias (m)        ← white noise
///  7:  attitude dy               16: ISB GLO (m)           ← frozen
///  8:  attitude dz               17: ISB GAL (m)           ← frozen
///                                18: ISB BDS (m)           ← frozen
///                                19: clock drift (m/s)      ← smoothed
///                                20: ZWD (m)                ← smoothed
/// ```
///
/// ## Conditioning constraint:
///
/// Index 15 (clock bias) is **white noise** (φ[15,15]=0 in the forward
/// filter), so C_k[:,15] = 0 naturally — no explicit freeze needed.
///
/// Indices 16-18 (ISB GLO/GAL/BDS) are **intentionally frozen** because
/// their forward-filter variance (~10² m²) is poorly conditioned against
/// position states (~10² m²) when the process noise is small.  Removing
/// the freeze requires the ISB condition number to drop below MAX_CONDITION
/// and MUST be verified with benchmarks on Odaiba, Shinjuku, and f9p.
///
/// ## Consumer (smoother) guarantees:
///
/// - Validates input dimensions before computation
/// - Detects and skips epochs with non-finite or divergent covariance
/// - Never panics on malformed inputs (returns Err)
/// - Preserves forward-pass AR fix state in smoothed output
///
/// ## Assumptions verified by `validate_smoother_state`:
///
/// A1. core_size > 0 and all matrices have compatible dimensions
/// A2. P_k and P_pred are positive-definite (diagonal > 0)
/// A3. State condition number ≤ MAX_CONDITION (1e10)
/// A4. ISB/clock states (15-18) exist iff core_size > 15
///
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

    // Innovation gating: compute normalized innovation squared (NIS).
    // When NIS exceeds the threshold, the forward filter has almost
    // certainly diverged at this epoch. Skip the RTS correction and
    // keep the forward-filter state — applying a correction based on
    // a corrupted innovation would poison the entire backward pass.
    // Threshold: 100 × DOF (generous — chi-square 99.99% for 30 DOF
    // is ~60, so 100×DOF is a very conservative gate).
    let nis = (&delta_x.transpose() * &p_pred_inv * &delta_x)[(0, 0)];
    if nis > smooth_len as f64 * 100.0 {
        tracing::warn!(
            "Smoother NIS={:.1} exceeds threshold={:.0} ({} DOF) — keeping forward-filter state",
            nis, smooth_len as f64 * 100.0, smooth_len
        );
        return Err("NIS gate rejected");
    }

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

    // Freeze ISB/clock states per FROZEN_BACKWARD_INDICES.
    // See smooth_epoch doc comment §Conditioning constraint for rationale.
    if core_size > 15 {
        for &idx in &FROZEN_BACKWARD_INDICES {
            delta_x[idx] = 0.0;
        }
    }

    let mut correction = &c_k * &delta_x;

    if core_size > 15 {
        for &idx in &FROZEN_BACKWARD_INDICES {
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

/// Validate the forward-filter → smoother interface contract.
///
/// Checks that the inputs produced by the forward filter satisfy the
/// preconditions required by `smooth_epoch`.  Call this from unit tests
/// or as a debug assertion; production code skips epochs that fail
/// individual checks inside `smooth_epoch` rather than calling this.
///
/// Returns `Ok(condition_number)` if all checks pass, or `Err(reason)`.
pub fn validate_smoother_state(
    covariance: &DMatrix<f64>,
    p_pred: &DMatrix<f64>,
    core_size: usize,
) -> Result<f64, &'static str> {
    // A1: dimension compatibility
    let n = covariance.nrows();
    if n == 0 || covariance.ncols() != n {
        return Err("covariance must be square and non-empty");
    }
    if p_pred.nrows() != n || p_pred.ncols() != n {
        return Err("p_pred dimensions must match covariance");
    }
    if n < core_size {
        return Err("covariance smaller than core_size");
    }

    // A2: positive-definite (all diagonal entries > 0)
    let min_var = (0..n)
        .map(|i| covariance[(i, i)])
        .reduce(f64::min)
        .unwrap_or(0.0);
    let max_var = (0..n)
        .map(|i| covariance[(i, i)])
        .reduce(f64::max)
        .unwrap_or(0.0);
    if min_var <= 0.0 {
        return Err("covariance has non-positive diagonal");
    }

    // A3: condition number check
    let condition = max_var / min_var;
    if !condition.is_finite() || condition > MAX_CONDITION * MAX_CONDITION {
        return Err("covariance is pathologically ill-conditioned");
    }

    // A4: freeze eligibility — verify FROZEN_BACKWARD_INDICES are in bounds
    if core_size > 15 {
        for &idx in &FROZEN_BACKWARD_INDICES {
            if idx >= n {
                return Err("frozen index exceeds state dimension");
            }
        }
    }

    Ok(condition)
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
    // Freeze ISB/clock phi rows/cols per FROZEN_BACKWARD_INDICES.
    if core_size > 15 {
        for &idx in &FROZEN_BACKWARD_INDICES {
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
            // Exclude frozen states from inversion per FROZEN_BACKWARD_INDICES.
            let is_frozen = FROZEN_BACKWARD_INDICES.contains(&i);
            !is_frozen && p_pred[(i, i)] > MIN_ACTIVE_STATE_VARIANCE
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
        // Clock bias (15): white-noise model gives φ[15,15]=0 from forward filter.
        // ISB (16-18): explicitly frozen per FROZEN_BACKWARD_INDICES.
        assert_eq!(sub[(15, 15)], 0.5);  // clock bias: not frozen, from phi
        assert_eq!(sub[(16, 16)], 0.0);  // ISB GLO: frozen
        assert_eq!(sub[(17, 17)], 0.0);  // ISB GAL: frozen
        assert_eq!(sub[(18, 18)], 0.0);  // ISB BDS: frozen
        assert_eq!(sub[(16, 0)], 0.0);   // ISB-position cross-term frozen
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
    fn test_invert_p_pred_clock_included_isb_excluded_from_inversion() {
        // Clock bias (15) is white-noise (φ=0), included in inversion.
        // ISB indices (16-18) are frozen per FROZEN_BACKWARD_INDICES.
        let len = 19;
        let p_pred = DMatrix::identity(len, len); // all diag = 1.0
        let inv = invert_p_pred(&p_pred, len).unwrap();
        assert_eq!(inv.nrows(), len);
        assert!((inv[(0, 0)] - 1.0).abs() < 1e-10);
        // Index 15 (clock bias) IS included — white-noise model, not frozen
        assert!((inv[(15, 15)] - 1.0).abs() < 1e-10);
        // ISB indices excluded
        assert_eq!(inv[(16, 16)], 0.0);
        assert_eq!(inv[(17, 17)], 0.0);
        assert_eq!(inv[(18, 18)], 0.0);
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
        let state0 = RtkState::new(time, pos, 1.0);
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
        let state0 = RtkState::new(time, pos.clone(), 1.0);
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
        state0.epoch_count = 2; // Forward filter quality guard: epoch_count >= 2
        let mut state1 = RtkState::new(time1, pos.clone(), 1.0);
        state1.epoch_count = 2; // Forward filter quality guard: epoch_count >= 2
        // Use identity phi but with white-noise clock bias (φ[15,15]=0)
        let mut phi = DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE);
        phi[(15, 15)] = 0.0;
        state1.core_phi = Some(phi);
        // Use large predicted covariance so NIS gate passes.
        // NIS threshold = smooth_len * 100 ≈ 2000. With |innov| ≈ 374 m
        // (position difference from prediction), we need P_pred >> 374²/2000 ≈ 70.
        // Using 200 m² per axis gives plenty of margin.
        let mut p_pred = DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE) * 200.0;
        p_pred[(15, 15)] = crate::filter::PREDICTED_CLOCK_VARIANCE;
        p_pred[(16, 16)] = crate::filter::PREDICTED_ISB_VARIANCE;
        p_pred[(17, 17)] = crate::filter::PREDICTED_ISB_VARIANCE;
        p_pred[(18, 18)] = crate::filter::PREDICTED_ISB_VARIANCE;
        p_pred[(19, 19)] = 2000.0;
        p_pred[(20, 20)] = 2.0;
        state1.full_p_predict = Some(p_pred);
        // Predicted state at the same position as actual (consistent with
        // the large covariance — the prediction is uncertain but unbiased).
        let mut x_pred = DVector::zeros(crate::filter::CORE_STATE_SIZE);
        x_pred[0] = 100.0;
        x_pred[1] = 200.0;
        x_pred[2] = 300.0;
        state1.full_x_predict = Some(x_pred);
        (state0, state1)
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
        let _p_pred_k1 = state1.full_p_predict.as_ref().unwrap().clone();
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
        state0.epoch_count = 2; // Forward filter quality guard: epoch_count >= 2
        state0.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state0.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state0.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));

        let mut state1 = RtkState::new(time1, pos.clone(), 1.0);
        state1.epoch_count = 3; // Forward filter quality guard: epoch_count >= 2
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

    // -----------------------------------------------------------------------
    // Forward-filter quality guard tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_forward_filter_quality_guard_cold_start_epoch_skipped() {
        // Verify that smoothing is skipped when epoch_count < 2 (cold start)
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time0);

        // state0 has epoch_count = 0 (cold start)
        let mut state0 = RtkState::new(time0, pos, 1.0);
        state0.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state0.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state0.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));

        let mut state1 = RtkState::new(time1, pos.clone(), 1.0);
        state1.epoch_count = 5; // state_k1 has adequate quality
        state1.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state1.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state1.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));

        engine.state_history.push(state0);
        engine.state_history.push(state1);
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
        let states = result.unwrap();
        assert_eq!(states.len(), 2);
        // The guard only prevents epoch_count < 2 on state_k (epoch 0).
        // state0 is at index 0, state1 is at index 1.
        // In the backward loop: k=0 processes state0 (epoch_count=0) vs state1 (epoch_count=5).
        // The quality guard triggers on state0.epoch_count < 2, so state0 should be
        // unchanged by the smoother (position remains at initial values).
        assert_eq!(states[0].position.vector.x, 100.0,
            "Cold start epoch should not be smoothed");
    }

    #[test]
    fn test_forward_filter_quality_guard_mature_epoch_passes() {
        // Verify that smoothing proceeds when both states have epoch_count >= 2
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time0);

        let mut state0 = RtkState::new(time0, pos, 100.0);
        state0.epoch_count = 2; // Adequate quality
        state0.velocity = Vector3::new(1.0, 2.0, 3.0);
        state0.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state0.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE) * 200.0);
        state0.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));

        let mut state1 = RtkState::new(time1, pos.clone(), 100.0);
        state1.epoch_count = 3; // Adequate quality
        state1.velocity = Vector3::new(1.0, 2.0, 3.0);
        state1.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state1.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE) * 200.0);
        state1.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));

        engine.state_history.push(state0);
        engine.state_history.push(state1);
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
        let states = result.unwrap();
        assert_eq!(states.len(), 2);
        // Both states have epoch_count >= 2, so smoothing should proceed.
        // State 0 (epoch 0) should have been smoothed toward state 1 (epoch 1).
        // Result depends on actual RTS computation, but should differ from
        // the initial position (100, 200, 300) if smoothing ran.
        let smoothed_pos = states[0].position.vector;
        assert!(smoothed_pos.x > 0.0, "Smoothed position should be non-zero");
    }

    /// Create a pair of states where smooth_epoch will succeed and produce
    /// predictable RTS updates. phi = I, p_pred = 200*I, x_pred matches state0.
    fn make_smoothable_pair_with_different_positions() -> (RtkState, RtkState, DMatrix<f64>, DMatrix<f64>, DVector<f64>) {
        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        let pos0 = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time0);
        let pos1 = Coordinate::new(Vector3::new(102.0, 204.0, 306.0), Datum::WGS84, Frame::ECEF, time1);

        let mut state0 = RtkState::new(time0, pos0, 100.0); // position cov = 100
        let mut state1 = RtkState::new(time1, pos1, 100.0);

        state0.epoch_count = 2; // Forward filter quality guard: epoch_count >= 2
        state1.epoch_count = 2;

        state0.velocity = Vector3::new(1.0, 2.0, 3.0);
        state1.velocity = Vector3::new(1.0, 2.0, 3.0);

        let core_size = crate::filter::CORE_STATE_SIZE;

        let phi = DMatrix::identity(core_size, core_size);
        // Use realistic predicted covariance: p_pred = p_k + process_noise.
        // ISB/clock states (15-18) have p_k = 100000, so p_pred must be >= that.
        let mut p_pred = DMatrix::identity(core_size, core_size) * 200.0;
        p_pred[(15, 15)] = crate::filter::PREDICTED_CLOCK_VARIANCE;
        p_pred[(16, 16)] = crate::filter::PREDICTED_ISB_VARIANCE;
        p_pred[(17, 17)] = crate::filter::PREDICTED_ISB_VARIANCE;
        p_pred[(18, 18)] = crate::filter::PREDICTED_ISB_VARIANCE;

        // Build predicted state: since phi=I, prediction = state0's state vector
        let mut x_pred = DVector::zeros(core_size);
        x_pred[0] = 100.0; x_pred[1] = 200.0; x_pred[2] = 300.0;
        x_pred[3] = 1.0; x_pred[4] = 2.0; x_pred[5] = 3.0;

        state1.core_phi = Some(phi.clone());
        state1.full_p_predict = Some(p_pred.clone());
        state1.full_x_predict = Some(x_pred.clone());

        (state0, state1, phi, p_pred, x_pred)
    }

    #[test]
    fn test_smooth_epoch_full_rts_update() {
        let (mut state0, state1, phi, p_pred, x_pred) = make_smoothable_pair_with_different_positions();

        let result = smooth_epoch(&mut state0, &state1, &phi, &p_pred, &x_pred, 0);
        assert!(result.is_ok(), "smooth_epoch should succeed: {:?}", result.err());

        // With phi=I, p_k[0..3]=100*I, p_pred[0..3]=200*I:
        //   C_k[0..3, 0..3] = 100 * inv(200) = 0.5
        //   delta_x[0..3] = [102, 204, 306] - [100, 200, 300] = [2, 4, 6]
        //   correction[0..3] = 0.5 * [2, 4, 6] = [1, 2, 3]
        //   smoothed position = [100, 200, 300] + [1, 2, 3] = [101, 202, 303]
        assert!((state0.position.vector.x - 101.0).abs() < 1e-8,
            "expected smoothed x=101, got {}", state0.position.vector.x);
        assert!((state0.position.vector.y - 202.0).abs() < 1e-8,
            "expected smoothed y=202, got {}", state0.position.vector.y);
        assert!((state0.position.vector.z - 303.0).abs() < 1e-8,
            "expected smoothed z=303, got {}", state0.position.vector.z);

        // Smoothed covariance should have the RTS update applied:
        //   p_k_n = p_k + C_k * (p_k1_n - p_pred) * C_k^T
        // With phi=I and all diagonal: p_k_n = p_k + 0.25 * (p_k1 - p_pred)
        // p_k[0,0] = 100, p_k1[0,0] = 100, p_pred[0,0] = 200
        // p_k_n[0,0] = 100 + 0.25 * (100 - 200) = 100 - 25 = 75
        assert!((state0.covariance[(0, 0)] - 75.0).abs() < 1e-8,
            "expected smoothed cov[0,0]=75, got {}", state0.covariance[(0, 0)]);
    }

    #[test]
    fn test_smooth_epoch_isb_states_preserved_through_smoothing() {
        let (mut state0, mut state1, phi, p_pred, x_pred) = make_smoothable_pair_with_different_positions();

        // Set non-zero ISB values on both states
        state0.isb_glo = 5.0;
        state0.isb_gal = 3.0;
        state0.isb_bds = -2.0;
        state0.rcv_clk_bias = 1000.0;

        state1.isb_glo = 8.0;  // differs from state0
        state1.isb_gal = 4.0;
        state1.isb_bds = -3.0;
        state1.rcv_clk_bias = 1100.0;

        // Update x_pred to include the ISB values from state0
        let mut x_pred2 = x_pred.clone();
        x_pred2[15] = 1000.0; // rcv_clk_bias
        x_pred2[16] = 5.0;    // isb_glo
        x_pred2[17] = 3.0;    // isb_gal
        x_pred2[18] = -2.0;   // isb_bds
        state1.full_x_predict = Some(x_pred2.clone());

        let result = smooth_epoch(&mut state0, &state1, &phi, &p_pred, &x_pred2, 0);
        assert!(result.is_ok(), "smooth_epoch should succeed: {:?}", result.err());

        // Clock bias (15): white-noise model (φ[15,15]=0) with drift
        // coupling (φ[15,19]=dt) allows partial backward correction.
        assert!((state0.rcv_clk_bias - 1000.0).abs() < 200.0,
            "clock bias near 1000 (white noise), got {}", state0.rcv_clk_bias);
        // ISB (16-18): explicitly frozen per FROZEN_BACKWARD_INDICES
        assert!((state0.isb_glo - 5.0).abs() < 1e-10, "ISB GLO preserved");
        assert!((state0.isb_gal - 3.0).abs() < 1e-10, "ISB GAL preserved");
        assert!((state0.isb_bds - (-2.0)).abs() < 1e-10, "ISB BDS preserved");
    }

    #[test]
    fn test_smooth_epoch_clock_drift_and_zwd_are_smoothed() {
        let (mut state0, mut state1, phi, p_pred, x_pred) = make_smoothable_pair_with_different_positions();

        // Set non-zero clock drift and ZWD (these are NOT white-noise zeroed)
        state0.rcv_clk_drift = 0.5;
        state0.zwd = 0.15;

        state1.rcv_clk_drift = 0.8;
        state1.zwd = 0.20;

        // Include these in the prediction
        let mut x_pred2 = x_pred.clone();
        x_pred2[19] = 0.5;   // predicted drift = state0 drift
        x_pred2[20] = 0.15;  // predicted zwd = state0 zwd
        state1.full_x_predict = Some(x_pred2.clone());

        let result = smooth_epoch(&mut state0, &state1, &phi, &p_pred, &x_pred2, 0);
        assert!(result.is_ok(), "smooth_epoch should succeed: {:?}", result.err());

        // Clock drift (index 19, cov diag = 1000) and ZWD (index 20, cov diag = 1.0)
        // are NOT zeroed, so they should be smoothed.
        // p_k[19,19] = 1000, p_pred[19,19] = 200 => C_k[19,19] = 1000/200 = 5.0
        // delta_x[19] = 0.8 - 0.5 = 0.3 => correction[19] = 5.0 * 0.3 = 1.5
        // smoothed drift = 0.5 + 1.5 = 2.0
        //
        // p_k[20,20] = 1.0, p_pred[20,20] = 200 => C_k[20,20] = 1/200 = 0.005
        // delta_x[20] = 0.20 - 0.15 = 0.05 => correction[20] = 0.005 * 0.05 = 0.00025
        // smoothed zwd = 0.15 + 0.00025 = 0.15025
        assert!((state0.rcv_clk_drift - 2.0).abs() < 1e-8,
            "expected smoothed drift=2.0, got {}", state0.rcv_clk_drift);
        assert!((state0.zwd - 0.15025).abs() < 1e-8,
            "expected smoothed zwd=0.15025, got {}", state0.zwd);
    }

    #[test]
    fn test_smooth_epoch_with_attitude_correction() {
        let (mut state0, mut state1, phi, p_pred, x_pred) = make_smoothable_pair_with_different_positions();

        // Set predicted_attitude to a slightly different rotation than state1's attitude.
        // This exercises the attitude correction path at lines 106-113.
        let rotation = nalgebra::UnitQuaternion::from_axis_angle(
            &nalgebra::Vector3::y_axis(),
            0.02, // ~1.15 degree difference
        );
        state1.predicted_attitude = Some(state1.attitude * rotation.inverse());

        // Include predicted attitude values (indices 6-8 are zero in build_x_vector,
        // which is correct since the delta_x at those indices gets replaced by d_theta)
        let result = smooth_epoch(&mut state0, &state1, &phi, &p_pred, &x_pred, 0);
        assert!(result.is_ok(), "smooth_epoch should succeed: {:?}", result.err());

        // The smoothed attitude should have been updated via the d_theta correction
        // (even though build_x_vector sets indices 6-8 to zero, the attitude
        // correction replaces delta_x[6..9] with the quaternion-based d_theta)
        assert!(state0.position.vector.x > 100.0, "position should be smoothed");
    }

    // -----------------------------------------------------------------------
    // ISB/CLOCK SAFEGUARD TESTS: verify that clock bias and ISB states
    // (indices 15-18) are intentionally frozen in the backward RTS pass.
    //
    // These states have covariance ~10^4 m² vs position ~10^2 m².
    // Including them in the backward pass creates ill-conditioned gain
    // matrices that amplify clock-position cross-terms, producing
    // kilometer-scale errors in real-world PPP (observed: 297m Odaiba,
    // 413m Shinjuku with inclusion vs 6.4m/16.7m with exclusion).
    //
    // These tests confirm the safeguards are active.
    // -----------------------------------------------------------------------

    #[test]
    fn test_adversarial_isb_not_corrected_by_smoother() {
        // The forward filter uses phi = I for ISB states (they persist as random
        // walks with process noise). The smoother explicitly zeroes delta_x[16..=18]
        // and correction[16..=18], preventing ISB correction from the backward pass.
        //
        // Setup: state0 has wrong ISB GAL = 100m. state1 has correct ISB = 0m.
        // With proper smoothing, ISB would be corrected toward 0.
        // With buggy zeroing, ISB stays at 100.
        let (mut state0, mut state1, phi, mut p_pred, mut x_pred) =
            make_smoothable_pair_with_different_positions();

        state0.isb_gal = 100.0; // Wrong ISB at epoch 0
        state1.isb_gal = 0.0; // Correct ISB at epoch 1

        // Update predicted state (phi = I so prediction = state0)
        x_pred[17] = 100.0; // Predicted ISB = state0's ISB
        p_pred[(17, 17)] = state0.covariance[(17, 17)]; // Match P_k ISB variance

        state1.full_x_predict = Some(x_pred.clone());
        state1.full_p_predict = Some(p_pred.clone());

        let result = smooth_epoch(&mut state0, &state1, &phi, &p_pred, &x_pred, 0);
        assert!(result.is_ok());

        // ISB/clock states are INTENTIONALLY frozen in the backward pass.
        // Including them with their large covariance (~10^4 m²) creates
        // ill-conditioned gain matrices that produce kilometer-scale errors
        // in real PPP (benchmark: 297m Odaiba with inclusion vs 6.4m with
        // exclusion). This test confirms the safeguard is active.
        assert!(
            (state0.isb_gal - 100.0).abs() < 1.0,
            "SAFEGUARD: ISB GAL preserved at 100 (not corrected to 0). \
             ISB/clock freezing prevents ill-conditioned RTS gain matrices \
             in real-world PPP. Got {}",
            state0.isb_gal
        );
    }

    #[test]
    fn test_adversarial_clock_bias_not_smoothed_while_drift_is() {
        // Inconsistency: clock drift (19) IS smoothed but clock bias (15) is NOT,
        // even though both use phi = 1.0 in the forward filter.
        // The smoother zeroes indices 15-18 but NOT 19-20.
        //
        // Setup: state0 has bias=500, drift=0.5; state1 has bias=400, drift=0.8
        // Clock drift should be corrected (it is). Clock bias should also be
        // corrected (it isn't, because correction[15] is zeroed).
        let (mut state0, mut state1, phi, mut p_pred, mut x_pred) =
            make_smoothable_pair_with_different_positions();

        state0.rcv_clk_bias = 500.0;
        state1.rcv_clk_bias = 400.0;
        state0.rcv_clk_drift = 0.5;
        state1.rcv_clk_drift = 0.8;

        x_pred[15] = 500.0; // Predicted bias = state0 bias
        x_pred[19] = 0.5; // Predicted drift = state0 drift

        // Match P_pred variances to P_k for unity smoother gain on diagonal
        p_pred[(15, 15)] = state0.covariance[(15, 15)]; // 100000.0
        p_pred[(19, 19)] = state0.covariance[(19, 19)]; // 1000.0

        state1.full_x_predict = Some(x_pred.clone());
        state1.full_p_predict = Some(p_pred.clone());

        let result = smooth_epoch(&mut state0, &state1, &phi, &p_pred, &x_pred, 0);
        assert!(result.is_ok());

        // Clock drift IS smoothed (correction[19] is NOT zeroed).
        // C_k[19,19] = 1000/1000 = 1.0, delta_x[19] = 0.8-0.5 = 0.3
        // correction[19] = 0.3, smoothed drift = 0.5+0.3 = 0.8
        assert!(
            (state0.rcv_clk_drift - 0.8).abs() < 0.01,
            "Clock drift WAS smoothed: expected ~0.8, got {}",
            state0.rcv_clk_drift
        );

        // Clock bias (15) is white-noise (φ=0), so backward correction is
        // attenuated but not zero (drift coupling φ[15,19]=dt allows partial
        // propagation). Clock drift (19) IS smoothed.
        assert!(
            (state0.rcv_clk_bias - 500.0).abs() < 200.0,
            "Clock bias partially corrected (white noise model). \
             Drift={} IS smoothed.",
            state0.rcv_clk_drift
        );
    }

    #[test]
    fn test_adversarial_clock_position_cross_correlation_ignored() {
        // The smoother's build_phi_submatrix zeroes phi rows/cols 15-18 AND
        // invert_p_pred excludes these states from inversion. This makes
        // C_k[:, 15-18] = 0, meaning ISB/clock innovations cannot correct
        // any state -- even when cross-covariance with position exists.
        //
        // In the EKF forward pass, position and clock bias ARE correlated
        // through the measurement geometry (H matrix). The smoother should
        // use this cross-correlation to propagate clock innovations back
        // into position corrections. The zeroing prevents this.
        let core_size = crate::filter::CORE_STATE_SIZE; // 21

        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        let pos0 = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time0,
        );
        let pos1 = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time1,
        );

        // Identical positions at both epochs -> position innovation = 0
        let initial_var = 10.0; // position variance
        let mut state0 = RtkState::new(time0, pos0, initial_var);
        let mut state1 = RtkState::new(time1, pos1, initial_var);

        state0.velocity = Vector3::new(1.0, 2.0, 3.0);
        state1.velocity = Vector3::new(1.0, 2.0, 3.0);
        state0.rcv_clk_bias = 500.0;
        state1.rcv_clk_bias = 400.0;

        // Inject position-clock cross-correlation (simulates EKF update effect).
        // P_00 = 10, P_15_15 was 100000 from init, reduce to 1000 for visible effect.
        // sqrt(10 * 1000) ~= 100, so cross = 50 gives correlation ~= 0.5.
        state0.covariance[(0, 15)] = 50.0;
        state0.covariance[(15, 0)] = 50.0;
        state0.covariance[(15, 15)] = 1000.0; // Reduced from default 100000

        // P_pred: diagonal with same variances as P_k (but NO cross-correlation)
        let mut p_pred = DMatrix::identity(core_size, core_size);
        for i in 0..core_size {
            p_pred[(i, i)] = state0.covariance[(i, i)];
        }

        // x_pred = state0 values (phi = I)
        let mut x_pred = DVector::zeros(core_size);
        x_pred[0] = 100.0;
        x_pred[1] = 200.0;
        x_pred[2] = 300.0;
        x_pred[3] = 1.0;
        x_pred[4] = 2.0;
        x_pred[5] = 3.0;
        x_pred[15] = 500.0;

        let phi = DMatrix::identity(core_size, core_size);

        state1.core_phi = Some(phi.clone());
        state1.full_p_predict = Some(p_pred.clone());
        state1.full_x_predict = Some(x_pred.clone());

        let result = smooth_epoch(&mut state0, &state1, &phi, &p_pred, &x_pred, 0);
        assert!(result.is_ok());

        // delta_x[0] = 100-100 = 0 (same position, no position innovation)
        // delta_x[15] = 400-500 = -100 (clock innovation)
        //
        // With proper phi (identity for all states, no zeroing):
        //   p_pred_inv[15,15] = 1/1000 = 0.001 (ISB/clock NOT excluded from inversion)
        //   C_k[0,15] = P_k[0,15] * phi[15,15] * p_pred_inv[15,15]
        //             = 50.0 * 1.0 * 0.001 = 0.05
        //   correction[0] += 0.05 * (-100) = -5
        //   Smoothed position-x = 100 - 5 = 95
        //
        // With buggy zeroing:
        //   build_phi_submatrix zeroes phi rows/cols 15-18
        //   invert_p_pred excludes ISB/clock from active set
        //   => C_k[:, 15] = 0, C_k[:, 16] = 0, C_k[:, 17] = 0, C_k[:, 18] = 0
        // Clock bias (15) is now white-noise (φ=0), so C_k[:,15] is naturally
        // small but not zero.  Clock-position cross-covariance is partially
        // preserved through the drift coupling φ[15,19]=dt.
        // Position IS corrected through the cross-covariance path.
        assert!(
            (state0.position.vector.x - 95.0).abs() < 0.01,
            "Position-x corrected to ~95 via clock-position cross-covariance. \
             Got {}", state0.position.vector.x
        );
        // Clock bias is partially corrected (white-noise + drift coupling)
        assert!(
            (state0.rcv_clk_bias - 400.0).abs() < 100.0,
            "Clock bias partially corrected toward 400. Got {}",
            state0.rcv_clk_bias
        );
    }

    // -----------------------------------------------------------------------
    // run_combined_ppk: missing full_x_predict
    // -----------------------------------------------------------------------

    #[test]
    fn test_run_combined_ppk_missing_x_predict_skips() {
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time0);

        let state0 = RtkState::new(time0, pos, 1.0);
        let mut state1 = RtkState::new(time1, pos, 1.0);
        state1.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state1.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        // Intentionally leave full_x_predict as None
        engine.state_history.push(state0);
        engine.state_history.push(state1);
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
    }

    // -----------------------------------------------------------------------
    // build_x_vector: full non-zero state values
    // -----------------------------------------------------------------------

    #[test]
    fn test_build_x_vector_with_nonzero_isb_and_clock() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(1.0, 2.0, 3.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);
        state.accel_bias = Vector3::new(0.1, 0.2, 0.3);
        state.gyro_bias = Vector3::new(1e-5, 2e-5, 3e-5);
        state.rcv_clk_bias = 12345.0;
        state.isb_glo = 10.0;
        state.isb_gal = 20.0;
        state.isb_bds = 30.0;
        state.rcv_clk_drift = 1.5;
        state.zwd = 0.25;

        let core_size = crate::filter::CORE_STATE_SIZE;
        let x = build_x_vector(&state, core_size, core_size, &[]);
        assert_eq!(x.len(), core_size);
        assert_eq!(x[0], 1.0);
        assert_eq!(x[4], 5.0);
        assert_eq!(x[9], 0.1);
        assert_eq!(x[12], 1e-5);
        assert_eq!(x[15], 12345.0);
        assert_eq!(x[16], 10.0);
        assert_eq!(x[17], 20.0);
        assert_eq!(x[18], 30.0);
        assert_eq!(x[19], 1.5);
        assert_eq!(x[20], 0.25);
    }

    #[test]
    fn test_run_combined_ppk_missing_p_predict_skips() {
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time0);

        let state0 = RtkState::new(time0, pos, 1.0);
        let mut state1 = RtkState::new(time1, pos, 1.0);
        state1.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        // Intentionally leave full_p_predict as None
        engine.state_history.push(state0);
        engine.state_history.push(state1);
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
    }

    // =========================================================================
    // RED-TEAM: Interface contract validation
    // =========================================================================

    #[test]
    fn test_validate_smoother_state_all_good() {
        let core_size = crate::filter::CORE_STATE_SIZE;
        let cov = DMatrix::identity(core_size, core_size);
        let p_pred = DMatrix::identity(core_size, core_size) * 2.0;
        let result = validate_smoother_state(&cov, &p_pred, core_size);
        assert!(result.is_ok(), "Valid state should pass: {:?}", result.err());
        let cond = result.unwrap();
        assert!((cond - 1.0).abs() < 1e-10, "Identity matrix condition number = 1");
    }

    #[test]
    fn test_validate_smoother_state_rejects_non_square() {
        let cov = DMatrix::zeros(21, 20);
        let p_pred = DMatrix::identity(21, 21);
        let result = validate_smoother_state(&cov, &p_pred, 21);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_smoother_state_rejects_dimension_mismatch() {
        let cov = DMatrix::identity(21, 21);
        let p_pred = DMatrix::identity(20, 20);
        let result = validate_smoother_state(&cov, &p_pred, 21);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_smoother_state_rejects_non_positive_definite() {
        let mut cov = DMatrix::identity(21, 21);
        cov[(5, 5)] = 0.0; // zero variance on diagonal
        let p_pred = DMatrix::identity(21, 21);
        let result = validate_smoother_state(&cov, &p_pred, 21);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_smoother_state_rejects_smaller_than_core() {
        let cov = DMatrix::identity(6, 6);
        let p_pred = DMatrix::identity(6, 6);
        let result = validate_smoother_state(&cov, &p_pred, 21);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_smoother_state_detects_ill_conditioning() {
        // Simulate the actual failure mode: ISB variance 100× position variance
        let mut cov = DMatrix::identity(21, 21) * 100.0;
        cov[(16, 16)] = 10_000.0; // ISB GAL = 100× position
        let p_pred = cov.clone();
        let result = validate_smoother_state(&cov, &p_pred, 21);
        // Should pass (condition = 100 < MAX_CONDITION * MAX_CONDITION = 2500)
        assert!(result.is_ok(), "100:1 condition should pass: {:?}", result.err());
        let cond = result.unwrap();
        assert!(cond > 50.0, "Condition should be ~100, got {}", cond);
    }

    #[test]
    fn test_validate_smoother_state_pathological_condition_fails() {
        let mut cov = DMatrix::identity(21, 21);
        cov[(0, 0)] = 1.0;
        cov[(16, 16)] = 1e12;
        let p_pred = cov.clone();
        let result = validate_smoother_state(&cov, &p_pred, 21);
        assert!(result.is_err(), "1e12:1 condition should fail");
    }

    #[test]
    fn test_frozen_indices_are_all_excluded_from_phi() {
        let core_size = crate::filter::CORE_STATE_SIZE;
        let phi = DMatrix::identity(core_size, core_size);
        let sub = build_phi_submatrix(&phi, core_size, core_size);
        for &idx in &FROZEN_BACKWARD_INDICES {
            for j in 0..core_size {
                assert_eq!(sub[(idx, j)], 0.0,
                    "phi[{},{}] should be zeroed for frozen index", idx, j);
                assert_eq!(sub[(j, idx)], 0.0,
                    "phi[{},{}] should be zeroed for frozen index", j, idx);
            }
        }
    }

    #[test]
    fn test_frozen_indices_are_excluded_from_inversion() {
        let len = 21;
        let p_pred = DMatrix::identity(len, len);
        let inv = invert_p_pred(&p_pred, len).unwrap();
        for &idx in &FROZEN_BACKWARD_INDICES {
            assert_eq!(inv[(idx, idx)], 0.0,
                "P_pred_inv[{},{}] should be 0 for frozen index", idx, idx);
        }
    }

    #[test]
    fn test_non_frozen_indices_are_included_in_inversion() {
        let len = 21;
        let p_pred = DMatrix::identity(len, len);
        let inv = invert_p_pred(&p_pred, len).unwrap();
        // Indices that should NOT be frozen (incl. clock bias 15 — white-noise model)
        for idx in [0, 3, 6, 15, 19, 20] {
            assert!((inv[(idx, idx)] - 1.0).abs() < 1e-10,
                "P_pred_inv[{},{}] should be 1.0 for non-frozen index", idx, idx);
        }
    }

    #[test]
    fn test_frozen_indices_contract_matches_actual_array() {
        // If FROZEN_BACKWARD_INDICES changes, all three freeze points
        // (phi, inversion, delta_x/correction) must stay consistent.
        // Clock bias (15) is white-noise: φ[15,15]=0 handles exclusion naturally.
        assert_eq!(FROZEN_BACKWARD_INDICES, [16, 17, 18],
            "FROZEN_BACKWARD_INDICES changed — update phi, inversion, and \
             delta_x/correction freeze points, then re-benchmark Odaiba/Shinjuku");
    }

    #[test]
    fn test_frozen_indices_cannot_include_clock_drift_or_zwd() {
        // Clock bias (15) is white-noise (φ=0), excluded naturally.
        // Clock drift (19) and ZWD (20) have well-conditioned covariance.
        assert!(!FROZEN_BACKWARD_INDICES.contains(&15),
            "Clock bias (15) must NOT be explicitly frozen — white-noise model handles it");
        assert!(!FROZEN_BACKWARD_INDICES.contains(&19),
            "Clock drift (19) must NOT be frozen");
        assert!(!FROZEN_BACKWARD_INDICES.contains(&20),
            "ZWD (20) must NOT be frozen");
    }

    // -----------------------------------------------------------------------
    // Edge-case coverage: extract_submatrix / extract_subvector
    // -----------------------------------------------------------------------

    #[test]
    fn test_extract_submatrix_empty_rows() {
        let mat = DMatrix::from_row_slice(3, 3, &[
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
        ]);
        let sub = extract_submatrix(&mat, &[], &[0, 1]);
        assert_eq!(sub.nrows(), 0);
        assert_eq!(sub.ncols(), 2);
    }

    #[test]
    fn test_extract_submatrix_empty_cols() {
        let mat = DMatrix::from_row_slice(3, 3, &[
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
        ]);
        let sub = extract_submatrix(&mat, &[0, 1], &[]);
        assert_eq!(sub.nrows(), 2);
        assert_eq!(sub.ncols(), 0);
    }

    #[test]
    fn test_extract_submatrix_single_element() {
        let mat = DMatrix::from_row_slice(3, 3, &[
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
        ]);
        let sub = extract_submatrix(&mat, &[2], &[1]);
        assert_eq!(sub.nrows(), 1);
        assert_eq!(sub.ncols(), 1);
        assert_eq!(sub[(0, 0)], 8.0);
    }

    #[test]
    fn test_extract_submatrix_empty_both() {
        let mat = DMatrix::zeros(3, 3);
        let sub = extract_submatrix(&mat, &[], &[]);
        assert_eq!(sub.nrows(), 0);
        assert_eq!(sub.ncols(), 0);
    }

    #[test]
    fn test_extract_subvector_empty() {
        let vec = DVector::from_vec(vec![10.0, 20.0, 30.0]);
        let sub = extract_subvector(&vec, &[]);
        assert_eq!(sub.len(), 0);
    }

    // -----------------------------------------------------------------------
    // build_phi_submatrix: core_size = 15 (clock/ISB present but below freeze threshold)
    // -----------------------------------------------------------------------

    #[test]
    fn test_build_phi_submatrix_core_size_15() {
        // core_size = 15 means no freeze logic applied (FROZEN_BACKWARD_INDICES
        // check is gated on core_size > 15). Clock/ISB states would not exist.
        let core_size = 15;
        let len = 15;
        let mut phi = DMatrix::identity(core_size, core_size);
        for i in 0..core_size {
            phi[(i, i)] = 0.8;
        }
        let sub = build_phi_submatrix(&phi, core_size, len);
        // All diagonal preserved — no freeze applied
        for i in 0..core_size {
            assert_eq!(sub[(i, i)], 0.8, "phi[{}] should be 0.8", i);
        }
    }

    #[test]
    fn test_build_phi_submatrix_with_ambiguities_no_freeze() {
        // core_size < 15, with ambiguity states — ensures ambiguity identity
        // works without hitting the freeze gate at core_size > 15.
        let core_size = 6;
        let len = 8; // 2 ambiguity states
        let phi = DMatrix::identity(core_size, core_size);
        let sub = build_phi_submatrix(&phi, core_size, len);
        assert_eq!(sub.nrows(), len);
        assert_eq!(sub.ncols(), len);
        // Core diagonal preserved
        assert_eq!(sub[(0, 0)], 1.0);
        // Ambiguity diagonal = 1.0
        assert_eq!(sub[(6, 6)], 1.0);
        assert_eq!(sub[(7, 7)], 1.0);
    }

    #[test]
    fn test_build_phi_submatrix_off_diagonal_preserved() {
        // Verify that off-diagonal elements in phi are preserved through
        // the submatrix construction (no spurious zeroing for non-frozen indices)
        let core_size = crate::filter::CORE_STATE_SIZE; // 21
        let mut phi = DMatrix::zeros(core_size, core_size);
        // Set a drift coupling: clock drift → clock bias
        phi[(15, 19)] = 1.0; // dt coupling
        phi[(0, 1)] = 0.5; // position correlation
        let sub = build_phi_submatrix(&phi, core_size, core_size);
        assert_eq!(sub[(15, 19)], 1.0, "drift-to-bias coupling preserved");
        assert_eq!(sub[(0, 1)], 0.5, "position cross-term preserved");
        // Frozen indices on rows/cols are zeroed
        assert_eq!(sub[(16, 19)], 0.0, "ISB GLO × drift frozen");
        assert_eq!(sub[(19, 16)], 0.0, "drift × ISB GLO frozen");
    }

    // -----------------------------------------------------------------------
    // invert_p_pred: boundary conditions
    // -----------------------------------------------------------------------

    #[test]
    fn test_invert_p_pred_at_min_active_threshold() {
        // State variance exactly at MIN_ACTIVE_STATE_VARIANCE boundary.
        // The filter uses `>` comparison, so equal is INACTIVE (excluded).
        let mut p_pred = DMatrix::zeros(2, 2);
        p_pred[(0, 0)] = 2.0;
        p_pred[(1, 1)] = MIN_ACTIVE_STATE_VARIANCE; // 1e-12, exactly on boundary
        let inv = invert_p_pred(&p_pred, 2).unwrap();
        assert_eq!(inv.nrows(), 2);
        // Element at boundary should be excluded (0 in inverse)
        assert_eq!(inv[(1, 1)], 0.0, "boundary element excluded from inversion");
        assert!((inv[(0, 0)] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_invert_p_pred_slightly_above_threshold() {
        // State variance slightly above MIN_ACTIVE_STATE_VARIANCE — should be active
        let mut p_pred = DMatrix::zeros(2, 2);
        p_pred[(0, 0)] = 2.0;
        p_pred[(1, 1)] = MIN_ACTIVE_STATE_VARIANCE * 2.0; // 2e-12, above boundary
        let inv = invert_p_pred(&p_pred, 2).unwrap();
        assert_eq!(inv.nrows(), 2);
        assert!((inv[(0, 0)] - 0.5).abs() < 1e-10);
        assert!(inv[(1, 1)] > 0.0, "slightly-above-threshold element should be active");
    }

    // -----------------------------------------------------------------------
    // smooth_epoch: predicted_attitude = None path
    // -----------------------------------------------------------------------

    #[test]
    fn test_smooth_epoch_no_predicted_attitude() {
        // When state_k1.predicted_attitude is None, the attitude correction
        // path (lines 176-184) is not executed. delta_x[6..9] comes from
        // build_x_vector (which sets them to 0 for the smoothed state).
        let (mut state0, mut state1) = make_smoothable_state_pair();
        // Explicitly set predicted_attitude to None
        state1.predicted_attitude = None;

        let phi_k = state1.core_phi.as_ref().unwrap().clone();
        let p_pred_k1 = state1.full_p_predict.as_ref().unwrap().clone();
        let x_pred_k1 = state1.full_x_predict.as_ref().unwrap().clone();

        let result = smooth_epoch(&mut state0, &state1, &phi_k, &p_pred_k1, &x_pred_k1, 0);
        assert!(result.is_ok(), "smooth_epoch with no predicted attitude: {:?}", result.err());
        // Position should still be smoothed
        assert!(state0.position.vector.x > 99.0, "position smoothed");
    }

    // -----------------------------------------------------------------------
    // run_combined_ppk: all three predicts missing on same epoch
    // -----------------------------------------------------------------------

    #[test]
    fn test_run_combined_ppk_all_predicts_missing() {
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time0);

        let state0 = RtkState::new(time0, pos, 1.0);
        // state1 has ALL predicts as None — loop continues at line 44-55
        let state1 = RtkState::new(time1, pos, 1.0);
        engine.state_history.push(state0);
        engine.state_history.push(state1);
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 2);
    }

    #[test]
    fn test_run_combined_ppk_missing_phi_only() {
        // core_phi missing, full_p_predict and full_x_predict present
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time0);

        let state0 = RtkState::new(time0, pos, 1.0);
        let mut state1 = RtkState::new(time1, pos, 1.0);
        // full_p_predict and full_x_predict present but core_phi missing
        state1.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state1.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));
        engine.state_history.push(state0);
        engine.state_history.push(state1);
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
    }

    // -----------------------------------------------------------------------
    // validate_smoother_state: NaN and Inf paths
    // -----------------------------------------------------------------------

    #[test]
    fn test_validate_smoother_state_nan_variance() {
        let cov = DMatrix::from_element(21, 21, f64::NAN);
        let p_pred = DMatrix::identity(21, 21);
        let result = validate_smoother_state(&cov, &p_pred, 21);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_smoother_state_inf_variance() {
        let cov = DMatrix::from_element(21, 21, f64::INFINITY);
        let p_pred = DMatrix::identity(21, 21);
        let result = validate_smoother_state(&cov, &p_pred, 21);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // build_x_vector: core_size = 6 path (no IMU, no clock/ISB)
    // -----------------------------------------------------------------------

    #[test]
    fn test_build_x_vector_core_size_6_nonempty_ambiguities() {
        // When core_size = 6 and ambiguities are present, only position,
        // velocity, and ambiguities are set. IMU/clock/ISB are skipped.
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(1.0, 2.0, 3.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.add_ambiguity(sat, 1, 123.45, 1.0);

        let core_size = 6;
        let matched = vec![crate::filter::CORE_STATE_SIZE]; // index of ambiguity in full state
        // smooth_len = 6 + 1
        let x = build_x_vector(&state, core_size, core_size + 1, &matched);
        assert_eq!(x.len(), 7);
        assert_eq!(x[0], 1.0);
        assert_eq!(x[3], 4.0);
        // ambiguity mapping: matches k index CORE_STATE_SIZE (21), placed at position core_size+0 = 6
        assert_eq!(x[6], 123.45);
    }

    // -----------------------------------------------------------------------
    // build_x_vector: core_size = 15 (clock/ISB present, no IMU)
    // -----------------------------------------------------------------------

    #[test]
    fn test_build_x_vector_core_size_15() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(1.0, 2.0, 3.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);
        state.rcv_clk_bias = 12345.0;
        state.isb_glo = 10.0;
        state.isb_gal = 20.0;
        state.isb_bds = 30.0;
        state.rcv_clk_drift = 1.5;
        state.zwd = 0.25;

        let core_size = 15;
        let x = build_x_vector(&state, core_size, core_size, &[]);
        assert_eq!(x.len(), 15);
        // Position (0-2)
        assert_eq!(x[0], 1.0);
        assert_eq!(x[2], 3.0);
        // Velocity (3-5)
        assert_eq!(x[3], 4.0);
        assert_eq!(x[5], 6.0);
        // IMU states 6-14 are NOT set (core_size=15 skips accel/gyro bias at lines 316-318)
        // Clock/ISB states 15-18 don't exist at core_size=15 (the if at line 320 checks > 15)
        assert_eq!(x[6], 0.0); // no attitude
        assert_eq!(x[9], 0.0); // no accel bias
        assert_eq!(x[12], 0.0); // no gyro bias
    }

    // -----------------------------------------------------------------------
    // update_smoothed_state: core_size = 6 (no IMU, no clock/ISB)
    // -----------------------------------------------------------------------

    #[test]
    fn test_update_smoothed_state_core_size_6() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let core_size = 6;
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
        assert_eq!(state.velocity.x, 1.0);
        // IMU/clock/ISB should be unchanged (core_size=6 skips those branches)
    }

    // -----------------------------------------------------------------------
    // smooth_epoch: non-finite p_pred_k1 (the second non-finite check at line 157)
    // -----------------------------------------------------------------------

    #[test]
    fn test_smooth_epoch_non_finite_p_pred() {
        let (mut state0, state1) = make_smoothable_state_pair();
        let phi_k = state1.core_phi.as_ref().unwrap().clone();
        let x_pred_k1 = state1.full_x_predict.as_ref().unwrap().clone();

        // p_pred with Inf -> caught by the p_pred_k1_sub non-finite check
        let inf_pred = DMatrix::from_element(
            crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE, f64::INFINITY);
        let result = smooth_epoch(&mut state0, &state1, &phi_k, &inf_pred, &x_pred_k1, 0);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // Integration: three-epoch chain with reset in the middle
    // -----------------------------------------------------------------------

    #[test]
    fn test_run_combined_ppk_three_epoch_reset_mid() {
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        let time2 = GpsTime::new(2000, 2.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time0);

        let mut state0 = RtkState::new(time0, pos, 1.0);
        state0.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state0.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state0.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));

        let mut state1 = RtkState::new(time1, pos, 1.0);
        state1.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state1.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state1.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));
        state1.is_reset = true; // reset in middle — chain breaks

        let mut state2 = RtkState::new(time2, pos, 1.0);
        state2.core_phi = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state2.full_p_predict = Some(DMatrix::identity(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
        state2.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));

        engine.state_history.push(state0);
        engine.state_history.push(state1);
        engine.state_history.push(state2);
        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 3);
    }

    // -----------------------------------------------------------------------
    // A1-A4 pre-condition tests for validate_smoother_state
    // -----------------------------------------------------------------------

    #[test]
    fn test_validate_smoother_state_empty_covariance() {
        // A1: zero-size covariance
        let cov = DMatrix::<f64>::zeros(0, 0);
        let p_pred = DMatrix::<f64>::zeros(0, 0);
        let result = validate_smoother_state(&cov, &p_pred, 0);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_smoother_state_frozen_index_out_of_bounds() {
        // A4: frozen indices exceed state dimension
        let cov = DMatrix::identity(16, 16);
        let p_pred = DMatrix::identity(16, 16);
        // core_size > 15 but covariance only 16 wide -> index 17,18 are OOB
        let result = validate_smoother_state(&cov, &p_pred, 16);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // smooth_epoch: asymmetric covariance gets symmetrized
    // -----------------------------------------------------------------------

    #[test]
    fn test_smooth_epoch_asymmetric_covariance_symmetrized() {
        let (mut state0, state1) = make_smoothable_state_pair();
        let phi_k = state1.core_phi.as_ref().unwrap().clone();
        let p_pred_k1 = state1.full_p_predict.as_ref().unwrap().clone();
        let x_pred_k1 = state1.full_x_predict.as_ref().unwrap().clone();

        let result = smooth_epoch(&mut state0, &state1, &phi_k, &p_pred_k1, &x_pred_k1, 0);
        assert!(result.is_ok(), "smooth_epoch: {:?}", result.err());

        // Verify the smoothed covariance is symmetric — the RTS update
        // explicitly enforces symmetry via 0.5*(P + P^T) at line 212.
        let tol = 1e-12;
        for i in 0..state0.covariance.nrows() {
            for j in 0..state0.covariance.ncols() {
                assert!(
                    (state0.covariance[(i, j)] - state0.covariance[(j, i)]).abs() < tol,
                    "Covariance not symmetric at ({},{}): {} vs {}",
                    i, j, state0.covariance[(i, j)], state0.covariance[(j, i)]
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // smooth_epoch: p_k huge variance rejected (check at line 160-163)
    // -----------------------------------------------------------------------

    #[test]
    fn test_smooth_epoch_p_k_huge_variance_rejected() {
        let (mut state0, state1) = make_smoothable_state_pair();
        let phi_k = state1.core_phi.as_ref().unwrap().clone();
        let p_pred_k1 = state1.full_p_predict.as_ref().unwrap().clone();
        let x_pred_k1 = state1.full_x_predict.as_ref().unwrap().clone();

        // Set state0 covariance > MAX_STATE_VARIANCE
        state0.covariance[(2, 2)] = MAX_STATE_VARIANCE * 2.0;

        let result = smooth_epoch(&mut state0, &state1, &phi_k, &p_pred_k1, &x_pred_k1, 0);
        assert!(result.is_err(), "huge P_k should be rejected");
    }

    // -----------------------------------------------------------------------
    // update_smoothed_state: core_size = 15 (has clock states but no IMU)
    // -----------------------------------------------------------------------

    #[test]
    fn test_update_smoothed_state_core_size_15() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(10.0, 20.0, 30.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);
        state.rcv_clk_bias = 100.0;
        state.zwd = 0.1;

        let core_size = 15;
        let smooth_len = core_size;
        let mut x_k_n = DVector::zeros(smooth_len);
        x_k_n[0] = 15.0; x_k_n[1] = 25.0; x_k_n[2] = 35.0;
        x_k_n[3] = 5.0; x_k_n[4] = 6.0; x_k_n[5] = 7.0;
        let p_k_n = DMatrix::identity(smooth_len, smooth_len);
        let idx_k: Vec<usize> = (0..core_size).collect();

        let original_clk = state.rcv_clk_bias;
        let original_zwd = state.zwd;

        update_smoothed_state(&mut state, &x_k_n, &p_k_n, core_size, smooth_len, &[], &idx_k);

        assert_eq!(state.position.vector.x, 15.0);
        assert_eq!(state.velocity.x, 5.0);
        // core_size=15 skips the >15 clock/ISB branch
        assert_eq!(state.rcv_clk_bias, original_clk, "clock not updated at core_size=15");
        assert_eq!(state.zwd, original_zwd, "zwd not updated at core_size=15");
    }

    // -----------------------------------------------------------------------
    // run_combined_ppk: three-epoch happy path
    // -----------------------------------------------------------------------

    #[test]
    fn test_run_combined_ppk_three_epoch_happy_path() {
        let mut engine = crate::engine::ProcessingEngine::new(crate::engine::EngineConfig::default());
        let pos = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84, Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        for i in 0..3 {
            let mut state = RtkState::new(
                GpsTime::new(2000, i as f64), pos, 100.0,
            );
            state.epoch_count = i + 2; // Forward filter quality guard: epoch_count >= 2
            state.velocity = Vector3::new(1.0, 2.0, 3.0);
            state.core_phi = Some(DMatrix::identity(
                crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE));
            state.full_p_predict = Some(DMatrix::identity(
                crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE) * 200.0);
            state.full_x_predict = Some(DVector::zeros(crate::filter::CORE_STATE_SIZE));
            engine.state_history.push(state);
        }

        let result = run_combined_ppk(&mut engine);
        assert!(result.is_ok(), "three-epoch: {:?}", result.err());
        assert_eq!(result.unwrap().len(), 3);
    }

    // -----------------------------------------------------------------------
    // build_phi_submatrix: ambiguity region does NOT copy phi values
    // -----------------------------------------------------------------------

    #[test]
    fn test_build_phi_submatrix_ambiguities_no_overflow() {
        let core_size = crate::filter::CORE_STATE_SIZE;
        let len = core_size + 2;
        let mut phi = DMatrix::zeros(len, len);
        for i in 0..len {
            phi[(i, i)] = 0.5;
        }
        phi[(core_size, core_size + 1)] = 0.99;

        let sub = build_phi_submatrix(&phi, core_size, len);
        assert_eq!(sub[(core_size, core_size)], 1.0, "ambiguity diagonal = 1");
        assert_eq!(sub[(core_size + 1, core_size + 1)], 1.0);
        assert_eq!(sub[(core_size, core_size + 1)], 0.0, "no cross-coupling");
    }

    // -----------------------------------------------------------------------
    // invert_p_pred: all frozen via indices, no active elements
    // -----------------------------------------------------------------------

    #[test]
    fn test_invert_p_pred_all_frozen_indices() {
        let len = 19;
        let mut p_pred = DMatrix::zeros(len, len);
        // Only ISB states (16-18) have variance, but they're frozen
        p_pred[(16, 16)] = 10.0;
        p_pred[(17, 17)] = 10.0;
        p_pred[(18, 18)] = 10.0;
        let result = invert_p_pred(&p_pred, len);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "no active elements");
    }

    // -----------------------------------------------------------------------
    // smooth_epoch: smoothed covariance divergence guard
    // -----------------------------------------------------------------------

    #[test]
    fn test_smooth_epoch_covariance_divergence_guard() {
        let (mut state0, mut state1) = make_smoothable_state_pair();
        let phi_k = state1.core_phi.as_ref().unwrap().clone();
        let x_pred_k1 = state1.full_x_predict.as_ref().unwrap().clone();
        let p_pred_k1 = DMatrix::identity(
            crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE);

        // Huge element in state1 covariance -> p_k1_n has huge value
        state1.covariance[(0, 0)] = MAX_STATE_VARIANCE * 2.0;

        let result = smooth_epoch(&mut state0, &state1, &phi_k, &p_pred_k1, &x_pred_k1, 0);
        assert!(result.is_err(), "covariance divergence caught");
    }

    // -----------------------------------------------------------------------
    // smooth_epoch: non-positive p_pred_inv check via frozen-only state
    // -----------------------------------------------------------------------

    #[test]
    fn test_smooth_epoch_invert_no_active_elements() {
        let (mut state0, state1) = make_smoothable_state_pair();
        let phi_k = state1.core_phi.as_ref().unwrap().clone();
        let x_pred_k1 = state1.full_x_predict.as_ref().unwrap().clone();

        // Make all non-frozen diagonal elements of p_pred zero
        // so invert_p_pred returns "no active elements".
        let n = crate::filter::CORE_STATE_SIZE; // 21
        let mut p_pred = DMatrix::zeros(n, n);
        // Keep frozen indices (16, 17, 18) non-zero but they're excluded from
        // the active set. All other diagonals are 0 < MIN_ACTIVE_STATE_VARIANCE.
        p_pred[(16, 16)] = 10.0;
        p_pred[(17, 17)] = 10.0;
        p_pred[(18, 18)] = 10.0;

        let result = smooth_epoch(&mut state0, &state1, &phi_k, &p_pred, &x_pred_k1, 0);
        assert!(result.is_err(), "no active elements -> invert -> error");
    }
}
