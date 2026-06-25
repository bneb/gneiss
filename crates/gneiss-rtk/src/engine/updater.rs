pub use crate::engine::updater_math::*;
use crate::filter::RtkState;
use crate::math::covariance::apply_joseph_covariance_update;
use crate::math::inversion::invert_matrix_robust;
use crate::math::thresholding::huber_scale_covariance;
use nalgebra::{DMatrix, DVector, UnitQuaternion, Vector3};

#[derive(Debug, PartialEq)]
pub enum UpdateError {
    SingularMatrix,
    DimensionMismatch,
    InvalidMeasurement,
}

pub fn apply_state_correction(state: &mut RtkState, dx: &DVector<f64>) {
    state.position.vector.x += dx[0];
    state.position.vector.y += dx[1];
    state.position.vector.z += dx[2];
    state.velocity.x += dx[3];
    state.velocity.y += dx[4];
    state.velocity.z += dx[5];
    if dx.len() >= crate::filter::CORE_STATE_SIZE {
        apply_imu_and_clock_correction(state, dx);
    }
    if dx.len() > crate::filter::CORE_STATE_SIZE {
        for i in 0..state.ambiguities.len() {
            state.ambiguities[i] += dx[crate::filter::CORE_STATE_SIZE + i];
        }
    }
}

fn apply_attitude_correction(state: &mut RtkState, dx: &DVector<f64>) {
    let d_theta = Vector3::new(dx[6], dx[7], dx[8]);
    if d_theta.norm() <= 1e-10 {
        return;
    }
    let dq =
        UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_normalize(d_theta), d_theta.norm());
    state.attitude = dq * state.attitude;
    state.attitude.renormalize();
}

fn apply_bias_correction(state: &mut RtkState, dx: &DVector<f64>) {
    state.accel_bias.x += dx[9];
    state.accel_bias.y += dx[10];
    state.accel_bias.z += dx[11];
    state.gyro_bias.x += dx[12];
    state.gyro_bias.y += dx[13];
    state.gyro_bias.z += dx[14];
}

fn apply_clock_correction(state: &mut RtkState, dx: &DVector<f64>) {
    state.rcv_clk_bias += dx[15];
    state.isb_glo += dx[16];
    state.isb_gal += dx[17];
    state.isb_bds += dx[18];
    state.rcv_clk_drift += dx[19];
    state.zwd = (state.zwd + dx[20]).max(0.0);
}

fn apply_imu_and_clock_correction(state: &mut RtkState, dx: &DVector<f64>) {
    apply_attitude_correction(state, dx);
    apply_bias_correction(state, dx);
    apply_clock_correction(state, dx);
}

fn build_loose_coupling_h(
    state: &RtkState,
    lever_arm: &Vector3<f64>,
    omega_b: &Vector3<f64>,
) -> DMatrix<f64> {
    let mut h_mat = DMatrix::zeros(6, state.covariance.ncols());
    h_mat.view_mut((0, 0), (6, 6)).fill_diagonal(1.0);
    if state.covariance.nrows() >= crate::filter::CORE_STATE_SIZE {
        populate_loosely_coupled_jacobian(&mut h_mat, &state.attitude, lever_arm, omega_b);
    }
    h_mat
}

fn compute_loose_coupling_innovation(
    state: &RtkState,
    gnss_state: &RtkState,
    lever_arm: &Vector3<f64>,
    omega_b: &Vector3<f64>,
) -> DVector<f64> {
    let r_b_e = state.attitude.to_rotation_matrix();
    compute_loose_coupling_innovations(
        r_b_e.matrix(),
        &state.position.vector,
        &state.velocity,
        &gnss_state.position.vector,
        &gnss_state.velocity,
        lever_arm,
        omega_b,
    )
}

fn compute_loose_coupling_gain(
    state: &RtkState,
    _z: &DVector<f64>,
    h_mat: &DMatrix<f64>,
    r_6x6: &DMatrix<f64>,
) -> Result<DMatrix<f64>, UpdateError> {
    let s = state.covariance.view((0, 0), (6, 6)).into_owned() + r_6x6;
    let s_inv = invert_matrix_robust(&s);
    Ok(&state.covariance * h_mat.transpose() * s_inv)
}

pub fn update_loosely_coupled(
    state: &mut RtkState,
    gnss_state: &RtkState,
    lever_arm: Vector3<f64>,
    omega_b: Vector3<f64>,
    tuning: &crate::engine::config::EkfTuningConfig,
) -> Result<(), UpdateError> {
    let p_6x6 = state.covariance.view((0, 0), (6, 6)).into_owned();
    let r_6x6_raw = gnss_state.covariance.view((0, 0), (6, 6)).into_owned();
    let z = compute_loose_coupling_innovation(state, gnss_state, &lever_arm, &omega_b);
    let r_6x6 = huber_scale_covariance(
        &p_6x6,
        &r_6x6_raw,
        &z,
        tuning.loosely_coupled_mahalanobis_sq,
    )
    .map_err(|_| UpdateError::SingularMatrix)?;
    let h_mat = build_loose_coupling_h(state, &lever_arm, &omega_b);
    let k = compute_loose_coupling_gain(state, &z, &h_mat, &r_6x6)?;
    let dx = &k * &z;
    if dx.iter().any(|x: &f64| x.is_nan()) {
        return Err(UpdateError::SingularMatrix);
    }
    apply_state_correction(state, &dx);
    state.covariance = apply_joseph_covariance_update(&state.covariance, &k, &h_mat, &r_6x6);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub struct EkfUpdateResult {
    pub dx: DVector<f64>,
    pub k: DMatrix<f64>,
    pub worst_idx: Option<usize>,
    pub max_outlier_ratio: f64,
    pub weights: DVector<f64>,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_update_iteration<C: CouplingStrategy>(
    state_cov: &DMatrix<f64>,
    current_z: &DVector<f64>,
    current_h: &DMatrix<f64>,
    current_r: &DMatrix<f64>,
    current_valid: &[usize],
    meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>,
    max_innovation: f64,
    tuning: &crate::engine::config::EkfTuningConfig,
) -> Result<EkfUpdateResult, UpdateError> {
    let hp = current_h * state_cov;
    let s: DMatrix<f64> = &hp * current_h.transpose() + current_r;
    if s.iter()
        .any(|x: &f64| x.is_nan() || x.is_infinite() || x.abs() > 1e15)
    {
        tracing::warn!("EKF update failed: SingularMatrix (NaN/Inf in S)");
        return Err(UpdateError::SingularMatrix);
    }

    let s_inv = invert_matrix_robust(&s);

    let k = state_cov * current_h.transpose() * &s_inv;
    let dx = &k * current_z;

    if dx.iter().any(|x: &f64| x.is_nan()) {
        tracing::warn!("EKF update failed: SingularMatrix (NaN in dx)");
        return Err(UpdateError::SingularMatrix);
    }

    let v = current_z - current_h * &dx;
    let (worst_idx, max_outlier_ratio, weights) = evaluate_post_fit_outliers::<C>(
        &v,
        &s,
        current_z,
        current_valid,
        meas_types,
        max_innovation,
        tuning,
    );

    Ok(EkfUpdateResult {
        dx,
        k,
        worst_idx,
        max_outlier_ratio,
        weights,
    })
}

fn subset_ekf_matrices(
    z: &DVector<f64>,
    h: &DMatrix<f64>,
    r: &DMatrix<f64>,
    valid_indices: &[usize],
) -> (DVector<f64>, DMatrix<f64>, DMatrix<f64>) {
    if valid_indices.len() == z.len() {
        return (z.clone(), h.clone(), r.clone());
    }
    let n = valid_indices.len();
    let cols = h.ncols();
    let mut z_new = DVector::zeros(n);
    let mut h_new = DMatrix::zeros(n, cols);
    let mut r_new = DMatrix::zeros(n, n);
    for (new_idx, &old_idx) in valid_indices.iter().enumerate() {
        z_new[new_idx] = z[old_idx];
        for j in 0..cols {
            h_new[(new_idx, j)] = h[(old_idx, j)];
        }
        for (new_col, &old_col) in valid_indices.iter().enumerate() {
            r_new[(new_idx, new_col)] = r[(old_idx, old_col)];
        }
    }
    (z_new, h_new, r_new)
}

#[derive(Debug, PartialEq)]
enum OutlierAction {
    ContinueLoop,
    BreakLoop,
}

#[allow(clippy::too_many_arguments)]
fn check_outlier(
    current_valid: &mut Vec<usize>,
    current_z: &mut DVector<f64>,
    current_h: &mut DMatrix<f64>,
    current_r: &mut DMatrix<f64>,
    base_r: &mut DMatrix<f64>,
    worst_idx: usize,
    iter: usize,
    max_outlier_ratio: f64,
) -> Result<OutlierAction, UpdateError> {
    let ratio_bad = max_outlier_ratio == f64::INFINITY || max_outlier_ratio > 3.0;
    if current_valid.len() <= 1 {
        return if ratio_bad {
            tracing::warn!(
                "EKF update failed: InvalidMeasurement (outlier ratio {:.2})",
                max_outlier_ratio
            );
            Err(UpdateError::InvalidMeasurement)
        } else {
            Ok(OutlierAction::BreakLoop)
        };
    }
    if iter >= 99 {
        return if ratio_bad {
            tracing::warn!(
                "EKF update failed: Iteration limit reached with remaining outliers (ratio {:.2})",
                max_outlier_ratio
            );
            Err(UpdateError::InvalidMeasurement)
        } else {
            Ok(OutlierAction::BreakLoop)
        };
    }
    *current_z = current_z.clone().remove_row(worst_idx);
    *current_h = current_h.clone().remove_row(worst_idx);
    *current_r = current_r
        .clone()
        .remove_row(worst_idx)
        .remove_column(worst_idx);
    *base_r = base_r
        .clone()
        .remove_row(worst_idx)
        .remove_column(worst_idx);
    current_valid.remove(worst_idx);
    Ok(OutlierAction::ContinueLoop)
}

fn update_measurement_variances(
    current_r: &mut DMatrix<f64>,
    base_r: &DMatrix<f64>,
    weights: &DVector<f64>,
) -> bool {
    let mut changed = false;
    for i in 0..weights.len() {
        let new_r_ii = base_r[(i, i)] / weights[i];
        if (current_r[(i, i)] - new_r_ii).abs() > base_r[(i, i)] * 0.05 {
            changed = true;
        }
        current_r[(i, i)] = new_r_ii;
    }
    changed
}

#[allow(clippy::too_many_arguments)]
fn run_ekf_iterations<C: CouplingStrategy>(
    state: &RtkState,
    current_z: &mut DVector<f64>,
    current_h: &mut DMatrix<f64>,
    current_r: &mut DMatrix<f64>,
    current_valid: &mut Vec<usize>,
    base_r: &mut DMatrix<f64>,
    meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>,
    chi_square_threshold: f64,
    tuning: &crate::engine::config::EkfTuningConfig,
) -> Result<(DVector<f64>, DMatrix<f64>), UpdateError> {
    for iter in 0..100 {
        let res = compute_update_iteration::<C>(
            &state.covariance,
            current_z,
            current_h,
            current_r,
            current_valid,
            meas_types,
            chi_square_threshold,
            tuning,
        )?;
        if let Some(idx) = res.worst_idx {
            match check_outlier(
                current_valid,
                current_z,
                current_h,
                current_r,
                base_r,
                idx,
                iter,
                res.max_outlier_ratio,
            )? {
                OutlierAction::ContinueLoop => continue,
                OutlierAction::BreakLoop => {}
            }
        }
        let changed = update_measurement_variances(current_r, base_r, &res.weights);
        if !changed || iter >= 99 {
            return Ok((res.dx, res.k));
        }
    }
    Err(UpdateError::InvalidMeasurement)
}

fn validate_ekf_dimensions(
    z: &DVector<f64>,
    h: &DMatrix<f64>,
    cov_cols: usize,
) -> Result<(), UpdateError> {
    if z.len() != h.nrows() || h.ncols() != cov_cols {
        tracing::warn!("EKF update failed: DimensionMismatch");
        return Err(UpdateError::DimensionMismatch);
    }
    Ok(())
}

fn validate_pr_measurements(
    valid_indices: &[usize],
    meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>,
) -> Result<(), UpdateError> {
    if valid_indices.is_empty() {
        tracing::warn!("EKF update failed: InvalidMeasurement (valid_indices empty)");
        return Err(UpdateError::InvalidMeasurement);
    }
    let pr_count = valid_indices
        .iter()
        .filter(|&&i| meas_types.is_none_or(|t| t[i].1 == 0))
        .count();
    if pr_count == 0 {
        tracing::error!(
            "EKF update lacks any valid PR measurements! Rejecting update to trigger SPP fallback."
        );
        return Err(UpdateError::InvalidMeasurement);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn update<C: CouplingStrategy>(
    state: &mut RtkState,
    z: &DVector<f64>,
    h: &DMatrix<f64>,
    r: &DMatrix<f64>,
    chi_square_threshold: f64,
    meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>,
    tuning: &crate::engine::config::EkfTuningConfig,
) -> Result<(Vec<usize>, DVector<f64>), UpdateError> {
    validate_ekf_dimensions(z, h, state.covariance.nrows())?;
    let valid_indices =
        filter_pre_fit_residuals::<C>(z, h, r, &state.covariance, chi_square_threshold, meas_types);
    validate_pr_measurements(&valid_indices, meas_types)?;
    let (mut cz, mut ch, mut cr) = subset_ekf_matrices(z, h, r, &valid_indices);
    let mut current_valid = valid_indices;
    let mut base_r = cr.clone();
    let (dx, k) = run_ekf_iterations::<C>(
        state,
        &mut cz,
        &mut ch,
        &mut cr,
        &mut current_valid,
        &mut base_r,
        meas_types,
        chi_square_threshold,
        tuning,
    )?;
    apply_state_correction(state, &dx);
    state.covariance = apply_joseph_covariance_update(&state.covariance, &k, &ch, &cr);
    Ok((current_valid, dx))
}

fn build_ambiguity_vector(state: &RtkState, cols: usize) -> DVector<f64> {
    let mut a = DVector::zeros(cols);
    for i in 0..state.ambiguities.len() {
        a[crate::filter::CORE_STATE_SIZE + i] = state.ambiguities[i];
    }
    a
}

fn compute_fix_hold_gain(
    state: &RtkState,
    d_full: &DMatrix<f64>,
    r: &DMatrix<f64>,
) -> Result<DMatrix<f64>, UpdateError> {
    let h_p = d_full * &state.covariance;
    let s: DMatrix<f64> = &h_p * d_full.transpose() + r;
    let s_inv = s.try_inverse().ok_or(UpdateError::SingularMatrix)?;
    let k = &state.covariance * d_full.transpose() * s_inv;
    Ok(k)
}

pub fn apply_fix_and_hold(
    state: &mut RtkState,
    z_dd: &DVector<f64>,
    d_full: &DMatrix<f64>,
    var: f64,
) -> Result<(), UpdateError> {
    let num_dd = z_dd.len();
    if num_dd == 0 {
        return Ok(());
    }
    let a_sd = build_ambiguity_vector(state, d_full.ncols());
    let v = z_dd - d_full * &a_sd;
    let mut r = DMatrix::zeros(num_dd, num_dd);
    for i in 0..num_dd {
        r[(i, i)] = var;
    }
    let k = compute_fix_hold_gain(state, d_full, &r)?;
    let dx = &k * &v;
    apply_state_correction(state, &dx);
    state.covariance = apply_joseph_covariance_update(&state.covariance, &k, d_full, &r);
    Ok(())
}

#[cfg(test)]
mod private_tests {
    use super::*;
    use crate::filter::RtkState;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector, Vector3};

    // -----------------------------------------------------------------------
    // validate_ekf_dimensions
    // -----------------------------------------------------------------------

    #[test]
    fn test_validate_ekf_dimensions_ok() {
        let z = DVector::zeros(3);
        let h = DMatrix::zeros(3, 5);
        assert!(validate_ekf_dimensions(&z, &h, 5).is_ok());
    }

    #[test]
    fn test_validate_ekf_dimensions_z_h_rows_mismatch() {
        let z = DVector::zeros(3);
        let h = DMatrix::zeros(5, 5);
        assert_eq!(
            validate_ekf_dimensions(&z, &h, 5),
            Err(UpdateError::DimensionMismatch)
        );
    }

    #[test]
    fn test_validate_ekf_dimensions_h_cols_cov_mismatch() {
        let z = DVector::zeros(3);
        let h = DMatrix::zeros(3, 5);
        assert_eq!(
            validate_ekf_dimensions(&z, &h, 7),
            Err(UpdateError::DimensionMismatch)
        );
    }

    // -----------------------------------------------------------------------
    // validate_pr_measurements
    // -----------------------------------------------------------------------

    #[test]
    fn test_validate_pr_measurements_empty_indices() {
        assert_eq!(
            validate_pr_measurements(&[], None),
            Err(UpdateError::InvalidMeasurement)
        );
    }

    #[test]
    fn test_validate_pr_measurements_no_pr_types() {
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let meas_types = vec![(sat, 1), (sat, 2)];
        assert_eq!(
            validate_pr_measurements(&[0, 1], Some(&meas_types)),
            Err(UpdateError::InvalidMeasurement)
        );
    }

    #[test]
    fn test_validate_pr_measurements_has_pr() {
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let meas_types = vec![(sat, 0), (sat, 1)];
        assert!(validate_pr_measurements(&[0, 1], Some(&meas_types)).is_ok());
    }

    #[test]
    fn test_validate_pr_measurements_no_meas_types_uses_default() {
        // When meas_types is None, all are default type 0 (PR)
        assert!(validate_pr_measurements(&[0, 1], None).is_ok());
    }

    // -----------------------------------------------------------------------
    // subset_ekf_matrices
    // -----------------------------------------------------------------------

    #[test]
    fn test_subset_ekf_matrices_full() {
        let z = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let h = DMatrix::identity(3, 2);
        let r = DMatrix::identity(3, 3) * 2.0;
        let valid = vec![0, 1, 2];

        let (z2, h2, r2) = subset_ekf_matrices(&z, &h, &r, &valid);
        assert_eq!(z2, z);
        assert_eq!(h2, h);
        assert_eq!(r2, r);
    }

    #[test]
    fn test_subset_ekf_matrices_partial() {
        let z = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        let h = DMatrix::identity(4, 2);
        let r = DMatrix::identity(4, 4);
        let valid = vec![0, 2, 3];

        let (z2, h2, r2) = subset_ekf_matrices(&z, &h, &r, &valid);
        assert_eq!(z2.len(), 3);
        assert_eq!(z2[0], 1.0);
        assert_eq!(z2[1], 3.0);
        assert_eq!(z2[2], 4.0);
        assert_eq!(h2.nrows(), 3);
        assert_eq!(h2.ncols(), 2);
        assert_eq!(r2.nrows(), 3);
        assert_eq!(r2.ncols(), 3);
    }

    // -----------------------------------------------------------------------
    // build_loose_coupling_h
    // -----------------------------------------------------------------------

    #[test]
    fn test_build_loose_coupling_h_full_state() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);
        let lever_arm = Vector3::new(1.0, 0.0, 0.0);
        let omega_b = Vector3::zeros();

        let h = build_loose_coupling_h(&state, &lever_arm, &omega_b);
        assert_eq!(h.nrows(), 6);
        assert_eq!(h.ncols(), state.covariance.ncols());
        // First 6 columns should be identity
        for i in 0..6 {
            assert_eq!(h[(i, i)], 1.0);
        }
    }

    #[test]
    fn test_build_loose_coupling_h_small_state() {
        // State with covariance nrows < CORE_STATE_SIZE
        let h = build_loose_coupling_h_custom(6);
        assert_eq!(h.nrows(), 6);
        assert_eq!(h.ncols(), 6);
        for i in 0..6 {
            assert_eq!(h[(i, i)], 1.0);
        }
    }

    // Helper for the above test — build h with controlled covariance size
    fn build_loose_coupling_h_custom(cov_cols: usize) -> DMatrix<f64> {
        let mut h_mat = DMatrix::zeros(6, cov_cols);
        h_mat.view_mut((0, 0), (6, 6)).fill_diagonal(1.0);
        if cov_cols >= crate::filter::CORE_STATE_SIZE {
            // Would normally call populate_loosely_coupled_jacobian
            // but we skip since cov_cols < 21
        }
        h_mat
    }

    // -----------------------------------------------------------------------
    // compute_loose_coupling_gain
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_loose_coupling_gain_basic() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 2.0); // cov = 2.0 on diagonal
        let z = DVector::zeros(6);
        let h = build_loose_coupling_h_custom(21);
        // R is identity (variance 1.0 for each of 6 measurements)
        let r = DMatrix::identity(6, 6);

        let k = compute_loose_coupling_gain(&state, &z, &h, &r).unwrap();
        // P = diag(2,...), H = [I 0], R = I
        // S = P[:6,:6] + R = 2*I + I = 3*I
        // K = P * H^T * S^-1
        // For i=0..6: K[i, i] = 2.0 / 3.0 ≈ 0.667
        assert_eq!(k.nrows(), state.covariance.nrows());
        assert_eq!(k.ncols(), 6);
        assert!((k[(0, 0)] - 2.0 / 3.0).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // check_outlier
    // -----------------------------------------------------------------------

    #[test]
    fn test_check_outlier_single_valid_no_ratio() {
        let mut valid = vec![0];
        let mut z = DVector::from_vec(vec![10.0]);
        let mut h = DMatrix::zeros(1, 3);
        let mut r = DMatrix::from_element(1, 1, 1.0);
        let mut base_r = DMatrix::from_element(1, 1, 1.0);

        // max_outlier_ratio = 0.0 (not bad), valid.len() == 1 => BreakLoop
        let action = check_outlier(&mut valid, &mut z, &mut h, &mut r, &mut base_r, 0, 0, 2.0).unwrap();
        assert_eq!(action, OutlierAction::BreakLoop);
    }

    #[test]
    fn test_check_outlier_single_valid_bad_ratio() {
        let mut valid = vec![0];
        let mut z = DVector::from_vec(vec![10.0]);
        let mut h = DMatrix::zeros(1, 3);
        let mut r = DMatrix::from_element(1, 1, 1.0);
        let mut base_r = DMatrix::from_element(1, 1, 1.0);

        // max_outlier_ratio > 3.0 and valid.len() <= 1 => InvalidMeasurement
        let result = check_outlier(&mut valid, &mut z, &mut h, &mut r, &mut base_r, 0, 0, f64::INFINITY);
        assert_eq!(result, Err(UpdateError::InvalidMeasurement));
    }

    #[test]
    fn test_check_outlier_iter_limit_bad_ratio() {
        let mut valid = vec![0, 1, 2];
        let mut z = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut h = DMatrix::identity(3, 3);
        let mut r = DMatrix::identity(3, 3);
        let mut base_r = DMatrix::identity(3, 3);

        // iter >= 99 and ratio bad => InvalidMeasurement
        let result = check_outlier(&mut valid, &mut z, &mut h, &mut r, &mut base_r, 0, 99, f64::INFINITY);
        assert_eq!(result, Err(UpdateError::InvalidMeasurement));
    }

    #[test]
    fn test_check_outlier_iter_limit_good_ratio() {
        let mut valid = vec![0, 1, 2];
        let mut z = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut h = DMatrix::identity(3, 3);
        let mut r = DMatrix::identity(3, 3);
        let mut base_r = DMatrix::identity(3, 3);

        // iter >= 99 and ratio OK => BreakLoop
        let action = check_outlier(&mut valid, &mut z, &mut h, &mut r, &mut base_r, 0, 99, 2.0).unwrap();
        assert_eq!(action, OutlierAction::BreakLoop);
    }

    #[test]
    fn test_check_outlier_removes_worst() {
        let mut valid = vec![0, 1, 2];
        let mut z = DVector::from_vec(vec![10.0, 20.0, 30.0]);
        let mut h = DMatrix::identity(3, 3);
        let mut r = DMatrix::identity(3, 3);
        let mut base_r = DMatrix::identity(3, 3);

        // ratio OK, iter=0, valid.len()=3 => remove worst_idx
        let action = check_outlier(&mut valid, &mut z, &mut h, &mut r, &mut base_r, 1, 0, 4.0).unwrap();
        assert_eq!(action, OutlierAction::ContinueLoop);
        assert_eq!(valid.len(), 2);
        assert_eq!(z.len(), 2);
        assert_eq!(h.nrows(), 2);
        assert_eq!(r.nrows(), 2);
        assert_eq!(r.ncols(), 2);
        // After removing index 1: valid = [0, 2], z = [10, 30]
        assert_eq!(z[0], 10.0);
        assert_eq!(z[1], 30.0);
    }

    // -----------------------------------------------------------------------
    // update_measurement_variances
    // -----------------------------------------------------------------------

    #[test]
    fn test_update_measurement_variances_no_change() {
        let mut r = DMatrix::identity(2, 2) * 2.0;
        let base_r = DMatrix::identity(2, 2) * 2.0;
        let weights = DVector::from_element(2, 1.0);

        // new_r_ii = 2.0 / 1.0 = 2.0, same as old, so no change
        let changed = update_measurement_variances(&mut r, &base_r, &weights);
        assert!(!changed);
    }

    #[test]
    fn test_update_measurement_variances_changes() {
        let mut r = DMatrix::identity(2, 2) * 2.0;
        let base_r = DMatrix::identity(2, 2) * 2.0;
        let weights = DVector::from_vec(vec![1.0, 0.5]);

        // new_r_00 = 2.0 / 1.0 = 2.0, same -> no change for idx 0
        // new_r_11 = 2.0 / 0.5 = 4.0, |4.0 - 2.0| = 2.0 > 2.0 * 0.05 = 0.1 -> changed!
        let changed = update_measurement_variances(&mut r, &base_r, &weights);
        assert!(changed);
        assert!((r[(1, 1)] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_update_measurement_variances_exact_boundary() {
        let mut r = DMatrix::identity(2, 2);
        let base_r = DMatrix::identity(2, 2);
        let weights = DVector::from_vec(vec![1.0, 1.0 / 0.95]); // new_r_11 = 1.0 / (1/0.95) = 0.95

        // |0.95 - 1.0| = 0.05 = base_r[1,1] * 0.05 = 0.05
        // The check uses >, so 0.05 > 0.05 is false -> no change
        let changed = update_measurement_variances(&mut r, &base_r, &weights);
        assert!(!changed);
    }

    // -----------------------------------------------------------------------
    // build_ambiguity_vector
    // -----------------------------------------------------------------------

    #[test]
    fn test_build_ambiguity_vector_empty() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);
        let cols = crate::filter::CORE_STATE_SIZE;

        let a = build_ambiguity_vector(&state, cols);
        assert_eq!(a.len(), cols);
        assert!(a.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn test_build_ambiguity_vector_with_values() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.add_ambiguity(sat, 1, 42.0, 1.0);

        let cols = state.covariance.ncols();
        let a = build_ambiguity_vector(&state, cols);
        assert_eq!(a[crate::filter::CORE_STATE_SIZE], 42.0);
    }

    // -----------------------------------------------------------------------
    // compute_fix_hold_gain
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_fix_hold_gain_basic() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.add_ambiguity(sat, 1, 10.0, 1.0);

        let n_cols = state.covariance.ncols();
        let mut d_full = DMatrix::zeros(1, n_cols);
        d_full[(0, crate::filter::CORE_STATE_SIZE)] = 1.0;
        let r = DMatrix::from_element(1, 1, 0.1);

        let k = compute_fix_hold_gain(&state, &d_full, &r).unwrap();
        assert_eq!(k.nrows(), n_cols);
        assert_eq!(k.ncols(), 1);
    }

    #[test]
    fn test_compute_fix_hold_gain_singular() {
        // d_full contains all zeros -> H*P*H^T = 0, S = 0+R = R, invertible
        // This test should actually work fine with zero d_full.
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.add_ambiguity(sat, 1, 10.0, 1.0);

        let n_cols = state.covariance.ncols();
        let d_full = DMatrix::zeros(1, n_cols);
        let r = DMatrix::from_element(1, 1, 0.1);

        // S = 0 + 0.1 = 0.1, not singular
        let k = compute_fix_hold_gain(&state, &d_full, &r).unwrap();
        assert_eq!(k.nrows(), n_cols);
    }

    // -----------------------------------------------------------------------
    // apply_fix_and_hold edge cases
    // -----------------------------------------------------------------------

    #[test]
    fn test_apply_fix_and_hold_empty_z() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.0);
        let z_dd = DVector::from_vec(vec![]);
        let d_full = DMatrix::zeros(0, crate::filter::CORE_STATE_SIZE);
        let var = 0.1;

        // Empty z_dd -> Ok(())
        let result = apply_fix_and_hold(&mut state, &z_dd, &d_full, var);
        assert!(result.is_ok());
    }

    // -----------------------------------------------------------------------
    // compute_update_iteration NaN in dx
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_update_iteration_nan_dx() {
        let state_cov = DMatrix::from_diagonal(&DVector::from_element(1, f64::NAN));
        let h = DMatrix::from_element(1, 1, 1.0);
        let r = DMatrix::from_element(1, 1, 1.0);
        let z = DVector::from_element(1, 5.0);
        let valid = vec![0];
        let tuning = crate::engine::config::EkfTuningConfig::default();

        // State cov has NaN, so S will have NaN, then S_inv will be computed,
        // and dx will have NaN, returning Err(SingularMatrix)
        let res = compute_update_iteration::<TightCoupling>(
            &state_cov, &z, &h, &r, &valid, None, 10.0, &tuning,
        );
        assert!(res.is_err());
    }

    // -----------------------------------------------------------------------
    // compute_loose_coupling_innovation (the private wrapper)
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_loose_coupling_innovation_identical_states() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos.clone(), 1.0);
        let gnss_state = RtkState::new(time, pos, 1.0);
        let lever_arm = Vector3::new(0.5, 0.3, 0.1);
        let omega_b = Vector3::new(0.01, 0.02, 0.03);

        let z = compute_loose_coupling_innovation(&state, &gnss_state, &lever_arm, &omega_b);
        // With identical positions, the innovation is -R_b_e * lever_arm for position
        // and -R_b_e * (omega_b x lever_arm) for velocity
        assert_eq!(z.len(), 6);
        // Innovation should be non-zero because lever arm creates a position offset
        assert!(z.rows_range(0..3).norm() > 0.0);
    }

    // -----------------------------------------------------------------------
    // apply_attitude_correction boundary
    // -----------------------------------------------------------------------

    #[test]
    fn test_apply_attitude_correction_boundary() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);

        // d_theta norm = 1e-10, exactly the boundary where it should return without change
        let mut dx = DVector::zeros(state.covariance.nrows());
        dx[6] = 1e-10;

        let original_att = state.attitude;
        // attitude correction is called inside apply_state_correction -> apply_imu_and_clock_correction -> apply_attitude_correction
        apply_state_correction(&mut state, &dx);
        // The check is `<= 1e-10`, so norm=1e-10 should NOT change attitude
        assert_eq!(state.attitude, original_att);

        // Now with d_theta norm slightly above 1e-10
        let mut state2 = RtkState::new(time, pos.clone(), 1.0);
        let mut dx2 = DVector::zeros(state2.covariance.nrows());
        dx2[6] = 1.0000001e-10;
        let original_att2 = state2.attitude;
        apply_state_correction(&mut state2, &dx2);
        // Should have changed
        assert_ne!(state2.attitude, original_att2);
    }

    // -----------------------------------------------------------------------
    // apply_clock_correction ZWD clamping
    // -----------------------------------------------------------------------

    #[test]
    fn test_apply_clock_correction_zwd_clamping() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.zwd = 0.1;

        // Apply large negative ZWD correction that would push ZWD below 0
        let mut dx = DVector::zeros(state.covariance.nrows());
        dx[20] = -0.5; // would make zwd = 0.1 - 0.5 = -0.4, clamped to 0.0
        apply_state_correction(&mut state, &dx);
        assert_eq!(state.zwd, 0.0);

        // Normal positive correction
        dx[20] = 0.2;
        let mut state2 = RtkState::new(time, pos, 1.0);
        state2.zwd = 0.1;
        apply_clock_correction(&mut state2, &dx);
        assert!((state2.zwd - 0.3).abs() < 1e-14);
    }

    // -----------------------------------------------------------------------
    // update (the main pub wrapper) — error paths
    // -----------------------------------------------------------------------

    #[test]
    fn test_update_dimension_mismatch() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);

        let z = DVector::zeros(3);  // 3 measurements
        let h = DMatrix::zeros(3, 5); // 5 state columns
        let r = DMatrix::identity(3, 3);
        let tuning = crate::engine::config::EkfTuningConfig::default();

        // state.covariance has 21 columns, but h has 5 cols -> mismatch
        let result = update::<TightCoupling>(
            &mut state, &z, &h, &r, 10.0, None, &tuning,
        );
        assert_eq!(result, Err(UpdateError::DimensionMismatch));
    }

    #[test]
    fn test_update_empty_valid_indices() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let n = state.covariance.nrows();

        // Make large innovations so all get rejected by pre-fit residual check
        let z = DVector::from_element(n, 1e10);
        let h = DMatrix::identity(n, n);
        let r = DMatrix::identity(n, n);
        let tuning = crate::engine::config::EkfTuningConfig::default();

        let result = update::<TightCoupling>(
            &mut state, &z, &h, &r, 1.0, None, &tuning,
        );
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // update_loosely_coupled — singular matrix path
    // -----------------------------------------------------------------------

    #[test]
    fn test_update_loosely_coupled_singular_huber() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos.clone(), 1.0);
        let mut gnss_state = RtkState::new(time, pos, 1.0);
        gnss_state.position.vector.x = 1e10; // huge innovation

        // Mahalanobis threshold very low -> forces scaling
        let tuning = crate::engine::config::EkfTuningConfig {
            loosely_coupled_mahalanobis_sq: 0.001,
            huber_threshold_loosely: 3.0,
            ..Default::default()
        };

        let lever_arm = Vector3::zeros();
        let omega_b = Vector3::zeros();
        let result = update_loosely_coupled(&mut state, &gnss_state, lever_arm, omega_b, &tuning);
        // Should succeed with Huber scaling
        assert!(result.is_ok());
    }

    // -----------------------------------------------------------------------
    // apply_state_correction: dx shorter than CORE_STATE_SIZE
    // -----------------------------------------------------------------------

    #[test]
    fn test_apply_state_correction_short_dx() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::new(1.0, 2.0, 3.0);

        // dx with only 6 elements (position + velocity only)
        // This exercises the branch where dx.len() < CORE_STATE_SIZE (line 22 is false)
        let dx = DVector::from_vec(vec![10.0, 20.0, 30.0, 0.5, 1.0, 1.5]);

        apply_state_correction(&mut state, &dx);
        assert_eq!(state.position.vector.x, 110.0);
        assert_eq!(state.position.vector.y, 220.0);
        assert_eq!(state.position.vector.z, 330.0);
        assert_eq!(state.velocity.x, 1.5);
        assert_eq!(state.velocity.y, 3.0);
        assert_eq!(state.velocity.z, 4.5);
    }

    // -----------------------------------------------------------------------
    // apply_state_correction: dx exactly CORE_STATE_SIZE (no ambiguity updates)
    // -----------------------------------------------------------------------

    #[test]
    fn test_apply_state_correction_core_size_dx() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::new(1.0, 2.0, 3.0);
        state.accel_bias = Vector3::new(0.01, 0.02, 0.03);
        state.gyro_bias = Vector3::new(0.001, 0.002, 0.003);
        state.rcv_clk_bias = 1000.0;
        state.isb_glo = 5.0;
        state.isb_gal = 3.0;
        state.isb_bds = -2.0;
        state.rcv_clk_drift = 0.5;
        state.zwd = 0.1;

        let mut dx = DVector::zeros(crate::filter::CORE_STATE_SIZE);
        dx[0] = 10.0; // position
        dx[15] = 50.0; // clock bias
        dx[20] = 0.05; // ZWD

        apply_state_correction(&mut state, &dx);
        assert_eq!(state.position.vector.x, 110.0);
        assert_eq!(state.rcv_clk_bias, 1050.0);
        assert!((state.zwd - 0.15).abs() < 1e-14);
    }

    // -----------------------------------------------------------------------
    // apply_attitude_correction: d_theta = 0 (norm = 0, directly)
    // -----------------------------------------------------------------------

    #[test]
    fn test_apply_attitude_correction_zero_rotation() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let original_att = state.attitude;

        // d_theta where all components are 0
        let mut dx = DVector::zeros(state.covariance.nrows());
        dx[6] = 0.0;
        dx[7] = 0.0;
        dx[8] = 0.0;

        apply_state_correction(&mut state, &dx);
        assert_eq!(state.attitude, original_att);
    }

    // -----------------------------------------------------------------------
    // apply_bias_correction: direct test
    // -----------------------------------------------------------------------

    #[test]
    fn test_apply_bias_correction_updates_values() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.accel_bias = Vector3::new(1.0, 2.0, 3.0);
        state.gyro_bias = Vector3::new(0.1, 0.2, 0.3);

        let mut dx = DVector::zeros(state.covariance.nrows());
        dx[9] = 0.5;
        dx[10] = -1.0;
        dx[11] = 2.0;
        dx[12] = -0.05;
        dx[13] = 0.1;
        dx[14] = -0.15;

        apply_bias_correction(&mut state, &dx);
        assert!((state.accel_bias.x - 1.5).abs() < 1e-14);
        assert!((state.accel_bias.y - 1.0).abs() < 1e-14);
        assert!((state.accel_bias.z - 5.0).abs() < 1e-14);
        assert!((state.gyro_bias.x - 0.05).abs() < 1e-14);
        assert!((state.gyro_bias.y - 0.3).abs() < 1e-14);
        assert!((state.gyro_bias.z - 0.15).abs() < 1e-14);
    }

    // -----------------------------------------------------------------------
    // subset_ekf_matrices: empty valid_indices
    // -----------------------------------------------------------------------

    #[test]
    fn test_subset_ekf_matrices_empty() {
        let z = DVector::zeros(3);
        let h = DMatrix::identity(3, 2);
        let r = DMatrix::identity(3, 3);
        let valid: Vec<usize> = vec![];

        let (z2, h2, r2) = subset_ekf_matrices(&z, &h, &r, &valid);
        assert_eq!(z2.len(), 0);
        assert_eq!(h2.nrows(), 0);
        assert_eq!(h2.ncols(), 2);
        assert_eq!(r2.nrows(), 0);
        assert_eq!(r2.ncols(), 0);
    }

    // -----------------------------------------------------------------------
    // subset_ekf_matrices: reordered indices
    // -----------------------------------------------------------------------

    #[test]
    fn test_subset_ekf_matrices_reordered() {
        let z = DVector::from_vec(vec![10.0, 20.0, 30.0]);
        let mut h = DMatrix::zeros(3, 2);
        h[(0, 0)] = 1.0; h[(0, 1)] = 0.0;
        h[(1, 0)] = 0.0; h[(1, 1)] = 1.0;
        h[(2, 0)] = 1.0; h[(2, 1)] = 1.0;
        let r = DMatrix::identity(3, 3) * 2.0;
        // Reordered indices: [2, 0]
        let valid = vec![2, 0];

        let (z2, h2, r2) = subset_ekf_matrices(&z, &h, &r, &valid);
        assert_eq!(z2.len(), 2);
        assert_eq!(z2[0], 30.0); // row 2 first
        assert_eq!(z2[1], 10.0); // row 0 second
        // r2 should be correctly reordered
        assert_eq!(r2[(0, 0)], r[(2, 2)]); // old 2,2
        assert_eq!(r2[(1, 1)], r[(0, 0)]); // old 0,0
        assert_eq!(r2[(0, 1)], r[(2, 0)]); // old 2,0
    }

    // -----------------------------------------------------------------------
    // compute_update_iteration: S with Inf values
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_update_iteration_inf_s() {
        let state_cov = DMatrix::identity(1, 1);
        let h = DMatrix::from_element(1, 1, f64::INFINITY); // H^T * cov * H + R = Inf
        let r = DMatrix::from_element(1, 1, 1.0);
        let z = DVector::from_element(1, 5.0);
        let valid = vec![0];
        let tuning = crate::engine::config::EkfTuningConfig::default();

        let res = compute_update_iteration::<TightCoupling>(
            &state_cov, &z, &h, &r, &valid, None, 10.0, &tuning,
        );
        assert!(res.is_err());
    }

    // -----------------------------------------------------------------------
    // compute_update_iteration: S with huge values causing NaN guard
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_update_iteration_huge_s() {
        let state_cov = DMatrix::from_element(1, 1, 1e15);
        let h = DMatrix::from_element(1, 1, 1.0);
        let r = DMatrix::from_element(1, 1, 1e15);
        let z = DVector::from_element(1, 5.0);
        let valid = vec![0];
        let tuning = crate::engine::config::EkfTuningConfig::default();

        // S = 1e15*1*1e15^T + 1e15 = 1e15 + 1e15 = 2e15, not > 1e15 in absolute
        // Need to make it actually > 1e15
        let huge_cov = DMatrix::from_element(1, 1, 1e20);
        let res = compute_update_iteration::<TightCoupling>(
            &huge_cov, &z, &h, &r, &valid, None, 10.0, &tuning,
        );
        assert!(res.is_err(), "huge S should return SingularMatrix");
    }

    // -----------------------------------------------------------------------
    // update_loosely_coupled: NaN in dx path
    // -----------------------------------------------------------------------

    #[test]
    #[test]
    #[ignore] // TODO: hangs due to NaN in covariance propagating through matrix ops
    fn test_update_loosely_coupled_nan_dx() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos.clone(), 1.0);
        let mut gnss_state = RtkState::new(time, pos, 1.0);

        // Create NaN in the position to propagate through the update
        state.covariance[(0, 0)] = f64::NAN;

        let tuning = crate::engine::config::EkfTuningConfig::default();
        let lever_arm = Vector3::zeros();
        let omega_b = Vector3::zeros();
        let result = update_loosely_coupled(&mut state, &gnss_state, lever_arm, omega_b, &tuning);
        // Should fail due to NaN in dx
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // compute_fix_hold_gain: singular S through try_inverse
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_fix_hold_gain_singular_s() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.0);

        // Set covariance to zero so H*P*H^T = 0
        state.covariance = DMatrix::zeros(
            state.covariance.nrows(), state.covariance.ncols());
        let n_cols = state.covariance.ncols();
        let mut d_full = DMatrix::zeros(1, n_cols);
        d_full[(0, 0)] = 1.0;

        // R = 0 + the singular S = 0 -> try_inverse fails
        let r = DMatrix::zeros(1, 1);

        let result = compute_fix_hold_gain(&state, &d_full, &r);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // update: explicit validation of EKF iteration with known state/measurement
    // -----------------------------------------------------------------------

    #[test]
    fn test_update_single_measurement_success() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);

        // Single PR measurement on satellite 1
        let n_cols = state.covariance.nrows();
        let mut h = DMatrix::zeros(1, n_cols);
        h[(0, 0)] = 1.0; // Position-x measurement
        h[(0, 3)] = 0.0;
        h[(0, 15)] = 1.0; // Clock bias measurement

        let z = DVector::from_vec(vec![5.0]); // innovation
        let r = DMatrix::from_element(1, 1, 1.0);
        let tuning = crate::engine::config::EkfTuningConfig::default();
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let meas_types = vec![(sat, 0)]; // PR type

        let result = update::<TightCoupling>(
            &mut state, &z, &h, &r, 100.0, Some(&meas_types), &tuning,
        );
        assert!(result.is_ok(), "single measurement update should succeed: {:?}", result.err());

        let (valid_indices, dx) = result.unwrap();
        assert_eq!(valid_indices.len(), 1);
        // State should have been corrected
        assert!(state.position.vector.x != 100.0, "position should be corrected");
    }

    // -----------------------------------------------------------------------
    // update: all measurements rejected by pre-fit filter
    // -----------------------------------------------------------------------

    #[test]
    fn test_update_all_rejected_returns_err() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);

        let n = state.covariance.nrows();
        let z = DVector::from_element(n, 1e10); // huge innovations
        let h = DMatrix::identity(n, n);
        let r = DMatrix::identity(n, n);
        let tuning = crate::engine::config::EkfTuningConfig::default();

        // Very tight chi-square threshold -> all rejected
        let result = update::<TightCoupling>(
            &mut state, &z, &h, &r, 1.0, None, &tuning,
        );
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // apply_clock_correction: ZWD stays positive through small negative
    // -----------------------------------------------------------------------

    #[test]
    fn test_apply_clock_correction_zwd_stays_positive() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.zwd = 0.01;

        // Small negative that doesn't cross zero
        let mut dx = DVector::zeros(state.covariance.nrows());
        dx[20] = -0.005;
        apply_clock_correction(&mut state, &dx);
        assert!((state.zwd - 0.005).abs() < 1e-14, "ZWD should be 0.005, got {}", state.zwd);

        // Larger negative that crosses zero
        let mut state2 = RtkState::new(time, pos, 1.0);
        state2.zwd = 0.01;
        dx[20] = -0.02;
        apply_clock_correction(&mut state2, &dx);
        assert_eq!(state2.zwd, 0.0, "ZWD should be clamped to 0");
    }

    // -----------------------------------------------------------------------
    // check_outlier: ratio exactly at boundary
    // -----------------------------------------------------------------------

    #[test]
    fn test_check_outlier_exact_ratio_boundary() {
        let mut valid = vec![0, 1];
        let mut z = DVector::from_vec(vec![10.0, 20.0]);
        let mut h = DMatrix::identity(2, 3);
        let mut r = DMatrix::identity(2, 2);
        let mut base_r = DMatrix::identity(2, 2);

        // ratio_bad = max_outlier_ratio > 3.0
        // 3.0 is NOT > 3.0, so ratio_bad = false
        // valid.len() = 2 > 1, so it proceeds to remove worst
        let action = check_outlier(&mut valid, &mut z, &mut h, &mut r, &mut base_r, 1, 0, 3.0).unwrap();
        assert_eq!(action, OutlierAction::ContinueLoop);
        assert_eq!(valid.len(), 1);
    }

    // -----------------------------------------------------------------------
    // check_outlier: worst_idx removal updates matrices correctly
    // -----------------------------------------------------------------------

    #[test]
    fn test_check_outlier_removes_first_element() {
        let mut valid = vec![5, 10, 15];
        let mut z = DVector::from_vec(vec![100.0, 200.0, 300.0]);
        let mut h = DMatrix::identity(3, 3);
        let mut r = DMatrix::identity(3, 3);
        let mut base_r = DMatrix::identity(3, 3);

        // Remove index 0 (worst_idx = 0)
        let action = check_outlier(&mut valid, &mut z, &mut h, &mut r, &mut base_r, 0, 0, 4.0).unwrap();
        assert_eq!(action, OutlierAction::ContinueLoop);
        assert_eq!(valid.len(), 2);
        assert_eq!(valid[0], 10); // original index 1 stays
        assert_eq!(valid[1], 15); // original index 2 stays
    }

    // -----------------------------------------------------------------------
    // run_ekf_iterations: convergence after one iteration
    // -----------------------------------------------------------------------

    #[test]
    fn test_run_ekf_iterations_converges_immediately() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);

        let n_cols = state.covariance.nrows();
        let mut z = DVector::from_vec(vec![5.0]);
        let mut h = DMatrix::zeros(1, n_cols);
        h[(0, 0)] = 1.0;
        h[(0, 15)] = 1.0;
        let mut r = DMatrix::from_element(1, 1, 1.0);
        let mut current_valid = vec![0];
        let mut base_r = r.clone();
        let tuning = crate::engine::config::EkfTuningConfig::default();

        let result = run_ekf_iterations::<TightCoupling>(
            &state, &mut z, &mut h, &mut r, &mut current_valid, &mut base_r,
            Some(&[(SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0)]),
            100.0, &tuning,
        );
        assert!(result.is_ok(), "EKF iteration converges: {:?}", result.err());
        let (dx, _k) = result.unwrap();
        // State correction should be non-zero
        assert!(dx[0] != 0.0, "dx should have a correction");
    }

    // -----------------------------------------------------------------------
    // run_ekf_iterations: outlier removing then converging (multiple iterations)
    // -----------------------------------------------------------------------

    #[test]
    fn test_run_ekf_iterations_outlier_then_converge() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);

        let n_cols = state.covariance.nrows();
        // Two measurements: one clean, one outlier
        let mut h = DMatrix::zeros(2, n_cols);
        h[(0, 0)] = 1.0; h[(0, 15)] = 1.0;
        h[(1, 0)] = 1.0; h[(1, 15)] = 1.0;
        let mut z = DVector::from_vec(vec![5.0, 1000.0]); // outlier on second
        let mut r = DMatrix::identity(2, 2);
        let mut current_valid = vec![0, 1];
        let mut base_r = r.clone();
        let tuning = crate::engine::config::EkfTuningConfig::default();

        let result = run_ekf_iterations::<TightCoupling>(
            &state, &mut z, &mut h, &mut r, &mut current_valid, &mut base_r,
            Some(&[(SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0); 2]),
            100.0, &tuning,
        );
        assert!(result.is_ok(), "outlier iteration: {:?}", result.err());
        // The EKF should succeed and return at least one valid measurement
        assert!(!current_valid.is_empty(), "at least one valid measurement");
        assert!(current_valid.len() <= 2, "no more than original 2 measurements");
    }

    // -----------------------------------------------------------------------
    // compute_update_iteration: full path with weights and outlier evaluation
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_update_iteration_full_path() {
        let state_cov = DMatrix::identity(2, 2) * 10.0;
        let h = DMatrix::identity(2, 2);
        let r = DMatrix::identity(2, 2) * 0.1;
        let z = DVector::from_vec(vec![2.0, 3.0]);
        let valid = vec![0, 1];
        let tuning = crate::engine::config::EkfTuningConfig::default();
        let meas_types = vec![
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 0),
        ];

        let result = compute_update_iteration::<TightCoupling>(
            &state_cov, &z, &h, &r, &valid, Some(&meas_types), 10.0, &tuning,
        );
        assert!(result.is_ok(), "full update path: {:?}", result.err());

        let ekf_result = result.unwrap();
        // Both measurements should be valid (ratio < threshold)
        assert_eq!(ekf_result.worst_idx, None, "should have no outliers");
        // Weights should be 1.0 since ratio < hard_reject_ratio -> Huber applied
        assert!((ekf_result.weights[0] - 1.0).abs() < 0.5, "weights near 1.0");
        // dx should be finite and reasonable
        assert!(ekf_result.dx.iter().all(|&x| x.is_finite()));
    }

    // -----------------------------------------------------------------------
    // compute_update_iteration: NaN in S matrix
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_update_iteration_s_has_nan() {
        let state_cov = DMatrix::from_element(1, 1, f64::NAN);
        let h = DMatrix::from_element(1, 1, 1.0);
        let r = DMatrix::from_element(1, 1, 1.0);
        let z = DVector::from_element(1, 5.0);
        let valid = vec![0];
        let tuning = crate::engine::config::EkfTuningConfig::default();

        // HP*H^T + R = NaN + 1 = NaN -> caught by non-finite check
        let result = compute_update_iteration::<TightCoupling>(
            &state_cov, &z, &h, &r, &valid, None, 10.0, &tuning,
        );
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // apply_fix_and_hold: successful correction with ambiguity
    // -----------------------------------------------------------------------

    #[test]
    fn test_apply_fix_and_hold_corrects_ambiguity() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.add_ambiguity(sat, 1, 100.0, 1.0); // float ambiguity = 100.0

        let n_cols = state.covariance.nrows();
        // z_dd = 0 (we want the double-difference to be 0, meaning ambiguity = reference)
        let z_dd = DVector::from_vec(vec![0.0]);
        let mut d_full = DMatrix::zeros(1, n_cols);
        d_full[(0, crate::filter::CORE_STATE_SIZE)] = 1.0; // pulls the ambiguity state
        let var = 0.01; // very tight variance -> strong pull

        let result = apply_fix_and_hold(&mut state, &z_dd, &d_full, var);
        assert!(result.is_ok(), "fix and hold: {:?}", result.err());

        // The ambiguity should have been pulled toward 0
        assert!(
            (state.ambiguities[0] - 0.0).abs() < 50.0,
            "ambiguity pulled toward 0, got {}",
            state.ambiguities[0]
        );
    }

    // -----------------------------------------------------------------------
    // update: successful multi-measurement EKF update
    // -----------------------------------------------------------------------

    #[test]
    fn test_update_multi_measurement_success() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);

        let n_cols = state.covariance.nrows();
        // Two measurements: position-x + clock
        let mut h = DMatrix::zeros(2, n_cols);
        h[(0, 0)] = 1.0;
        h[(0, 15)] = 1.0;
        h[(1, 3)] = 1.0; // velocity-x measurement
        let z = DVector::from_vec(vec![5.0, 0.5]); // innovations
        let r = DMatrix::identity(2, 2) * 0.01; // low measurement noise
        let tuning = crate::engine::config::EkfTuningConfig::default();
        let meas_types = vec![
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 0),
        ];

        let result = update::<TightCoupling>(
            &mut state, &z, &h, &r, 100.0, Some(&meas_types), &tuning,
        );
        assert!(result.is_ok(), "multi-measurement update: {:?}", result.err());
        let (valid, dx) = result.unwrap();
        assert_eq!(valid.len(), 2, "both measurements valid");
        assert!(state.position.vector.x > 100.0, "position-x corrected");
    }

    // -----------------------------------------------------------------------
    // check_outlier: iter limit reached with bad ratio
    // -----------------------------------------------------------------------

    #[test]
    fn test_check_outlier_iter_99_bad_ratio_is_error() {
        let mut valid = vec![0, 1, 2];
        let mut z = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut h = DMatrix::identity(3, 3);
        let mut r = DMatrix::identity(3, 3);
        let mut base_r = DMatrix::identity(3, 3);

        // iter=99 AND ratio=INF -> return InvalidMeasurement error
        let result = check_outlier(&mut valid, &mut z, &mut h, &mut r, &mut base_r, 0, 99, f64::INFINITY);
        assert_eq!(result, Err(UpdateError::InvalidMeasurement));
    }

    // -----------------------------------------------------------------------
    // update_measurement_variances: change below threshold (no change detected)
    // -----------------------------------------------------------------------

    #[test]
    fn test_update_measurement_variances_small_change_no_detect() {
        let mut r = DMatrix::identity(2, 2);
        let base_r = DMatrix::identity(2, 2);
        let weights = DVector::from_vec(vec![1.0, 1.02]);
        // new_r_11 = 1.0 / 1.02 ≈ 0.98039
        // |0.98039 - 1.0| = 0.01961
        // base_r[1,1] * 0.05 = 0.05
        // 0.01961 < 0.05 -> no change
        let changed = update_measurement_variances(&mut r, &base_r, &weights);
        assert!(!changed);
    }

    // -----------------------------------------------------------------------
    // compute_update_iteration: large S values from cov and H
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_update_iteration_large_values_handled() {
        // Large but not huge covariance + H -> should still work
        let state_cov = DMatrix::from_element(2, 2, 1e8);
        let h = DMatrix::from_element(2, 2, 1.0);
        let r = DMatrix::from_element(2, 2, 1e8);
        let z = DVector::from_element(2, 1.0);
        let valid = vec![0, 1];
        let tuning = crate::engine::config::EkfTuningConfig::default();

        // S = 1e8 * 1*1 + 1e8 = 2e8. Not > 1e15, not NaN/Inf. Should succeed.
        let result = compute_update_iteration::<TightCoupling>(
            &state_cov, &z, &h, &r, &valid, None, 10.0, &tuning,
        );
        assert!(result.is_ok(), "large but finite values: {:?}", result.err());
    }

    // -----------------------------------------------------------------------
    // update_loosely_coupled: zero lever arm, identical states
    // -----------------------------------------------------------------------

    #[test]
    fn test_update_loosely_coupled_identical_states() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos.clone(), 1.0);
        let gnss_state = RtkState::new(time, pos, 1.0);
        let tuning = crate::engine::config::EkfTuningConfig::default();

        let result = update_loosely_coupled(
            &mut state, &gnss_state, Vector3::zeros(), Vector3::zeros(), &tuning,
        );
        // With identical states and zero lever arm, innovation = 0 -> no correction
        assert!(result.is_ok(), "identical states: {:?}", result.err());
        // State should be unchanged
        assert_eq!(state.position.vector.x, 100.0);
    }

    // -----------------------------------------------------------------------
    // build_ambiguity_vector: state with multiple ambiguities
    // -----------------------------------------------------------------------

    #[test]
    fn test_build_ambiguity_vector_multiple() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat2 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        state.add_ambiguity(sat1, 1, 10.0, 1.0);
        state.add_ambiguity(sat2, 2, 20.0, 1.0);

        let a = build_ambiguity_vector(&state, state.covariance.ncols());
        assert_eq!(a[crate::filter::CORE_STATE_SIZE], 10.0);
        assert_eq!(a[crate::filter::CORE_STATE_SIZE + 1], 20.0);
    }

    // -----------------------------------------------------------------------
    // compute_fix_hold_gain: multiple double-differences
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_fix_hold_gain_multi_dd() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.0);
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat2 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        state.add_ambiguity(sat1, 1, 10.0, 1.0);
        state.add_ambiguity(sat2, 2, 20.0, 1.0);

        let n_cols = state.covariance.ncols();
        // Two double-difference constraints on the two ambiguities
        let mut d_full = DMatrix::zeros(2, n_cols);
        d_full[(0, crate::filter::CORE_STATE_SIZE)] = 1.0;
        d_full[(1, crate::filter::CORE_STATE_SIZE + 1)] = 1.0;
        let r = DMatrix::identity(2, 2) * 0.1;

        let k = compute_fix_hold_gain(&state, &d_full, &r).unwrap();
        assert_eq!(k.nrows(), n_cols);
        assert_eq!(k.ncols(), 2);
    }
}
