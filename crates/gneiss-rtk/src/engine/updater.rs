pub use crate::engine::updater_math::*;
use crate::filter::RtkState;
use crate::math::covariance::apply_joseph_covariance_update;
use crate::math::inversion::invert_matrix_robust;
use crate::math::thresholding::huber_scale_covariance;
use nalgebra::{DMatrix, DVector, UnitQuaternion, Vector3};

#[derive(Debug)]
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

#[derive(PartialEq)]
enum OutlierAction {
    ContinueLoop,
    BreakLoop,
}

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
