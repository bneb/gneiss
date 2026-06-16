use crate::filter::RtkState;
use nalgebra::{DMatrix, DVector, Vector3, UnitQuaternion};
pub use crate::engine::updater_math::*;
use crate::math::inversion::invert_matrix_robust;
use crate::math::covariance::apply_joseph_covariance_update;
use crate::math::thresholding::huber_scale_covariance;

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
    if dx.len() >= crate::filter::CORE_STATE_SIZE { apply_imu_and_clock_correction(state, dx); }
    if dx.len() > crate::filter::CORE_STATE_SIZE {
        for i in 0..state.ambiguities.len() { state.ambiguities[i] += dx[crate::filter::CORE_STATE_SIZE + i]; }
    }
}

fn apply_imu_and_clock_correction(state: &mut RtkState, dx: &DVector<f64>) {
    let d_theta = Vector3::new(dx[6], dx[7], dx[8]);
    if d_theta.norm() > 1e-10 {
        let dq = UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_normalize(d_theta), d_theta.norm());
        state.attitude = dq * state.attitude;
        state.attitude.renormalize();
    }
    state.accel_bias.x += dx[9]; state.accel_bias.y += dx[10]; state.accel_bias.z += dx[11];
    state.gyro_bias.x += dx[12]; state.gyro_bias.y += dx[13]; state.gyro_bias.z += dx[14];
    if crate::filter::CORE_STATE_SIZE > 15 {
        state.rcv_clk_bias += dx[15];
        state.isb_glo += dx[16];
        state.isb_gal += dx[17];
        state.isb_bds += dx[18];
        state.rcv_clk_drift += dx[19];
        state.zwd = (state.zwd + dx[20]).max(0.0);
    }
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
    
    let r_b_e = state.attitude.to_rotation_matrix();
    let z = compute_loose_coupling_innovations(r_b_e.matrix(), &state.position.vector, &state.velocity, &gnss_state.position.vector, &gnss_state.velocity, &lever_arm, &omega_b);

    // Compute raw Mahalanobis distance and apply Huber scaling if needed
    let r_6x6 = huber_scale_covariance(&p_6x6, &r_6x6_raw, &z, tuning.loosely_coupled_mahalanobis_sq).map_err(|_| UpdateError::SingularMatrix)?;

    let s = &p_6x6 + &r_6x6;
    let s_inv = invert_matrix_robust(&s);
    
    let mut h_mat = DMatrix::zeros(6, state.covariance.ncols());
    h_mat.view_mut((0, 0), (6, 6)).fill_diagonal(1.0);

    if state.covariance.nrows() >= crate::filter::CORE_STATE_SIZE {
        populate_loosely_coupled_jacobian(&mut h_mat, &state.attitude, &lever_arm, &omega_b);
    }
    
    let k = &state.covariance * h_mat.transpose() * s_inv;
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

fn compute_update_iteration<C: CouplingStrategy>(
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
    if s.iter().any(|x: &f64| x.is_nan() || x.is_infinite() || x.abs() > 1e15) {
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
    let (worst_idx, max_outlier_ratio, weights) = evaluate_post_fit_outliers::<C>(&v, &s, current_z, current_valid, meas_types, max_innovation, tuning);
    
    Ok(EkfUpdateResult { dx, k, worst_idx, max_outlier_ratio, weights })
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
    if z.len() != h.nrows() || h.ncols() != state.covariance.nrows() {
        tracing::warn!("EKF update failed: DimensionMismatch");
        return Err(UpdateError::DimensionMismatch);
    }
    
    let valid_indices = filter_pre_fit_residuals::<C>(z, h, r, &state.covariance, chi_square_threshold, meas_types);

    let pr_valid_count = valid_indices.iter().filter(|&&i| meas_types.is_none_or(|t| t[i].1 == 0)).count();
    let _cp_valid_count = valid_indices.iter().filter(|&&i| meas_types.is_none_or(|t| t[i].1 == 1 || t[i].1 == 2)).count();
    
    if pr_valid_count == 0 {
        tracing::error!("EKF update lacks any valid PR measurements! Rejecting update to trigger SPP fallback.");
        return Err(UpdateError::InvalidMeasurement);
    }

    if valid_indices.is_empty() {
        tracing::warn!("EKF update failed: InvalidMeasurement (valid_indices empty)");
        return Err(UpdateError::InvalidMeasurement);
    }
    
    let (mut current_z, mut current_h, mut current_r) = if valid_indices.len() == z.len() {
        (z.clone(), h.clone(), r.clone())
    } else {
        let mut z_new = DVector::zeros(valid_indices.len());
        let mut h_new = DMatrix::zeros(valid_indices.len(), h.ncols());
        let mut r_new = DMatrix::zeros(valid_indices.len(), valid_indices.len());
        
        for (new_idx, &old_idx) in valid_indices.iter().enumerate() {
            z_new[new_idx] = z[old_idx];
            for j in 0..h.ncols() { h_new[(new_idx, j)] = h[(old_idx, j)]; }
            for (new_col, &old_col) in valid_indices.iter().enumerate() {
                r_new[(new_idx, new_col)] = r[(old_idx, old_col)];
            }
        }
        (z_new, h_new, r_new)
    };

    let mut current_valid = valid_indices;
    let mut dx = DVector::zeros(state.covariance.nrows());
    let mut k = DMatrix::zeros(state.covariance.nrows(), current_z.len());
    let mut base_r = current_r.clone();
    
    for _iter in 0..100 {
        let iter_res = compute_update_iteration::<C>(
            &state.covariance, &current_z, &current_h, &current_r, &current_valid, 
            meas_types, chi_square_threshold, tuning
        )?;
        let (dx_iter, k_iter, worst_idx, max_outlier_ratio, weights) = (iter_res.dx, iter_res.k, iter_res.worst_idx, iter_res.max_outlier_ratio, iter_res.weights);
        
        dx = dx_iter;
        k = k_iter;
        
        if let Some(idx) = worst_idx {
            if current_valid.len() > 1 {
                if _iter >= 99 { 
                    if max_outlier_ratio == f64::INFINITY || max_outlier_ratio > 3.0 {
                        tracing::warn!("EKF update failed: Iteration limit reached with remaining outliers (ratio {:.2})", max_outlier_ratio);
                        return Err(UpdateError::InvalidMeasurement);
                    }
                    break; 
                }
                current_z = current_z.remove_row(idx);
                current_h = current_h.remove_row(idx);
                current_r = current_r.remove_row(idx).remove_column(idx);
                base_r = base_r.remove_row(idx).remove_column(idx);
                current_valid.remove(idx);
                continue;
            } else if max_outlier_ratio == f64::INFINITY || max_outlier_ratio > 3.0 {
                tracing::warn!("EKF update failed: InvalidMeasurement (outlier ratio {:.2})", max_outlier_ratio);
                return Err(UpdateError::InvalidMeasurement);
            }
        }

        let mut weights_changed = false;
        for i in 0..weights.len() {
            let new_r_ii = base_r[(i, i)] / weights[i];
            if (current_r[(i, i)] - new_r_ii).abs() > base_r[(i, i)] * 0.05 {
                weights_changed = true;
            }
            current_r[(i, i)] = new_r_ii;
        }

        if !weights_changed || _iter >= 99 {
            break;
        }
    }

    apply_state_correction(state, &dx);
    state.covariance = apply_joseph_covariance_update(&state.covariance, &k, &current_h, &current_r);
    
    Ok((current_valid, dx))
}

pub fn apply_fix_and_hold(state: &mut RtkState, z_dd: &DVector<f64>, d_full: &DMatrix<f64>, var: f64) -> Result<(), UpdateError> {
    let num_dd = z_dd.len();
    if num_dd == 0 { return Ok(()); }
    
    let mut a_sd_full = DVector::zeros(d_full.ncols());
    for i in 0..state.ambiguities.len() {
        a_sd_full[crate::filter::CORE_STATE_SIZE + i] = state.ambiguities[i];
    }
    
    let v = z_dd - d_full * &a_sd_full;
    let mut r = DMatrix::zeros(num_dd, num_dd);
    for i in 0..num_dd { r[(i, i)] = var; }
    
    let h_p = d_full * &state.covariance;
    let s: DMatrix<f64> = &h_p * d_full.transpose() + &r;
    
    let s_inv: DMatrix<f64> = match s.try_inverse() {
        Some(inv) => inv,
        None => return Err(UpdateError::SingularMatrix),
    };
    
    let mut k = &state.covariance * d_full.transpose() * s_inv;
    // Dampen (don't zero) attitude and IMU bias gain rows to prevent
    // aggressive feedback while still allowing gradual coupling
    const FIX_HOLD_IMU_GAIN_DAMPING: f64 = 0.1;
    if state.covariance.nrows() > 15 {
        for i in 6..15 {
            for j in 0..k.ncols() { k[(i, j)] *= FIX_HOLD_IMU_GAIN_DAMPING; }
        }
    }
    let dx = &k * &v;
    apply_state_correction(state, &dx);
    state.covariance = apply_joseph_covariance_update(&state.covariance, &k, d_full, &r);
    
    Ok(())
}
