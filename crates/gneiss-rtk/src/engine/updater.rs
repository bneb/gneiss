use crate::filter::RtkState;
use nalgebra::{DMatrix, DVector, Vector3, UnitQuaternion};

/// Pre-fit chi-squared threshold for carrier phase (normalized innovation squared).
/// A CP innovation of ~0.5m with sigma ~0.01m gives chi2 ~2500. Reject only extreme outliers.
const CP_PRE_FIT_CHI2_THRESHOLD: f64 = 100.0;

/// Pre-fit chi-squared threshold for Doppler measurements.
/// A Doppler innovation of 1 m/s with sigma ~0.3 m/s gives chi2 ~11. Reject above 50.
const DOPPLER_PRE_FIT_CHI2_THRESHOLD: f64 = 50.0;

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

pub fn apply_joseph_covariance_update(
    p: &DMatrix<f64>,
    k: &DMatrix<f64>,
    h: &DMatrix<f64>,
    r: &DMatrix<f64>,
) -> DMatrix<f64> {
    let identity = DMatrix::identity(p.nrows(), p.ncols());
    let i_kh = identity - k * h;
    let mut p_new = &i_kh * p * i_kh.transpose() + k * r * k.transpose();
    
    for r_idx in 0..p_new.nrows() {
        for c_idx in 0..r_idx {
            let avg = (p_new[(r_idx, c_idx)] + p_new[(c_idx, r_idx)]) * 0.5;
            p_new[(r_idx, c_idx)] = avg;
            p_new[(c_idx, r_idx)] = avg;
        }
    }
    p_new
}

pub fn filter_pre_fit_residuals(
    z: &DVector<f64>,
    h: &DMatrix<f64>,
    r: &DMatrix<f64>,
    p: &DMatrix<f64>,
    max_innovation: f64,
    meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>,
) -> Vec<usize> {
    let mut valid_indices = Vec::new();
    for i in 0..z.len() {
        let nu = z[i];
        let h_row = h.row(i);
        let s_ii = (h_row * p * h_row.transpose())[(0, 0)] + r[(i, i)];
        
        let meas_type = meas_types.map_or(0, |m| m[i].1);
        let threshold = match meas_type {
            1 | 2 => CP_PRE_FIT_CHI2_THRESHOLD,
            3 => DOPPLER_PRE_FIT_CHI2_THRESHOLD,
            _ => max_innovation * max_innovation,
        };
        
        if nu * nu / s_ii < threshold {
            valid_indices.push(i);
        } else {
            let r_ii = r[(i, i)];
            if meas_type != 1 && meas_type != 2 && nu.abs() < 1000.0 && r_ii < 1.0 {
                tracing::debug!("EKF rejected Doppler/PR measurement! type={}, nu={:.2}, s_ii={:.2}, r_ii={:.4}", meas_type, nu, s_ii, r_ii);
            } else {
                tracing::debug!("EKF rejected meas: type={}, nu={:.2}, s_ii={:.2}, r_ii={:.4}", meas_type, nu, s_ii, r_ii);
            }
        }
    }
    valid_indices
}

#[allow(clippy::too_many_arguments)]
pub fn evaluate_post_fit_outliers(
    v: &DVector<f64>,
    s: &DMatrix<f64>,
    current_z: &DVector<f64>,
    current_valid: &[usize],
    meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>,
    max_innovation: f64,
    is_tightly_coupled: bool,
    tuning: &crate::engine::config::EkfTuningConfig,
) -> (Option<usize>, f64) {
    let mut max_outlier_ratio = 0.0;
    let mut worst_idx = None;
    
    for i in 0..v.len() {
        let orig_idx = current_valid[i];
        let meas_type = meas_types.map_or(0, |m| m[orig_idx].1);
        let s_ii = s[(i, i)];
        let ratio = v[i].abs() / s_ii.sqrt();
        
        let thresh = match meas_type {
            0 => max_innovation, 
            1 | 2 => tuning.phase_outlier_ratio_thresh,
            3 => max_innovation * tuning.doppler_outlier_ratio_mult,
            _ => max_innovation,
        };
        
        let abs_thresh = match meas_type {
            0 => tuning.pr_abs_thresh,
            1 | 2 => tuning.cp_abs_thresh,
            3 => tuning.dop_abs_thresh,
            _ => 40.0,
        };
        
        if (v[i].abs() > thresh && ratio > max_outlier_ratio) || (is_tightly_coupled && current_z[i].abs() > abs_thresh) {
            if v[i].abs() > thresh && ratio > max_outlier_ratio {
                max_outlier_ratio = ratio;
                worst_idx = Some(i);
            } else if is_tightly_coupled && current_z[i].abs() > abs_thresh {
                worst_idx = Some(i);
                max_outlier_ratio = f64::INFINITY;
            }
        }
    }
    (worst_idx, max_outlier_ratio)
}

fn huber_scale_covariance(
    p: &DMatrix<f64>,
    r: &DMatrix<f64>,
    z: &DVector<f64>,
    tuning: &crate::engine::config::EkfTuningConfig,
) -> Result<DMatrix<f64>, UpdateError> {
    let s_raw = p + r;
    let s_raw_inv = match s_raw.clone().cholesky() {
        Some(chol) => chol.inverse(),
        None => return Err(UpdateError::SingularMatrix),
    };
    let mahal_sq = (&z.transpose() * &s_raw_inv * z)[(0, 0)];
    let huber_sq = tuning.huber_threshold_loosely.powi(2);

    if mahal_sq <= huber_sq {
        if mahal_sq > tuning.loosely_coupled_mahalanobis_sq {
            return Err(UpdateError::InvalidMeasurement);
        }
        return Ok(r.clone());
    }
    // Inflate R so effective Mahalanobis clamps to huber_sq
    let scale = mahal_sq / huber_sq;
    let r_scaled = r * scale;
    if huber_sq > tuning.loosely_coupled_mahalanobis_sq {
        return Err(UpdateError::InvalidMeasurement);
    }
    Ok(r_scaled)
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
    let l_e = r_b_e * lever_arm;
    let pos_apc = state.position.vector + l_e;
    let v_apc = state.velocity + r_b_e * omega_b.cross(&lever_arm);

    let mut z = DVector::zeros(6);
    z.rows_mut(0, 3).copy_from(&(gnss_state.position.vector - pos_apc));
    z.rows_mut(3, 3).copy_from(&(gnss_state.velocity - v_apc));

    // Compute raw Mahalanobis distance and apply Huber scaling if needed
    let r_6x6 = huber_scale_covariance(&p_6x6, &r_6x6_raw, &z, tuning)?;

    let s = &p_6x6 + &r_6x6;
    let s_inv = match s.clone().cholesky() {
        Some(chol) => chol.inverse(),
        None => match (s + DMatrix::identity(6, 6) * 1e-6).cholesky() {
            Some(chol) => chol.inverse(),
            None => return Err(UpdateError::SingularMatrix),
        }
    };
    
    let mut h_mat = DMatrix::zeros(6, state.covariance.ncols());
    h_mat.view_mut((0, 0), (6, 6)).fill_diagonal(1.0);

    if state.covariance.nrows() >= crate::filter::CORE_STATE_SIZE {
        let h_pos_att = -l_e.cross_matrix();
        let a_e = r_b_e * omega_b.cross(&lever_arm);
        let h_vel_att = -a_e.cross_matrix();
        let h_vel_bg = r_b_e.matrix() * lever_arm.cross_matrix();
        
        for i in 0..3 {
            for j in 0..3 {
                h_mat[(i, 6 + j)] = h_pos_att[(i, j)];
                h_mat[(3 + i, 6 + j)] = h_vel_att[(i, j)];
                h_mat[(3 + i, 12 + j)] = h_vel_bg[(i, j)];
            }
        }
    }
    
    let k = &state.covariance * h_mat.transpose() * s_inv;
    let dx = &k * &z;
    
    if dx.iter().any(|x| x.is_nan()) {
        return Err(UpdateError::SingularMatrix);
    }
    
    apply_state_correction(state, &dx);
    state.covariance = apply_joseph_covariance_update(&state.covariance, &k, &h_mat, &r_6x6);
    
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn update(state: &mut RtkState, z: &DVector<f64>, h: &DMatrix<f64>, r: &DMatrix<f64>, max_innovation: f64, meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>, is_tightly_coupled: bool, tuning: &crate::engine::config::EkfTuningConfig) -> Result<Vec<usize>, UpdateError> {
    if z.len() != h.nrows() || h.ncols() != state.covariance.nrows() {
        return Err(UpdateError::DimensionMismatch);
    }
    
    let valid_indices = filter_pre_fit_residuals(z, h, r, &state.covariance, max_innovation, meas_types);

    let pr_valid_count = valid_indices.iter().filter(|&&i| meas_types.map_or(true, |t| t[i].1 == 0)).count();
    let total_pr = meas_types.map_or(0, |t| t.iter().filter(|&&type_| type_.1 == 0).count());
    
    if total_pr > 0 && pr_valid_count == 0 {
        tracing::error!("EKF rejected ALL {} pseudoranges! Force reset to SPP.", total_pr);
        return Err(UpdateError::InvalidMeasurement);
    }

    if valid_indices.is_empty() {
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
    
    for _iter in 0..21 {
        let hp = &current_h * &state.covariance;
        let s = &hp * current_h.transpose() + &current_r;
        if s.iter().any(|x| x.is_nan() || x.is_infinite() || x.abs() > 1e15) {
            return Err(UpdateError::SingularMatrix);
        }
        
        let s_inv = match s.clone().cholesky() {
            Some(chol) => chol.inverse(),
            None => {
                let regularized = s.clone() + DMatrix::identity(s.nrows(), s.ncols()) * 1e-6;
                match regularized.try_inverse() {
                    Some(inv) => inv,
                    None => return Err(UpdateError::SingularMatrix),
                }
            }
        };
        
        k = &state.covariance * current_h.transpose() * &s_inv;
        dx = &k * &current_z;
        
        if dx.iter().any(|x| x.is_nan()) {
            return Err(UpdateError::SingularMatrix);
        }

        let v = &current_z - &current_h * &dx;
        let (worst_idx, max_outlier_ratio) = evaluate_post_fit_outliers(&v, &s, &current_z, &current_valid, meas_types, max_innovation, is_tightly_coupled, tuning);
        
        if let Some(idx) = worst_idx {
            if current_valid.len() > 4 {
                if _iter >= 20 { break; }
                current_z = current_z.remove_row(idx);
                current_h = current_h.remove_row(idx);
                current_r = current_r.remove_row(idx).remove_column(idx);
                current_valid.remove(idx);
                continue;
            } else if max_outlier_ratio == f64::INFINITY || max_outlier_ratio > 3.0 {
                return Err(UpdateError::InvalidMeasurement);
            }
        }
        break;
    }

    apply_state_correction(state, &dx);
    state.covariance = apply_joseph_covariance_update(&state.covariance, &k, &current_h, &current_r);
    
    Ok(current_valid)
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
    let s = &h_p * d_full.transpose() + &r;
    
    let s_inv = match s.try_inverse() {
        Some(inv) => inv,
        None => return Err(UpdateError::SingularMatrix),
    };
    
    let mut k = &state.covariance * d_full.transpose() * s_inv;
    if state.covariance.nrows() > 15 {
        for i in 6..15 {
            for j in 0..k.ncols() { k[(i, j)] = 0.0; }
        }
    }
    
    let dx = &k * &v;
    apply_state_correction(state, &dx);
    state.covariance = apply_joseph_covariance_update(&state.covariance, &k, d_full, &r);
    
    Ok(())
}
