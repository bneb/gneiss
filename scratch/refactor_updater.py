import re

with open("crates/gneiss-rtk/src/engine/updater.rs", "r") as f:
    content = f.read()

replacement = """use crate::filter::RtkState;
use nalgebra::{DMatrix, DVector, Vector3, UnitQuaternion};

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
        let d_theta = Vector3::new(dx[6], dx[7], dx[8]);
        if d_theta.norm() > 1e-10 {
            let dq = UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_normalize(d_theta), d_theta.norm());
            state.attitude = state.attitude * dq;
            state.attitude.renormalize();
        }
        
        state.accel_bias.x += dx[9];
        state.accel_bias.y += dx[10];
        state.accel_bias.z += dx[11];
        state.gyro_bias.x += dx[12];
        state.gyro_bias.y += dx[13];
        state.gyro_bias.z += dx[14];
        
        if crate::filter::CORE_STATE_SIZE > 15 {
            state.rcv_clk_bias += dx[15];
            state.rcv_clk_drift += dx[16];
            state.zwd += dx[17];
            state.zwd = state.zwd.max(0.0);
        }
    }
    
    if dx.len() > crate::filter::CORE_STATE_SIZE {
        for i in 0..state.ambiguities.len() {
            state.ambiguities[i] += dx[crate::filter::CORE_STATE_SIZE + i];
        }
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
            1 | 2 => max_innovation * 10000.0,
            3 => max_innovation * 1000.0,
            _ => max_innovation * max_innovation,
        };
        
        if nu * nu / s_ii < threshold {
            valid_indices.push(i);
        } else {
            tracing::debug!("EKF pre-fit rejected meas: type={}, nu={:.2}, s_ii={:.2}", meas_type, nu, s_ii);
        }
    }
    valid_indices
}

fn apply_dynamic_covariance_scaling(chi2: f64, threshold: f64, s_inv: &mut DMatrix<f64>) {
    if chi2 > threshold {
        let s_dcs = ((2.0 * threshold) / (threshold + chi2)).powi(2);
        *s_inv *= s_dcs;
        tracing::debug!("DCS Applied: chi2={:.2}, threshold={:.2}, scaling={:.4}", chi2, threshold, s_dcs);
    }
}

pub fn update_loosely_coupled(
    state: &mut RtkState,
    gnss_state: &RtkState,
    lever_arm: Vector3<f64>,
    omega_b: Vector3<f64>,
    tuning: &crate::engine::config::EkfTuningConfig,
) -> Result<(), UpdateError> {
    let p_6x6 = state.covariance.view((0, 0), (6, 6));
    let r_6x6 = gnss_state.covariance.view((0, 0), (6, 6));
    let s = p_6x6 + r_6x6;
    let mut s_inv = s.clone().cholesky().map(|c| c.inverse()).unwrap_or_else(|| (s + DMatrix::identity(6, 6) * 1e-6).try_inverse().unwrap());
    
    let r_b_e = state.attitude.to_rotation_matrix();
    let l_e = r_b_e * lever_arm;
    let pos_apc = state.position.vector + l_e;
    let v_apc = state.velocity + r_b_e * omega_b.cross(&lever_arm);

    let mut z = DVector::zeros(6);
    z.rows_mut(0, 3).copy_from(&(gnss_state.position.vector - pos_apc));
    z.rows_mut(3, 3).copy_from(&(gnss_state.velocity - v_apc));
    
    let chi2 = (&z.transpose() * &s_inv * &z)[(0, 0)];
    apply_dynamic_covariance_scaling(chi2, tuning.loosely_coupled_mahalanobis_sq, &mut s_inv);
    
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
    if dx.iter().any(|x| x.is_nan()) { return Err(UpdateError::SingularMatrix); }
    
    apply_state_correction(state, &dx);
    state.covariance = apply_joseph_covariance_update(&state.covariance, &k, &h_mat, &gnss_state.covariance.view((0, 0), (6, 6)).into_owned());
    Ok(())
}

fn compute_huber_weights(v: &DVector<f64>, s: &DMatrix<f64>, current_valid: &[usize], meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>, tuning: &crate::engine::config::EkfTuningConfig) -> (DVector<f64>, Vec<usize>) {
    let mut weights = DVector::from_element(v.len(), 1.0);
    let mut to_drop = Vec::new();
    
    for i in 0..v.len() {
        let orig_idx = current_valid[i];
        let meas_type = meas_types.map_or(0, |m| m[orig_idx].1);
        let s_ii = s[(i, i)];
        let u = v[i].abs() / s_ii.sqrt();
        
        let k_huber = match meas_type {
            1 | 2 => tuning.phase_outlier_ratio_thresh,
            3 => 3.0 * tuning.doppler_outlier_ratio_mult,
            _ => 3.0,
        };
        
        if u > k_huber * 5.0 {
            to_drop.push(i);
        } else if u > k_huber {
            weights[i] = k_huber / u;
        }
    }
    (weights, to_drop)
}

#[allow(clippy::too_many_arguments)]
pub fn update(state: &mut RtkState, z: &DVector<f64>, h: &DMatrix<f64>, r: &DMatrix<f64>, max_innovation: f64, meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>, _is_tightly_coupled: bool, tuning: &crate::engine::config::EkfTuningConfig) -> Result<Vec<usize>, UpdateError> {
    if z.len() != h.nrows() || h.ncols() != state.covariance.nrows() { return Err(UpdateError::DimensionMismatch); }
    
    let mut valid_indices = filter_pre_fit_residuals(z, h, r, &state.covariance, max_innovation, meas_types);
    if valid_indices.is_empty() { return Err(UpdateError::InvalidMeasurement); }
    
    let mut dx = DVector::zeros(state.covariance.nrows());
    let mut final_k = DMatrix::zeros(state.covariance.nrows(), valid_indices.len());
    let mut final_h = DMatrix::zeros(valid_indices.len(), h.ncols());
    let mut final_r = DMatrix::zeros(valid_indices.len(), valid_indices.len());
    
    for _iter in 0..5 {
        let (mut current_z, mut current_h, mut current_r) = extract_matrices(z, h, r, &valid_indices);
        
        let hp = &current_h * &state.covariance;
        let s = &hp * current_h.transpose() + &current_r;
        let s_inv = match s.clone().cholesky() {
            Some(chol) => chol.inverse(),
            None => (s.clone() + DMatrix::identity(s.nrows(), s.ncols()) * 1e-6).try_inverse().unwrap_or_else(|| DMatrix::identity(s.nrows(), s.ncols()))
        };
        
        let k = &state.covariance * current_h.transpose() * &s_inv;
        dx = &k * &current_z;
        if dx.iter().any(|x| x.is_nan()) { return Err(UpdateError::SingularMatrix); }

        let v = &current_z - &current_h * &dx;
        let (weights, to_drop) = compute_huber_weights(&v, &s, &valid_indices, meas_types, tuning);
        
        if !to_drop.is_empty() && valid_indices.len() > 4 {
            for &idx in to_drop.iter().rev() { valid_indices.remove(idx); }
            continue;
        }
        
        let mut converged = true;
        for i in 0..weights.len() {
            if weights[i] < 0.99 {
                r[(valid_indices[i], valid_indices[i])] /= weights[i];
                converged = false;
            }
        }
        
        final_k = k; final_h = current_h; final_r = current_r;
        if converged { break; }
    }
    
    apply_state_correction(state, &dx);
    state.covariance = apply_joseph_covariance_update(&state.covariance, &final_k, &final_h, &final_r);
    Ok(valid_indices)
}

fn extract_matrices(z: &DVector<f64>, h: &DMatrix<f64>, r: &DMatrix<f64>, valid_indices: &[usize]) -> (DVector<f64>, DMatrix<f64>, DMatrix<f64>) {
    let mut z_new = DVector::zeros(valid_indices.len());
    let mut h_new = DMatrix::zeros(valid_indices.len(), h.ncols());
    let mut r_new = DMatrix::zeros(valid_indices.len(), valid_indices.len());
    for (new_idx, &old_idx) in valid_indices.iter().enumerate() {
        z_new[new_idx] = z[old_idx];
        for j in 0..h.ncols() { h_new[(new_idx, j)] = h[(old_idx, j)]; }
        r_new[(new_idx, new_idx)] = r[(old_idx, old_idx)];
    }
    (z_new, h_new, r_new)
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
    let s_inv = s.try_inverse().ok_or(UpdateError::SingularMatrix)?;
    
    let mut k = &state.covariance * d_full.transpose() * s_inv;
    if state.covariance.nrows() > 15 {
        for i in 6..15 { for j in 0..k.ncols() { k[(i, j)] = 0.0; } }
    }
    
    let dx = &k * &v;
    apply_state_correction(state, &dx);
    state.covariance = apply_joseph_covariance_update(&state.covariance, &k, d_full, &r);
    Ok(())
}
"""

with open("crates/gneiss-rtk/src/engine/updater.rs", "w") as f:
    f.write(replacement)

print("Refactored updater.rs")
