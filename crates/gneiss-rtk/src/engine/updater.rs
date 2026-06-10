use crate::filter::RtkState;
use nalgebra::{DMatrix, DVector, Vector3, UnitQuaternion};

#[derive(Debug)]
pub enum UpdateError {
    SingularMatrix,
    DimensionMismatch,
    InvalidMeasurement,
}

pub fn update_loosely_coupled(
    state: &mut RtkState,
    gnss_state: &RtkState,
    lever_arm: Vector3<f64>,
    omega_b: Vector3<f64>,
) -> Result<(), UpdateError> {
    let p_6x6 = state.covariance.view((0, 0), (6, 6));
    let r_6x6 = gnss_state.covariance.view((0, 0), (6, 6));
    
    let s = p_6x6 + r_6x6;
    let s_inv = match s.clone().cholesky() {
        Some(chol) => chol.inverse(),
        None => {
            let reg = DMatrix::identity(6, 6) * 1e-6;
            match (s.clone() + reg).cholesky() {
                Some(chol) => chol.inverse(),
                None => return Err(UpdateError::SingularMatrix),
            }
        }
    };
    
    let r_b_e = state.attitude.to_rotation_matrix();
    let l_e = r_b_e * lever_arm;
    let pos_apc = state.position.vector + l_e;
    let v_apc = state.velocity + r_b_e * omega_b.cross(&lever_arm);

    let mut z = DVector::zeros(6);
    z.rows_mut(0, 3).copy_from(&(gnss_state.position.vector - pos_apc));
    z.rows_mut(3, 3).copy_from(&(gnss_state.velocity - v_apc));
    
    let mahalanobis_sq = (&z.transpose() * &s_inv * &z)[(0, 0)];
    
    // 6 DOF chi-square 99.9% is 22.46. 
    // We use a generous threshold because SPP covariance can be overly optimistic
    if mahalanobis_sq > 250.0 {
        return Err(UpdateError::InvalidMeasurement); // Trigger reset
    }
    
    let mut h_mat = DMatrix::zeros(6, state.covariance.ncols());
    h_mat.view_mut((0, 0), (6, 6)).fill_diagonal(1.0);

    if state.covariance.nrows() >= crate::filter::CORE_STATE_SIZE {
        let h_pos_att = -l_e.cross_matrix();
        for i in 0..3 {
            for j in 0..3 {
                h_mat[(i, 6 + j)] = h_pos_att[(i, j)];
            }
        }
        
        let a_b = omega_b.cross(&lever_arm);
        let a_e = r_b_e * a_b;
        let h_vel_att = -a_e.cross_matrix();
        for i in 0..3 {
            for j in 0..3 {
                h_mat[(3 + i, 6 + j)] = h_vel_att[(i, j)];
            }
        }
        
        let h_vel_bg = r_b_e.matrix() * lever_arm.cross_matrix();
        for i in 0..3 {
            for j in 0..3 {
                h_mat[(3 + i, 12 + j)] = h_vel_bg[(i, j)];
            }
        }
    }
    
    let k = &state.covariance * h_mat.transpose() * s_inv;
    let dx = &k * &z;
    
    if dx.iter().any(|x| x.is_nan()) {
        return Err(UpdateError::SingularMatrix);
    }
    
    state.position.vector += dx.rows(0, 3);
    state.velocity += dx.rows(3, 3);
    
    if state.covariance.nrows() >= crate::filter::CORE_STATE_SIZE {
        let d_theta = dx.rows(6, 3);
        state.accel_bias += dx.rows(9, 3);
        state.gyro_bias += dx.rows(12, 3);
        let angle = d_theta.norm();
        if angle > 1e-8 {
            let axis_vec3 = nalgebra::Vector3::new(d_theta[0], d_theta[1], d_theta[2]);
            let axis = nalgebra::Unit::new_normalize(axis_vec3);
            let d_quat = nalgebra::UnitQuaternion::from_axis_angle(&axis, angle);
            state.attitude = d_quat * state.attitude;
        }
    }
    
    let id = DMatrix::identity(state.covariance.nrows(), state.covariance.ncols());
    let ikh = id - &k * h_mat;
    state.covariance = &ikh * &state.covariance * ikh.transpose() + &k * r_6x6 * k.transpose();
    
    Ok(())
}

pub fn update(state: &mut RtkState, z: &DVector<f64>, h: &DMatrix<f64>, r: &DMatrix<f64>, max_innovation: f64, meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>) -> Result<Vec<usize>, UpdateError> {
    if z.len() != h.nrows() || h.ncols() != state.covariance.nrows() {
        return Err(UpdateError::DimensionMismatch);
    }
    
    let mut valid_indices = Vec::new();

    for i in 0..z.len() {
        let nu = z[i];
        let h_row = h.row(i);
        let p_ht = &state.covariance * h_row.transpose();
        let s_ii = (h_row * p_ht)[(0, 0)] + r[(i, i)];
        
        let meas_type = meas_types.map_or(0, |m| m[i].1);
        
        let threshold = match meas_type {
            1 | 2 => max_innovation * 10000.0, // CP
            3 => max_innovation * 1000.0,        // Doppler (allow large innovations to correct INS velocity drift)
            _ => max_innovation * max_innovation,      // PR
        };
        
        let is_valid = nu * nu / s_ii < threshold;

        
        if is_valid {
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

    // Check how many PR measurements were valid
    let pr_valid_count = valid_indices.iter().filter(|&&i| meas_types.as_ref().map(|t| t[i].1 == 0).unwrap_or(true)).count();
    
    // If no pseudoranges are valid but we had some in the input, we are totally lost.
    // We should reject everything so the filter resets to SPP!
    let total_pr = meas_types.as_ref().map(|t| t.iter().filter(|&&type_| type_.1 == 0).count()).unwrap_or(0);
    if total_pr > 0 && pr_valid_count == 0 {
        tracing::error!("EKF rejected ALL {} pseudoranges! Force reset to SPP.", total_pr);
        return Err(UpdateError::InvalidMeasurement);
    }

    if valid_indices.is_empty() {
        return Err(UpdateError::InvalidMeasurement);
    }
    
    let (z_filt, h_filt, r_filt) = if valid_indices.len() == z.len() {
        (z.clone(), h.clone(), r.clone())
    } else {
        let mut z_new = DVector::zeros(valid_indices.len());
        let mut h_new = DMatrix::zeros(valid_indices.len(), h.ncols());
        let mut r_new = DMatrix::zeros(valid_indices.len(), valid_indices.len());
        
        for (new_idx, &old_idx) in valid_indices.iter().enumerate() {
            z_new[new_idx] = z[old_idx];
            for j in 0..h.ncols() {
                h_new[(new_idx, j)] = h[(old_idx, j)];
            }
            // Assume R is diagonal
            r_new[(new_idx, new_idx)] = r[(old_idx, old_idx)];
        }
        (z_new, h_new, r_new)
    };

    let mut current_z = z_filt;
    let mut current_h = h_filt;
    let mut current_r = r_filt;
    let mut current_valid = valid_indices;

    let mut dx = DVector::zeros(state.covariance.nrows());
    let mut k = DMatrix::zeros(state.covariance.nrows(), current_z.len());
    
    for _iter in 0..21 {
        let hp = &current_h * &state.covariance;
        let s = &hp * current_h.transpose() + &current_r;
        if s.iter().any(|x| x.is_nan() || x.is_infinite() || x.abs() > 1e15) {
            tracing::error!("EKF error: Non-finite values in innovation covariance matrix S.");
            return Err(UpdateError::SingularMatrix);
        }
        let s_inv = match s.clone().cholesky() {
            Some(chol) => chol.inverse(),
            None => {
                let reg = DMatrix::identity(s.nrows(), s.ncols()) * 1e-6;
                (s.clone() + reg).try_inverse().unwrap_or_else(|| DMatrix::identity(s.nrows(), s.ncols()))
            }
        };
        k = &state.covariance * current_h.transpose() * &s_inv;
        dx = &k * &current_z;
        
        if dx.iter().any(|x| x.is_nan()) {
            tracing::error!("EKF error: NaN detected in state correction (dx). Cov max: {}", state.covariance.abs().max());
            return Err(UpdateError::SingularMatrix);
        }

        // Post-fit residual check
        let v = &current_z - &current_h * &dx;
        
        let mut max_outlier_ratio = 0.0;
        let mut worst_idx = None;
        
        for i in 0..v.len() {
            let orig_idx = current_valid[i];
            let meas_type = meas_types.map_or(0, |m| m[orig_idx].1);
            
            let s_ii = s[(i, i)];
            let ratio = v[i].abs() / s_ii.sqrt();
            
            let thresh = match meas_type {
                0 => max_innovation, 
                1 | 2 => 5.0, // 5.0m for Phase to allow tracking vehicle dynamics over outages
                3 => max_innovation * 2.0, // Doppler max innovation (typically 15-30m/s)
                _ => max_innovation,
            };
            
            let abs_thresh = match meas_type {
                0 => 40.0, // Pseudorange max 40m error
                1 | 2 => 1.0, // Phase max 1m error
                3 => 15.0, // Doppler max 15m/s error
                _ => 40.0,
            };
            
            let is_tightly_coupled = state.covariance.nrows() > 6;
            
            if (v[i].abs() > thresh && ratio > max_outlier_ratio) || (is_tightly_coupled && current_z[i].abs() > abs_thresh) {
                if v[i].abs() > thresh && ratio > max_outlier_ratio {
                    max_outlier_ratio = ratio;
                    worst_idx = Some(i);
                } else if is_tightly_coupled && current_z[i].abs() > abs_thresh {
                    // Force rejection if absolute value is too high (only for INS to protect attitude)
                    worst_idx = Some(i);
                    max_outlier_ratio = f64::INFINITY;
                }
            }
        }
        
        if let Some(idx) = worst_idx {
            if current_valid.len() > 4 {
                if _iter >= 20 { 
                    tracing::debug!("Max outlier rejection iterations reached. Not dropping meas type {} with post-fit residual {:.2}", meas_types.map_or(0, |m| m[current_valid[idx]].1), v[idx]);
                    break;
                }
                current_z = current_z.remove_row(idx);
                current_h = current_h.remove_row(idx);
                current_r = current_r.remove_row(idx).remove_column(idx);
                current_valid.remove(idx);
                continue;
            } else if max_outlier_ratio == f64::INFINITY || max_outlier_ratio > 3.0 {
                // If we hit the minimum measurement limit but the remaining ones are still garbage
                // (e.g. max_outlier_ratio == INFINITY from our absolute cap, or ratio > 3.0)
                tracing::warn!("EKF update rejected: Remaining measurements are extreme outliers (ratio={:.1}, z={:.1}).", max_outlier_ratio, current_z[idx]);
                return Err(UpdateError::InvalidMeasurement);
            } else {
                break;
            }
        }
        
        break;
    }

    let h_filt = current_h;
    let r_filt = current_r;
    let valid_indices = current_valid;

    let pos_change = (dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]).sqrt();
    if pos_change > 250.0 {
        tracing::debug!("EKF large state correction: dx_pos = {:.1}m, z_max = {:.1}. Allowing update for convergence.", pos_change, z.abs().max());
    }

    if state.epoch_count < 5 || state.epoch_count.is_multiple_of(100) {
        tracing::info!("EKF Update dx: pos=[{:.3}, {:.3}, {:.3}] vel=[{:.3}, {:.3}, {:.3}] att=[{:.5}, {:.5}, {:.5}] | z len: {}", dx[0], dx[1], dx[2], dx[3], dx[4], dx[5], dx[6], dx[7], dx[8], z.len());
    }

    // Apply state correction
    state.position.vector.x += dx[0];
    state.position.vector.y += dx[1];
    state.position.vector.z += dx[2];
    state.velocity.x += dx[3];
    state.velocity.y += dx[4];
    state.velocity.z += dx[5];
    
    // Apply Attitude correction (Small-angle rotation error in ECEF)
    let d_psi = Vector3::new(dx[6], dx[7], dx[8]);
    let angle = d_psi.norm();
    if angle > 1e-12 {
        let dq = UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_unchecked(d_psi / angle), angle);
        // R_fixed = exp([d_psi x]) * R_nominal
        state.attitude = dq * state.attitude;
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
        // Physical constraint: Tropospheric wet delay must be non-negative
        state.zwd = state.zwd.max(0.0);
    }
    
    for i in 0..state.ambiguities.len() {
        state.ambiguities[i] += dx[crate::filter::CORE_STATE_SIZE + i];
    }
    
    // Covariance update (Joseph Form)
    let identity = DMatrix::identity(state.covariance.nrows(), state.covariance.ncols());
    let i_kh = identity - &k * h_filt;
    let mut p_new = &i_kh * &state.covariance * i_kh.transpose() + &k * r_filt * k.transpose();
    
    // Force symmetry to prevent numerical divergence
    for r in 0..p_new.nrows() {
        for c in 0..r {
            let avg = (p_new[(r, c)] + p_new[(c, r)]) * 0.5;
            p_new[(r, c)] = avg;
            p_new[(c, r)] = avg;
        }
    }
    state.covariance = p_new;
    
    Ok(valid_indices)
}

pub fn apply_fix_and_hold(state: &mut RtkState, z_dd: &DVector<f64>, d_full: &DMatrix<f64>, var: f64) -> Result<(), UpdateError> {
    let num_dd = z_dd.len();
    if num_dd == 0 {
        return Ok(());
    }
    let state_size = d_full.ncols();
    
    let mut a_sd_full = DVector::zeros(state_size);
    for i in 0..state.ambiguities.len() {
        a_sd_full[crate::filter::CORE_STATE_SIZE + i] = state.ambiguities[i];
    }
    let a_dd = d_full * &a_sd_full;
    let v = z_dd - a_dd;
    
    let mut r = DMatrix::zeros(num_dd, num_dd);
    for i in 0..num_dd {
        r[(i, i)] = var;
    }
    
    let p = &state.covariance;
    let h_p = d_full * p;
    let s = &h_p * d_full.transpose() + &r;
    
    let s_inv = match s.try_inverse() {
        Some(inv) => inv,
        None => return Err(UpdateError::SingularMatrix),
    };
    
    let mut k = p * d_full.transpose() * s_inv;
    
    // DECOUPLE IMU STATES FROM AR SNAPS:
    if state.covariance.nrows() > 15 {
        for i in 6..15 {
            for j in 0..k.ncols() {
                k[(i, j)] = 0.0;
            }
        }
    }
    
    let dx = &k * &v;
    
    state.position.vector.x += dx[0];
    state.position.vector.y += dx[1];
    state.position.vector.z += dx[2];
    state.velocity.x += dx[3];
    state.velocity.y += dx[4];
    state.velocity.z += dx[5];
    
    let d_psi = Vector3::new(dx[6], dx[7], dx[8]);
    let angle = d_psi.norm();
    if angle > 1e-12 {
        let dq = nalgebra::UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_unchecked(d_psi / angle), angle);
        state.attitude = dq * state.attitude;
    }
    
    state.accel_bias.x += dx[9];
    state.accel_bias.y += dx[10];
    state.accel_bias.z += dx[11];
    state.gyro_bias.x += dx[12];
    state.gyro_bias.y += dx[13];
    state.gyro_bias.z += dx[14];
    
    state.rcv_clk_bias += dx[15];
    state.rcv_clk_drift += dx[16];
    state.zwd += dx[17];
    state.zwd = state.zwd.max(0.0);
    
    for i in 0..state.ambiguities.len() {
        state.ambiguities[i] += dx[crate::filter::CORE_STATE_SIZE + i];
    }
    
    let identity = DMatrix::identity(p.nrows(), p.ncols());
    let j = identity - &k * d_full;
    let mut p_new = &j * p * j.transpose() + &k * r * k.transpose();
    
    for r_idx in 0..p_new.nrows() {
        for c_idx in 0..r_idx {
            let avg = (p_new[(r_idx, c_idx)] + p_new[(c_idx, r_idx)]) * 0.5;
            p_new[(r_idx, c_idx)] = avg;
            p_new[(c_idx, r_idx)] = avg;
        }
    }
    state.covariance = p_new;
    
    Ok(())
}
