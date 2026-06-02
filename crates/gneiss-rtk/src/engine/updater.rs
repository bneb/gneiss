use crate::filter::RtkState;
use nalgebra::{DMatrix, DVector, Vector3, UnitQuaternion};

#[derive(Debug)]
pub enum UpdateError {
    SingularMatrix,
    DimensionMismatch,
    InvalidMeasurement,
}

pub fn update(state: &mut RtkState, z: &DVector<f64>, h: &DMatrix<f64>, r: &DMatrix<f64>) -> Result<(), UpdateError> {
    if z.len() != h.nrows() || h.ncols() != state.covariance.nrows() {
        return Err(UpdateError::DimensionMismatch);
    }
    
    let h_t_full = h.transpose();
    let s_full = h * &state.covariance * &h_t_full + r;
    
    let mut valid_indices = Vec::new();
    let gamma = 1e8; // Very loose validation gate (10000-sigma) to prevent 1600km divergence but allow normal GNSS corrections

    for i in 0..z.len() {
        if s_full[(i, i)] < 0.0 {
            tracing::error!("NEGATIVE VARIANCE! s_ii: {:.4}", s_full[(i, i)]);
        }
        // Enforce minimum variance bound on S_ii to prevent gain collapse
        let s_ii = s_full[(i, i)].max(1e-4);
        let nu = z[i];
        
        if nu * nu / s_ii < gamma {
            valid_indices.push(i);
        } else {
            tracing::warn!("EKF rejected measurement! nu={:.2}, s_ii={:.2}", nu, s_ii);
        }
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

    let hp = &h_filt * &state.covariance;
    let s = &hp * h_filt.transpose() + &r_filt;
    let s_inv = nalgebra::linalg::SVD::new(s.clone(), true, true).pseudo_inverse(1e-12)
        .unwrap_or_else(|_| s.clone().try_inverse().unwrap_or_else(|| DMatrix::identity(s.nrows(), s.ncols())));
    let k = &state.covariance * h_filt.transpose() * &s_inv;
    let dx = &k * &z_filt;
    
    if dx.iter().any(|x| x.is_nan()) {
        tracing::error!("EKF NaN EXPLOSION! Cov max: {}", state.covariance.abs().max());
        return Err(UpdateError::SingularMatrix);
    }

    let pos_change = (dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]).sqrt();
    if pos_change > 100.0 {
        tracing::warn!("EKF LARGE JUMP! dx_pos: {:.1}m, z max: {:.1}", pos_change, z.abs().max());
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
    
    // Force perfect symmetry to prevent numerical divergence
    for r in 0..p_new.nrows() {
        for c in 0..r {
            let avg = (p_new[(r, c)] + p_new[(c, r)]) * 0.5;
            p_new[(r, c)] = avg;
            p_new[(c, r)] = avg;
        }
    }
    state.covariance = p_new;
    
    Ok(())
}
