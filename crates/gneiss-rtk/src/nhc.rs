use crate::filter::RtkState;
use crate::engine::updater;
use nalgebra::{DMatrix, DVector, Matrix3};

/// Applies Non-Holonomic Constraints (NHC) to the EKF state.
/// This assumes the vehicle's lateral and vertical velocity in the body frame is zero.
pub fn apply_nhc(state: &mut RtkState, sigma_lateral: f64, sigma_vertical: f64) -> Result<(), &'static str> {
    let r_e_b = state.attitude.inverse().to_rotation_matrix();
    let v_b = r_e_b * state.velocity;
    
    // Measurement z: we want v_y and v_z to be 0
    let z = DVector::from_column_slice(&[-v_b.y, -v_b.z]);
    
    let n = state.covariance.nrows();
    let mut h = DMatrix::<f64>::zeros(2, n);
    
    // Jacobian w.r.t velocity in ECEF
    // v_b = R_e_b * v_e => d(v_b)/dv_e = R_e_b
    let dr_dv = r_e_b.matrix();
    h[(0, 3)] = dr_dv[(1, 0)]; h[(0, 4)] = dr_dv[(1, 1)]; h[(0, 5)] = dr_dv[(1, 2)];
    h[(1, 3)] = dr_dv[(2, 0)]; h[(1, 4)] = dr_dv[(2, 1)]; h[(1, 5)] = dr_dv[(2, 2)];
    
    // Jacobian w.r.t attitude
    // d(v_b)/d(psi) = [v_b x] * R_e^b
    let v_b_skew = Matrix3::new(
        0.0, -v_b.z,  v_b.y,
        v_b.z,  0.0, -v_b.x,
       -v_b.y,  v_b.x,  0.0
    );
    let dr_dpsi = v_b_skew * r_e_b.matrix();
    h[(0, 6)] = dr_dpsi[(1, 0)]; h[(0, 7)] = dr_dpsi[(1, 1)]; h[(0, 8)] = dr_dpsi[(1, 2)];
    h[(1, 6)] = dr_dpsi[(2, 0)]; h[(1, 7)] = dr_dpsi[(2, 1)]; h[(1, 8)] = dr_dpsi[(2, 2)];

    let r = DMatrix::from_diagonal(&DVector::from_column_slice(&[
        sigma_lateral * sigma_lateral,
        sigma_vertical * sigma_vertical
    ]));

    updater::update(state, &z, &h, &r, 1e9, None).map_err(|_| "NHC update failed")?;
    Ok(())
}

pub fn apply_zupt(state: &mut RtkState, sigma: f64) -> Result<(), &'static str> {
    let z = DVector::from_column_slice(&[-state.velocity.x, -state.velocity.y, -state.velocity.z]);
    let n = state.covariance.nrows();
    let mut h = DMatrix::<f64>::zeros(3, n);
    for i in 0..3 { h[(i, 3 + i)] = 1.0; }
    
    let r = DMatrix::from_diagonal(&DVector::from_element(3, sigma * sigma));
    updater::update(state, &z, &h, &r, 1e9, None).map_err(|_| "ZUPT update failed")?;
    Ok(())}
