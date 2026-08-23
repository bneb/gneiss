//! Iterated Extended Kalman Filter (IEKF) Double-Difference Measurement Update.

use nalgebra::{DMatrix, DVector, Vector3};
use super::state::{DoubleDiffKey, RtkState};

/// Double-difference observation for a satellite pair on a single frequency band.
#[derive(Debug, Clone)]
pub struct DoubleDiffMeasurement {
    pub key: DoubleDiffKey,
    pub dd_pr_m: f64,
    pub dd_cp_cycles: Option<f64>,
    pub sat_pos: Vector3<f64>,
    pub ref_pos: Vector3<f64>,
    pub base_pos: Vector3<f64>,
    pub lambda: f64,
    pub pr_var_m2: f64,
    pub cp_var_cycles2: f64,
}

/// Perform Iterated Extended Kalman Filter (IEKF) update with double-differenced measurements.
pub fn iekf_update(state: &mut RtkState, measurements: &[DoubleDiffMeasurement]) -> Result<f64, String> {
    if measurements.is_empty() {
        return Ok(0.0);
    }

    let p0 = state.cov.clone();
    let x0 = state.to_dvector();
    let mut x = x0.clone();
    let mut last_rms = 0.0;

    for _iter in 0..4 {
        let (h, y, r) = build_measurement_system(state, &x, measurements);
        if y.is_empty() {
            return Ok(0.0);
        }
        let delta_x = match compute_iekf_step(&p0, &h, &y, &r, &x, &x0) {
            Some(dx) => dx,
            None => return Err("Innovation covariance matrix inversion failed".to_string()),
        };
        x = &x0 + &delta_x;
        last_rms = delta_x.rows(0, 3).norm();
        if last_rms < 1e-4 {
            break;
        }
    }

    apply_joseph_update(state, &p0, &x, measurements);
    state.update_from_dvector(&x);
    Ok(last_rms)
}

fn compute_iekf_step(
    p0: &DMatrix<f64>,
    h: &DMatrix<f64>,
    y: &DVector<f64>,
    r: &DMatrix<f64>,
    x: &DVector<f64>,
    x0: &DVector<f64>,
) -> Option<DVector<f64>> {
    let s = h * p0 * h.transpose() + r;
    let s_inv = s.try_inverse()?;
    let k = p0 * h.transpose() * s_inv;
    let dx_iter = x - x0;
    Some(&k * (y + h * dx_iter))
}

fn apply_joseph_update(
    state: &mut RtkState,
    p0: &DMatrix<f64>,
    x: &DVector<f64>,
    measurements: &[DoubleDiffMeasurement],
) {
    let (h_final, y_final, r_final) = build_measurement_system(state, x, measurements);
    if y_final.is_empty() {
        return;
    }
    let s = &h_final * p0 * h_final.transpose() + &r_final;
    if let Some(s_inv) = s.try_inverse() {
        let k = p0 * h_final.transpose() * s_inv;
        let i_kh = DMatrix::identity(state.dim(), state.dim()) - &k * &h_final;
        state.cov = &i_kh * p0 * i_kh.transpose() + &k * &r_final * k.transpose();
    }
}

/// Build measurement Jacobians, innovation residuals, and covariance matrix R.
fn build_measurement_system(
    state: &RtkState,
    x_current: &DVector<f64>,
    measurements: &[DoubleDiffMeasurement],
) -> (DMatrix<f64>, DVector<f64>, DMatrix<f64>) {
    let mut h_rows = Vec::new();
    let mut y_vals = Vec::new();
    let mut r_diag = Vec::new();

    let cur_pos = Vector3::new(x_current[0], x_current[1], x_current[2]);
    let state_dim = state.dim();

    for m in measurements {
        append_dd_meas_rows(m, cur_pos, state, x_current, state_dim, &mut h_rows, &mut y_vals, &mut r_diag);
    }

    assemble_matrices(h_rows, y_vals, r_diag, state_dim)
}

/// Double-difference troposphere delay (Saastamoinen, RTKLIB coefficients).
pub(crate) fn compute_tropo_dd(
    sat_pos: Vector3<f64>,
    ref_pos: Vector3<f64>,
    base_pos: Vector3<f64>,
    rx_pos: Vector3<f64>,
) -> f64 {
    let rx_llh = gneiss_core::coords::ecef_to_llh(rx_pos);
    let bs_llh = gneiss_core::coords::ecef_to_llh(base_pos);
    let params = gneiss_core::atmosphere::TropoParams::default();
    let t_rs = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&params, rx_llh, gneiss_core::coords::az_el(rx_llh, rx_pos, sat_pos).1);
    let t_rr = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&params, rx_llh, gneiss_core::coords::az_el(rx_llh, rx_pos, ref_pos).1);
    let t_bs = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&params, bs_llh, gneiss_core::coords::az_el(bs_llh, base_pos, sat_pos).1);
    let t_br = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&params, bs_llh, gneiss_core::coords::az_el(bs_llh, base_pos, ref_pos).1);
    (t_rs - t_rr) - (t_bs - t_br)
}

#[allow(clippy::too_many_arguments)]
fn append_dd_meas_rows(
    m: &DoubleDiffMeasurement,
    cur_pos: Vector3<f64>,
    state: &RtkState,
    x_current: &DVector<f64>,
    state_dim: usize,
    h_rows: &mut Vec<DVector<f64>>,
    y_vals: &mut Vec<f64>,
    r_diag: &mut Vec<f64>,
) {
    let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
    let r_sat = (m.sat_pos - cur_pos).norm();
    let r_ref = (m.ref_pos - cur_pos).norm();
    let trop_dd = compute_tropo_dd(m.sat_pos, m.ref_pos, m.base_pos, cur_pos);

    let geom_dd = (r_sat - r_ref) - base_dd + trop_dd;

    let los_sat = (m.sat_pos - cur_pos) / r_sat.max(1e-3);
    let los_ref = (m.ref_pos - cur_pos) / r_ref.max(1e-3);
    let d_geom_dpos = los_ref - los_sat;

    let mut pr_h = DVector::zeros(state_dim);
    pr_h[0] = d_geom_dpos.x;
    pr_h[1] = d_geom_dpos.y;
    pr_h[2] = d_geom_dpos.z;
    h_rows.push(pr_h);
    y_vals.push(m.dd_pr_m - geom_dd);
    r_diag.push(m.pr_var_m2.max(0.01));

    if let (Some(cp_obs), Some(amb_idx)) = (m.dd_cp_cycles, state.get_amb_idx(&m.key)) {
        let amb_val = x_current[amb_idx];
        let pred_cp = geom_dd / m.lambda + amb_val;
        let mut cp_h = DVector::zeros(state_dim);
        cp_h[0] = d_geom_dpos.x / m.lambda;
        cp_h[1] = d_geom_dpos.y / m.lambda;
        cp_h[2] = d_geom_dpos.z / m.lambda;
        cp_h[amb_idx] = 1.0;
        h_rows.push(cp_h);
        y_vals.push(cp_obs - pred_cp);
        r_diag.push(m.cp_var_cycles2.max(1e-4));
    }
}

fn assemble_matrices(
    h_rows: Vec<DVector<f64>>,
    y_vals: Vec<f64>,
    r_diag: Vec<f64>,
    state_dim: usize,
) -> (DMatrix<f64>, DVector<f64>, DMatrix<f64>) {
    let n_meas = h_rows.len();
    let mut h = DMatrix::zeros(n_meas, state_dim);
    let mut y = DVector::zeros(n_meas);
    let mut r = DMatrix::zeros(n_meas, n_meas);

    for (i, row) in h_rows.into_iter().enumerate() {
        for c in 0..state_dim {
            h[(i, c)] = row[c];
        }
        y[i] = y_vals[i];
        r[(i, i)] = r_diag[i];
    }
    (h, y, r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::time::GpsTime;

    #[test]
    fn test_iekf_update_reduces_position_error() {
        let true_pos = Vector3::new(100.0, 200.0, 300.0);
        let base_pos = Vector3::new(0.0, 0.0, 0.0);
        let time = GpsTime::new(2000, 100.0);

        let mut state = RtkState::new(true_pos + Vector3::new(2.0, -2.0, 1.0), time);
        let key = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        state.ensure_ambiguity(key, 0.0, 100.0);

        let sat_pos = Vector3::new(10_000.0, 20_000.0, 20_000.0);
        let ref_pos = Vector3::new(5_000.0, 25_000.0, 20_000.0);
        let base_dd = (sat_pos - base_pos).norm() - (ref_pos - base_pos).norm();
        let true_dd = (sat_pos - true_pos).norm() - (ref_pos - true_pos).norm() - base_dd;

        let meas = vec![DoubleDiffMeasurement {
            key,
            dd_pr_m: true_dd,
            dd_cp_cycles: Some(true_dd / 0.190),
            sat_pos,
            ref_pos,
            base_pos,
            lambda: 0.190,
            pr_var_m2: 0.04,
            cp_var_cycles2: 0.0001,
        }];

        let res = iekf_update(&mut state, &meas);
        assert!(res.is_ok());
        assert!((state.pos_ecef.x - true_pos.x).abs() < 2.0);
    }
}
