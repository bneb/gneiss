//! Measurement system matrix construction for double-differenced observations.

use nalgebra::{DMatrix, DVector, Vector3};
use crate::estimators::rtk_iekf::state::RtkState;
use super::robust::robust_inflate;
use super::{pcv_corrected_cp, DoubleDiffMeasurement};

#[derive(Clone, Copy, PartialEq, Eq)]
enum ObsKind {
    Code,
    Phase,
}

struct RowMeta {
    constellation_id: u8,
    ref_sat: u16,
    freq_band: u8,
    kind: ObsKind,
    ref_var: f64,
}

/// Build measurement Jacobians, innovation residuals, and covariance matrix R.
pub fn build_measurement_system(
    state: &RtkState,
    x_current: &DVector<f64>,
    measurements: &[DoubleDiffMeasurement],
    innov_gate_scale: f64,
) -> (DMatrix<f64>, DVector<f64>, DMatrix<f64>) {
    let mut h_rows = Vec::new();
    let mut y_vals = Vec::new();
    let mut r_diag = Vec::new();
    let mut row_metas = Vec::new();

    let cur_pos = Vector3::new(x_current[0], x_current[1], x_current[2]);
    let state_dim = state.dim();

    for m in measurements {
        append_dd_meas_rows(
            m, cur_pos, state, x_current, state_dim,
            &mut h_rows, &mut y_vals, &mut r_diag, &mut row_metas, innov_gate_scale,
        );
    }

    assemble_matrices(h_rows, y_vals, r_diag, &row_metas, state_dim)
}

/// Double-difference troposphere delay (Saastamoinen, RTKLIB coefficients).
pub fn compute_tropo_dd(
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
    row_metas: &mut Vec<RowMeta>,
    innov_gate_scale: f64,
) {
    let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
    let r_sat = (m.sat_pos - cur_pos).norm();
    let r_ref = (m.ref_pos - cur_pos).norm();
    let trop_dd = compute_tropo_dd(m.sat_pos, m.ref_pos, m.base_pos, cur_pos);
    let geom_dd = (r_sat - r_ref) - base_dd + trop_dd + m.tide_dd_m;

    let los_sat = (m.sat_pos - cur_pos) / r_sat.max(1e-3);
    let los_ref = (m.ref_pos - cur_pos) / r_ref.max(1e-3);
    let d_geom_dpos = los_ref - los_sat;

    let zwd_idx = state.zwd_idx();
    let zwd_val = zwd_idx.map(|i| x_current[i]).unwrap_or(0.0);
    let (grad_idx, grad_n_val, grad_e_val) = match state.grad_idx() {
        Some((gn, ge)) => (Some((gn, ge)), x_current[gn], x_current[ge]),
        None => (None, 0.0, 0.0),
    };
    let iono_val = state.get_iono_idx(&m.key)
        .map(|ii| x_current[ii])
        .unwrap_or(0.0);

    let si_is = state.get_sat_iono_key_idx(m.key.constellation_id, m.key.sat)
        .or_else(|| state.get_sat_iono_idx(m.key.sat));
    let si_ir = state.get_sat_iono_key_idx(m.key.constellation_id, m.key.ref_sat)
        .or_else(|| state.get_sat_iono_idx(m.key.ref_sat));
    let sat_iono_dd = match (si_is, si_ir) {
        (Some(a), Some(b)) => x_current[a] - x_current[b],
        _ => 0.0,
    };

    let mut pr_h = DVector::zeros(state_dim);
    pr_h[0] = d_geom_dpos.x;
    pr_h[1] = d_geom_dpos.y;
    pr_h[2] = d_geom_dpos.z;
    if let Some(zi) = zwd_idx {
        pr_h[zi] = m.dm_wet_rov;
    }
    if let Some((gn, ge)) = grad_idx {
        pr_h[gn] = m.dgrad_n_rov;
        pr_h[ge] = m.dgrad_e_rov;
    }
    if let Some(ii) = state.get_iono_idx(&m.key) {
        pr_h[ii] = 1.0;
    }
    h_rows.push(pr_h);
    let grad_pr = m.dgrad_n_rov * grad_n_val + m.dgrad_e_rov * grad_e_val;
    let pr_y = m.dd_pr_m - geom_dd - m.dm_wet_rov * zwd_val - grad_pr - iono_val - sat_iono_dd;
    let pr_r = m.pr_var_m2.max(0.01);
    y_vals.push(pr_y);
    r_diag.push(robust_inflate(pr_y, pr_r, innov_gate_scale));
    row_metas.push(RowMeta {
        constellation_id: m.key.constellation_id,
        ref_sat: m.key.ref_sat,
        freq_band: m.key.freq_band,
        kind: ObsKind::Code,
        ref_var: m.pr_ref_var_m2,
    });

    if let (Some(cp_obs), Some(amb_idx)) = (pcv_corrected_cp(m), state.get_amb_idx(&m.key)) {
        let amb_val = x_current[amb_idx];
        let pred_cp = geom_dd / m.lambda + amb_val
            + (m.dm_wet_rov * zwd_val + grad_pr) / m.lambda
            - iono_val / m.lambda
            - sat_iono_dd / m.lambda;
        let mut cp_h = DVector::zeros(state_dim);
        cp_h[0] = d_geom_dpos.x / m.lambda;
        cp_h[1] = d_geom_dpos.y / m.lambda;
        cp_h[2] = d_geom_dpos.z / m.lambda;
        cp_h[amb_idx] = 1.0;
        if let Some(ii) = state.get_iono_idx(&m.key) {
            cp_h[ii] = -1.0 / m.lambda;
        }
        if let (Some(a), Some(b)) = (si_is, si_ir) {
            cp_h[a] -= 1.0 / m.lambda;
            cp_h[b] += 1.0 / m.lambda;
        }
        if let Some(zi) = zwd_idx {
            cp_h[zi] = m.dm_wet_rov / m.lambda;
        }
        if let Some((gn, ge)) = grad_idx {
            cp_h[gn] = m.dgrad_n_rov / m.lambda;
            cp_h[ge] = m.dgrad_e_rov / m.lambda;
        }
        h_rows.push(cp_h);
        let cp_y = cp_obs - pred_cp;
        let cp_r = m.cp_var_cycles2.max(1e-4);
        y_vals.push(cp_y);
        r_diag.push(robust_inflate(cp_y, cp_r, innov_gate_scale));
        row_metas.push(RowMeta {
            constellation_id: m.key.constellation_id,
            ref_sat: m.key.ref_sat,
            freq_band: m.key.freq_band,
            kind: ObsKind::Phase,
            ref_var: m.cp_ref_var_cycles2,
        });
    }
}

fn assemble_matrices(
    h_rows: Vec<DVector<f64>>,
    y_vals: Vec<f64>,
    r_diag: Vec<f64>,
    metas: &[RowMeta],
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
    fill_dd_covariances(&mut r, metas);
    (h, y, r)
}

fn fill_dd_covariances(r: &mut DMatrix<f64>, metas: &[RowMeta]) {
    for (i, mi) in metas.iter().enumerate() {
        for (j, mj) in metas.iter().enumerate().skip(i + 1) {
            let same_group = mi.kind == mj.kind
                && mi.constellation_id == mj.constellation_id
                && mi.ref_sat == mj.ref_sat
                && mi.freq_band == mj.freq_band;
            if same_group {
                let cov = mi.ref_var.min(mj.ref_var);
                r[(i, j)] = cov;
                r[(j, i)] = cov;
            }
        }
    }
}
