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
struct MeasGeom {
    geom_dd: f64,
    d_geom_dpos: Vector3<f64>,
    zwd_val: f64,
    grad_pr: f64,
    iono_val: f64,
    sat_iono_dd: f64,
    zwd_idx: Option<usize>,
    grad_idx: Option<(usize, usize)>,
    sat_iono_idxs: (Option<usize>, Option<usize>),
}

fn compute_meas_geom(
    m: &DoubleDiffMeasurement,
    cur_pos: Vector3<f64>,
    state: &RtkState,
    x: &DVector<f64>,
) -> MeasGeom {
    let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
    let r_sat = (m.sat_pos - cur_pos).norm();
    let r_ref = (m.ref_pos - cur_pos).norm();
    let trop_dd = compute_tropo_dd(m.sat_pos, m.ref_pos, m.base_pos, cur_pos);
    let geom_dd = (r_sat - r_ref) - base_dd + trop_dd + m.tide_dd_m;
    let d_geom_dpos = (m.ref_pos - cur_pos) / r_ref.max(1e-3) - (m.sat_pos - cur_pos) / r_sat.max(1e-3);
    let zwd_idx = state.zwd_idx();
    let zwd_val = zwd_idx.map(|i| x[i]).unwrap_or(0.0);
    let (grad_idx, gn_val, ge_val) = match state.grad_idx() {
        Some((gn, ge)) => (Some((gn, ge)), x[gn], x[ge]),
        None => (None, 0.0, 0.0),
    };
    let iono_val = state.get_iono_idx(&m.key).map(|ii| x[ii]).unwrap_or(0.0);
    let si_is = state.get_sat_iono_key_idx(m.key.constellation_id, m.key.sat).or_else(|| state.get_sat_iono_idx(m.key.sat));
    let si_ir = state.get_sat_iono_key_idx(m.key.constellation_id, m.key.ref_sat).or_else(|| state.get_sat_iono_idx(m.key.ref_sat));
    let sat_iono_dd = match (si_is, si_ir) {
        (Some(a), Some(b)) => x[a] - x[b],
        _ => 0.0,
    };
    let grad_pr = m.dgrad_n_rov * gn_val + m.dgrad_e_rov * ge_val;
    MeasGeom { geom_dd, d_geom_dpos, zwd_val, grad_pr, iono_val, sat_iono_dd, zwd_idx, grad_idx, sat_iono_idxs: (si_is, si_ir) }
}

#[allow(clippy::too_many_arguments)]
fn append_dd_code_row(
    m: &DoubleDiffMeasurement,
    g: &MeasGeom,
    state: &RtkState,
    state_dim: usize,
    h_rows: &mut Vec<DVector<f64>>,
    y_vals: &mut Vec<f64>,
    r_diag: &mut Vec<f64>,
    row_metas: &mut Vec<RowMeta>,
    gate_scale: f64,
) {
    let pr_y = m.dd_pr_m - g.geom_dd - m.dm_wet_rov * g.zwd_val - g.grad_pr - g.iono_val - g.sat_iono_dd;
    let pr_r = m.pr_var_m2.max(0.01);
    if pr_y.abs() > 30.0 && (pr_y * pr_y / pr_r) > 100.0 {
        return; // RAIM: exclude gross pseudorange multipath blunders
    }
    let mut pr_h = DVector::zeros(state_dim);
    pr_h[0] = g.d_geom_dpos.x;
    pr_h[1] = g.d_geom_dpos.y;
    pr_h[2] = g.d_geom_dpos.z;
    if let Some(zi) = g.zwd_idx { pr_h[zi] = m.dm_wet_rov; }
    if let Some((gn, ge)) = g.grad_idx { pr_h[gn] = m.dgrad_n_rov; pr_h[ge] = m.dgrad_e_rov; }
    if let Some(ii) = state.get_iono_idx(&m.key) { pr_h[ii] = 1.0; }
    h_rows.push(pr_h);
    y_vals.push(pr_y);
    r_diag.push(robust_inflate(pr_y, pr_r, gate_scale));
    row_metas.push(RowMeta {
        constellation_id: m.key.constellation_id,
        ref_sat: m.key.ref_sat,
        freq_band: m.key.freq_band,
        kind: ObsKind::Code,
        ref_var: m.pr_ref_var_m2,
    });
}

#[allow(clippy::too_many_arguments)]
fn append_dd_phase_row(
    m: &DoubleDiffMeasurement,
    g: &MeasGeom,
    state: &RtkState,
    x: &DVector<f64>,
    state_dim: usize,
    h_rows: &mut Vec<DVector<f64>>,
    y_vals: &mut Vec<f64>,
    r_diag: &mut Vec<f64>,
    row_metas: &mut Vec<RowMeta>,
    gate_scale: f64,
) {
    let (Some(cp_obs), Some(amb_idx)) = (pcv_corrected_cp(m), state.get_amb_idx(&m.key)) else { return };
    let pred_cp = g.geom_dd / m.lambda + x[amb_idx]
        + (m.dm_wet_rov * g.zwd_val + g.grad_pr) / m.lambda
        - g.iono_val / m.lambda
        - g.sat_iono_dd / m.lambda;
    let mut cp_h = DVector::zeros(state_dim);
    cp_h[0] = g.d_geom_dpos.x / m.lambda;
    cp_h[1] = g.d_geom_dpos.y / m.lambda;
    cp_h[2] = g.d_geom_dpos.z / m.lambda;
    cp_h[amb_idx] = 1.0;
    if let Some(ii) = state.get_iono_idx(&m.key) { cp_h[ii] = -1.0 / m.lambda; }
    if let (Some(a), Some(b)) = g.sat_iono_idxs { cp_h[a] -= 1.0 / m.lambda; cp_h[b] += 1.0 / m.lambda; }
    if let Some(zi) = g.zwd_idx { cp_h[zi] = m.dm_wet_rov / m.lambda; }
    if let Some((gn, ge)) = g.grad_idx { cp_h[gn] = m.dgrad_n_rov / m.lambda; cp_h[ge] = m.dgrad_e_rov / m.lambda; }
    h_rows.push(cp_h);
    let cp_y = cp_obs - pred_cp;
    let cp_r = m.cp_var_cycles2.max(1e-4);
    y_vals.push(cp_y);
    r_diag.push(robust_inflate(cp_y, cp_r, gate_scale));
    row_metas.push(RowMeta {
        constellation_id: m.key.constellation_id,
        ref_sat: m.key.ref_sat,
        freq_band: m.key.freq_band,
        kind: ObsKind::Phase,
        ref_var: m.cp_ref_var_cycles2,
    });
}

#[allow(clippy::too_many_arguments)]
fn append_dd_meas_rows(
    m: &DoubleDiffMeasurement,
    cur_pos: Vector3<f64>,
    state: &RtkState,
    x: &DVector<f64>,
    state_dim: usize,
    h_rows: &mut Vec<DVector<f64>>,
    y_vals: &mut Vec<f64>,
    r_diag: &mut Vec<f64>,
    row_metas: &mut Vec<RowMeta>,
    gate_scale: f64,
) {
    let geom = compute_meas_geom(m, cur_pos, state, x);
    append_dd_code_row(m, &geom, state, state_dim, h_rows, y_vals, r_diag, row_metas, gate_scale);
    append_dd_phase_row(m, &geom, state, x, state_dim, h_rows, y_vals, r_diag, row_metas, gate_scale);
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
