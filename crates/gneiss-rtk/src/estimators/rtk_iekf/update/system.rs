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

    let geoms: Vec<MeasGeom> = measurements
        .iter()
        .map(|m| compute_meas_geom(m, cur_pos, state, x_current))
        .collect();
    let rel = code_robust_centre_scale(measurements, &geoms);

    for (m, geom) in measurements.iter().zip(geoms.iter()) {
        append_dd_code_row(m, geom, state, state_dim, &mut h_rows, &mut y_vals, &mut r_diag, &mut row_metas, innov_gate_scale, rel);
        append_dd_phase_row(m, geom, state, x_current, state_dim, &mut h_rows, &mut y_vals, &mut r_diag, &mut row_metas, innov_gate_scale);
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

#[derive(Clone, Copy)]
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

/// Canonical double-difference pseudorange innovation (metres).
///
/// Single definition shared by the measurement row and the per-epoch robust
/// gate, so the gate can never score a different quantity than the row does.
fn code_innovation(m: &DoubleDiffMeasurement, g: &MeasGeom) -> f64 {
    m.dd_pr_m - g.geom_dd - m.dm_wet_rov * g.zwd_val - g.grad_pr - g.iono_val - g.sat_iono_dd
}

/// Minimum absolute code-innovation deviation (m) before the per-epoch
/// relative gate may reject a pair. NLOS bias below this is not worth
/// discarding a usable observation over.
pub const CODE_REL_MIN_DEV_M: f64 = 1.0;

/// Robust-sigma multiplier (k) for the per-epoch relative code gate.
pub const CODE_REL_K: f64 = 3.0;

/// Floor (metres) on the robust scale, so a degenerate MAD -- when most pairs
/// share an identical residual -- cannot silently disable the gate. Kept far
/// below [`CODE_REL_MIN_DEV_M`] so it never changes behaviour on real data.
pub const CODE_REL_MIN_SCALE_M: f64 = 0.05;

/// Median of a slice (0.0 when empty).
fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut v = values.to_vec();
    v.sort_by(|a, b| a.total_cmp(b));
    let mid = v.len() / 2;
    if v.len().is_multiple_of(2) { 0.5 * (v[mid - 1] + v[mid]) } else { v[mid] }
}

/// Per-epoch robust centre and scale of the code innovations.
///
/// NLOS bias is *relative* to the epoch: most pairs stay clean while a few
/// satellites return tens of metres of reflected signal. A median/MAD scale
/// therefore detects the corrupted minority without punishing epochs where the
/// filter as a whole is displaced (where every innovation shifts together and
/// the median moves with them).
fn code_robust_centre_scale(
    measurements: &[DoubleDiffMeasurement],
    geoms: &[MeasGeom],
) -> (f64, f64) {
    let ys: Vec<f64> = measurements.iter().zip(geoms).map(|(m, g)| code_innovation(m, g)).collect();
    if ys.len() < MIN_CODE_ROWS_FOR_REL_GATE { return (0.0, 0.0); }
    let med = median(&ys);
    let devs: Vec<f64> = ys.iter().map(|y| (y - med).abs()).collect();
    (med, (1.4826 * median(&devs)).max(CODE_REL_MIN_SCALE_M))
}

/// Smallest epoch population for which a relative gate is meaningful.
const MIN_CODE_ROWS_FOR_REL_GATE: usize = 6;

/// Reject a code pair that is an outlier against its own epoch.
fn is_code_epoch_outlier(pr_y: f64, centre: f64, scale: f64) -> bool {
    if scale <= 0.0 { return false; }
    let dev = (pr_y - centre).abs();
    dev > (CODE_REL_K * scale).max(CODE_REL_MIN_DEV_M)
}

fn is_code_blunder(m: &DoubleDiffMeasurement, pr_y: f64, pr_r: f64) -> bool {
    let nis = pr_y * pr_y / pr_r;
    if m.dd_cp_cycles.is_none() {
        pr_y.abs() > 15.0 && nis > 36.0
    } else {
        pr_y.abs() > 10.0 && nis > 25.0
    }
}

fn effective_code_variance(m: &DoubleDiffMeasurement, pr_y: f64) -> f64 {
    let base_r = m.pr_var_m2.max(0.01);
    if m.dd_cp_cycles.is_some() && pr_y.abs() > 3.0 {
        let excess = pr_y.abs() - 3.0;
        base_r * (1.0 + excess * excess)
    } else {
        base_r
    }
}

#[allow(clippy::too_many_arguments)]
fn append_dd_code_row(
    m: &DoubleDiffMeasurement, g: &MeasGeom, state: &RtkState, state_dim: usize,
    h_rows: &mut Vec<DVector<f64>>, y_vals: &mut Vec<f64>, r_diag: &mut Vec<f64>,
    row_metas: &mut Vec<RowMeta>, gate_scale: f64, rel: (f64, f64),
) {
    let pr_y = code_innovation(m, g);
    // NLOS-corrupted code is biased relative to the rest of the epoch, not just
    // large in absolute terms. Measured on Whampoa Survey, mean code innovation
    // reaches 5.0 m at the worst epochs while the absolute RAIM gate (10-15 m)
    // rejects under 1 pair; the relative gate is what actually fires.
    let pr_r = effective_code_variance(m, pr_y);
    if is_code_epoch_outlier(pr_y, rel.0, rel.1) || is_code_blunder(m, pr_y, pr_r) {
        return; // RAIM: exclude gross pseudorange multipath / NLOS blunders
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::estimators::rtk_iekf::state::DoubleDiffKey;
    use crate::estimators::rtk_iekf::update::DoubleDiffMeasurement;
    use nalgebra::Vector3;

    /// A measurement whose only relevant field is its pseudorange value.
    fn meas(dd_pr_m: f64) -> DoubleDiffMeasurement {
        DoubleDiffMeasurement {
            key: DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 },
            dd_pr_m,
            dd_cp_cycles: None,
            sat_pos: Vector3::new(2.0e7, 1.0e7, 1.5e7),
            ref_pos: Vector3::new(-2.0e7, 1.5e7, 1.0e7),
            base_pos: Vector3::zeros(),
            lambda: 0.19,
            pr_var_m2: 0.05,
            cp_var_cycles2: 0.0,
            pr_ref_var_m2: 0.02,
            cp_ref_var_cycles2: 0.0,
            dm_wet_rov: 0.0,
            dgrad_n_rov: 0.0,
            dgrad_e_rov: 0.0,
            tide_dd_m: 0.0,
            dd_pcv_m: 0.0,
        }
    }

    /// Geometry with `geom_dd = 0`, so innovation equals `dd_pr_m`.
    fn geom() -> MeasGeom {
        MeasGeom {
            geom_dd: 0.0,
            d_geom_dpos: Vector3::new(0.1, 0.2, 0.3),
            zwd_val: 0.0,
            grad_pr: 0.0,
            iono_val: 0.0,
            sat_iono_dd: 0.0,
            zwd_idx: None,
            grad_idx: None,
            sat_iono_idxs: (None, None),
        }
    }

    #[test]
    fn median_odd_and_even_counts() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]), 2.5);
        assert_eq!(median(&[]), 0.0);
    }

    #[test]
    fn code_innovation_is_the_canonical_residual() {
        // The row builder and the gate must score the same quantity.
        assert_eq!(code_innovation(&meas(2.5), &geom()), 2.5);
    }

    #[test]
    fn gate_ignores_small_epochs() {
        let ms = vec![meas(0.1), meas(0.2)];
        let gs = vec![geom(), geom()];
        assert_eq!(code_robust_centre_scale(&ms, &gs), (0.0, 0.0));
        assert!(!is_code_epoch_outlier(5.0, 0.0, 0.0));
    }

    #[test]
    fn gate_flags_a_single_nlos_outlier_among_clean_pairs() {
        // 11 clean pairs near zero, one NLOS pair 12 m out.
        let mut ms: Vec<DoubleDiffMeasurement> = (0..11).map(|_| meas(0.1)).collect();
        ms.push(meas(12.0));
        let gs = vec![geom(); ms.len()];
        let (centre, scale) = code_robust_centre_scale(&ms, &gs);
        assert!((centre - 0.1).abs() < 1e-9, "median must ignore the outlier");
        assert!(scale > 0.0);
        assert!(is_code_epoch_outlier(12.0, centre, scale), "NLOS pair must be flagged");
        assert!(!is_code_epoch_outlier(0.1, centre, scale), "clean pair must survive");
    }

    #[test]
    fn gate_does_not_flag_a_whole_epoch_displacement() {
        // Every pair is 20 m from the model (filter displaced, data clean).
        // The median follows the displacement, so nothing is flagged.
        let ms: Vec<DoubleDiffMeasurement> = (0..10).map(|_| meas(20.0)).collect();
        let gs = vec![geom(); ms.len()];
        let (centre, scale) = code_robust_centre_scale(&ms, &gs);
        assert!((centre - 20.0).abs() < 1e-9);
        assert!(!is_code_epoch_outlier(20.0, centre, scale));
    }

    #[test]
    fn gate_respects_minimum_absolute_deviation() {
        // Tiny spread with a sub-metre difference must not trigger rejection.
        assert!(!is_code_epoch_outlier(0.4, 0.0, 0.01));
        assert!(is_code_epoch_outlier(3.0, 0.0, 0.01));
    }
}
