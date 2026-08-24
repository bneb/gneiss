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
    /// Wet-mapping difference (satellite − reference) at the rover:
    /// sensitivity of this DD to the rover ZWD residual state.
    pub dm_wet_rov: f64,
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

    // Non-dispersive wet residual: identical additive delay on code and
    // phase, mapped by the satellite/reference elevation difference.
    let zwd_idx = state.zwd_idx();
    let zwd_val = zwd_idx.map(|i| x_current[i]).unwrap_or(0.0);

    let mut pr_h = DVector::zeros(state_dim);
    pr_h[0] = d_geom_dpos.x;
    pr_h[1] = d_geom_dpos.y;
    pr_h[2] = d_geom_dpos.z;
    if let Some(zi) = zwd_idx {
        pr_h[zi] = m.dm_wet_rov;
    }
    h_rows.push(pr_h);
    y_vals.push(m.dd_pr_m - geom_dd - m.dm_wet_rov * zwd_val);
    r_diag.push(m.pr_var_m2.max(0.01));

    if let (Some(cp_obs), Some(amb_idx)) = (m.dd_cp_cycles, state.get_amb_idx(&m.key)) {
        let amb_val = x_current[amb_idx];
        let pred_cp = geom_dd / m.lambda + amb_val + m.dm_wet_rov * zwd_val / m.lambda;
        let mut cp_h = DVector::zeros(state_dim);
        cp_h[0] = d_geom_dpos.x / m.lambda;
        cp_h[1] = d_geom_dpos.y / m.lambda;
        cp_h[2] = d_geom_dpos.z / m.lambda;
        cp_h[amb_idx] = 1.0;
        if let Some(zi) = zwd_idx {
            cp_h[zi] = m.dm_wet_rov / m.lambda;
        }
        h_rows.push(cp_h);
        y_vals.push(cp_obs - pred_cp);
        r_diag.push(m.cp_var_cycles2.max(1e-4));
    }
}

/// Seed variance of the rover ZWD residual state (m^2): ~15 cm zenith.
pub const ZWD_INIT_VAR_M2: f64 = 0.0225;
/// Random-walk variance rate of the rover ZWD residual (m^2/s).
pub const ZWD_RW_M2_PER_S: f64 = 3e-7;
/// Per-epoch saturation bound on the ZWD correction (m). Wet zenith delay
/// moves at millimetres per second; an unconstrained step means the
/// innovation carried some other model error (slip, multipath burst).
pub const MAX_ZWD_STEP_M: f64 = 0.05;
// MEASURED OUTCOME (CORS day set): unbounded steps diverge (SLAC vertical
// RMS 15 m); with 0.05 m/epoch saturation vertical is STILL net-negative
// everywhere (P222 sm 0.75->1.43 m, SLAC sm 0.78->2.19 m) while horizontal
// improves marginally. The innovations carry non-tropo model error that a
// single rover-side state faithfully absorbs into the wrong bucket.
// Dormant until two-station tropo estimation exists.

/// Scalar random-walk Kalman update for the rover ZWD residual with a
/// per-epoch step saturation. `pairs` are `(h, y, r)` per phase
/// observation, where `h = dm_wet/lambda` is the sensitivity of predicted
/// phase (cycles) to zenith wet delay and `y` the innovation.
pub fn update_zwd_scalar(
    zwd: f64,
    var: f64,
    rw_m2_per_s: f64,
    dt_s: f64,
    pairs: &[(f64, f64, f64)],
) -> (f64, f64) {
    let prior_var = var + rw_m2_per_s * dt_s.abs().max(1e-3);
    let mut denom = 1.0 / prior_var;
    let mut num = 0.0;
    for &(h, y, r) in pairs {
        if r <= 0.0 {
            continue;
        }
        denom += h * h / r;
        num += h * y / r;
    }
    let post_var = 1.0 / denom;
    let delta = (post_var * num).clamp(-MAX_ZWD_STEP_M, MAX_ZWD_STEP_M);
    (zwd + delta, post_var)
}
// TUNING NOTE (CORS day set, full-day sweep): rates 3e-7..2.8e-6 monotonically
// improve SLAC vertical with tighter values, but even at 3e-7 SLAC smoothed
// vertical RMS stays ~2.2 m vs 0.78 m WITHOUT the state. A single rover-side
// ZWD cannot represent DD wet-delay error that includes a BASE-side residual
// (long baselines): both station residuals are needed, or none. State left
// dormant until two-station estimation is implemented.

/// Phase-innovation cycle-slip gate (cycles). A genuine DD phase residual
/// beyond this cannot come from orbit/model error; it means the arc's
/// integer assumption broke (undetected slip) and the pair must be
/// re-seeded instead of being fitted into the state.
pub const PHASE_INNOVATION_GATE_CYCLES: f64 = 500.0;

/// Keys whose DD phase innovation exceeds the slip gate under the current
/// state (evaluated pre-update, so the prediction is the prior one).
pub fn phase_innovation_outliers(
    state: &RtkState,
    measurements: &[DoubleDiffMeasurement],
    max_cycles: f64,
) -> Vec<DoubleDiffKey> {
    let mut out = Vec::new();
    for m in measurements {
        let Some(cp_obs) = m.dd_cp_cycles else { continue };
        let Some(amb_idx) = state.get_amb_idx(&m.key) else { continue };
        let cur_pos = state.pos_ecef;
        let r_sat = (m.sat_pos - cur_pos).norm();
        let r_ref = (m.ref_pos - cur_pos).norm();
        let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
        let geom_dd = (r_sat - r_ref) - base_dd
            + compute_tropo_dd(m.sat_pos, m.ref_pos, m.base_pos, cur_pos);
        let pred = geom_dd / m.lambda + state.to_dvector()[amb_idx];
        if (cp_obs - pred).abs() > max_cycles {
            out.push(m.key);
        }
    }
    out
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
            dm_wet_rov: 0.0,
        }];

        let res = iekf_update(&mut state, &meas);
        assert!(res.is_ok());
        assert!((state.pos_ecef.x - true_pos.x).abs() < 2.0);
    }

    #[test]
    fn test_zwd_state_tracks_wet_delay_ramp() {
        use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        let f1 = 1575.42e6_f64;
        let l1 = SPEED_OF_LIGHT_M_S / f1;
        let true_pos = Vector3::new(100.0, 200.0, 300.0);
        let base_pos = Vector3::new(0.0, 0.0, 0.0);
        let time = GpsTime::new(2000, 100.0);

        let mut state = RtkState::new(true_pos + Vector3::new(0.3, -0.2, 0.1), time);
        state.enable_zwd(ZWD_INIT_VAR_M2);

        // Four satellites at distinct elevations -> distinct wet mappings.
        let dirs = [
            (Vector3::new(20_000_000.0, 5_000_000.0, 8_000_000.0), 0.9),   // low elev
            (Vector3::new(-6_000_000.0, 18_000_000.0, 19_000_000.0), 0.55),
            (Vector3::new(9_000_000.0, -14_000_000.0, 21_000_000.0), 0.35),
            (Vector3::new(2_000_000.0, 7_000_000.0, 22_500_000.0), 0.15),  // near zenith
        ];
        let ref_pos = Vector3::new(5_000.0, 25_000.0, 20_000.0);
        let mut meas = Vec::new();
        for (p, (dir, dm)) in dirs.iter().enumerate() {
            let key = DoubleDiffKey { constellation_id: 0, sat: 2 + p as u16, ref_sat: 1, freq_band: 1 };
            state.ensure_ambiguity(key, 50.0 + p as f64, 100.0);
            meas.push(DoubleDiffMeasurement {
                key,
                dd_pr_m: 0.0,
                dd_cp_cycles: Some(0.0),
                sat_pos: true_pos + *dir * 20.0,
                ref_pos,
                base_pos,
                lambda: l1,
                pr_var_m2: 0.04,
                cp_var_cycles2: 1e-4,
                dm_wet_rov: *dm,
            });
        }

        let n = meas.len();
        for k in 0..40 {
            let zwd_true = 0.0025 * k as f64; // ramps 0 -> ~10 cm
            for (i, m) in meas.iter_mut().enumerate() {
                let pos = true_pos;
                let geom = (m.sat_pos - pos).norm() - (ref_pos - pos).norm()
                    - (m.sat_pos - base_pos).norm() + (ref_pos - base_pos).norm();
                m.dd_pr_m = geom + m.dm_wet_rov * zwd_true;
                m.dd_cp_cycles = Some(geom / l1 + 50.0 + i as f64 + m.dm_wet_rov * zwd_true / l1);
            }
            // Real flow: predict_state adds the random-walk process noise
            // each epoch (30 s sampling here).
            let zi = state.zwd_idx().unwrap();
            state.cov[(zi, zi)] += ZWD_RW_M2_PER_S * 30.0;
            iekf_update(&mut state, &meas).expect("update ok");
        }

        assert!((state.zwd_m - 0.0975).abs() < 0.04,
            "zwd must track the wet ramp, got {:.4} vs {}", state.zwd_m, 0.0975);
        let pos_err = (state.pos_ecef - true_pos).norm();
        assert!(pos_err < 0.2, "position must stay pinned while zwd absorbs drift: {:.3}", pos_err);
    }

}
