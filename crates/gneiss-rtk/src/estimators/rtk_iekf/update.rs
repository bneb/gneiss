//! Iterated Extended Kalman Filter (IEKF) Double-Difference Measurement Update.

use std::collections::HashMap;

use nalgebra::{DMatrix, DVector, Vector3};
use super::state::{DoubleDiffKey, RtkState};

/// Robust estimation: inflate measurement variance when its normalized
/// innovation squared exceeds this threshold (~3-sigma point). Below it,
/// measurements keep nominal weight; above it, weight decays as 1/innov².
/// This limits multipath and unmodeled-atmosphere corruption without the
/// information loss of hard exclusion.
pub const ROBUST_INNOVATION_THRESHOLD: f64 = 9.0;

/// Gradient state seed variance (m^2): ~2 mm of horizon slant delay each.
pub const GRAD_INIT_VAR_M2: f64 = 4.0e-6;
/// Gradient random-walk rate (m^2/s): the field is mm-scale and drifts
/// slowly; this keeps session-long sigma near the seed so gradients track
/// weather, never noise.
pub const GRAD_RW_M2_PER_S: f64 = 5.0e-11;
/// Elevation floor (rad) inside the cot(el) gradient mapping, bounding
/// low-elevation leverage like the variance model's sin floor does.
pub(crate) const GRAD_MIN_SIN_EL: f64 = 0.17;

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
    /// Gradient mapping differences (satellite − reference) at the rover:
    /// [north, east] sensitivities cot(el)·{cos az, sin az} differenced
    /// between the pair members. Zero-mean under no gradient.
    pub dgrad_n_rov: f64,
    pub dgrad_e_rov: f64,
    /// Solid Earth tide DD correction (metres, LOS-projected).
    pub tide_dd_m: f64,
    /// Differential receiver-antenna PCV embedded in this DD pair
    /// (metres): `[PCV_rov(z_sat) - PCV_rov(z_ref)] -
    /// [PCV_base(z_sat) - PCV_base(z_ref)]`. Subtracted from the carrier
    /// phase by [`pcv_corrected_cp`] before any model comparison; zero
    /// unless receiver calibrations are loaded and `GNEISS_RECV_PCV=1`.
    pub dd_pcv_m: f64,
}

/// DD phase observation with the embedded differential receiver-PCV
/// signature removed (`cp - dd_pcv_m / lambda` cycles), matching the
/// phase-windup convention. All consumers of [`DoubleDiffMeasurement::
/// dd_cp_cycles`] must go through here so the correction is applied
/// exactly once per comparison.
pub(crate) fn pcv_corrected_cp(m: &DoubleDiffMeasurement) -> Option<f64> {
    m.dd_cp_cycles.map(|cp| cp - m.dd_pcv_m / m.lambda)
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
    let geom_dd = (r_sat - r_ref) - base_dd + trop_dd + m.tide_dd_m;

    let los_sat = (m.sat_pos - cur_pos) / r_sat.max(1e-3);
    let los_ref = (m.ref_pos - cur_pos) / r_ref.max(1e-3);
    let d_geom_dpos = los_ref - los_sat;

    // Non-dispersive wet residual: identical additive delay on code and
    // phase, mapped by the satellite/reference elevation difference.
    let zwd_idx = state.zwd_idx();
    let zwd_val = zwd_idx.map(|i| x_current[i]).unwrap_or(0.0);
    let (grad_idx, grad_n_val, grad_e_val) = match state.grad_idx() {
        Some((gn, ge)) => (Some((gn, ge)), x_current[gn], x_current[ge]),
        None => (None, 0.0, 0.0),
    };
    let iono_val = state.get_iono_idx(&m.key)
        .map(|ii| x_current[ii])
        .unwrap_or(0.0);

    let mut pr_h = DVector::zeros(state_dim);
    pr_h[0] = d_geom_dpos.x;
    pr_h[1] = d_geom_dpos.y;
    pr_h[2] = d_geom_dpos.z;
    if let Some(zi) = zwd_idx {
        pr_h[zi] = m.dm_wet_rov;
    }
    if let Some((gn, ge)) = grad_idx {
        // Slant delay = cot(el) * (gN cos az + gE sin az), differenced
        // satellite-minus-reference like the zenith mapping.
        pr_h[gn] = m.dgrad_n_rov;
        pr_h[ge] = m.dgrad_e_rov;
    }
    // Iono state: DD code includes +I_DD
    if let Some(ii) = state.get_iono_idx(&m.key) {
        pr_h[ii] = 1.0;
    }
    h_rows.push(pr_h);
    let grad_pr = m.dgrad_n_rov * grad_n_val + m.dgrad_e_rov * grad_e_val;
    let pr_y = m.dd_pr_m - geom_dd - m.dm_wet_rov * zwd_val - grad_pr - iono_val;
    let pr_r = m.pr_var_m2.max(0.01);
    y_vals.push(pr_y);
    r_diag.push(robust_inflate(pr_y, pr_r));

    if let (Some(cp_obs), Some(amb_idx)) = (pcv_corrected_cp(m), state.get_amb_idx(&m.key)) {
        let amb_val = x_current[amb_idx];
        let pred_cp = geom_dd / m.lambda + amb_val
            + (m.dm_wet_rov * zwd_val + grad_pr) / m.lambda
            - iono_val / m.lambda;
        let mut cp_h = DVector::zeros(state_dim);
        cp_h[0] = d_geom_dpos.x / m.lambda;
        cp_h[1] = d_geom_dpos.y / m.lambda;
        cp_h[2] = d_geom_dpos.z / m.lambda;
        cp_h[amb_idx] = 1.0;
        // Iono state: DD phase includes -I_DD/λ
        if let Some(ii) = state.get_iono_idx(&m.key) {
            cp_h[ii] = -1.0 / m.lambda;
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
        r_diag.push(robust_inflate(cp_y, cp_r));
    }
}

/// Huber-style variance inflation: keep nominal R below the normalized
/// innovation threshold, decay weight smoothly above it.
fn robust_inflate(innovation: f64, variance: f64) -> f64 {
    let nis = innovation * innovation / variance;
    if nis > ROBUST_INNOVATION_THRESHOLD {
        variance * (nis / ROBUST_INNOVATION_THRESHOLD)
    } else {
        variance
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
        let Some(cp_obs) = pcv_corrected_cp(m) else { continue };
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
            dgrad_n_rov: 0.0,
            dgrad_e_rov: 0.0,
                tide_dd_m: 0.0,
            dm_wet_rov: 0.0,
            dd_pcv_m: 0.0,
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
                dgrad_n_rov: 0.0,
                dgrad_e_rov: 0.0,
                tide_dd_m: 0.0,
                dd_pcv_m: 0.0,
            });
        }

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

    // --- IF residual outlier screen (TDD) --------------------------------
    #[allow(clippy::type_complexity)]
    fn if_meas_fixture(
        slip_cycles: Option<f64>,
    ) -> (Vec<DoubleDiffMeasurement>, Vector3<f64>, HashMap<DoubleDiffKey, f64>, HashMap<DoubleDiffKey, f64>) {
        use gneiss_core::time::GpsTime;
        use super::super::state::DoubleDiffKey;
        let true_pos = Vector3::new(0.0, 0.0, 0.0);
        let base_pos = Vector3::new(4_000_000.0, 500_000.0, 1_000_000.0);
        let dirs = [
            Vector3::new(-20_000.0, 15_000.0, 18_000.0),
            Vector3::new(5_000.0, 24_000.0, -14_000.0),
            Vector3::new(16_000.0, -18_000.0, 22_000.0),
            Vector3::new(-9_000.0, -6_000.0, -26_000.0),
            Vector3::new(22_000.0, 9_000.0, -19_000.0),
            Vector3::new(-15_000.0, 21_000.0, 8_000.0),
        ];
        let ref_dir = dirs[0];
        let f1 = 1575.42e6_f64;
        let f2 = 1227.60e6_f64;
        let l1 = 299792458.0 / f1;
        let l2 = 299792458.0 / f2;
        let mut meas = Vec::new();
        for (i, dir) in dirs.iter().enumerate().skip(1) {
            let key = DoubleDiffKey { constellation_id: 0, sat: 1 + i as u16, ref_sat: 1, freq_band: 1 };
            let key2 = DoubleDiffKey { freq_band: 2, ..key };
            let sat_p = *dir * 20_000_000.0;
            let ref_p = ref_dir * 20_000_000.0;
            let rho = |p: Vector3<f64>| math_dist(true_pos, p);
            let geom_s = rho(sat_p); let geom_r = rho(ref_p);
            let base_dd = math_dist(base_pos, sat_p) - math_dist(base_pos, ref_p);
            for (k, lam, slip_here) in [(key, l1, false), (key2, l2, true)] {
                let n = 10.0 + i as f64;
                let mut cp = (geom_s - geom_r) / lam - base_dd / lam + n;
                // Inject on the FIRST pair only: equal slips on every pair
                // become common-mode and the median cancels them.
                if slip_here && i == 1 {
                    if let Some(sl) = slip_cycles { cp += sl; }
                }
                meas.push(DoubleDiffMeasurement {
                    key: k,
                    dd_pr_m: (geom_s - geom_r) - base_dd,
                    dd_cp_cycles: Some(cp),
                    sat_pos: sat_p,
                    ref_pos: ref_p,
                    base_pos,
                    lambda: lam,
                    pr_var_m2: 0.04,
                    cp_var_cycles2: 1e-4,
                    dm_wet_rov: 0.0,
                    dgrad_n_rov: 0.0,
                    dgrad_e_rov: 0.0,
                tide_dd_m: 0.0,
                dd_pcv_m: 0.0,
                });
            }
        }
        // Truth integers for every pair on both bands.
        let mut n1 = HashMap::new();
        let mut n2 = HashMap::new();
        for (i, dir) in dirs.iter().enumerate().skip(1) {
            let key = DoubleDiffKey { constellation_id: 0, sat: 1 + i as u16, ref_sat: 1, freq_band: 1 };
            let key2 = DoubleDiffKey { freq_band: 2, ..key };
            let n = 10.0 + i as f64;
            n1.insert(key, n);
            n2.insert(key2, n);
        }
        (meas, true_pos, n1, n2)
    }
    // local dist helper to avoid name clash
    #[allow(non_snake_case)]
    fn math_dist(a: Vector3<f64>, b: Vector3<f64>) -> f64 {
        (a - b).norm()
    }

    #[test]
    fn test_if_screen_clean_data_has_no_outliers() {
        let (meas, pos, n1, n2) = if_meas_fixture(None);
        let out = if_residual_outliers(pos, &meas, &n1, &n2);
        assert!(out.is_empty(), "clean data flagged: {:?}", out);
    }

    #[test]
    fn test_if_screen_flags_single_slipped_pair() {
        // One pair slipped +2 cycles on BOTH bands: geometry-free blind,
        // but its post-fix IF residual jumps ~2 * lambda_IF.
        let (meas, pos, n1, n2) = if_meas_fixture(Some(2.0));
        let out = if_residual_outliers(pos, &meas, &n1, &n2);
        assert_eq!(out.len(), 1, "expected exactly the slipped pair: {:?}", out);
        assert_eq!(out[0].sat, 2, "slipped pair is sat=2 vs ref=1");
    }

    #[test]
    fn test_if_screen_common_mode_position_error_not_flagged() {
        // A biased position shifts every pair coherently; the median
        // cancels it and no pair may be flagged.
        let (meas, _, n1, n2) = if_meas_fixture(None);
        let biased = Vector3::new(0.02, -0.012, 0.008);
        let out = if_residual_outliers(biased, &meas, &n1, &n2);
        assert!(out.is_empty(), "common-mode error flagged pairs: {:?}", out);
    }

    // --- Receiver-antenna PCV correction (TDD) ---------------------------

    /// Minimal single-pair fixture: state parked on the truth position,
    /// one ambiguity initialised, PCV signature already baked into the
    /// observed phase.
    fn pcv_fixture(dd_pcv_m: f64) -> (RtkState, DoubleDiffMeasurement) {
        let truth = Vector3::new(100.0, 200.0, 300.0);
        let mut state = RtkState::new(truth, GpsTime::new(2000, 100.0));
        let key = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        state.ensure_ambiguity(key, 7.25, 100.0);
        let m = DoubleDiffMeasurement {
            key,
            dd_pr_m: 10.0,
            dd_cp_cycles: Some(55.0),
            sat_pos: truth + Vector3::new(1.2e7, 0.4e7, 1.8e7),
            ref_pos: truth + Vector3::new(-0.6e7, 1.9e7, 1.1e7),
            base_pos: truth,
            lambda: 0.190,
            pr_var_m2: 0.04,
            cp_var_cycles2: 1e-4,
            dm_wet_rov: 0.0,
            dgrad_n_rov: 0.0,
            dgrad_e_rov: 0.0,
            tide_dd_m: 0.0,
            dd_pcv_m,
        };
        (state, m)
    }

    #[test]
    fn pcv_corrected_cp_subtracts_pcv_over_lambda() {
        let (_, m) = pcv_fixture(0.008);
        let corrected = pcv_corrected_cp(&m).expect("phase present");
        assert!(
            (corrected - (55.0 - 0.008 / 0.190)).abs() < 1e-12,
            "corrected={corrected}"
        );
    }

    #[test]
    fn pcv_corrected_cp_zero_pcv_is_identity() {
        let (_, m) = pcv_fixture(0.0);
        assert_eq!(pcv_corrected_cp(&m), Some(55.0));
    }

    #[test]
    fn pcv_corrected_cp_passes_none_through() {
        let (_, mut m) = pcv_fixture(0.008);
        m.dd_cp_cycles = None;
        assert_eq!(pcv_corrected_cp(&m), None);
    }

    /// The phase innovation must decrease by exactly `dd_pcv_m / lambda`
    /// when the correction is enabled: the embedded antenna signature is
    /// removed before model comparison.
    #[test]
    fn phase_innovation_shifts_exactly_by_dd_pcv_over_lambda() {
        let pcv_m = 0.008_f64;
        let (state_off, m_off) = pcv_fixture(0.0);
        let (state_on, m_on) = pcv_fixture(pcv_m);
        let x = |s: &RtkState| s.to_dvector();
        let (_, y_off, _) = build_measurement_system(&state_off, &x(&state_off), &[m_off]);
        let (_, y_on, _) = build_measurement_system(&state_on, &x(&state_on), &[m_on]);
        assert_eq!(y_off.len(), 2, "code + phase rows expected");
        assert_eq!(y_on.len(), 2);
        let shift = y_off[1] - y_on[1];
        let expected = pcv_m / 0.190;
        assert!(
            (shift - expected).abs() < 1e-12,
            "shift={shift} expected={expected}"
        );
    }

}
/// Per-pair iono-free residual outlier screen over a fixed ambiguity set.
/// For each band-1 pair with band-2 phases and fixed integers on both
/// bands, computes the DD iono-free float ambiguity against geometry at
/// `pos`. A same-cycle dual-frequency slip leaves the geometry-free
/// combination unchanged but shifts that pair's IF ambiguity by exactly
/// the slip count in cycles of lambda_IF (~10.7 cm GPS), so deviating
/// pairs stand out once the cross-pair MEDIAN cancels common-mode
/// position/model error. Returns keys deviating more than MAX_DEV metres.
pub fn if_residual_outliers(
    pos: Vector3<f64>,
    measurements: &[DoubleDiffMeasurement],
    fixed_n1: &HashMap<DoubleDiffKey, f64>,
    fixed_n2: &HashMap<DoubleDiffKey, f64>,
) -> Vec<DoubleDiffKey> {
    /// ~4x the ~1 cm iono-free observation noise, safely below half of
    /// lambda_IF (~10.7 cm GPS): a single whole-cycle slip cannot hide.
    const MAX_DEV_M: f64 = 0.05;

    let mut b2: HashMap<DoubleDiffKey, &DoubleDiffMeasurement> = HashMap::new();
    for m in measurements {
        if m.key.freq_band == 2 && m.dd_cp_cycles.is_some() {
            b2.insert(DoubleDiffKey { freq_band: 1, ..m.key }, m);
        }
    }

    let mut items: Vec<(DoubleDiffKey, f64)> = Vec::new();
    for m in measurements {
        if m.key.freq_band != 1 {
            continue;
        }
        let Some(cp1) = pcv_corrected_cp(m) else { continue };
        let Some(m2) = b2.get(&m.key) else { continue };
        let Some(cp2) = pcv_corrected_cp(m2) else { continue };
        let k2 = DoubleDiffKey { freq_band: 2, ..m.key };
        let (Some(n1), Some(n2)) = (fixed_n1.get(&m.key), fixed_n2.get(&k2)) else {
            continue;
        };

        let r_sat = (m.sat_pos - pos).norm();
        let r_ref = (m.ref_pos - pos).norm();
        let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
        let geom_m = (r_sat - r_ref) - base_dd;
        // Per-band post-fix range residuals (metres): observation minus
        // geometry minus the committed integer's range contribution.
        // Zero when integers are correct; a same-cycle dual-frequency slip
        // delta leaves d1 = delta*lambda1, d2 = delta*lambda2.
        let d1 = cp1 * m.lambda - geom_m - n1 * m.lambda;
        let d2 = cp2 * m2.lambda - geom_m - n2 * m2.lambda;
        // Ionosphere-free combination of the two range residuals:
        // (f1^2 d1 - f2^2 d2)/(f1^2 - f2^2), expressed via lambda ratios.
        // A same-cycle slip maps to delta * lambda_IF (~10.7 cm GPS).
        let res = (d1 / (m.lambda * m.lambda) - d2 / (m2.lambda * m2.lambda))
            / (1.0 / (m.lambda * m.lambda) - 1.0 / (m2.lambda * m2.lambda));
        items.push((m.key, res));
    }
    if items.len() < 3 {
        return Vec::new();
    }
    let mut rs: Vec<f64> = items.iter().map(|(_, v)| *v).collect();
    rs.sort_by(|a, b| a.total_cmp(b));
    let med = rs[rs.len() / 2];
    items.into_iter()
        .filter(|(_, v)| (v - med).abs() > MAX_DEV_M)
        .map(|(k, _)| k)
        .collect()
}
