//! Ionosphere-free (L1/L2) fixed-position update.
//!
//! The per-band fixed solution carries the double-difference ionospheric
//! residual, which grows with baseline (~20 cm at 50 km). Combining the L1/L2
//! carrier phases iono-free cancels first-order ionosphere exactly, so once
//! the per-band integers N1, N2 are fixed, the combined ambiguity
//! N_IF = (f1*N1 - f2*N2)/(f1 - f2) is known and the position can be
//! re-estimated from the iono-free phase without the bias.

use std::collections::HashMap;

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};

use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::obs::SatObs;
use gneiss_core::sat::SatelliteId;

use super::ar::ArResult;
use super::state::{DoubleDiffKey, RtkState};
use super::update::compute_tropo_dd;

/// Iono-free phase combination: phi_IF = (f1*phi1 - f2*phi2)/(f1 - f2).
pub fn combine_iono_free(f1_hz: f64, f2_hz: f64, phi1_cycles: f64, phi2_cycles: f64) -> f64 {
    (f1_hz * phi1_cycles - f2_hz * phi2_cycles) / (f1_hz - f2_hz)
}

/// Iono-free wavelength: lambda_IF = c/(f1 + f2).
pub fn lambda_iono_free(f1_hz: f64, f2_hz: f64) -> f64 {
    SPEED_OF_LIGHT_M_S / (f1_hz + f2_hz)
}

/// Iono-free double-difference phase for a satellite pair.
#[derive(Debug, Clone)]
pub struct IonoFreeMeasurement {
    /// Band-1 key identifying the pair.
    pub key: DoubleDiffKey,
    pub sat_pos: Vector3<f64>,
    pub ref_pos: Vector3<f64>,
    pub base_pos: Vector3<f64>,
    pub lambda_if: f64,
    pub f1_hz: f64,
    pub f2_hz: f64,
    pub dd_phase_if_cycles: f64,
    pub variance_cycles2: f64,
    /// Approximate rover-base distance (m): scales the atmospheric
    /// residual variance of this observation.
    pub baseline_m: f64,
    /// Gradient mapping differences (satellite − reference) at the rover:
    /// [north, east], cot(el)·{cos az, sin az} form. Consumed as a known
    /// correction when the fixed position is re-estimated.
    pub dgrad_n_rov: f64,
    pub dgrad_e_rov: f64,
}

/// Form the iono-free DD phase for a pair when both L1 and L2 are observed.
#[allow(clippy::too_many_arguments)] // same obs bundle as build_single_dd_pair
pub fn form_iono_free_dd(
    sat_id: SatelliteId,
    rov_s: &SatObs,
    bas_s: &SatObs,
    rov_ref: &SatObs,
    bas_ref: &SatObs,
    sat_pos: Vector3<f64>,
    ref_pos: Vector3<f64>,
    base_pos: Vector3<f64>,
    rover_pos: Vector3<f64>,
    key: DoubleDiffKey,
    glo_k: i8,
) -> Option<IonoFreeMeasurement> {
    let f1 = gneiss_core::frequencies::track_c_frequency(sat_id.constellation, 1, glo_k);
    // Secondary band: L2 for GPS/GLONASS; E5a (band 5) for Galileo
    // exports that carry no L2 slot. All four stations must have it.
    let b2 = if rov_s.get_observable_phase(2).is_some()
        && bas_s.get_observable_phase(2).is_some()
        && rov_ref.get_observable_phase(2).is_some()
        && bas_ref.get_observable_phase(2).is_some()
    {
        2
    } else if rov_s.get_observable_phase(7).is_some()
        && bas_s.get_observable_phase(7).is_some()
        && rov_ref.get_observable_phase(7).is_some()
        && bas_ref.get_observable_phase(7).is_some()
    {
        7
    } else {
        5
    };
    let f2 = gneiss_core::frequencies::track_c_frequency(sat_id.constellation, b2, glo_k);

    if f1 <= 0.0 || f2 <= 0.0 || (f1 - f2).abs() < 1e6 {
        return None;
    }
    let (p1_rs, p2_rs) = (rov_s.get_observable_phase(1)?, rov_s.get_observable_phase(b2)?);
    let (p1_rr, p2_rr) = (rov_ref.get_observable_phase(1)?, rov_ref.get_observable_phase(b2)?);
    let (p1_bs, p2_bs) = (bas_s.get_observable_phase(1)?, bas_s.get_observable_phase(b2)?);
    let (p1_br, p2_br) = (bas_ref.get_observable_phase(1)?, bas_ref.get_observable_phase(b2)?);

    let dd_if = (combine_iono_free(f1, f2, p1_rs, p2_rs) - combine_iono_free(f1, f2, p1_rr, p2_rr))
        - (combine_iono_free(f1, f2, p1_bs, p2_bs) - combine_iono_free(f1, f2, p1_br, p2_br));

    let lambda_if = lambda_iono_free(f1, f2);
    // Gradient mapping differences at the rover for this pair.
    let rx_llh = gneiss_core::coords::ecef_to_llh(rover_pos);
    let (az_s, el_s) = gneiss_core::coords::az_el(rx_llh, rover_pos, sat_pos);
    let (az_r, el_r) = gneiss_core::coords::az_el(rx_llh, rover_pos, ref_pos);
    let gterm = |az: f64, el: f64| -> (f64, f64) {
        let se = el.sin().max(0.17);
        let m = el.cos().max(0.0) / se;
        (m * az.cos(), m * az.sin())
    };
    let (gn_s, ge_s) = gterm(az_s, el_s);
    let (gn_r, ge_r) = gterm(az_r, el_r);

    Some(IonoFreeMeasurement {
        key,
        sat_pos,
        ref_pos,
        base_pos,
        lambda_if,
        f1_hz: f1,
        f2_hz: f2,
        dd_phase_if_cycles: dd_if,
        // IF noise: a 3 mm floor plus a wet-delay residual term growing
        // with baseline (0.15 mm/km class). Without the length term the
        // acceptance gate assumes short-baseline noise everywhere and
        // rejects nearly every long-baseline solution.
        variance_cycles2: {
            let baseline_km = (base_pos - rover_pos).norm() / 1000.0;
            let sigma_m = 0.003 + 0.00015 * baseline_km;
            2.0 * (sigma_m / lambda_if).powi(2)
        },
        baseline_m: (base_pos - rover_pos).norm(),
        dgrad_n_rov: gn_s - gn_r,
        dgrad_e_rov: ge_s - ge_r,
    })
}

/// Outcome of the fixed iono-free re-estimation.
#[derive(Debug, Clone)]
pub enum IonoFreeOutcome {
    /// Fewer than six both-band fixed pairs: stage not applicable.
    NotEngaged,
    /// Engaged but the post-fit residual / displacement gates rejected the
    /// integers — strong evidence the fixed set is wrong (e.g. a common-mode
    /// dual-band slip, which widelane consistency cannot see).
    Rejected,
    /// Iono-free position solution.
    Solution(Vector3<f64>, Matrix3<f64>),
}

/// Re-estimate the rover position from iono-free phase with known integers.
///
/// Returns the position and covariance only if a residual gate and a loose
/// physical bound pass (a wrong integer fix produces large post-fit
/// residuals and is rejected).
pub fn apply_fixed_iono_free(
    state: &RtkState,
    meas: &[IonoFreeMeasurement],
    ar: &ArResult,
) -> IonoFreeOutcome {
    let n1: HashMap<DoubleDiffKey, f64> = ar.fixed_ambiguities.iter()
        .filter(|(k, _)| k.freq_band == 1)
        .map(|(k, v)| (*k, *v))
        .collect();
    let n2: HashMap<DoubleDiffKey, f64> = ar.fixed_ambiguities.iter()
        .filter(|(k, _)| k.freq_band == 2)
        .map(|(k, v)| (*k, *v))
        .collect();
    // Current gradient estimates act as a known correction here (like the
    // troposphere model): the LSQ re-estimates position only.
    let (grad_n, grad_e) = if state.grad_enabled {
        (state.grad_n_m, state.grad_e_m)
    } else {
        (0.0, 0.0)
    };

    let mut h_rows = Vec::new();
    let mut y_vals = Vec::new();
    let mut r_diag = Vec::new();
    let cur_pos = state.pos_ecef;

    for m in meas {
        let key2 = DoubleDiffKey { freq_band: 2, ..m.key };
        let (Some(n1v), Some(n2v)) = (n1.get(&m.key), n2.get(&key2)) else { continue };
        let n_if = (m.f1_hz * n1v - m.f2_hz * n2v) / (m.f1_hz - m.f2_hz);

        let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
        let r_sat = (m.sat_pos - cur_pos).norm();
        let r_ref = (m.ref_pos - cur_pos).norm();
        let geom_dd = (r_sat - r_ref) - base_dd + compute_tropo_dd(m.sat_pos, m.ref_pos, m.base_pos, cur_pos);
        let grad_slant_m = m.dgrad_n_rov * grad_n + m.dgrad_e_rov * grad_e;

        let los_sat = (m.sat_pos - cur_pos) / r_sat.max(1e-3);
        let los_ref = (m.ref_pos - cur_pos) / r_ref.max(1e-3);
        let d = (los_ref - los_sat) / m.lambda_if;
        h_rows.push(d);
        y_vals.push(m.dd_phase_if_cycles - geom_dd / m.lambda_if - n_if
                    - grad_slant_m / m.lambda_if);
        r_diag.push(m.variance_cycles2.max(1e-4));
    }
    // Require a well-conditioned geometry: with fewer pairs the prior pulls
    // the estimate toward the float position and can regress the already
    // strong conditional per-band fix. Activates once the AR fixes both
    // bands on most pairs (better float ambiguity quality).
    if h_rows.len() < 6 {
        tracing::debug!("if-outcome: pairs={} NOT_ENGAGED", h_rows.len());
        return IonoFreeOutcome::NotEngaged;
    }
    // Expose the raw residual scale so the acceptance gate can be
    // calibrated against measured noise rather than assumptions.
    let mut sq = 0.0_f64;
    for &y in &y_vals {
        sq += y * y;
    }
    let prefit_rms = (sq / y_vals.len().max(1) as f64).sqrt();
    match solve_position_lsq(cur_pos, &h_rows, &y_vals, &r_diag, &state.extract_pos_cov()) {
        Some((pos, cov)) => {
            tracing::debug!("if-outcome: SOLUTION pairs={} prefit_rms={:.3}", h_rows.len(), prefit_rms);
            IonoFreeOutcome::Solution(pos, cov)
        }
        None => {
            tracing::debug!("if-outcome: REJECTED pairs={} prefit_rms={:.3}", h_rows.len(), prefit_rms);
            IonoFreeOutcome::Rejected
        }
    }
}

/// Solve the fixed iono-free position as a MAP estimate: the float position
/// covariance regularizes weak geometry (few both-bands-fixed pairs), the
/// same principle as the conditional per-band correction.
fn solve_position_lsq(
    cur_pos: Vector3<f64>,
    h_rows: &[Vector3<f64>],
    y_vals: &[f64],
    r_diag: &[f64],
    prior_cov: &Matrix3<f64>,
) -> Option<(Vector3<f64>, Matrix3<f64>)> {
    let n = h_rows.len();
    let mut h = DMatrix::zeros(n, 3);
    let mut r = DMatrix::zeros(n, n);
    let mut z = DVector::zeros(n);
    for i in 0..n {
        h[(i, 0)] = h_rows[i].x;
        h[(i, 1)] = h_rows[i].y;
        h[(i, 2)] = h_rows[i].z;
        z[i] = y_vals[i];
        r[(i, i)] = r_diag[i];
    }
    let r_inv = r.try_inverse()?;
    let mut prior_inv = Matrix3::zeros();
    let prior = prior_cov.map(|v| v.max(1e-4));
    for row in 0..3 {
        for col in 0..3 {
            prior_inv[(row, col)] = prior[(row, col)];
        }
    }
    let prior_inv = prior_inv.try_inverse()?;
    let n_plus = h.transpose() * &r_inv * &h + prior_inv;
    let n_inv = n_plus.try_inverse()?;
    let dx = n_inv * (h.transpose() * &r_inv * &z);

    // Reject wrong fixes: post-fit residual RMS must stay near measurement
    // noise. Log both sides so the gate stays calibrated against reality.
    let v = &z - &h * dx;
    let rms = (v.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();
    let exp_rms = (r_diag.iter().sum::<f64>() / n as f64).sqrt();
    tracing::debug!(
        "if-gate: postfit_rms={rms:.3} exp={exp_rms:.3} dx_norm={:.3}",
        dx.norm()
    );
    if rms > 4.0 * exp_rms || dx.norm() > 1.5 {
        return None;
    }
    let mut cov = Matrix3::zeros();
    for row in 0..3 {
        for col in 0..3 {
            cov[(row, col)] = n_inv[(row, col)];
        }
    }
    Some((cur_pos + Vector3::new(dx[0], dx[1], dx[2]), cov))
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::time::GpsTime;

    const F1: f64 = 1575.42e6;
    const F2: f64 = 1227.60e6;

    fn make_key(sat: u16) -> DoubleDiffKey {
        DoubleDiffKey { constellation_id: 0, sat, ref_sat: 1, freq_band: 1 }
    }

    #[test]
    fn test_iono_free_cancels_ionosphere_exactly() {
        let range = 23_000_000.0;
        let n1 = 5.0;
        let n2 = -3.0;
        let iono_m = 0.25; // L1 path iono delay (m)
        // Phase iono is an advance: -iono/λ on each band, L2 scaled by (f1/f2)^2.
        let phi1 = range * F1 / SPEED_OF_LIGHT_M_S + n1 - iono_m * F1 / SPEED_OF_LIGHT_M_S;
        let phi2 = range * F2 / SPEED_OF_LIGHT_M_S + n2 - iono_m * (F1 / F2).powi(2) * F2 / SPEED_OF_LIGHT_M_S;

        let phi_if = combine_iono_free(F1, F2, phi1, phi2);
        let n_if = (F1 * n1 - F2 * n2) / (F1 - F2);
        let expect = range / lambda_iono_free(F1, F2) + n_if;
        assert!((phi_if - expect).abs() < 1e-9, "IF combination must cancel iono, got {:.12e}", phi_if - expect);
    }

    #[test]
    fn test_fixed_iono_free_position_recovers_truth_with_iono() {
        let true_pos = Vector3::new(100.0, 200.0, 300.0);
        let base_pos = Vector3::new(0.0, 0.0, 0.0);
        let state = RtkState::new(true_pos + Vector3::new(0.8, -0.4, 0.2), GpsTime::new(2000, 100.0));

        let sats: [(Vector3<f64>, f64, f64, f64); 6] = [
            (Vector3::new(10_000.0, 20_000.0, 20_000.0), 1.0, 2.0, 0.30),
            (Vector3::new(25_000.0, 5_000.0, 18_000.0), 3.0, -1.0, 0.22),
            (Vector3::new(8_000.0, 30_000.0, 12_000.0), -2.0, 4.0, 0.18),
            (Vector3::new(20_000.0, 12_000.0, 25_000.0), 2.0, 0.0, 0.26),
            (Vector3::new(15_000.0, 22_000.0, 15_000.0), 0.0, -2.0, 0.20),
            (Vector3::new(22_000.0, 18_000.0, 22_000.0), 4.0, 1.0, 0.24),
        ];
        let ref_pos = Vector3::new(30_000.0, 5_000.0, 10_000.0);

        let mut meas = Vec::new();
        for (i, (sat_pos, n1, n2, iono_m)) in sats.iter().enumerate() {
            let base_dd = (sat_pos - base_pos).norm() - (ref_pos - base_pos).norm();
            let geom = (sat_pos - true_pos).norm() - (ref_pos - true_pos).norm() - base_dd;
            let lambda1 = SPEED_OF_LIGHT_M_S / F1;
            let lambda2 = SPEED_OF_LIGHT_M_S / F2;
            let iono_l2 = iono_m * (F1 / F2).powi(2);
            let dd1 = geom / lambda1 + n1 - iono_m / lambda1;
            let dd2 = geom / lambda2 + n2 - iono_l2 / lambda2;
            meas.push(IonoFreeMeasurement {
                key: make_key(2 + i as u16),
                sat_pos: *sat_pos,
                ref_pos,
                base_pos,
                lambda_if: lambda_iono_free(F1, F2),
                f1_hz: F1,
                f2_hz: F2,
                dd_phase_if_cycles: combine_iono_free(F1, F2, dd1, dd2),
                variance_cycles2: 1e-4,
                baseline_m: (base_pos - true_pos).norm() + 1000.0,
                dgrad_n_rov: 0.0,
                dgrad_e_rov: 0.0,
            });
        }
        let fixed: Vec<(DoubleDiffKey, f64)> = sats.iter().enumerate()
            .flat_map(|(i, (_, n1, n2, _))| vec![(make_key(2 + i as u16), *n1), (DoubleDiffKey { freq_band: 2, ..make_key(2 + i as u16) }, *n2)])
            .collect();
        let ar = ArResult {
            position_ecef: true_pos,
            cov_position: Matrix3::identity(),
            ratio: 3.0,
            is_fixed: true,
            num_ambiguities: 10,
            fixed_ambiguities: fixed,
        };

        let (pos, _) = match apply_fixed_iono_free(&state, &meas, &ar) {
            IonoFreeOutcome::Solution(p, c) => (p, c),
            _ => panic!("IF stage should produce a solution"),
        };
        let err = (pos - true_pos).norm();
        assert!(err < 0.01, "IF fixed position should recover truth under iono, got {:.4}m", err);
    }

    #[test]
    fn test_rejected_when_integers_shifted_by_common_mode_slip() {
        // Same fixture as above but every N1 and N2 shifted +1: the common
        // mode cancels in N1-N2 (widelane-blind) yet biases the iono-free
        // combination by a full cycle, which the residual gate must reject.
        let true_pos = Vector3::new(100.0, 200.0, 300.0);
        let base_pos = Vector3::new(0.0, 0.0, 0.0);
        let state = RtkState::new(true_pos, GpsTime::new(2000, 100.0));

        let sats: [(Vector3<f64>, f64, f64); 6] = [
            (Vector3::new(10_000.0, 20_000.0, 20_000.0), 1.0, 2.0),
            (Vector3::new(25_000.0, 5_000.0, 18_000.0), 3.0, -1.0),
            (Vector3::new(8_000.0, 30_000.0, 12_000.0), -2.0, 4.0),
            (Vector3::new(20_000.0, 12_000.0, 25_000.0), 2.0, 0.0),
            (Vector3::new(15_000.0, 22_000.0, 15_000.0), 0.0, -2.0),
            (Vector3::new(22_000.0, 18_000.0, 22_000.0), 4.0, 1.0),
        ];
        let ref_pos = Vector3::new(30_000.0, 5_000.0, 10_000.0);
        let lambda1 = SPEED_OF_LIGHT_M_S / F1;
        let lambda2 = SPEED_OF_LIGHT_M_S / F2;

        let mut meas = Vec::new();
        for (i, (sat_pos, n1, n2)) in sats.iter().enumerate() {
            let geom = (*sat_pos - true_pos).norm() - (ref_pos - true_pos).norm();
            let dd1 = geom / lambda1; // integer-free phases: truth-consistent
            let dd2 = geom / lambda2;
            meas.push((make_key(2 + i as u16), *sat_pos, combine_iono_free(F1, F2, dd1, dd2), *n1 + 1.0, *n2 + 1.0));
        }
        let if_meas: Vec<IonoFreeMeasurement> = meas.iter().map(|(k, sat_pos, phase_if, _, _)| {
            IonoFreeMeasurement {
                key: *k,
                sat_pos: *sat_pos,
                ref_pos,
                base_pos,
                lambda_if: lambda_iono_free(F1, F2),
                f1_hz: F1,
                f2_hz: F2,
                dd_phase_if_cycles: *phase_if,
                variance_cycles2: 1e-4,
                dgrad_n_rov: 0.0,
                dgrad_e_rov: 0.0,
                baseline_m: 1000.0,
            }
        }).collect();
        let fixed: Vec<(DoubleDiffKey, f64)> = meas.iter()
            .flat_map(|(k, _, _, n1, n2)| vec![(*k, *n1), (DoubleDiffKey { freq_band: 2, ..*k }, *n2)])
            .collect();
        let ar = ArResult {
            position_ecef: true_pos,
            cov_position: Matrix3::identity(),
            ratio: 3.0,
            is_fixed: true,
            num_ambiguities: 12,
            fixed_ambiguities: fixed,
        };
        assert!(matches!(
            apply_fixed_iono_free(&state, &if_meas, &ar),
            IonoFreeOutcome::Rejected
        ));
    }
}
