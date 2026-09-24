//! Tightly-coupled integer ambiguity conditioning inside the 15-state ESKF.
//!
//! Conditioning equations (Teunissen 1995, Groves 2013 Eq 14.65-14.66):
//! - State error injection: dx = -P_xa Q_aa^-1 (a_float - a_fixed)
//! - Covariance collapse: P_check = 0.5 * ((P_xx - P_xa Q_aa^-1 P_ax) + sym^T)
//! - Positive definiteness: diagonal variance floored at 10^-10

use nalgebra::{DMatrix, DVector, SymmetricEigen};

use super::types::{EngineError, EskfState, Matrix15, Vector15};
use super::update::apply_error_injection;

const DIAG_FLOOR: f64 = 1e-10;
const MAX_POS_JUMP_M: f64 = 2.0;
const MIN_POS_GATE_M: f64 = 0.50;

/// Result summary of integer ambiguity conditioning on ESKF.
#[derive(Debug, Clone, PartialEq)]
pub struct ConditionSummary {
    pub dx: Vector15<f64>,
    pub trace_reduction: f64,
    pub min_eigenvalue: f64,
    pub applied: bool,
    pub accepted: bool,
}

fn check_dimensions(
    p_xa: &DMatrix<f64>,
    q_aa: &DMatrix<f64>,
    a_float: &DVector<f64>,
    a_fixed: &DVector<f64>,
) -> Result<usize, EngineError> {
    let n = q_aa.nrows();
    if q_aa.ncols() != n || p_xa.nrows() != 15 || p_xa.ncols() != n {
        return Err(EngineError::InvalidMeasurement("Matrix dimension mismatch".into()));
    }
    if a_float.len() != n || a_fixed.len() != n {
        return Err(EngineError::InvalidMeasurement("Vector dimension mismatch".into()));
    }
    Ok(n)
}

fn compute_dx_and_cov_reduction(
    p_xa: &DMatrix<f64>,
    q_inv: &DMatrix<f64>,
    a_float: &DVector<f64>,
    a_fixed: &DVector<f64>,
) -> (Vector15<f64>, Matrix15<f64>) {
    let da = a_float - a_fixed;
    let dx_inj = -(p_xa * (q_inv * da));
    let cov_red = p_xa * q_inv * p_xa.transpose();
    let mut dx = Vector15::zeros();
    let mut delta_p = Matrix15::zeros();
    for i in 0..15 {
        dx[i] = dx_inj[i];
        for j in 0..15 {
            delta_p[(i, j)] = cov_red[(i, j)];
        }
    }
    (dx, delta_p)
}

fn check_jump_gate(dx: &Vector15<f64>, cov: &Matrix15<f64>) -> bool {
    let dx_pos = dx.fixed_rows::<3>(0);
    let dx_norm = dx_pos.norm();
    let sigma_3d = (cov[(0, 0)] + cov[(1, 1)] + cov[(2, 2)]).sqrt();
    let gate = (3.0 * sigma_3d).max(MIN_POS_GATE_M);
    dx_norm <= gate && dx_norm <= MAX_POS_JUMP_M
}

fn collapse_covariance(cov: &Matrix15<f64>, delta_p: &Matrix15<f64>) -> Matrix15<f64> {
    let p_raw = cov - delta_p;
    let mut p_check = 0.5 * (p_raw + p_raw.transpose());
    for i in 0..15 {
        if p_check[(i, i)] < DIAG_FLOOR {
            p_check[(i, i)] = DIAG_FLOOR;
        }
    }
    p_check
}

fn make_empty_summary(cov: &Matrix15<f64>) -> ConditionSummary {
    let min_eig = SymmetricEigen::new(*cov).eigenvalues.min();
    ConditionSummary {
        dx: Vector15::zeros(),
        trace_reduction: 0.0,
        min_eigenvalue: min_eig,
        applied: true,
        accepted: true,
    }
}

fn make_rejected_summary(dx: Vector15<f64>, cov: &Matrix15<f64>) -> ConditionSummary {
    let min_eig = SymmetricEigen::new(*cov).eigenvalues.min();
    ConditionSummary {
        dx,
        trace_reduction: 0.0,
        min_eigenvalue: min_eig,
        applied: false,
        accepted: false,
    }
}

fn apply_conditioned_state_and_cov(
    state: &mut EskfState,
    dx: &Vector15<f64>,
    delta_p: &Matrix15<f64>,
) -> ConditionSummary {
    let p_check = collapse_covariance(&state.cov, delta_p);
    let min_eig = SymmetricEigen::new(p_check).eigenvalues.min();
    if min_eig < DIAG_FLOOR {
        return make_rejected_summary(*dx, &state.cov);
    }
    let trace_reduction = state.cov.trace() - p_check.trace();
    apply_error_injection(state, dx);
    state.cov = p_check;
    ConditionSummary {
        dx: *dx,
        trace_reduction,
        min_eigenvalue: min_eig,
        applied: true,
        accepted: true,
    }
}

/// Applies integer ambiguity conditioning to the 15-state ESKF.
pub fn apply_integer_conditioning(
    state: &mut EskfState,
    p_xa: &DMatrix<f64>,
    q_aa: &DMatrix<f64>,
    a_float: &DVector<f64>,
    a_fixed: &DVector<f64>,
) -> Result<ConditionSummary, EngineError> {
    let n = check_dimensions(p_xa, q_aa, a_float, a_fixed)?;
    if n == 0 {
        return Ok(make_empty_summary(&state.cov));
    }
    let q_inv = q_aa.clone().try_inverse().ok_or(EngineError::InversionError)?;
    let (dx, delta_p) = compute_dx_and_cov_reduction(p_xa, &q_inv, a_float, a_fixed);
    if !check_jump_gate(&dx, &state.cov) {
        return Ok(make_rejected_summary(dx, &state.cov));
    }
    Ok(apply_conditioned_state_and_cov(state, &dx, &delta_p))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::estimators::eskf::dd_update::{compute_dd_jacobian_15, DdSatGeometry};
    use nalgebra::{Matrix3, UnitQuaternion, Vector3};

    fn make_test_state(pos: Vector3<f64>, att: UnitQuaternion<f64>, diag: &[f64; 15]) -> EskfState {
        let mut cov = Matrix15::zeros();
        for i in 0..15 {
            cov[(i, i)] = diag[i];
        }
        EskfState::with_cov(pos, Vector3::zeros(), att, Vector3::zeros(), Vector3::zeros(), cov)
    }

    const BASE_DIAG: [f64; 15] = [
        1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01, 0.04, 0.04, 0.04, 1e-4, 1e-4, 1e-4,
    ];

    #[test]
    fn test_tier1_golden_1d_integer_conditioning() {
        let mut st = make_test_state(Vector3::zeros(), UnitQuaternion::identity(), &BASE_DIAG);
        let mut p_xa = DMatrix::zeros(15, 1);
        p_xa[(0, 0)] = 0.02;
        p_xa[(1, 0)] = -0.01;
        p_xa[(2, 0)] = 0.03;
        p_xa[(6, 0)] = 0.001;
        let q_aa = DMatrix::from_element(1, 1, 0.04);
        let a_flt = DVector::from_element(1, 2.15);
        let a_fix = DVector::from_element(1, 2.00);

        let res = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix)
            .expect("conditioning should succeed");
        assert!(res.applied && res.accepted);
        assert!((res.dx[0] - (-0.075)).abs() < 1e-12);
        assert!((res.dx[1] - 0.0375).abs() < 1e-12);
        assert!((res.dx[2] - (-0.1125)).abs() < 1e-12);
        assert!((res.dx[6] - (-0.00375)).abs() < 1e-12);
        assert_eq!(res.dx.fixed_rows::<3>(3).norm(), 0.0);
        assert_eq!(res.dx.fixed_rows::<6>(9).norm(), 0.0);
        assert!((st.cov[(0, 0)] - 0.9900).abs() < 1e-12);
        assert!((st.cov[(1, 1)] - 0.9975).abs() < 1e-12);
        assert!((st.cov[(2, 2)] - 0.9775).abs() < 1e-12);
        assert!((st.cov[(6, 6)] - 0.009975).abs() < 1e-12);
        assert!((st.cov[(0, 1)] - 0.0050).abs() < 1e-12);
        assert!((st.cov[(0, 2)] - (-0.0150)).abs() < 1e-12);
        assert!((st.cov[(1, 2)] - 0.0075).abs() < 1e-12);
    }

    #[test]
    fn test_tier1_golden_2d_correlated_integer_conditioning() {
        let mut st = make_test_state(Vector3::zeros(), UnitQuaternion::identity(), &BASE_DIAG);
        let q_aa = DMatrix::from_row_slice(2, 2, &[0.08, 0.04, 0.04, 0.08]);
        let mut p_xa = DMatrix::zeros(15, 2);
        p_xa[(0, 0)] = 0.01;
        p_xa[(0, 1)] = 0.02;
        p_xa[(1, 0)] = -0.02;
        p_xa[(1, 1)] = 0.01;
        p_xa[(2, 0)] = 0.03;
        p_xa[(2, 1)] = -0.01;
        p_xa[(6, 0)] = 0.001;
        p_xa[(6, 1)] = 0.002;
        let a_flt = DVector::from_column_slice(&[1.12, -2.06]);
        let a_fix = DVector::from_column_slice(&[1.0, -2.0]);

        let res = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix)
            .expect("2D conditioning succeeds");
        assert!(res.applied);
        assert!((res.dx[0] - 0.0150).abs() < 1e-12);
        assert!((res.dx[1] - 0.0700).abs() < 1e-12);
        assert!((res.dx[2] - (-0.0950)).abs() < 1e-12);
        assert!((res.dx[6] - 0.0015).abs() < 1e-12);
        assert!((st.cov[(0, 0)] - 0.9950).abs() < 1e-12);
        assert!((st.cov[(0, 1)] - (-0.0025)).abs() < 1e-12);
        assert!((st.cov[(1, 1)] - (593.0 / 600.0)).abs() < 1e-12);
        assert!((st.cov[(2, 2)] - (587.0 / 600.0)).abs() < 1e-12);
    }

    #[test]
    fn test_tier1_groves_ch14_lever_arm_attitude_benchmark() {
        let (r, p, y) = (5.0_f64.to_radians(), (-3.0_f64).to_radians(), 45.0_f64.to_radians());
        let att = UnitQuaternion::from_euler_angles(r, p, y);
        let mut st = make_test_state(Vector3::new(10.0, 20.0, 30.0), att, &BASE_DIAG);
        let lever = Vector3::new(0.20, -0.10, 1.20);
        let (u_s, u_ref) = (Vector3::new(0.267261, 0.534522, 0.801784), Vector3::new(-0.408248, 0.408248, 0.816497));
        let geom = DdSatGeometry::new(u_s, u_ref, lever);
        let h = compute_dd_jacobian_15(&st, &geom);
        let p_xa = DMatrix::from_fn(15, 1, |row, _| 0.1 * st.cov[(row, row)] * h[(0, row)]);
        let q_aa = DMatrix::from_element(1, 1, 0.04);
        let (a_flt, a_fix) = (DVector::from_element(1, 2.1), DVector::from_element(1, 2.0));
        let res = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix).expect("ok");

        let delta_u = u_s - u_ref;
        let dx_pos = res.dx.fixed_rows::<3>(0);
        let cos_angle = dx_pos.dot(&delta_u) / (dx_pos.norm() * delta_u.norm());
        assert!((cos_angle - 1.0).abs() < 1e-10);
        let l_e = att.to_rotation_matrix() * lever;
        let dx_att = res.dx.fixed_rows::<3>(6);
        assert!(dx_att.dot(&l_e).abs() < 1e-10);
        assert!(dx_att.dot(&delta_u).abs() < 1e-10);
    }

    #[test]
    fn test_tier1_identity_fix_zero_injection_schur_contraction() {
        let mut st = make_test_state(Vector3::new(1.0, 2.0, 3.0), UnitQuaternion::identity(), &BASE_DIAG);
        let p_prior = st.cov;
        let p_xa = DMatrix::from_fn(15, 2, |r, c| if r < 3 { 0.02 * (c as f64 + 1.0) } else { 0.0 });
        let q_aa = DMatrix::from_diagonal(&DVector::from_column_slice(&[0.05, 0.05]));
        let a_same = DVector::from_column_slice(&[3.0, -1.0]);

        let res = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_same, &a_same)
            .expect("identity fix succeeds");
        assert_eq!(res.dx.norm(), 0.0);
        assert_eq!(st.pos_ecef, Vector3::new(1.0, 2.0, 3.0));
        assert!(st.cov.trace() < p_prior.trace());
        let diff = p_prior - st.cov;
        let eig = SymmetricEigen::new(diff);
        assert!(eig.eigenvalues.min() >= -1e-14);
    }

    #[test]
    fn test_tier2_conditioning_jacobian_finite_difference() {
        let st = make_test_state(Vector3::zeros(), UnitQuaternion::identity(), &BASE_DIAG);
        let p_xa = DMatrix::from_fn(15, 2, |r, c| 0.01 * (r as f64 + 1.0) * (c as f64 + 1.0));
        let q_aa = DMatrix::from_row_slice(2, 2, &[0.06, 0.02, 0.02, 0.05]);
        let a_fix = DVector::from_column_slice(&[2.0, 3.0]);
        let a_flt = DVector::from_column_slice(&[2.1, 3.05]);
        let q_inv = q_aa.clone().try_inverse().expect("invertible");
        let j_analytical = -(&p_xa * &q_inv);

        let eps = 1e-5;
        for col in 0..2 {
            let mut flt_plus = a_flt.clone();
            let mut flt_minus = a_flt.clone();
            flt_plus[col] += eps;
            flt_minus[col] -= eps;
            let mut st_p = st.clone();
            let mut st_m = st.clone();
            let res_p = apply_integer_conditioning(&mut st_p, &p_xa, &q_aa, &flt_plus, &a_fix).expect("ok");
            let res_m = apply_integer_conditioning(&mut st_m, &p_xa, &q_aa, &flt_minus, &a_fix).expect("ok");
            for r in 0..15 {
                let d_num = (res_p.dx[r] - res_m.dx[r]) / (2.0 * eps);
                assert!((d_num - j_analytical[(r, col)]).abs() < 1e-7);
            }
        }
    }

    #[test]
    fn test_tier2_dd_carrier_jacobian_15state_finite_difference() {
        let att = UnitQuaternion::from_euler_angles(0.08, -0.05, 0.2);
        let st = make_test_state(Vector3::new(100.0, -200.0, 300.0), att, &BASE_DIAG);
        let lever = Vector3::new(0.25, -0.1, 0.9);
        let r_s = Vector3::new(20_000_000.0, 5_000_000.0, 10_000_000.0);
        let r_ref = Vector3::new(-10_000_000.0, 15_000_000.0, 18_000_000.0);
        let p_ant = st.pos_ecef + att.to_rotation_matrix() * lever;
        let geom = DdSatGeometry::new((r_s - p_ant).normalize(), (r_ref - p_ant).normalize(), lever);
        let h_ana = compute_dd_jacobian_15(&st, &geom);

        let range_fn = |pos: &Vector3<f64>, q: &UnitQuaternion<f64>| {
            let ant = pos + q.to_rotation_matrix() * lever;
            (r_s - ant).norm() - (r_ref - ant).norm()
        };
        for i in 0..3 {
            let mut p_p = st.pos_ecef;
            let mut p_m = st.pos_ecef;
            p_p[i] += 0.05;
            p_m[i] -= 0.05;
            let d_num = (range_fn(&p_p, &att) - range_fn(&p_m, &att)) / 0.10;
            assert!((d_num - h_ana[(0, i)]).abs() < 1e-7);
        }
        for i in 0..3 {
            let mut v = Vector3::zeros();
            v[i] = 1e-3;
            let q_p = UnitQuaternion::from_scaled_axis(v) * att;
            let q_m = UnitQuaternion::from_scaled_axis(-v) * att;
            let d_num = (range_fn(&st.pos_ecef, &q_p) - range_fn(&st.pos_ecef, &q_m)) / 2e-3;
            assert!((d_num - h_ana[(0, 6 + i)]).abs() < 1e-5);
        }
    }

    #[test]
    fn test_tier2_taylor_quadratic_step_halving_convergence() {
        let att = UnitQuaternion::from_euler_angles(0.1, -0.1, 0.3);
        let pos = Vector3::new(50.0, -100.0, 150.0);
        let lever = Vector3::new(0.3, -0.2, 1.0);
        let r_s = Vector3::new(1000.0, 500.0, 800.0);
        let r_ref = Vector3::new(-800.0, 600.0, 700.0);
        let st = make_test_state(pos, att, &BASE_DIAG);
        let p_ant = pos + att.to_rotation_matrix() * lever;
        let geom = DdSatGeometry::new((r_s - p_ant).normalize(), (r_ref - p_ant).normalize(), lever);
        let h_ana = compute_dd_jacobian_15(&st, &geom);

        let eval_rem = |eps: f64| {
            let v = Vector3::new(eps, 0.0, 0.0);
            let q_p = UnitQuaternion::from_scaled_axis(v) * att;
            let ant_p = pos + q_p.to_rotation_matrix() * lever;
            let h_val = (r_s - ant_p).norm() - (r_ref - ant_p).norm();
            let h_0 = (r_s - p_ant).norm() - (r_ref - p_ant).norm();
            let lin_pred = h_0 + h_ana[(0, 6)] * eps;
            (h_val - lin_pred).abs()
        };
        let r1 = eval_rem(0.02);
        let r2 = eval_rem(0.01);
        let ratio = r1 / r2;
        assert!((ratio - 4.0).abs() < 0.5, "Ratio was {ratio}");
    }

    #[test]
    fn test_tier3_post_fix_covariance_strictly_positive_definite() {
        let mut st = make_test_state(Vector3::zeros(), UnitQuaternion::identity(), &BASE_DIAG);
        let q_aa = DMatrix::from_diagonal(&DVector::from_column_slice(&[0.05, 0.06, 0.04]));
        let p_xa = DMatrix::from_fn(15, 3, |r, c| {
            0.2 * (st.cov[(r, r)] * q_aa[(c, c)]).sqrt() / (c as f64 + 1.0)
        });
        let a_flt = DVector::from_column_slice(&[1.1, 2.05, -3.08]);
        let a_fix = DVector::from_column_slice(&[1.0, 2.0, -3.0]);

        let res = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix).expect("ok");
        assert!(res.min_eigenvalue >= 1e-10);
        let eig = SymmetricEigen::new(st.cov);
        assert!(eig.eigenvalues.min() >= 1e-10);
        assert!((st.cov - st.cov.transpose()).norm() < 1e-12);
    }

    #[test]
    fn test_tier3_monotonic_covariance_loewner_contraction() {
        let mut st = make_test_state(Vector3::zeros(), UnitQuaternion::identity(), &BASE_DIAG);
        let p_before = st.cov;
        let q_aa = DMatrix::from_diagonal(&DVector::from_column_slice(&[0.04, 0.04]));
        let p_xa = DMatrix::from_fn(15, 2, |r, c| {
            0.2 * (st.cov[(r, r)] * q_aa[(c, c)]).sqrt() / (c as f64 + 1.0)
        });
        let a_flt = DVector::from_column_slice(&[0.1, -0.05]);
        let a_fix = DVector::from_column_slice(&[0.0, 0.0]);

        let res = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix).expect("ok");
        assert!(res.trace_reduction > 0.0);
        assert!(st.cov.trace() < p_before.trace());
        let delta_p = p_before - st.cov;
        let eig = SymmetricEigen::new(delta_p);
        assert!(eig.eigenvalues.min() >= -1e-14);
        for i in 0..15 {
            assert!(st.cov[(i, i)] <= p_before[(i, i)] + 1e-14);
        }
    }

    #[test]
    fn test_tier3_lyapunov_closed_loop_error_convergence() {
        let true_pos = Vector3::new(10.0, 20.0, 30.0);
        let mut st = make_test_state(
            true_pos + Vector3::new(0.18, -0.22, 0.15),
            UnitQuaternion::identity(),
            &BASE_DIAG,
        );
        let prior_err = (st.pos_ecef - true_pos).norm();
        let h_pos = Matrix3::<f64>::identity();
        let mut p_xa = DMatrix::zeros(15, 3);
        for i in 0..3 {
            for j in 0..3 {
                p_xa[(i, j)] = 0.5 * h_pos[(i, j)];
            }
        }
        let q_aa = DMatrix::from_diagonal(&DVector::from_column_slice(&[0.5, 0.5, 0.5]));
        let a_flt = DVector::from_column_slice(&[0.18, -0.22, 0.15]);
        let a_fix = DVector::zeros(3);

        let res = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix).expect("ok");
        assert!(res.applied);
        let post_err = (st.pos_ecef - true_pos).norm();
        assert!(post_err < 0.03, "Post error was {post_err}");
        assert!(post_err < 0.10 * prior_err);
    }

    #[test]
    fn test_tier3_blunder_jump_gate_rejection_float_fallback() {
        let diag = [0.01; 15]; // sigma_3d = sqrt(0.03) = 0.173m -> gate = max(3*0.173, 0.50) = 0.52m
        let mut st = make_test_state(Vector3::new(5.0, 5.0, 5.0), UnitQuaternion::identity(), &diag);
        let prior_st = st.clone();
        let mut p_xa = DMatrix::zeros(15, 1);
        p_xa[(0, 0)] = 0.05;
        let q_aa = DMatrix::from_element(1, 1, 0.02);
        let a_flt = DVector::from_element(1, 1.35); // da = 1.35 -> dx_0 = -0.05/0.02 * 1.35 = -3.375m > 2.0m
        let a_fix = DVector::from_element(1, 0.0);

        let res = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix).expect("returns summary");
        assert!(!res.applied && !res.accepted);
        assert_eq!(st.pos_ecef, prior_st.pos_ecef);
        assert_eq!(st.cov, prior_st.cov);
    }

    #[test]
    fn test_tier3_rejects_negative_eigenvalue_indefinite_covariance() {
        let diag = [0.01; 15];
        let mut st = make_test_state(Vector3::zeros(), UnitQuaternion::identity(), &diag);
        let prior_st = st.clone();
        let mut p_xa = DMatrix::zeros(15, 1);
        p_xa[(0, 0)] = 0.20;
        p_xa[(1, 0)] = 0.20;
        let q_aa = DMatrix::from_element(1, 1, 0.01);
        let a_flt = DVector::from_element(1, 1.001);
        let a_fix = DVector::from_element(1, 1.0);

        let res = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix).expect("returns summary");
        assert!(!res.applied && !res.accepted);
        assert_eq!(st.pos_ecef, prior_st.pos_ecef);
        assert_eq!(st.cov, prior_st.cov);
    }
}
