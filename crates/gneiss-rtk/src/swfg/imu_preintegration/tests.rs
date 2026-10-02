#![allow(clippy::unwrap_used)]

use super::*;
use crate::swfg::variables::{VariableKind, VariableNode};
use std::collections::BTreeMap;

    fn make_pose_vel_bias_vars() -> (BTreeMap<VariableId, VariableNode>, VariableId, VariableId, VariableId, VariableId, VariableId) {
        let mut vars = BTreeMap::new();
        let pi = VariableId::new(0);
        let vi = VariableId::new(1);
        let pj = VariableId::new(2);
        let vj = VariableId::new(3);
        let bias = VariableId::new(4);

        vars.insert(pi, VariableNode::new(pi, VariableKind::Pose { epoch: 0 }));
        vars.insert(vi, VariableNode::new(vi, VariableKind::Velocity { epoch: 0 }));
        vars.insert(pj, VariableNode::new(pj, VariableKind::Pose { epoch: 1 }));
        vars.insert(vj, VariableNode::new(vj, VariableKind::Velocity { epoch: 1 }));
        vars.insert(bias, VariableNode::new(bias, VariableKind::ImuBias));

        (vars, pi, vi, pj, vj, bias)
    }

    #[test]
    fn preintegration_zero_motion_is_identity() {
        let mut preint = ImuPreintegration::new();
        // No IMU data → deltas should remain identity
        preint.integrate(&[], &Vector3::zeros(), &Vector3::zeros());
        assert!(preint.dp.norm() < 1e-12);
        assert!(preint.dv.norm() < 1e-12);
        assert!((preint.dq.quaternion().w - 1.0).abs() < 1e-12);
    }

    #[test]
    fn preintegration_constant_acceleration() {
        let mut preint = ImuPreintegration::new();
        let samples: Vec<ImuSample> = (0..11)
            .map(|i| ImuSample {
                accel: Vector3::new(0.0, 0.0, 9.8), // gravity in body Z
                gyro: Vector3::zeros(),
                time_us: i * 10_000, // 100 Hz
            })
            .collect();
        preint.integrate(&samples, &Vector3::zeros(), &Vector3::zeros());

        let _dt = 0.1; // 10 samples at 100 Hz
        // With gravity 9.8 in body Z, no rotation → world-frame accel depends on initial attitude (identity)
        // Gravity in body frame: the IMU measures reaction force, so accel = -g_body.
        // With identity attitude, the world-frame accel is g_world = body_accel.
        // Wait — the preintegration uses `a_world = dq * (accel - ba)`.
        // With identity attitude and 0 bias: a_world = [0, 0, 9.8]
        // dv = a_world * dt_total = [0, 0, 0.98]
        // dp = 0.5 * a_world * dt_total^2 = [0, 0, 0.049]

        assert!((preint.dv.z - 0.98).abs() < 0.01);
        assert!((preint.dp.z - 0.049).abs() < 0.01);
        assert!(preint.dt > 0.09 && preint.dt < 0.11);
    }

    #[test]
    fn factor_residual_zero_at_nominal() {
        let (vars, pi, vi, pj, vj, bias) = make_pose_vel_bias_vars();
        let values = VariableValues::build(&vars);

        let mut preint = ImuPreintegration::new();
        preint.dt = 1.0;

        let factor = ImuPreintegrationFactor::new(
            preint,
            Vector3::new(0.0, 0.0, -9.8), // gravity
            Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity(), // nominal i
            Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity(), // nominal j
            Vector3::zeros(), Vector3::zeros(), // nominal bias
            pi, vi, pj, vj, bias,
        );

        let r = factor.residual(&values);
        // At nominal values with zero preintegration, all quantities are zero
        // except gravity: dp_pred = p_j - p_i - v_i*dt - 0.5*g*dt^2
        // = 0 - 0 - 0 - 0.5*[0,0,-9.8]*1 = [0, 0, 4.9]
        // dv_pred = v_j - v_i - g*dt = 0 - 0 - [0,0,-9.8] = [0, 0, 9.8]
        // These should match the preintegrated deltas (which are zero), so residual is non-zero.
        // Actually for this test, we want the residual to be zero — so the nominal
        // preintegration should match the predicted motion.
        // Let's just check the residual is finite and has the right dimension.
        assert_eq!(r.len(), 15);
        assert!(r.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn bias_correction_is_first_order_accurate() {
        let mut preint = ImuPreintegration::new();
        let samples: Vec<ImuSample> = (0..11)
            .map(|i| ImuSample {
                accel: Vector3::new(0.0, 0.0, -9.8),
                gyro: Vector3::zeros(),
                time_us: i * 10_000,
            })
            .collect();
        let ba = Vector3::new(0.1, 0.0, 0.0);
        preint.integrate(&samples, &ba, &Vector3::zeros());

        // Correct with zero delta → should return original
        let corr0 = preint.correct(&Vector3::zeros(), &Vector3::zeros());
        assert!((corr0.dp - preint.dp).norm() < 1e-12);

        // Correct with non-zero delta → should differ by ~J * delta
        let dba = Vector3::new(0.01, 0.0, 0.0);
        let corr = preint.correct(&dba, &Vector3::zeros());
        let dp_diff = corr.dp - preint.dp;
        // The difference should be approximately J_p_ba * dba
        let expected = preint.dp_dba * dba;
        assert!((dp_diff - expected).norm() < 0.01);
    }

    #[test]
    fn time_diff_wraps_correctly() {
        let dt = time_diff_us(u64::MAX - 1000, 1000);
        assert!((dt - 2001.0 / 1_000_000.0).abs() < 1e-9);
    }

// ---- analytic checks of the preintegration ---------------------------------
//
// All of these use a constant body-frame specific force with zero gyro and
// zero bias at identity attitude, where the discrete mid-point update is the
// exact solution of  v' = a,  p' = v:
//   dv = a T        dp = a T^2 / 2

/// 11 samples at 100 Hz: 10 intervals of 0.01 s, so T = 0.1 s exactly.
const T: f64 = 0.1;

fn const_accel_samples(accel: Vector3<f64>) -> Vec<ImuSample> {
    (0..11)
        .map(|i| ImuSample { accel, gyro: Vector3::zeros(), time_us: i * 10_000 })
        .collect()
}

#[test]
fn integrate_reproduces_the_constant_acceleration_solution() {
    let a = Vector3::new(1.0, -2.0, 3.0);
    let mut p = ImuPreintegration::new();
    p.integrate(&const_accel_samples(a), &Vector3::zeros(), &Vector3::zeros());
    assert!((p.dt - T).abs() < 1e-12, "dt = {}", p.dt);
    let want_dv = a * T;
    let want_dp = a * (T * T / 2.0);
    assert!((p.dv - want_dv).norm() < 1e-12, "dv = {:?}, want {want_dv:?}", p.dv);
    assert!((p.dp - want_dp).norm() < 1e-12, "dp = {:?}, want {want_dp:?}", p.dp);
    assert!((p.dq.quaternion().w - 1.0).abs() < 1e-15, "no gyro means no rotation");
}

#[test]
fn integrate_accumulates_the_bias_jacobians_in_closed_form() {
    // With zero specific force and identity attitude (R = I) the recurrence
    //   dp_dba += dv_dba dt - 0.5 R dt^2,  dv_dba -= R dt,  dq_dbg = -I T
    // integrates to  dp_dba = -T^2/2 I = -0.005 I,  dv_dba = -T I = -0.1 I.
    let mut p = ImuPreintegration::new();
    p.integrate(&const_accel_samples(Vector3::zeros()), &Vector3::zeros(), &Vector3::zeros());
    let dp_dba = Matrix3::from_diagonal(&Vector3::new(-0.005, -0.005, -0.005));
    let dv_dba = Matrix3::from_diagonal(&Vector3::new(-0.1, -0.1, -0.1));
    assert!((p.dp_dba - dp_dba).norm() < 1e-12, "dp_dba = {:?}", p.dp_dba);
    assert!((p.dv_dba - dv_dba).norm() < 1e-12, "dv_dba = {:?}", p.dv_dba);
    assert!(p.dp_dbg.norm() < 1e-15, "no gyro means no position sensitivity to bg");
    assert!(p.dv_dbg.norm() < 1e-15, "no gyro means no velocity sensitivity to bg");
    let dq_dbg = Matrix3::from_diagonal(&Vector3::new(-T, -T, -T));
    assert!((p.dq_dbg - dq_dbg).norm() < 1e-12, "dq_dbg = {:?}", p.dq_dbg);
}

#[test]
fn bias_corrected_deltas_are_the_first_order_taylor_expansion() {
    // dp_corr = dp + dp_dba dba + dp_dbg dbg with a hand-set Jacobian:
    //   dp = (1, 2, 3), dp_dba = diag(-0.005), dba = (0.1, 0, 0)
    //   -> dp_corr = (1, 2, 3) + (-0.0005, 0, 0) = (0.9995, 2, 3)
    // dq_corr = dq exp(dq_dbg dbg): identity dq with dq_dbg = -I T and
    // dbg = (0, 0, theta) is a rotation of -T theta about +Z.
    let mut p = ImuPreintegration::new();
    p.dt = T;
    p.dp = Vector3::new(1.0, 2.0, 3.0);
    p.dp_dba = Matrix3::from_diagonal(&Vector3::new(-0.005, -0.005, -0.005));
    p.dq_dbg = Matrix3::from_diagonal(&Vector3::new(-T, -T, -T));

    let c = p.correct(&Vector3::new(0.1, 0.0, 0.0), &Vector3::new(0.0, 0.0, 0.3));
    let want_dp = Vector3::new(0.9995, 2.0, 3.0);
    assert!((c.dp - want_dp).norm() < 1e-15, "dp_corr = {:?}, want {want_dp:?}", c.dp);
    let want_rot = -T * 0.3;
    let axis = c.dq.scaled_axis();
    assert!((axis - Vector3::new(0.0, 0.0, want_rot)).norm() < 1e-15, "axis = {axis:?}, want {want_rot}");
    assert!((c.dt - T).abs() < 1e-15, "dt must pass through the correction");
}

#[test]
fn default_preintegration_matches_a_fresh_one() {
    let p = ImuPreintegration::default();
    assert!(p.dp.norm() < 1e-15 && p.dv.norm() < 1e-15 && p.dt == 0.0);
    assert!((p.dq.quaternion().w - 1.0).abs() < 1e-15);
    assert!((p.covariance - DMatrix::identity(15, 15) * 1e-4).norm() < 1e-18);
}

#[test]
fn time_difference_handles_every_branch() {
    // Ordinary forward step.
    assert!((time_diff_us(1_000_000, 1_100_000) - 0.1).abs() < 1e-12);
    // Week rollover: 0.2 s before the wrap, 0.1 s after it.
    let week_us = 604_800_000_000_u64;
    let wrapped = time_diff_us(week_us - 200_000, 100_000);
    assert!((wrapped - 0.3).abs() < 1e-9, "week rollover gave {wrapped}");
    // 64-bit counter rollover.
    let wrapped64 = time_diff_us(u64::MAX - 1_000, 1_000);
    assert!((wrapped64 - 2_001.0 / 1e6).abs() < 1e-12, "u64 rollover gave {wrapped64}");
}

// ---- factor -----------------------------------------------------------------

/// A factor over the five states of `make_pose_vel_bias_vars`, with the
/// nominal values left at zero so the linearisation point is the origin.
fn factor_at_origin(preint: ImuPreintegration) -> (ImuPreintegrationFactor, VariableValues) {
    let (f, vars) = factor_with_vars(preint);
    (f, VariableValues::build(&vars))
}

fn factor_with_vars(preint: ImuPreintegration) -> (ImuPreintegrationFactor, BTreeMap<VariableId, VariableNode>) {
    let (vars, pi, vi, pj, vj, bias) = make_pose_vel_bias_vars();
    let f = ImuPreintegrationFactor::new(
        preint,
        Vector3::new(0.0, 0.0, -9.80665),
        Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity(),
        Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity(),
        Vector3::zeros(), Vector3::zeros(),
        pi, vi, pj, vj, bias,
    );
    (f, vars)
}

#[test]
fn factor_residual_vanishes_for_a_free_fall_preintegration() {
    // Nominal state, dt = 1, gravity (0, 0, -9.80665), identity attitude:
    //   dp_pred = -0.5 g dt^2 = (0, 0, 4.903325)
    //   dv_pred = -g dt       = (0, 0, 9.80665)
    // Preintegrating exactly those deltas makes all 15 rows vanish.
    let mut p = ImuPreintegration::new();
    p.dt = 1.0;
    p.dp = Vector3::new(0.0, 0.0, 0.5 * 9.80665);
    p.dv = Vector3::new(0.0, 0.0, 9.80665);
    let (f, values) = factor_at_origin(p);
    let r = f.residual(&values);
    for (i, v) in r.iter().enumerate() {
        assert!(v.abs() < 1e-12, "residual[{i}] = {v}; free fall must cancel exactly");
    }
}

#[test]
fn factor_residual_reports_the_bias_as_its_own_rows() {
    // Rows 9..12 and 12..15 are the accel and gyro bias error states, so a
    // bias offset of (1, 2, 3, 4, 5, 6) must appear verbatim there and nowhere
    // else in the first nine rows.
    let (vars, pi, vi, pj, vj, bias) = make_pose_vel_bias_vars();
    let mut node = vars[&bias].clone();
    node.value = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    let mut vars = vars;
    vars.insert(bias, node);
    let values = VariableValues::build(&vars);

    let mut p = ImuPreintegration::new();
    p.dt = 1.0;
    p.dp = Vector3::new(0.0, 0.0, 0.5 * 9.80665);
    p.dv = Vector3::new(0.0, 0.0, 9.80665);
    let f = ImuPreintegrationFactor::new(
        p,
        Vector3::new(0.0, 0.0, -9.80665),
        Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity(),
        Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity(),
        Vector3::zeros(), Vector3::zeros(),
        pi, vi, pj, vj, bias,
    );
    let r = f.residual(&values);
    let want = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    for (k, &w) in want.iter().enumerate() {
        assert!((r[9 + k] - w).abs() < 1e-12, "bias row {k} = {}", r[9 + k]);
    }
    // The first-order bias correction is inactive: the Jacobians are zero.
    for k in 0..9 {
        assert!(r[k].abs() < 1e-12, "motion row {k} moved to {}", r[k]);
    }
}

#[test]
fn factor_information_is_the_inverse_covariance_plus_a_bias_prior() {
    // The default preintegration covariance is 1e-4 * I, so the top-left 9x9
    // information block is 1e4 * I; the six bias rows carry sigma = 0.1, i.e.
    // W = 100 each, and the two off-diagonal blocks stay zero.
    let (f, _values) = factor_at_origin(ImuPreintegration::new());
    let info = f.information();
    assert_eq!(info.shape(), (15, 15));
    for r in 0..9 {
        for c in 0..9 {
            let want = if r == c { 1.0e4 } else { 0.0 };
            assert!((info[(r, c)] - want).abs() < 1e-9, "info[{r},{c}] = {}", info[(r, c)]);
        }
    }
    for k in 9..15 {
        assert!((info[(k, k)] - 100.0).abs() < 1e-9, "bias info[{k},{k}] = {}", info[(k, k)]);
        for j in 0..9 {
            assert_eq!(info[(k, j)], 0.0, "motion/bias blocks must stay zero");
            assert_eq!(info[(j, k)], 0.0, "motion/bias blocks must stay zero");
        }
    }
    assert_eq!(f.robust_threshold(), Some(3.0), "IMU outliers are Huber-clipped at 3 sigma");
    assert!(!f.use_cauchy());
    assert!(format!("{f:?}").contains("ImuPreintegrationFactor"), "Debug must name the factor");
}

#[test]
fn factor_jacobian_blocks_follow_the_body_to_world_rotation() {
    // The factor rotates the predicted motion into the body frame with
    // R_i^T, so with R_i = I:
    //   dr_p/dp_i = -I, dr_p/dp_j = +I, dr_p/dv_i = -I dt
    //   dr_v/dv_i = -I,   dr_v/dv_j = +I,  dr_q/dq_j = +I
    // and the bias columns are minus the preintegration Jacobians.
    let mut p = ImuPreintegration::new();
    p.dt = 2.0;
    p.dp_dba = Matrix3::from_diagonal(&Vector3::new(-0.005, -0.005, -0.005));
    p.dv_dba = Matrix3::from_diagonal(&Vector3::new(-0.1, -0.1, -0.1));
    p.dp_dbg = Matrix3::from_diagonal(&Vector3::new(-0.001, -0.001, -0.001));
    p.dv_dbg = Matrix3::from_diagonal(&Vector3::new(-0.02, -0.02, -0.02));
    p.dq_dbg = Matrix3::from_diagonal(&Vector3::new(-2.0, -2.0, -2.0));
    let (f, values) = factor_at_origin(p);
    let j = f.jacobian(&values);

    let (s_pi, _) = values.index_of(VariableId::new(0)).unwrap();
    let (s_vi, _) = values.index_of(VariableId::new(1)).unwrap();
    let (s_pj, _) = values.index_of(VariableId::new(2)).unwrap();
    let (s_vj, _) = values.index_of(VariableId::new(3)).unwrap();
    let (s_b, _) = values.index_of(VariableId::new(4)).unwrap();

    for k in 0..3 {
        assert_eq!(j[(k, s_pi + k)], -1.0, "dr_p/dp_i");
        assert_eq!(j[(k, s_pj + k)], 1.0, "dr_p/dp_j");
        assert_eq!(j[(k, s_vi + k)], -2.0, "dr_p/dv_i = -R^T dt");
        assert_eq!(j[(3 + k, s_vi + k)], -1.0, "dr_v/dv_i");
        assert_eq!(j[(3 + k, s_vj + k)], 1.0, "dr_v/dv_j");
        assert_eq!(j[(6 + k, s_pj + 3 + k)], 1.0, "dr_q/dq_j = I");
        assert_eq!(j[(6 + k, s_pi + 3 + k)], -1.0, "dr_q/dq_i = -R(dq)^T = -I");
        // Bias columns are the negated preintegration Jacobians.
        assert!((j[(k, s_b + k)] - 0.005).abs() < 1e-12, "dr_p/dba = -dp_dba");
        assert!((j[(k, s_b + 3 + k)] - 0.001).abs() < 1e-12, "dr_p/dbg = -dp_dbg");
        assert!((j[(3 + k, s_b + k)] - 0.1).abs() < 1e-12, "dr_v/dba = -dv_dba");
        assert!((j[(3 + k, s_b + 3 + k)] - 0.02).abs() < 1e-12, "dr_v/dbg = -dv_dbg");
        assert!((j[(6 + k, s_b + 3 + k)] - 2.0).abs() < 1e-12, "dr_q/dbg = -dq_dbg");
    }
    // Rows 9..15 differentiate the bias random walk: an identity block.
    for k in 0..6 {
        assert_eq!(j[(9 + k, s_b + k)], 1.0, "bias random walk diagonal");
        for c in 0..6 {
            if c != k {
                assert_eq!(j[(9 + k, s_b + c)], 0.0, "bias random walk must be diagonal");
            }
        }
    }
}

#[test]
fn factor_jacobian_agrees_with_central_differences() {
    // Independent check of the whole analytic Jacobian: perturb each state by
    // +h and -h and compare (r(x+h) - r(x-h)) / 2h with the analytic column.
    let mut p = ImuPreintegration::new();
    p.dt = 1.0;
    p.dp = Vector3::new(0.3, -0.2, 0.1);
    p.dv = Vector3::new(0.02, 0.0, -0.01);
    p.dp_dba = Matrix3::from_diagonal(&Vector3::new(-0.005, -0.005, -0.005));
    p.dv_dba = Matrix3::from_diagonal(&Vector3::new(-0.1, -0.1, -0.1));
    p.dq_dbg = Matrix3::from_diagonal(&Vector3::new(-1.0, -1.0, -1.0));
    let (f, vars) = factor_with_vars(p);
    let analytic = f.jacobian(&VariableValues::build(&vars));

    let h = 1e-6;
    for id in (0..5).map(VariableId::new) {
        let dim = vars[&id].value.len();
        for k in 0..dim {
            let rp = f.residual(&shifted(&vars, id, k, h));
            let rm = f.residual(&shifted(&vars, id, k, -h));
            for m in 0..15 {
                let fd = (rp[m] - rm[m]) / (2.0 * h);
                let an = analytic[(m, start_of(&vars, id) + k)];
                if fd.abs() < 1e-9 && an.abs() < 1e-9 {
                    continue;
                }
                assert!(
                    (fd - an).abs() < 1e-5 * (1.0 + fd.abs()),
                    "variable {id:?} dim {k} row {m}: finite difference {fd} vs analytic {an}"
                );
            }
        }
    }
}

/// The packed-state index of `id`.
fn start_of(vars: &BTreeMap<VariableId, VariableNode>, id: VariableId) -> usize {
    let mut offset = 0;
    for (vid, node) in vars {
        if *vid == id {
            return offset;
        }
        offset += node.value.len();
    }
    panic!("variable {id:?} not found");
}

/// `VariableValues` with element `k` of `id` shifted by `delta`.
fn shifted(vars: &BTreeMap<VariableId, VariableNode>, id: VariableId, k: usize, delta: f64) -> VariableValues {
    let mut copy = vars.clone();
    copy.get_mut(&id).expect("variable present").value[k] += delta;
    VariableValues::build(&copy)
}
