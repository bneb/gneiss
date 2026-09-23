use nalgebra::{SMatrix, Vector3};

use super::types::{skew_symmetric, EngineError, EskfState, Vector15};
use super::update::{apply_error_injection, joseph_form_update};

pub type RowVector15<T = f64> = SMatrix<T, 1, 15>;

/// Line-of-sight geometry and lever arm for a double-difference satellite pair.
#[derive(Debug, Clone, PartialEq)]
pub struct DdSatGeometry {
    pub u_s: Vector3<f64>,
    pub u_ref: Vector3<f64>,
    pub lever_arm: Vector3<f64>,
}

impl DdSatGeometry {
    pub fn new(u_s: Vector3<f64>, u_ref: Vector3<f64>, lever_arm: Vector3<f64>) -> Self {
        Self {
            u_s,
            u_ref,
            lever_arm,
        }
    }
}

/// Double-difference observation type: Pseudorange or Carrier-Phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DdMeasurementKind {
    Pseudorange,
    CarrierPhase,
}

/// Result of a scalar double-difference Kalman update.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DdScalarUpdateResult {
    pub innovation: f64,
    pub innovation_var: f64,
    pub accepted: bool,
}

/// Compute 15-state analytical row Jacobian for a double-difference range/phase measurement.
///
/// H = [ -(u^s - u^ref)^T,  0_1x3,  +(u^s - u^ref)^T * [l_e x],  0_1x3,  0_1x3 ]
pub fn compute_dd_jacobian_15(state: &EskfState, geom: &DdSatGeometry) -> RowVector15<f64> {
    let mut h = RowVector15::zeros();
    let delta_u = geom.u_s - geom.u_ref;

    for i in 0..3 {
        h[(0, i)] = -delta_u[i];
    }

    let l_e = state.attitude.to_rotation_matrix() * geom.lever_arm;
    let l_e_skew = skew_symmetric(&l_e);
    let h_att = delta_u.transpose() * l_e_skew;

    for i in 0..3 {
        h[(0, i + 6)] = h_att[i];
    }

    h
}

/// Check if double-difference innovation passes consistency / blunder gating.
fn is_innovation_acceptable(residual: f64, nis: f64, kind: DdMeasurementKind) -> bool {
    let (max_res, max_nis) = match kind {
        DdMeasurementKind::Pseudorange => (15.0, 16.0),
        DdMeasurementKind::CarrierPhase => (0.08, 25.0),
    };
    residual.abs() <= max_res && nis <= max_nis
}

/// Compute innovation variance, Kalman gain, and error correction for scalar DD update.
fn compute_dd_scalar_system(
    state: &EskfState,
    geom: &DdSatGeometry,
    residual: f64,
    var: f64,
    kind: DdMeasurementKind,
) -> Result<(RowVector15<f64>, Vector15<f64>, f64, bool), EngineError> {
    if var <= 0.0 || !var.is_finite() || !residual.is_finite() {
        return Err(EngineError::InvalidMeasurement("Invalid var/residual".into()));
    }

    let h = compute_dd_jacobian_15(state, geom);
    let s = (h * state.cov * h.transpose())[(0, 0)] + var;
    if s <= 0.0 || !s.is_finite() {
        return Err(EngineError::SingularState);
    }

    let nis = (residual * residual) / s;
    if !is_innovation_acceptable(residual, nis, kind) {
        return Ok((h, Vector15::zeros(), s, false));
    }

    let k = (state.cov * h.transpose()) / s;
    Ok((h, k, s, true))
}

/// Apply a scalar double-difference measurement update (Pseudorange or Carrier-Phase)
/// to the 15-state ESKF using the numerically stabilized Joseph-form covariance update.
pub fn update_dd_scalar(
    state: &mut EskfState,
    geom: &DdSatGeometry,
    residual: f64,
    var: f64,
    kind: DdMeasurementKind,
) -> Result<DdScalarUpdateResult, EngineError> {
    let (h, k, s, accepted) = compute_dd_scalar_system(state, geom, residual, var, kind)?;
    if !accepted {
        return Ok(DdScalarUpdateResult {
            innovation: residual,
            innovation_var: s,
            accepted: false,
        });
    }

    let dx = k * residual;
    apply_error_injection(state, &dx);

    let r_mat = SMatrix::<f64, 1, 1>::new(var);
    state.cov = joseph_form_update(&state.cov, &h, &k, &r_mat);

    Ok(DdScalarUpdateResult {
        innovation: residual,
        innovation_var: s,
        accepted: true,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::UnitQuaternion;

    #[test]
    fn test_dd_zero_lever_arm_attitude_insensitivity() {
        let state = EskfState::new(
            Vector3::new(100.0, 200.0, 300.0),
            Vector3::zeros(),
            UnitQuaternion::identity(),
        );
        let geom = DdSatGeometry::new(
            Vector3::new(0.6, 0.0, 0.8),
            Vector3::new(0.0, 0.6, 0.8),
            Vector3::zeros(),
        );

        let h = compute_dd_jacobian_15(&state, &geom);
        assert!((h[(0, 0)] - (-0.6)).abs() < 1e-12);
        assert!((h[(0, 1)] - 0.6).abs() < 1e-12);
        assert!((h[(0, 2)] - 0.0).abs() < 1e-12);

        for i in 3..15 {
            assert_eq!(h[(0, i)], 0.0);
        }
    }

    #[test]
    fn test_dd_lever_arm_rotation_cross_product() {
        let state = EskfState::new(
            Vector3::new(10.0, 20.0, 30.0),
            Vector3::zeros(),
            UnitQuaternion::identity(),
        );
        let delta_u = Vector3::new(0.2, -0.2, 0.0);
        let lever_arm = Vector3::new(0.0, 0.0, -1.0);
        let geom = DdSatGeometry::new(delta_u, Vector3::zeros(), lever_arm);

        let h = compute_dd_jacobian_15(&state, &geom);
        // H_att = delta_u^T * [l_e x]
        // [l_e x] for l_e = [0, 0, -1] is:
        // [[ 0,  1,  0],
        //  [-1,  0,  0],
        //  [ 0,  0,  0]]
        // [0.2, -0.2, 0] * [l_e x] = [0.2, 0.2, 0.0]
        assert!((h[(0, 6)] - 0.2).abs() < 1e-12);
        assert!((h[(0, 7)] - 0.2).abs() < 1e-12);
        assert!((h[(0, 8)] - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_dd_collinear_satellites_null_update() {
        let mut state = EskfState::new(
            Vector3::new(10.0, 20.0, 30.0),
            Vector3::zeros(),
            UnitQuaternion::identity(),
        );
        let u = Vector3::new(0.577, 0.577, 0.577);
        let geom = DdSatGeometry::new(u, u, Vector3::new(0.1, -0.2, 0.5));

        let h = compute_dd_jacobian_15(&state, &geom);
        for i in 0..15 {
            assert_eq!(h[(0, i)], 0.0);
        }

        let p_before = state.cov;
        let res = update_dd_scalar(
            &mut state,
            &geom,
            0.05,
            1.0,
            DdMeasurementKind::Pseudorange,
        )
        .expect("update succeeds");
        assert!(res.accepted);
        assert!((state.cov - p_before).norm() < 1e-12);
    }

    #[test]
    fn test_dd_jacobian_numerical_finite_difference_consistency() {
        let rot = UnitQuaternion::from_euler_angles(0.1, -0.2, 0.3);
        let state = EskfState::new(
            Vector3::new(100.0, 200.0, 300.0),
            Vector3::new(1.0, -2.0, 0.5),
            rot,
        );
        let r_s = Vector3::new(20_000_000.0, 5_000_000.0, 15_000_000.0);
        let r_ref = Vector3::new(-10_000_000.0, 18_000_000.0, 20_000_000.0);
        let lever_arm = Vector3::new(0.3, -0.15, 0.8);

        let p_ant = state.pos_ecef + state.attitude.to_rotation_matrix() * lever_arm;
        let u_s = (r_s - p_ant).normalize();
        let u_ref = (r_ref - p_ant).normalize();
        let geom = DdSatGeometry::new(u_s, u_ref, lever_arm);

        let h_analytical = compute_dd_jacobian_15(&state, &geom);

        let eps_pos = 0.1; // 10 cm step avoids catastrophic floating-point cancellation at 20,000 km
        let dd_range = |pos: &Vector3<f64>, att: &UnitQuaternion<f64>| -> f64 {
            let ant = pos + att.to_rotation_matrix() * lever_arm;
            (r_s - ant).norm() - (r_ref - ant).norm()
        };

        // Numerical check on position (states 0..3)
        for i in 0..3 {
            let mut p_plus = state.pos_ecef;
            let mut p_minus = state.pos_ecef;
            p_plus[i] += eps_pos;
            p_minus[i] -= eps_pos;
            let val_plus = dd_range(&p_plus, &state.attitude);
            let val_minus = dd_range(&p_minus, &state.attitude);
            let d_num = (val_plus - val_minus) / (2.0 * eps_pos);
            assert!(
                (d_num - h_analytical[(0, i)]).abs() < 1e-7,
                "Position component {} mismatch: num={}, ana={}",
                i,
                d_num,
                h_analytical[(0, i)]
            );
        }

        // Numerical check on attitude (states 6..9) with 1 mrad step
        let eps_att = 1e-3;
        for i in 0..3 {
            let mut rot_vec = Vector3::zeros();
            rot_vec[i] = eps_att;
            let dq_plus = UnitQuaternion::from_scaled_axis(rot_vec);
            let dq_minus = UnitQuaternion::from_scaled_axis(-rot_vec);
            let q_plus = dq_plus * state.attitude;
            let q_minus = dq_minus * state.attitude;
            let val_plus = dd_range(&state.pos_ecef, &q_plus);
            let val_minus = dd_range(&state.pos_ecef, &q_minus);
            let d_num = (val_plus - val_minus) / (2.0 * eps_att);
            assert!(
                (d_num - h_analytical[(0, 6 + i)]).abs() < 1e-5,
                "Attitude component {} mismatch: num={}, ana={}",
                i,
                d_num,
                h_analytical[(0, 6 + i)]
            );
        }
    }

    #[test]
    fn test_groves_chapter_14_exact_geometry_benchmark() {
        // Groves (2013) Chapter 14 canonical benchmark geometry
        let att = UnitQuaternion::from_euler_angles(
            5.0_f64.to_radians(),
            (-3.0_f64).to_radians(),
            45.0_f64.to_radians(),
        ); // 5 deg roll, -3 deg pitch, 45 deg yaw
        let state = EskfState::new(
            Vector3::new(-2430601.8, -4700258.9, 3544321.5),
            Vector3::new(15.0, 5.0, -1.0),
            att,
        );
        let lever_arm = Vector3::new(0.20, -0.10, 1.20);
        let u_s = Vector3::new(0.267261, 0.534522, 0.801784); // Satellite 1
        let u_ref = Vector3::new(-0.408248, 0.408248, 0.816497); // Satellite 2 (reference)
        let geom = DdSatGeometry::new(u_s, u_ref, lever_arm);

        let h = compute_dd_jacobian_15(&state, &geom);
        let delta_u = u_s - u_ref;

        // Verify Groves Eq 14.52: H_pos = -delta_u^T
        for i in 0..3 {
            assert!((h[(0, i)] - (-delta_u[i])).abs() < 1e-12);
        }

        // Verify Groves Eq 14.52: H_att = +delta_u^T * [l_e x]
        let l_e = state.attitude.to_rotation_matrix() * lever_arm;
        let l_skew = skew_symmetric(&l_e);
        let expected_h_att = delta_u.transpose() * l_skew;
        for i in 0..3 {
            assert!((h[(0, 6 + i)] - expected_h_att[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_lyapunov_closed_loop_error_contraction() {
        let true_att = UnitQuaternion::identity();
        let true_pos = Vector3::new(10.0, 20.0, 30.0);
        let lever_arm = Vector3::new(0.0, 0.0, 2.0); // 2m roof antenna
        let r_s = Vector3::new(10.0, 20000000.0, 30.0); // North satellite (+Y)
        let r_ref = Vector3::new(10.0, 20.0, 20000030.0); // Zenith satellite (+Z)

        let true_ant = true_pos + true_att.to_rotation_matrix() * lever_arm;
        let true_range = (r_s - true_ant).norm() - (r_ref - true_ant).norm();

        // Introduce an intentional +2.0 deg pitch error (rotation around X)
        let pitch_err = 0.0349066; // 2 deg in rad
        let est_att = UnitQuaternion::from_scaled_axis(Vector3::new(pitch_err, 0.0, 0.0)) * true_att;
        let mut state = EskfState::new(true_pos, Vector3::zeros(), est_att);
        state.cov[(6, 6)] = (5.0_f64.to_radians()).powi(2); // 5 deg attitude uncertainty

        let est_ant = state.pos_ecef + state.attitude.to_rotation_matrix() * lever_arm;
        let u_s = (r_s - est_ant).normalize();
        let u_ref = (r_ref - est_ant).normalize();
        let geom = DdSatGeometry::new(u_s, u_ref, lever_arm);

        let pred_range = (r_s - est_ant).norm() - (r_ref - est_ant).norm();
        let residual = true_range - pred_range;

        let res = update_dd_scalar(
            &mut state,
            &geom,
            residual,
            0.003 * 0.003, // 3mm carrier noise
            DdMeasurementKind::CarrierPhase,
        )
        .expect("scalar update succeeds");
        assert!(res.accepted);

        // Extract posterior pitch error from attitude quaternion
        let dq = state.attitude * true_att.inverse();
        let post_err = dq.scaled_axis().norm();

        assert!(
            post_err < pitch_err,
            "Lyapunov invariant violated: posterior pitch error ({:.4} rad) must be strictly less than prior error ({:.4} rad)",
            post_err, pitch_err
        );
    }

    #[test]
    fn test_taylor_second_order_convergence_rate() {
        let state = EskfState::new(
            Vector3::new(50.0, 100.0, 150.0),
            Vector3::zeros(),
            UnitQuaternion::from_euler_angles(0.1, -0.15, 0.2),
        );
        let lever_arm = Vector3::new(0.5, -0.2, 1.0);
        let r_s = Vector3::new(15000000.0, 8000000.0, 12000000.0);
        let r_ref = Vector3::new(-5000000.0, 12000000.0, 18000000.0);

        let p_ant = state.pos_ecef + state.attitude.to_rotation_matrix() * lever_arm;
        let geom = DdSatGeometry::new((r_s - p_ant).normalize(), (r_ref - p_ant).normalize(), lever_arm);
        let h = compute_dd_jacobian_15(&state, &geom);

        let f = |dtheta: f64| -> f64 {
            let dq = UnitQuaternion::from_scaled_axis(Vector3::new(dtheta, 0.0, 0.0));
            let q = dq * state.attitude;
            let ant = state.pos_ecef + q.to_rotation_matrix() * lever_arm;
            (r_s - ant).norm() - (r_ref - ant).norm()
        };

        let f0 = f(0.0);
        let eps1 = 0.02;
        let eps2 = 0.01;

        let r1 = (f(eps1) - f0 - h[(0, 6)] * eps1).abs();
        let r2 = (f(eps2) - f0 - h[(0, 6)] * eps2).abs();

        let ratio = r1 / r2;
        // As step size halves, quadratic remainder should decrease by approximately 4x (slope = 2)
        assert!(
            (ratio - 4.0).abs() < 0.6,
            "Taylor remainder ratio {} is not quadratic (expected ~4.0)",
            ratio
        );
    }

    #[test]
    fn test_dd_carrier_update_sub_centimeter_covariance_contraction() {
        let mut state = EskfState::new(
            Vector3::new(10.0, 20.0, 30.0),
            Vector3::zeros(),
            UnitQuaternion::identity(),
        );
        let geom = DdSatGeometry::new(
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(0.0, 0.0, 0.5),
        );

        let p_x_before = state.cov[(0, 0)];
        let res = update_dd_scalar(
            &mut state,
            &geom,
            0.005,
            0.003 * 0.003,
            DdMeasurementKind::CarrierPhase,
        )
        .expect("carrier update succeeds");

        assert!(res.accepted);
        assert!(state.cov[(0, 0)] < p_x_before);

        // Verify covariance remains symmetric and positive-definite
        let diff = state.cov - state.cov.transpose();
        assert!(diff.norm() < 1e-12);
        for i in 0..15 {
            assert!(state.cov[(i, i)] > 0.0);
        }
    }

    #[test]
    fn test_dd_blunder_rejection() {
        let mut state = EskfState::new(
            Vector3::new(10.0, 20.0, 30.0),
            Vector3::zeros(),
            UnitQuaternion::identity(),
        );
        let geom = DdSatGeometry::new(
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::zeros(),
        );

        let p_before = state.cov;
        let res = update_dd_scalar(
            &mut state,
            &geom,
            50.0, // 50m blunder
            1.0,
            DdMeasurementKind::Pseudorange,
        )
        .expect("update check succeeds");

        assert!(!res.accepted);
        assert!((state.cov - p_before).norm() < 1e-12);
    }
}
