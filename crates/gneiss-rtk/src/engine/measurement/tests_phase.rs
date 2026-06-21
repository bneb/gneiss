use crate::engine::measurement_math::{doppler_attitude_jacobian, range_attitude_jacobian};
use crate::filter::{DdObservation, RtkState};
use gneiss_core::coords::{Coordinate, Datum, Frame};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use nalgebra::Vector3;

#[test]
fn test_phase_windup_correction_sign_rtk() {
    use super::windup::apply_windup_to_obs;

    let mut obs = DdObservation {
        sat: SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        },
        pr_l1: 0.0,
        pr_l2: None,
        cp_l1: Some(10.0),
        cp_l2: Some(20.0),
        doppler: 0.0,
        snr: 45.0,
        locktime: None,
    };

    let windup = 0.25; // 0.25 cycles of positive wind-up
    apply_windup_to_obs(&mut obs, windup);

    // Corrected carrier phase = raw_cp - windup
    assert_eq!(obs.cp_l1.unwrap(), 9.75);
    assert_eq!(obs.cp_l2.unwrap(), 19.75);
}

#[test]
fn test_compute_phase_windup() {
    use super::{DdContext, SatState};

    let time = GpsTime::new(2137, 422922.0);
    let state = RtkState::new(
        time,
        Coordinate::new(
            Vector3::new(1000.0, 2000.0, 3000.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ),
        10.0,
    );

    let sat1 = SatelliteId {
        constellation: Constellation::Gps,
        prn: 1,
    };
    let sat2 = SatelliteId {
        constellation: Constellation::Gps,
        prn: 2,
    };

    let mut rov_sat = DdObservation {
        sat: sat1,
        pr_l1: 0.0,
        pr_l2: None,
        cp_l1: Some(10.0),
        cp_l2: Some(20.0),
        doppler: 0.0,
        snr: 45.0,
        locktime: None,
    };
    let mut base_sat = DdObservation {
        sat: sat1,
        pr_l1: 0.0,
        pr_l2: None,
        cp_l1: Some(10.0),
        cp_l2: Some(20.0),
        doppler: 0.0,
        snr: 45.0,
        locktime: None,
    };
    let mut rov_ref = DdObservation {
        sat: sat2,
        pr_l1: 0.0,
        pr_l2: None,
        cp_l1: Some(10.0),
        cp_l2: Some(20.0),
        doppler: 0.0,
        snr: 45.0,
        locktime: None,
    };
    let mut ref_base = DdObservation {
        sat: sat2,
        pr_l1: 0.0,
        pr_l2: None,
        cp_l1: Some(10.0),
        cp_l2: Some(20.0),
        doppler: 0.0,
        snr: 45.0,
        locktime: None,
    };

    let sat_state = SatState {
        rov_pos: Vector3::new(20000000.0, 0.0, 0.0),
        rov_vel: Vector3::zeros(),
        bas_pos: Vector3::new(0.0, 20000000.0, 0.0),
        bas_vel: Vector3::zeros(),
        f1: 1575.42e6,
        f2: 1227.60e6,
    };
    let ref_state = SatState {
        rov_pos: Vector3::new(0.0, 20000000.0, 0.0),
        rov_vel: Vector3::zeros(),
        bas_pos: Vector3::new(20000000.0, 0.0, 0.0),
        bas_vel: Vector3::zeros(),
        f1: 1575.42e6,
        f2: 1227.60e6,
    };

    let ctx = DdContext {
        rov_sat: &mut rov_sat,
        base_sat: &mut base_sat,
        rov_ref: &mut rov_ref,
        ref_base: &mut ref_base,
        sat_state: &sat_state,
        ref_state: &ref_state,
    };

    let sun_pos = gneiss_core::sun::sun_position_ecef(state.time);
    let crate::engine::measurement_math::WindupUpdates {
        w_sat,
        w_ref,
        w_bas_sat,
        w_bas_ref,
    } = crate::engine::measurement_math::compute_phase_windup(
        Vector3::zeros(),
        Vector3::zeros(),
        sun_pos,
        ctx.sat_state.rov_pos,
        ctx.ref_state.rov_pos,
        ctx.sat_state.bas_pos,
        ctx.ref_state.bas_pos,
        0.0,
        0.0,
        0.0,
        0.0,
    );
    if let Some(cp) = &mut ctx.rov_sat.cp_l1 {
        *cp += w_sat;
    }
    if let Some(cp2) = &mut ctx.rov_sat.cp_l2 {
        *cp2 += w_sat;
    }
    if let Some(cp) = &mut ctx.rov_ref.cp_l1 {
        *cp += w_ref;
    }
    if let Some(cp2) = &mut ctx.rov_ref.cp_l2 {
        *cp2 += w_ref;
    }
    if let Some(cp) = &mut ctx.base_sat.cp_l1 {
        *cp += w_bas_sat;
    }
    if let Some(cp2) = &mut ctx.base_sat.cp_l2 {
        *cp2 += w_bas_sat;
    }
    if let Some(cp) = &mut ctx.ref_base.cp_l1 {
        *cp += w_bas_ref;
    }
    if let Some(cp2) = &mut ctx.ref_base.cp_l2 {
        *cp2 += w_bas_ref;
    }

    assert!(w_sat != 0.0, "windup correction should be non-zero");
    assert_eq!(
        ctx.rov_sat.cp_l1.unwrap(),
        10.0 + w_sat,
        "corrected L1 phase should differ from original by w_sat"
    );
    assert_eq!(
        ctx.rov_sat.cp_l2.unwrap(),
        20.0 + w_sat,
        "corrected L2 phase should differ from original by w_sat"
    );
}

// -----------------------------------------------------------------------
// Numerical Jacobian verification for attitude coupling
// -----------------------------------------------------------------------

/// Compute DD range given a rotation applied to base position + lever arm.
fn dd_range(
    pos: Vector3<f64>,
    rot: nalgebra::UnitQuaternion<f64>,
    lever_body: Vector3<f64>,
    sat_pos: Vector3<f64>,
    ref_pos: Vector3<f64>,
) -> f64 {
    let pos_apc = pos + rot * lever_body;
    (sat_pos - pos_apc).norm() - (ref_pos - pos_apc).norm()
}

/// Compute DD range-rate (rover portion only).
fn dd_range_rate(
    pos: Vector3<f64>,
    vel: Vector3<f64>,
    rot: nalgebra::UnitQuaternion<f64>,
    lever_body: Vector3<f64>,
    omega_b: Vector3<f64>,
    sat_pos: Vector3<f64>,
    ref_pos: Vector3<f64>,
    sat_vel: Vector3<f64>,
    ref_vel: Vector3<f64>,
) -> f64 {
    let pos_apc = pos + rot * lever_body;
    let v_apc = vel + rot * omega_b.cross(&lever_body);
    let e_sat = (sat_pos - pos_apc).normalize();
    let e_ref = (ref_pos - pos_apc).normalize();
    e_sat.dot(&(sat_vel - v_apc)) - e_ref.dot(&(ref_vel - v_apc))
}

/// Apply a small rotation perturbation (left-multiplicative).
fn perturb_attitude(
    rot: nalgebra::UnitQuaternion<f64>,
    d_theta: Vector3<f64>,
) -> nalgebra::UnitQuaternion<f64> {
    let angle = d_theta.norm();
    if angle < 1e-15 {
        return rot;
    }
    let dq = nalgebra::UnitQuaternion::from_axis_angle(
        &nalgebra::Unit::new_normalize(d_theta),
        angle,
    );
    dq * rot
}

#[test]
fn test_range_attitude_jacobian_numerical() {
    let pos = Vector3::new(4_000_000.0, 1_000_000.0, 4_500_000.0);
    let lever_body = Vector3::new(0.0, 0.0, 1.5);
    let rot = nalgebra::UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3);
    let lever_ecef = rot * lever_body;

    let sat_pos = Vector3::new(20_000_000.0, 10_000_000.0, 15_000_000.0);
    let ref_pos = Vector3::new(15_000_000.0, 20_000_000.0, 10_000_000.0);

    let pos_apc = pos + lever_ecef;
    let e_sat = (sat_pos - pos_apc).normalize();
    let e_ref = (ref_pos - pos_apc).normalize();
    let h_r = e_ref - e_sat;

    let j_analytical = range_attitude_jacobian(&lever_ecef, &h_r);

    let eps = 1e-7;
    let mut j_numerical = Vector3::zeros();
    for axis in 0..3 {
        let mut d_theta = Vector3::zeros();
        d_theta[axis] = eps;
        let dd_plus = dd_range(
            pos,
            perturb_attitude(rot, d_theta),
            lever_body,
            sat_pos,
            ref_pos,
        );
        let dd_minus = dd_range(
            pos,
            perturb_attitude(rot, -d_theta),
            lever_body,
            sat_pos,
            ref_pos,
        );
        j_numerical[axis] = (dd_plus - dd_minus) / (2.0 * eps);
    }

    let err = (j_analytical - j_numerical).norm();
    let scale = j_numerical.norm().max(1e-12);
    assert!(
        err / scale < 0.05,
        "Range attitude Jacobian sign/magnitude error!\n  analytical: {:?}\n  numerical:  {:?}\n  rel_error: {:.2e}",
        j_analytical, j_numerical, err / scale
    );
}

#[test]
fn test_range_attitude_jacobian_linearized_exact() {
    let lever_ecef = Vector3::new(0.3, -0.7, 1.2);
    let h_r = Vector3::new(0.4, -0.1, 0.6);
    let j = range_attitude_jacobian(&lever_ecef, &h_r);

    for d_theta in &[
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.3, -0.5, 0.8),
    ] {
        let lhs: f64 = h_r.dot(&d_theta.cross(&lever_ecef));
        let rhs: f64 = j.dot(d_theta);
        assert!(
            (lhs - rhs).abs() < 1e-14,
            "Linearized check failed: lhs={}, rhs={}, δθ={:?}",
            lhs,
            rhs,
            d_theta
        );
    }
}

#[test]
fn test_range_attitude_jacobian_axis_aligned() {
    let lever_ecef = Vector3::new(0.0, 0.0, 1.5);
    let h_r = Vector3::new(1.0, 0.0, 0.0);

    let j = range_attitude_jacobian(&lever_ecef, &h_r);
    assert!((j.x - 0.0).abs() < 1e-12);
    assert!((j.y - 1.5).abs() < 1e-12);
    assert!((j.z - 0.0).abs() < 1e-12);
}

#[test]
fn test_range_attitude_jacobian_zero_lever_arm() {
    let lever_ecef = Vector3::zeros();
    let h_r = Vector3::new(0.3, -0.5, 0.8);
    let j = range_attitude_jacobian(&lever_ecef, &h_r);
    assert!(j.norm() < 1e-15);
}

#[test]
fn test_doppler_attitude_jacobian_numerical() {
    let pos = Vector3::new(4_000_000.0, 1_000_000.0, 4_500_000.0);
    let vel = Vector3::new(1.0, 2.0, 3.0);
    let lever_body = Vector3::new(0.0, 0.0, 1.5);
    let omega_b = Vector3::new(0.01, -0.02, 0.005);
    let rot = nalgebra::UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3);
    let r_b_e = rot.to_rotation_matrix().into_inner();

    let sat_pos = Vector3::new(20_000_000.0, 10_000_000.0, 15_000_000.0);
    let ref_pos = Vector3::new(15_000_000.0, 20_000_000.0, 10_000_000.0);
    let sat_vel = Vector3::new(-500.0, 200.0, 3000.0);
    let ref_vel = Vector3::new(300.0, -800.0, 2500.0);

    let pos_apc = pos + rot * lever_body;
    let e_sat = (sat_pos - pos_apc).normalize();
    let e_ref = (ref_pos - pos_apc).normalize();
    let h_r = e_ref - e_sat;

    let j_analytical = doppler_attitude_jacobian(&r_b_e, &omega_b, &lever_body, &h_r);

    let eps = 1e-7;
    let mut j_numerical = Vector3::zeros();
    for axis in 0..3 {
        let mut d_theta = Vector3::zeros();
        d_theta[axis] = eps;
        let rr_plus = dd_range_rate(
            pos, vel, perturb_attitude(rot, d_theta), lever_body, omega_b,
            sat_pos, ref_pos, sat_vel, ref_vel,
        );
        let rr_minus = dd_range_rate(
            pos, vel, perturb_attitude(rot, -d_theta), lever_body, omega_b,
            sat_pos, ref_pos, sat_vel, ref_vel,
        );
        j_numerical[axis] = (rr_plus - rr_minus) / (2.0 * eps);
    }

    let err = (j_analytical - j_numerical).norm();
    let scale = j_numerical.norm().max(1e-12);
    assert!(
        err / scale < 0.05,
        "Doppler attitude Jacobian sign/magnitude error!\n  analytical: {:?}\n  numerical:  {:?}\n  rel_error: {:.2e}",
        j_analytical, j_numerical, err / scale
    );
}

#[test]
fn test_doppler_attitude_jacobian_linearized_exact() {
    let r_b_e = nalgebra::UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3)
        .to_rotation_matrix()
        .into_inner();
    let omega_b = Vector3::new(0.01, -0.02, 0.005);
    let lever_body = Vector3::new(0.0, 0.0, 1.5);
    let h_r = Vector3::new(0.4, -0.1, 0.6);

    let j = doppler_attitude_jacobian(&r_b_e, &omega_b, &lever_body, &h_r);
    let a = r_b_e * omega_b.cross(&lever_body);

    for d_theta in &[
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(-0.7, 0.3, 0.9),
    ] {
        let lhs = h_r.dot(&d_theta.cross(&a));
        let rhs = j.dot(d_theta);
        assert!(
            (lhs - rhs).abs() < 1e-14,
            "Doppler linearized check failed: lhs={}, rhs={}, δθ={:?}",
            lhs,
            rhs,
            d_theta
        );
    }
}

#[test]
fn test_doppler_attitude_jacobian_zero_omega() {
    let r_b_e = nalgebra::Matrix3::identity();
    let omega_b = Vector3::zeros();
    let lever_body = Vector3::new(0.0, 0.0, 1.5);
    let h_r = Vector3::new(0.3, -0.5, 0.8);
    let j = doppler_attitude_jacobian(&r_b_e, &omega_b, &lever_body, &h_r);
    assert!(j.norm() < 1e-15);
}

#[test]
fn test_range_jacobian_antisymmetry() {
    let lever_ecef = Vector3::new(0.5, -1.0, 1.5);
    let h_r = Vector3::new(0.3, 0.7, -0.2);
    let j1 = range_attitude_jacobian(&lever_ecef, &h_r);
    let j2 = range_attitude_jacobian(&h_r, &lever_ecef);
    assert!((j1 + j2).norm() < 1e-15, "Should be antisymmetric");
}
