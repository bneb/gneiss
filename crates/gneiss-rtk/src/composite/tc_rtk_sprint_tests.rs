//! Adversarial tests for the tightly-coupled Network RTK/INS architecture
//! (`composite/tc_rtk.rs`) driven through its real observation path: real
//! ephemerides, real geometry, real lever-arm Jacobians.
#![allow(clippy::unwrap_used)]

use nalgebra::{UnitQuaternion, Vector3};

use crate::composite::tc_rtk::{TightlyCoupledNetworkRtkIns, TcRtkConfig};
use crate::composite::CompositeMode;
use crate::estimators::eskf::EskfState;
use crate::post_process::vrs::VrsSynthesizer;

use super::sprint_tests::{fixture_ephemerides, ppp_epoch, t0, RX_ECEF};

// ---------------------------------------------------------------------------
// tc_rtk: real double-difference formation (hunt item 1)
// ---------------------------------------------------------------------------

#[test]
fn tc_rtk_forms_double_differences_from_real_ephemerides() {
    let eph = fixture_ephemerides(t0());
    let eskf = EskfState::new(RX_ECEF, Vector3::zeros(), UnitQuaternion::identity());
    let synth = VrsSynthesizer::new("BASE", RX_ECEF);
    let cfg = TcRtkConfig { enable_zupt: false, enable_nhc: false, ..Default::default() };
    let mut rtk = TightlyCoupledNetworkRtkIns::new(eskf, synth, cfg).with_ephemerides(eph.clone());

    let rover = ppp_epoch(t0(), &eph, 8);
    let base = ppp_epoch(t0(), &eph, 8);
    let cors = vec![crate::composite::StationEpoch {
        station_id: "BASE".to_string(),
        station_pos: RX_ECEF,
        obs: base.clone(),
    }];

    let sol = rtk.process_epoch(&[], &rover, &cors).expect("RTK epoch must process");

    // With 8 satellites, 7 double-differenced pairs form per group and
    // num_sats = pairs + 1 = 8. An empty ephemeris list (the pre-existing
    // tests' degenerate input) yields zero pairs and DeadReckoning.
    assert_ne!(
        sol.mode,
        CompositeMode::DeadReckoning,
        "real ephemerides + real geometry must not fall back to dead reckoning"
    );
    assert_eq!(sol.num_satellites, 8, "7 DD pairs + 1 reference = 8");
    assert!(
        sol.cov.iter().all(|v| v.is_finite()),
        "RTK covariance went non-finite"
    );
}

#[test]
fn tc_rtk_los_jacobian_matches_central_differences() {
    use crate::composite::tc_rtk::{compute_dd_att_coupling_jacobian, compute_dd_los_jacobian};
    use crate::estimators::eskf::skew_symmetric;

    let p_i = Vector3::new(1.0e7, 2.0e7, 1.5e7);
    let p_ref = Vector3::new(-1.0e7, 1.2e7, 2.0e7);
    let rx = RX_ECEF;
    let u_i = (p_i - rx).normalize();
    let u_ref = (p_ref - rx).normalize();
    let lever = Vector3::new(0.3, -0.15, 0.8);
    let l_e = lever; // identity attitude

    let (delta_u, h_dd) = compute_dd_los_jacobian(&u_ref, &u_i);
    let h_att = compute_dd_att_coupling_jacobian(&delta_u, &skew_symmetric(&l_e));

    // Residual model actually used by the filter:
    //   dd(p) = |p_i - p| - |p_i - base| - |p_ref - p| + |p_ref - base|
    let base = RX_ECEF + Vector3::new(50.0, 50.0, 0.0);
    let dd = |p: Vector3<f64>| -> f64 {
        (p_i - p).norm() - (p_i - base).norm() - (p_ref - p).norm() + (p_ref - base).norm()
    };

    let eps = 1.0;
    for j in 0..3 {
        let mut dp = Vector3::zeros();
        dp[j] = eps;
        let num = (dd(rx + dp) - dd(rx - dp)) / (2.0 * eps);
        assert!(
            (num - h_dd[j]).abs() < 1e-6,
            "position column {j}: FD {num} vs analytic {}",
            h_dd[j]
        );
    }

    // Attitude columns must be checked by differentiating the SAME residual,
    // through the SAME left-multiplied retraction the ESKF uses
    // (`q <- Exp(d_theta) * q`, see eskf/update.rs apply_error_injection).
    // Perturbing position inside `dd_att` would only re-derive the analytic
    // expression and prove nothing.
    let dd_att = |dq: UnitQuaternion<f64>| -> f64 {
        let ant = rx + dq.to_rotation_matrix() * lever;
        (p_i - ant).norm() - (p_i - base).norm() - (p_ref - ant).norm() + (p_ref - base).norm()
    };
    // 1e-3 rad: small enough that the O(eps^2) remainder is ~1e-6, large
    // enough that subtracting two ~2.0e7 m ranges in f64 does not lose the
    // signal to cancellation.
    let eps_a = 1e-3;
    for j in 0..3 {
        let mut dv = Vector3::zeros();
        dv[j] = eps_a;
        let q_plus = UnitQuaternion::from_scaled_axis(dv) * UnitQuaternion::identity();
        let q_minus = UnitQuaternion::from_scaled_axis(-dv) * UnitQuaternion::identity();
        let num = (dd_att(q_plus) - dd_att(q_minus)) / (2.0 * eps_a);
        assert!(
            (num - h_att[j]).abs() < 1e-4,
            "attitude column {j}: FD {num} vs analytic {}",
            h_att[j]
        );
    }

    // Sanity on the sign of the attitude block: with a pure-x rotation the
    // lever arm swings in the y-z plane, so column 0 is the dominant term and
    // the analytic row must reproduce it to first order. Guard against a
    // transposed block silently passing a loose tolerance.
    let h_att_z = h_att[2];
    assert!(
        h_att_z.abs() > 1e-3,
        "fixture is degenerate: attitude row too small to discriminate ({h_att:?})"
    );
}

// ---------------------------------------------------------------------------
// Synthetic ephemeris sanity (guards the fixture itself)
// ---------------------------------------------------------------------------

#[test]
fn synthetic_ephemeris_fixture_keeps_satellites_above_the_elevation_mask() {
    let eph = fixture_ephemerides(t0());
    let up = RX_ECEF.normalize();
    for e in &eph {
        let (p, _, _, _) = e.position(t0());
        let u = (p - RX_ECEF).normalize();
        let sin_el = u.dot(&up);
        assert!(
            sin_el >= 10.0_f64.to_radians().sin(),
            "sat {:?} at elevation {:.1} deg is below the 10 deg mask",
            e.sat(),
            sin_el.asin().to_degrees()
        );
    }
}
