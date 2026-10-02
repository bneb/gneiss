//! RED phase for the IMU preintegration frame-convention defect.
//!
//! Every existing IMU-residual test routes through
//! `swfg/engine/setup_tests.rs::preint_factor`, which deliberately forces
//! `UnitQuaternion::identity()` — and that helper's own comment (setup_tests.rs:
//! 231-236) says why: *"the residual then compares a body-frame prediction
//! against the world-frame vector `ImuPreintegration::integrate` accumulates"*, and
//! points at a test named `imu_residual_rotates_only_the_prediction` that was
//! never written. Identity attitude makes body and ECEF coincide, so the mixing
//! cancels and the defect is invisible.
//!
//! These tests use a REAL attitude so the defect becomes detectable. This is the
//! test the comment already names. It is expected to FAIL until the frame
//! convention is fixed — that failure is the point.
//!
//! Geometry, derived by hand:
//!   Receiver at the WGS84 equator/prime meridian. Gravity in ECEF is
//!   (0, 0, -9.80665) — along -Up. NED->ECEF at the equator is the exact
//!   permutation north->+Y, east->+Z, down->-X, i.e. ECEF->NED maps X->down.
//!   Hence body->ECEF at that site maps body X onto ECEF -X, a real 90 deg
//!   rotation, not the identity.
//!
//! One second of free fall from rest: dp = -0.5*g*dt^2 = (+4.903325, 0, 0) m and
//! dv = g*dt = (-9.80665, 0, 0) m/s, both vertical and therefore ECEF-aligned.
//! The factor predicts p_j - p_i - v_i*dt - 0.5*g*dt^2 in ECEF, and is fed the
//! same dp/dv, so a consistent-frame implementation gives residual 0. Rotating
//! the prediction into the body frame before subtracting the ECEF measurement
//! does not — that is the defect.

use gneiss_rtk::swfg::imu_preintegration::{ImuPreintegration, ImuPreintegrationFactor};
use gneiss_rtk::swfg::factor::Factor;
use gneiss_rtk::swfg::variables::{VariableId, VariableKind, VariableNode, VariableValues};
use std::collections::BTreeMap;
use nalgebra::{UnitQuaternion, Vector3};

const G: f64 = 9.80665;

/// Body->ECEF rotation at the equator/prime meridian.
/// ECEF->NED there is the permutation X->down(-X), Y->north(+Y), Z->east(+Z);
/// its transpose (NED->ECEF) has columns (north, east, down) expressed in ECEF:
///   north = +Y, east = +Z, down = -X
/// so R = [[0,0,-1],[1,0,0],[0,1,0]], which maps body X onto ECEF -X.
fn equator_attitude() -> UnitQuaternion<f64> {
    let r = nalgebra::Matrix3::new(0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    UnitQuaternion::from_rotation_matrix(&nalgebra::Rotation3::from_matrix_unchecked(r))
}

/// One second of free fall at `att`, returned as (preintegration, residual).
fn free_fall_residual() -> (f64, f64) {
    let dt = 1.0;
    let q = equator_attitude();
    let p_i = Vector3::new(6_378_137.0, 0.0, 0.0);
    let v_i = Vector3::zeros();
    // Free fall: a = -g = (0,0,+G), so v_j = v_i + a*dt and p_j = p_i + 0.5*a*dt^2.
    let p_j = p_i + Vector3::new(0.0, 0.0, 0.5 * G);
    let v_j = v_i + Vector3::new(0.0, 0.0, G);
    let gravity = Vector3::new(0.0, 0.0, -G);

    // The factor predicts  dp = p_j - p_i - v_i*dt - 0.5*g*dt^2  and  dv = v_j - v_i - g*dt.
    // With v_i = 0 and the free-fall endpoints above those evaluate to exactly
    // dp = (0,0,G) and dv = (0,0,2G). Feeding those in, in the SAME frame as the
    // prediction, is what makes a frame-consistent residual exactly zero.
    let mut preint = ImuPreintegration::new();
    preint.dp = Vector3::new(0.0, 0.0, G);
    preint.dv = Vector3::new(0.0, 0.0, 2.0 * G);
    preint.dt = dt;
    let bias = Vector3::zeros();

    let factor = ImuPreintegrationFactor::new(
        preint, gravity,
        p_i, v_i, q,
        p_j, v_j, q,
        bias, bias,
        VariableId::new(0), VariableId::new(1), VariableId::new(2),
        VariableId::new(3), VariableId::new(4),
    );

    // Build the graph values: pose_i, vel_i, pose_j, vel_j, bias.
    let mut nodes: BTreeMap<VariableId, VariableNode> = BTreeMap::new();
    for (k, p) in [(0u64, p_i), (2u64, p_j)] {
        let mut n = VariableNode::new(VariableId::new(k), VariableKind::Pose { epoch: 0 });
        n.value[0] = p.x; n.value[1] = p.y; n.value[2] = p.z;
        nodes.insert(VariableId::new(k), n);
    }
    for (k, v) in [(1u64, v_i), (3u64, v_j)] {
        let mut n = VariableNode::new(VariableId::new(k), VariableKind::Velocity { epoch: 0 });
        n.value[0] = v.x; n.value[1] = v.y; n.value[2] = v.z;
        nodes.insert(VariableId::new(k), n);
    }
    nodes.insert(VariableId::new(4), VariableNode::new(VariableId::new(4), VariableKind::ImuBias));
    let vv = VariableValues::build(&nodes);

    let r = factor.residual(&vv);
    let pos_norm = r.rows(0, 3).norm();
    let vel_norm = r.rows(3, 3).norm();
    (pos_norm, vel_norm)
}

#[test]
fn imu_residual_rotates_only_the_prediction() {
    // The test swfg/engine/setup_tests.rs:236 names but never wrote.
    // Free fall is a known-motion case: the prediction and the preintegrated
    // delta describe the SAME trajectory, so a frame-consistent residual is 0.
    let (pos_norm, vel_norm) = free_fall_residual();
    assert!(
        pos_norm < 1e-9,
        "free-fall position residual must vanish when prediction and preintegration \
         share a frame; got norm {pos_norm}. The residual is rotating only the \
         prediction into the body frame while ImuPreintegration::integrate \
         accumulates in the world frame."
    );
    assert!(
        vel_norm < 1e-9,
        "free-fall velocity residual must vanish when prediction and preintegration \
         share a frame; got norm {vel_norm}."
    );
}

#[test]
fn red_phase_guard_defect_is_reproducible() {
    // Non-vacuity guard: proves the test above is capable of failing, i.e. that
    // the defect is real at a non-identity attitude. DELETE together with the
    // fix, once free_fall residual becomes zero.
    let (pos_norm, vel_norm) = free_fall_residual();
    assert!(
        pos_norm > 1e-6,
        "expected the known body/ECEF mixing to be observable at a real attitude, \
         but the position residual was {pos_norm}"
    );
    assert!(
        vel_norm > 1e-6,
        "expected the known body/ECEF mixing to be observable in velocity too, \
         but the velocity residual was {vel_norm}"
    );
}