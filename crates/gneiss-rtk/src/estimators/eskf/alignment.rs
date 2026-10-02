use nalgebra::{Rotation3, UnitQuaternion, Vector3};

use crate::estimators::eskf::types::EskfState;
use crate::swfg::imu_preintegration::ImuSample;

/// Estimate stationary gyroscope bias by averaging static samples.
pub fn compute_gyro_bias(imu_samples: &[ImuSample], max_samples: usize) -> Vector3<f64> {
    let n = imu_samples.len().clamp(1, max_samples);
    let mut sum_g = Vector3::zeros();
    for s in &imu_samples[..n] {
        sum_g += s.gyro;
    }
    sum_g / n as f64
}

/// Compute leveling roll and pitch angles from specific force under gravity.
pub fn compute_leveling_angles(imu_samples: &[ImuSample], max_samples: usize) -> (f64, f64) {
    let n = imu_samples.len().clamp(1, max_samples);
    let mut sum_a = Vector3::zeros();
    for s in &imu_samples[..n] {
        sum_a += s.accel;
    }
    let mean_a = sum_a / n as f64;
    let pitch = (mean_a.x / 9.7803).clamp(-0.5, 0.5);
    let roll = (-mean_a.y / 9.7803).clamp(-0.5, 0.5);
    (roll, pitch)
}

/// Compute initial body-to-ECEF attitude quaternion from leveling angles, position, and heading.
pub fn compute_initial_attitude(
    imu_samples: &[ImuSample],
    init_pos: Vector3<f64>,
    heading_rad: f64,
    max_samples: usize,
) -> UnitQuaternion<f64> {
    let (roll, pitch) = compute_leveling_angles(imu_samples, max_samples);
    let llh = gneiss_core::coords::ecef_to_llh(init_pos);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(llh).transpose();
    let r_body = Rotation3::from_euler_angles(roll, pitch, heading_rad);
    let rot = Rotation3::from_matrix_unchecked(ned_to_ecef * r_body.matrix());
    UnitQuaternion::from_rotation_matrix(&rot)
}

/// Initialize an ESKF state with specified position, velocity, attitude, and gyro bias.
pub fn init_eskf_filter(
    init_pos: Vector3<f64>,
    init_vel: Vector3<f64>,
    init_att: UnitQuaternion<f64>,
    gyro_bias: Vector3<f64>,
) -> EskfState {
    let mut state = EskfState::new(init_pos, init_vel, init_att);
    state.gyro_bias = gyro_bias;
    state
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Matrix3;

    fn make_sample(accel: Vector3<f64>, gyro: Vector3<f64>) -> ImuSample {
        ImuSample {
            accel,
            gyro,
            time_us: 0,
        }
    }

    #[test]
    fn test_compute_gyro_bias_averaging() {
        let samples = vec![
            make_sample(Vector3::zeros(), Vector3::new(0.01, -0.02, 0.03)),
            make_sample(Vector3::zeros(), Vector3::new(0.03, -0.04, 0.05)),
        ];
        let bias = compute_gyro_bias(&samples, 10);
        assert!((bias.x - 0.02).abs() < 1e-12);
        assert!((bias.y - (-0.03)).abs() < 1e-12);
        assert!((bias.z - 0.04).abs() < 1e-12);
    }

    #[test]
    fn test_compute_leveling_angles_level() {
        // Specific force when upright and stationary: +9.7803 upwards in NED or z-axis
        let samples = vec![make_sample(Vector3::new(0.0, 0.0, 9.7803), Vector3::zeros())];
        let (roll, pitch) = compute_leveling_angles(&samples, 10);
        assert!(roll.abs() < 1e-6);
        assert!(pitch.abs() < 1e-6);
    }

    #[test]
    fn test_init_eskf_filter() {
        let p = Vector3::new(100.0, 200.0, 300.0);
        let v = Vector3::new(1.0, 2.0, 3.0);
        let q = UnitQuaternion::identity();
        let bg = Vector3::new(0.001, -0.002, 0.003);
        let state = init_eskf_filter(p, v, q, bg);
        assert_eq!(state.pos_ecef, p);
        assert_eq!(state.vel_ecef, v);
        assert_eq!(state.gyro_bias, bg);
    }

    /// The leveling constant the module divides by (equatorial normal gravity).
    const G_REF: f64 = 9.7803;

    /// One specific-force sample for a platform tilted by `pitch` and `roll`
    /// (radians), held stationary in NED.
    ///
    /// At rest the accelerometer measures specific force f = a - g, i.e. minus
    /// the gravity vector. Gravity in NED is (0, 0, +g), so a level platform
    /// reads (0, 0, -g). Rotating into a body frame tilted by roll `phi` about
    /// NED-East and pitch `theta` about NED-North,
    ///     f_body = ( g*sin(theta),  -g*cos(theta)*sin(phi),  -g*cos(theta)*cos(phi) )
    /// (the vertical is negated because specific force opposes gravity).
    fn tilted_stationary(accel: Vector3<f64>) -> ImuSample {
        ImuSample { accel, gyro: Vector3::zeros(), time_us: 0 }
    }

    #[test]
    fn leveling_pitch_is_positive_nose_up_and_recovers_tilt_to_first_order() {
        // theta = +3 deg nose-up, phi = 0, so f = (g sin th, 0, -g cos th).
        // The module inverts this LINEARLY: pitch = f_x / g = sin(theta), which
        // under-reads the true tilt by the classic small-angle remainder
        // theta - sin(theta) = theta^3/6 = (0.0523599)^3 / 6 = 2.3925e-5 rad.
        let th = 3.0_f64.to_radians();
        let g = G_REF;
        let f = Vector3::new(g * th.sin(), 0.0, -g * th.cos());
        let (roll, pitch) = compute_leveling_angles(&[tilted_stationary(f)], 10);
        assert!(pitch > 0.0, "nose-up must give POSITIVE pitch, got {pitch}");
        let remainder = th.powi(3) / 6.0;
        assert!(
            (pitch - th).abs() <= remainder * 1.001,
            "pitch {pitch} differs from theta {th} by more than theta^3/6 = {remainder}"
        );
        // An exact (non-linearised) inversion would return exactly theta; the
        // residual must be the under-read, i.e. pitch < theta.
        assert!(pitch < th, "expected the first-order shortfall, got {pitch} >= {th}");
        // The exact shortfall is theta - sin(theta) = theta^3/6 - theta^5/120 + ...,
        // so it must sit within one neglected theta^5/120 of theta^3/6:
        // theta^5/120 = (0.0523599)^5 / 120 = 3.279e-9 rad.
        let next = th.powi(5) / 120.0;
        assert!(((th - pitch) - remainder).abs() < next * 1.01,
            "shortfall {} differs from theta^3/6 = {remainder} by more than {next}", th - pitch);
        assert!(roll.abs() < 1e-15, "roll = {roll}");
    }

    #[test]
    fn leveling_roll_is_positive_right_wing_down_with_a_cos_pitch_attenuation() {
        // Tilted by pitch theta = 1 deg AND roll phi = 2 deg, the exact
        // specific force is f = (g sin th, -g cos th sin ph, -g cos th cos ph).
        // The module's linear inversion therefore returns
        //   pitch =  sin(theta)          and
        //   roll  =  cos(theta) * sin(phi) = sin(2deg) * cos(1deg).
        // The cos(theta) factor is the cross-coupling term: dropping it (or
        // flipping the minus sign on f_y) fails this assertion.
        let th = 1.0_f64.to_radians();
        let ph = 2.0_f64.to_radians();
        let g = G_REF;
        let f = Vector3::new(
            g * th.sin(),
            -g * th.cos() * ph.sin(),
            -g * th.cos() * ph.cos(),
        );
        let (roll, pitch) = compute_leveling_angles(&[tilted_stationary(f)], 10);
        assert!(roll > 0.0, "right-wing-down must give POSITIVE roll, got {roll}");
        assert!((pitch - th.sin()).abs() < 1e-15, "pitch = {pitch}, want {}", th.sin());
        assert!((roll - th.cos() * ph.sin()).abs() < 1e-15,
            "roll = {roll}, want {}", th.cos() * ph.sin());
        // Without the cos(theta) attenuation the answer would be sin(phi); the
        // two differ by sin(phi) * (1 - cos theta) = 0.0348995 * 1.523048e-4
        // = 5.3153e-6 rad here, so a tolerance tighter than that discriminates.
        assert!((roll - ph.sin()).abs() > 1e-6, "cos(theta) attenuation is missing");
    }

    #[test]
    fn leveling_clamps_both_angles_to_half_a_radian() {
        // |roll|, |pitch| <= 0.5 rad (~28.6 deg) is the documented clamp; an
        // accelerometer reading far past 90 deg of tilt is a fault, not attitude.
        let big = Vector3::new(1.0e6, -1.0e6, 0.0);
        let (roll, pitch) = compute_leveling_angles(&[tilted_stationary(big)], 10);
        assert!((roll - 0.5).abs() < 1e-15 && (pitch - 0.5).abs() < 1e-15, "roll {roll} pitch {pitch}");
    }

    #[test]
    fn leveling_averages_the_first_max_samples_only() {
        // Mean of the first two specific forces: ((1+3)/2, (2+4)/2, 0) = (2, 3, 0).
        let s = [
            tilted_stationary(Vector3::new(1.0, 2.0, 0.0)),
            tilted_stationary(Vector3::new(3.0, 4.0, 0.0)),
            tilted_stationary(Vector3::new(100.0, 100.0, 0.0)),
        ];
        let (roll, pitch) = compute_leveling_angles(&s, 2);
        assert!((pitch - 2.0 / G_REF).abs() < 1e-15, "pitch = {pitch}");
        assert!((roll + 3.0 / G_REF).abs() < 1e-15, "roll = {roll}");
    }

    #[test]
    fn gyro_bias_averages_the_first_max_samples_only() {
        // First two gyros only: ((1+3)/2, (2+4)/2, (3+5)/2) = (2, 3, 4).
        let s = [
            ImuSample { accel: Vector3::zeros(), gyro: Vector3::new(1.0, 2.0, 3.0), time_us: 0 },
            ImuSample { accel: Vector3::zeros(), gyro: Vector3::new(3.0, 4.0, 5.0), time_us: 0 },
            ImuSample { accel: Vector3::zeros(), gyro: Vector3::new(50.0, 50.0, 50.0), time_us: 0 },
        ];
        let bias = compute_gyro_bias(&s, 2);
        assert!((bias - Vector3::new(2.0, 3.0, 4.0)).norm() < 1e-15, "got {bias:?}");
    }

    #[test]
    fn initial_attitude_at_equator_prime_meridian_maps_body_axes_exactly() {
        // At (lat, lon) = (0, 0) the NED -> ECEF rotation is the exact
        // permutation [[0,0,-1],[0,1,0],[1,0,0]] (North=+x_ec, East=+y_ec,
        // Down=-z_ec). With a level IMU and heading 0 the body axes must be
        //   forward -> ECEF +z (Up), right -> ECEF +y (East), down -> ECEF -x.
        let level = ImuSample { accel: Vector3::new(0.0, 0.0, G_REF), gyro: Vector3::zeros(), time_us: 0 };
        let q = compute_initial_attitude(&[level], Vector3::new(6_378_137.0, 0.0, 0.0), 0.0, 10);
        let rot = q.to_rotation_matrix().into_inner();
        let fwd = rot * Vector3::new(1.0, 0.0, 0.0);
        let right = rot * Vector3::new(0.0, 1.0, 0.0);
        let down = rot * Vector3::new(0.0, 0.0, 1.0);
        assert!((fwd - Vector3::new(0.0, 0.0, 1.0)).norm() < 1e-12, "fwd = {fwd:?}");
        assert!((right - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12, "right = {right:?}");
        assert!((down - Vector3::new(-1.0, 0.0, 0.0)).norm() < 1e-12, "down = {down:?}");
        // The basis must be a proper rotation, not a scaled or mirrored one.
        assert!((rot * rot.transpose() - Matrix3::identity()).norm() < 1e-12);
        assert!((rot.determinant() - 1.0).abs() < 1e-12);
        assert!((q.norm() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn initial_attitude_heading_90_degrees_points_forward_east() {
        // Heading is measured clockwise from North, so 90 deg must face East
        // (ECEF +y at the equator / prime meridian). from_euler_angles is
        // intrinsic (R = Rz(yaw) Ry(pitch) Rx(roll)), so body forward
        // (1,0,0) maps to (cos yaw, sin yaw, 0) in NED = East for yaw = 90 deg.
        let level = ImuSample { accel: Vector3::new(0.0, 0.0, G_REF), gyro: Vector3::zeros(), time_us: 0 };
        let yaw = 90.0_f64.to_radians();
        let rot = compute_initial_attitude(&[level], Vector3::new(6_378_137.0, 0.0, 0.0), yaw, 10)
            .to_rotation_matrix();
        let fwd = rot * Vector3::new(1.0, 0.0, 0.0);
        assert!((fwd - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12, "fwd = {fwd:?}");
        // Facing East, the right wing points South. At this site (ECEF +x axis)
        // North = ECEF +z, East = ECEF +y, so South = ECEF -z.
        let right = rot * Vector3::new(0.0, 1.0, 0.0);
        assert!((right - Vector3::new(0.0, 0.0, -1.0)).norm() < 1e-12, "right = {right:?}");
    }

    #[test]
    fn initial_attitude_composes_the_measured_levelling_with_the_heading() {
        // A 3 deg nose-up tilt with heading 0 (North). At this site the NED ->
        // ECEF permutation maps (n, e, d) -> (-d, e, n), and a nose-up body
        // forward axis in NED is Ry(theta) (1,0,0) = (cos theta, 0, -sin theta),
        // so the exact ECEF forward axis is (sin theta, 0, cos theta).
        //
        // The levelling block is a first-order inversion, so the attitude it
        // produces is tilted by theta_hat = sin(theta), leaving a residual of at
        // most theta^3/6 = 2.3925e-5 rad. Pin both the bound and the SIGN of the
        // residual: the recovered axis must under-tilt (ECEF x below sin theta),
        // because theta_hat < theta.
        let th = 3.0_f64.to_radians();
        let f = Vector3::new(G_REF * th.sin(), 0.0, -G_REF * th.cos());
        let q = compute_initial_attitude(
            &[tilted_stationary(f)],
            Vector3::new(6_378_137.0, 0.0, 0.0),
            0.0,
            10,
        );
        let fwd = q.to_rotation_matrix() * Vector3::new(1.0, 0.0, 0.0);
        let exact = Vector3::new(th.sin(), 0.0, th.cos());
        let tol = 1.5 * th.powi(3) / 6.0;
        assert!((fwd - exact).norm() < tol, "fwd = {fwd:?}, exact = {exact:?}, tol {tol}");
        assert!(fwd.x < exact.x, "first-order levelling must under-tilt, got {fwd:?}");
        assert!((fwd.norm() - 1.0).abs() < 1e-15);
    }
}
