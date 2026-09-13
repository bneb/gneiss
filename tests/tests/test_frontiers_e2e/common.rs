//! Shared test utilities, mathematical reference oracles, and fixtures.

use nalgebra::{Matrix3, UnitQuaternion, Vector2, Vector3};

/// WGS-84 Earth semi-major axis (meters).
pub const WGS84_A: f64 = 6378137.0;
/// Earth rotation rate (rad/s).
pub const OMEGA_EARTH: f64 = 7.2921151467e-5;
/// Speed of light (m/s).
pub const SPEED_OF_LIGHT: f64 = 299792458.0;

/// Earth rotation vector in ECEF frame.
pub fn omega_ie_ecef() -> Vector3<f64> {
    Vector3::new(0.0, 0.0, OMEGA_EARTH)
}

/// Compute 3x3 skew-symmetric matrix [v×].
pub fn skew_symmetric(v: &Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(
        0.0, -v.z, v.y,
        v.z, 0.0, -v.x,
        -v.y, v.x, 0.0,
    )
}

/// Compute nominal rotation matrix from roll, pitch, yaw (radians).
pub fn r_b_e_from_rpy(roll: f64, pitch: f64, yaw: f64) -> UnitQuaternion<f64> {
    UnitQuaternion::from_euler_angles(roll, pitch, yaw)
}

/// Calculate barycentric coordinates (l1, l2, l3) of point p in triangle (a, b, c).
pub fn compute_barycentric(
    p: &Vector2<f64>,
    a: &Vector2<f64>,
    b: &Vector2<f64>,
    c: &Vector2<f64>,
) -> Option<Vector3<f64>> {
    let det = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
    if det.abs() < 1e-12 {
        return None;
    }
    let l1 = ((b.y - c.y) * (p.x - c.x) + (c.x - b.x) * (p.y - c.y)) / det;
    let l2 = ((c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y)) / det;
    let l3 = 1.0 - l1 - l2;
    Some(Vector3::new(l1, l2, l3))
}

/// Circumcircle test: returns true if point d lies inside circumcircle of triangle (a, b, c).
/// Points a, b, c must be oriented counter-clockwise.
pub fn in_circumcircle(
    a: &Vector2<f64>,
    b: &Vector2<f64>,
    c: &Vector2<f64>,
    d: &Vector2<f64>,
) -> bool {
    let adx = a.x - d.x;
    let ady = a.y - d.y;
    let bdx = b.x - d.x;
    let bdy = b.y - d.y;
    let cdx = c.x - d.x;
    let cdy = c.y - d.y;
    let abdet = adx * (bdy * (cdx * cdx + cdy * cdy) - cdy * (bdx * bdx + bdy * bdy));
    let bcdet = bdx * (cdy * (adx * adx + ady * ady) - ady * (cdx * cdx + cdy * cdy));
    let cadet = cdx * (ady * (bdx * bdx + bdy * bdy) - bdy * (adx * adx + ady * ady));
    (abdet + bcdet + cadet) > 1e-12
}

/// Compute Melbourne-Wübbena combination in cycles given dual-freq carrier & code.
#[allow(dead_code)]
pub fn melbourne_wubbena_cycles(
    cp1_m: f64,
    cp2_m: f64,
    pr1_m: f64,
    pr2_m: f64,
    f1: f64,
    f2: f64,
) -> f64 {
    let lam_wl = SPEED_OF_LIGHT / (f1 - f2);
    let phi_wl_m = (f1 * cp1_m - f2 * cp2_m) / (f1 - f2);
    let rho_nl_m = (f1 * pr1_m + f2 * pr2_m) / (f1 + f2);
    (phi_wl_m - rho_nl_m) / lam_wl
}

/// Form double difference between receivers (r1, r2) and satellites (s1, s2).
pub fn form_double_difference(
    obs_r1_s1: f64,
    obs_r1_s2: f64,
    obs_r2_s1: f64,
    obs_r2_s2: f64,
) -> f64 {
    (obs_r2_s2 - obs_r2_s1) - (obs_r1_s2 - obs_r1_s1)
}

/// Statistical evaluation metrics for 3D error trajectory.
#[allow(dead_code)]
pub struct TrajectoryMetrics {
    pub p50: f64,
    pub p68: f64,
    pub p95: f64,
    pub rms: f64,
}

/// Compute percentile from sorted slice.
fn compute_percentile(sorted: &[f64], pct: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() - 1) as f64 * pct).round() as usize;
    sorted[idx]
}

/// Compute trajectory error metrics from a slice of horizontal error magnitudes.
pub fn compute_trajectory_metrics(mut errors: Vec<f64>) -> TrajectoryMetrics {
    if errors.is_empty() {
        return TrajectoryMetrics { p50: 0.0, p68: 0.0, p95: 0.0, rms: 0.0 };
    }
    errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let p50 = compute_percentile(&errors, 0.50);
    let p68 = compute_percentile(&errors, 0.68);
    let p95 = compute_percentile(&errors, 0.95);
    let sum_sq: f64 = errors.iter().map(|e| e * e).sum();
    let rms = (sum_sq / errors.len() as f64).sqrt();
    TrajectoryMetrics { p50, p68, p95, rms }
}

/// Standard 5-station regional CORS network coordinates in ENU (km).
pub fn get_standard_cors_network_2d() -> Vec<Vector2<f64>> {
    vec![
        Vector2::new(0.0, 0.0),       // Master base
        Vector2::new(25.0, 10.0),     // Station East
        Vector2::new(12.0, 30.0),     // Station North
        Vector2::new(-20.0, 15.0),    // Station West
        Vector2::new(-5.0, -25.0),    // Station South
    ]
}
