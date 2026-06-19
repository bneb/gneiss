//! Factor definitions for the Factor Graph Optimization (FGO) backend.
//!
//! Includes IMU pre-integration, GNSS pseudorange, carrier phase, and doppler factors.

pub mod carrier_phase;
pub mod doppler;
pub mod imu;
pub mod pseudorange;

use nalgebra::{Matrix3, Vector3};

/// Computes the skew-symmetric matrix of a 3D vector.
pub fn skew_symmetric(v: &Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(0.0, -v.z, v.y, v.z, 0.0, -v.x, -v.y, v.x, 0.0)
}
