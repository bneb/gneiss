use nalgebra::Vector3;
use alloc::vec::Vec;

/// Statistical summary of a set of error values.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ErrorStats {
    /// Number of samples
    pub count: usize,
    /// Median (50th percentile)
    pub median: f64,
    /// 95th percentile
    pub p95: f64,
    /// 99th percentile
    pub p99: f64,
    /// Arithmetic mean
    pub mean: f64,
    /// Root Mean Square
    pub rms: f64,
    /// Maximum value
    pub max: f64,
}

/// Computes the horizontal (2D) error between two ECEF positions.
/// Projects the 3D difference into the local East-North plane at the truth position.
/// Returns error in meters.
pub fn horizontal_error(pos_ecef: Vector3<f64>, truth_ecef: Vector3<f64>) -> f64 {
    let truth_llh = crate::coords::ecef_to_llh(truth_ecef);
    let lat = truth_llh.x;
    let lon = truth_llh.y;

    let sin_lat = libm::sin(lat);
    let cos_lat = libm::cos(lat);
    let sin_lon = libm::sin(lon);
    let cos_lon = libm::cos(lon);

    let dx = pos_ecef.x - truth_ecef.x;
    let dy = pos_ecef.y - truth_ecef.y;
    let dz = pos_ecef.z - truth_ecef.z;

    // ECEF to ENU
    let e = -sin_lon * dx + cos_lon * dy;
    let n = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;

    libm::sqrt(e * e + n * n)
}

/// Computes the vertical error between two ECEF positions.
/// Projects the 3D difference onto the local Up axis at the truth position.
/// Returns signed error in meters (positive = above truth).
pub fn vertical_error(pos_ecef: Vector3<f64>, truth_ecef: Vector3<f64>) -> f64 {
    let truth_llh = crate::coords::ecef_to_llh(truth_ecef);
    let lat = truth_llh.x;
    let lon = truth_llh.y;

    let cos_lat = libm::cos(lat);
    let sin_lat = libm::sin(lat);
    let cos_lon = libm::cos(lon);
    let sin_lon = libm::sin(lon);

    let dx = pos_ecef.x - truth_ecef.x;
    let dy = pos_ecef.y - truth_ecef.y;
    let dz = pos_ecef.z - truth_ecef.z;

    // Up component in ENU
    cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz
}

/// Computes the 3D Euclidean distance error between two ECEF positions.
pub fn error_3d(pos_ecef: Vector3<f64>, truth_ecef: Vector3<f64>) -> f64 {
    (pos_ecef - truth_ecef).norm()
}

/// Computes statistical summary of a set of error values.
/// Returns None if the input slice is empty.
pub fn compute_statistics(errors: &[f64]) -> Option<ErrorStats> {
    if errors.is_empty() {
        return None;
    }

    let n = errors.len();
    let mut sorted = Vec::with_capacity(n);
    sorted.extend_from_slice(errors);
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));

    let sum: f64 = errors.iter().sum();
    let mean = sum / (n as f64);

    let sum_sq: f64 = errors.iter().map(|e| e * e).sum();
    let rms = libm::sqrt(sum_sq / (n as f64));

    let median = if (n & 1) == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    };

    let p95_idx = libm::ceil((n as f64) * 0.95) as usize;
    let p95 = sorted[p95_idx.min(n - 1)];

    let p99_idx = libm::ceil((n as f64) * 0.99) as usize;
    let p99 = sorted[p99_idx.min(n - 1)];

    let max = sorted[n - 1];

    Some(ErrorStats {
        count: n,
        median,
        p95,
        p99,
        mean,
        rms,
        max,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_horizontal_error_zero_at_same_point() {
        let pos = Vector3::new(6378137.0, 0.0, 0.0);
        assert!(horizontal_error(pos, pos) < 1e-10);
    }

    #[test]
    fn test_horizontal_error_north_displacement() {
        // At equator/prime meridian, a 1m north displacement is +1m in Z (approximately)
        let truth = Vector3::new(6378137.0, 0.0, 0.0);
        // Move 1m north at equator = move in +Z direction
        let pos = Vector3::new(6378137.0, 0.0, 1.0);
        let h_err = horizontal_error(pos, truth);
        // Should be approximately 1m
        assert!((h_err - 1.0).abs() < 0.01, "Expected ~1m horizontal error, got {}", h_err);
    }

    #[test]
    fn test_vertical_error_upward_displacement() {
        // At equator/prime meridian, radial displacement = vertical
        let truth = Vector3::new(6378137.0, 0.0, 0.0);
        let pos = Vector3::new(6378138.0, 0.0, 0.0); // 1m radially outward
        let v_err = vertical_error(pos, truth);
        assert!((v_err - 1.0).abs() < 0.01, "Expected ~1m vertical error, got {}", v_err);
    }

    #[test]
    fn test_vertical_error_signed() {
        let truth = Vector3::new(6378137.0, 0.0, 0.0);
        let above = Vector3::new(6378138.0, 0.0, 0.0);
        let below = Vector3::new(6378136.0, 0.0, 0.0);
        assert!(vertical_error(above, truth) > 0.0, "Above should be positive");
        assert!(vertical_error(below, truth) < 0.0, "Below should be negative");
    }

    #[test]
    fn test_error_3d() {
        let truth = Vector3::new(6378137.0, 0.0, 0.0);
        let pos = Vector3::new(6378138.0, 1.0, 1.0);
        let err = error_3d(pos, truth);
        let expected = libm::sqrt(1.0 + 1.0 + 1.0);
        assert!((err - expected).abs() < 1e-10);
    }

    #[test]
    fn test_compute_statistics_basic() {
        let errors = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let stats = compute_statistics(&errors).unwrap();
        
        assert_eq!(stats.count, 10);
        assert!((stats.median - 5.5).abs() < 1e-10, "Median should be 5.5, got {}", stats.median);
        assert!((stats.mean - 5.5).abs() < 1e-10, "Mean should be 5.5, got {}", stats.mean);
        assert!(stats.p95 >= 9.0, "P95 should be >= 9.0, got {}", stats.p95);
        assert!(stats.p99 >= 9.0, "P99 should be >= 9.0, got {}", stats.p99);
        assert!((stats.max - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_statistics_single_value() {
        let errors = [42.0];
        let stats = compute_statistics(&errors).unwrap();
        assert_eq!(stats.count, 1);
        assert!((stats.median - 42.0).abs() < 1e-10);
        assert!((stats.mean - 42.0).abs() < 1e-10);
        assert!((stats.rms - 42.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_statistics_empty() {
        assert!(compute_statistics(&[]).is_none());
    }

    #[test]
    fn test_rms_calculation() {
        // RMS of [3, 4] = sqrt((9 + 16) / 2) = sqrt(12.5)
        let errors = [3.0, 4.0];
        let stats = compute_statistics(&errors).unwrap();
        let expected_rms = libm::sqrt(12.5);
        assert!((stats.rms - expected_rms).abs() < 1e-10);
    }
}
