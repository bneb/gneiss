use crate::filter::RtkState;
use gneiss_core::coords::{ecef_to_llh, ecef_to_ned_matrix, Coordinate};
use statrs::distribution::ContinuousCDF;

const DEFAULT_P_FA: f64 = 1e-5;
const DEFAULT_P_MD: f64 = 1e-3;
const DEFAULT_P_HMI: f64 = 1e-7;

/// Advanced RAIM (ARAIM) / Solution Separation integrity monitor.
pub struct AraimMonitor {
    pub p_fa: f64,  // Probability of False Alert
    pub p_md: f64,  // Probability of Missed Detection
    pub p_hmi: f64, // Target Integrity Risk
}

impl Default for AraimMonitor {
    fn default() -> Self {
        Self {
            p_fa: DEFAULT_P_FA,
            p_md: DEFAULT_P_MD,
            p_hmi: DEFAULT_P_HMI,
        }
    }
}

/// The result of an ARAIM evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct IntegrityStatus {
    pub hpl: f64,    // Horizontal Protection Level (meters)
    pub vpl: f64,    // Vertical Protection Level (meters)
    pub alert: bool, // True if a fault was detected
}

impl AraimMonitor {
    pub fn new() -> Self {
        Self::default()
    }

    /// Evaluates the protection levels and detects anomalies using Solution Separation.
    pub fn evaluate_solution_separation(
        &self,
        full: &RtkState,
        subsets: &[RtkState],
    ) -> IntegrityStatus {
        if subsets.is_empty() {
            return IntegrityStatus {
                hpl: f64::INFINITY,
                vpl: f64::INFINITY,
                alert: false,
            };
        }

        let (k_fa, k_md) = compute_thresholds(self.p_fa, self.p_md, subsets.len());

        let mut status = IntegrityStatus {
            hpl: 0.0,
            vpl: 0.0,
            alert: false,
        };

        for subset in subsets {
            let sub_status = evaluate_single_subset(full, subset, k_fa, k_md);
            status.hpl = status.hpl.max(sub_status.hpl);
            status.vpl = status.vpl.max(sub_status.vpl);
            status.alert |= sub_status.alert;
        }

        status
    }
}

fn compute_thresholds(p_fa: f64, p_md: f64, num_subsets: usize) -> (f64, f64) {
    let dist = statrs::distribution::Normal::new(0.0, 1.0).unwrap();
    let k_fa = dist.inverse_cdf(1.0 - p_fa / (2.0 * num_subsets as f64));
    let k_md = dist.inverse_cdf(1.0 - p_md / 2.0);
    (k_fa, k_md)
}

fn evaluate_single_subset(
    full: &RtkState,
    subset: &RtkState,
    k_fa: f64,
    k_md: f64,
) -> IntegrityStatus {
    let dx = full.position.vector.x - subset.position.vector.x;
    let dy = full.position.vector.y - subset.position.vector.y;
    let dz = full.position.vector.z - subset.position.vector.z;

    let (dn, de, du) = ecef_to_enu_diff(full.position.clone(), dx, dy, dz);
    let d_horiz = (dn.powi(2) + de.powi(2)).sqrt();
    let d_vert = du.abs();

    let (sigma_diff_h, sigma_diff_v) = project_covariance_diff_to_ned(full, subset);

    let t_h = k_fa * sigma_diff_h;
    let t_v = k_fa * sigma_diff_v;
    let alert = d_horiz > t_h || d_vert > t_v;

    let p_sub = subset.covariance.view((0, 0), (3, 3));
    let var_sub_h = p_sub[(0, 0)] + p_sub[(1, 1)]; // Simplification
    let var_sub_v = p_sub[(2, 2)];

    let hpl = d_horiz + k_md * var_sub_h.sqrt();
    let vpl = d_vert + k_md * var_sub_v.sqrt();

    IntegrityStatus { hpl, vpl, alert }
}

fn project_covariance_diff_to_ned(full: &RtkState, subset: &RtkState) -> (f64, f64) {
    let p_full = full.covariance.view((0, 0), (3, 3));
    let p_sub = subset.covariance.view((0, 0), (3, 3));
    let dp = &p_sub - &p_full;

    let llh = ecef_to_llh(full.position.vector);
    let rot = ecef_to_ned_matrix(llh);
    let dp_ned = rot * dp * rot.transpose();

    let var_diff_h = dp_ned[(0, 0)].max(0.0) + dp_ned[(1, 1)].max(0.0);
    let var_diff_v = dp_ned[(2, 2)].max(0.0);

    (var_diff_h.sqrt(), var_diff_v.sqrt())
}

fn ecef_to_enu_diff(coord: Coordinate, dx: f64, dy: f64, dz: f64) -> (f64, f64, f64) {
    let llh = ecef_to_llh(coord.vector);
    let lat = llh.x;
    let lon = llh.y;

    let slon = lon.sin();
    let clon = lon.cos();
    let slat = lat.sin();
    let clat = lat.cos();

    let de = -slon * dx + clon * dy;
    let dn = -slat * clon * dx - slat * slon * dy + clat * dz;
    let du = clat * clon * dx + clat * slon * dy + slat * dz;

    (dn, de, du)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::RtkState;
    use gneiss_core::coords::Coordinate;
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, Vector3};

    fn mock_state(x: f64, y: f64, z: f64, cov_scale: f64) -> RtkState {
        use gneiss_core::coords::{Datum, Frame};
        let epoch = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(x, y, z), Datum::WGS84, Frame::ECEF, epoch);
        let mut state = RtkState::new(epoch, pos, 30.0);
        let mut cov = DMatrix::zeros(30, 30);
        for i in 0..3 {
            cov[(i, i)] = cov_scale;
        }
        state.covariance = cov;
        state
    }

    #[test]
    fn test_araim_no_fault() {
        let monitor = AraimMonitor::new();
        let full = mock_state(1000.0, 1000.0, 1000.0, 0.1);
        let sub1 = mock_state(1000.01, 1000.01, 1000.01, 0.2); // slight noise, larger cov
        let sub2 = mock_state(999.99, 999.99, 999.99, 0.2);

        let status = monitor.evaluate_solution_separation(&full, &[sub1, sub2]);
        assert!(!status.alert); // Should not alert for minor differences
        assert!(status.hpl > 0.0);
        assert!(status.vpl > 0.0);
    }

    #[test]
    fn test_araim_fault_detected() {
        let monitor = AraimMonitor::new();
        let full = mock_state(1000.0, 1000.0, 1000.0, 0.1);
        let sub1 = mock_state(1000.01, 1000.01, 1000.01, 0.2);
        let sub2 = mock_state(1005.0, 1000.0, 1000.0, 0.2); // 5 meter fault

        let status = monitor.evaluate_solution_separation(&full, &[sub1, sub2]);
        assert!(status.alert); // Should alert due to the 5 meter discrepancy
    }

    #[test]
    fn test_compute_thresholds() {
        let (k_fa, k_md) = compute_thresholds(1e-5, 1e-3, 5);
        assert!(k_fa > 4.0); // k_fa for p_fa=1e-5 / 10 should be ~4.7
        assert!(k_md > 3.0); // k_md for p_md=1e-3 / 2 should be ~3.29
    }
}
