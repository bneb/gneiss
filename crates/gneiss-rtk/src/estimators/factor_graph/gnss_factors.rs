use super::Factor;
use nalgebra::{DMatrix, DVector, Vector3};

macro_rules! compute_clock_delta {
    ($self:ident, $delta:ident) => {{
        let mut dt = $self.nominal_dt + $delta[$self.index_dt];
        match $self.sat_id.constellation {
            gneiss_core::sat::Constellation::Galileo => {
                if let Some(idx) = $self.index_dt_gal {
                    dt += $self.nominal_dt_gal + $delta[idx];
                }
            }
            gneiss_core::sat::Constellation::Beidou => {
                if let Some(idx) = $self.index_dt_bds {
                    dt += $self.nominal_dt_bds + $delta[idx];
                }
            }
            gneiss_core::sat::Constellation::Glonass => {
                if let Some(idx) = $self.index_dt_glo {
                    dt += $self.nominal_dt_glo + $delta[idx];
                }
            }
            _ => {}
        }
        dt
    }};
}

macro_rules! apply_clock_jacobian {
    ($self:ident, $jac:ident) => {
        $jac[(0, $self.index_dt)] = -1.0;
        match $self.sat_id.constellation {
            gneiss_core::sat::Constellation::Galileo => {
                if let Some(idx) = $self.index_dt_gal {
                    $jac[(0, idx)] = -1.0;
                }
            }
            gneiss_core::sat::Constellation::Beidou => {
                if let Some(idx) = $self.index_dt_bds {
                    $jac[(0, idx)] = -1.0;
                }
            }
            gneiss_core::sat::Constellation::Glonass => {
                if let Some(idx) = $self.index_dt_glo {
                    $jac[(0, idx)] = -1.0;
                }
            }
            _ => {}
        }
    };
}

/// Factor for a Pseudorange measurement.
pub struct PseudorangeFactor {
    pub sat_pos: Vector3<f64>,
    pub measured_pr: f64,
    pub variance: f64,
    pub sat_clock_bias: f64,
    pub tropo_dry_delay: f64,
    pub map_wet: f64,
    pub index_x: usize,           // index of x in state
    pub index_y: usize,           // index of y
    pub index_z: usize,           // index of z
    pub index_dt: usize,          // index of receiver clock bias
    pub index_zwd: Option<usize>, // index of zenith wet delay
    pub robust_threshold: f64,
}

impl Factor for PseudorangeFactor {
    fn residual(&self, state: &DVector<f64>) -> DVector<f64> {
        let rx = state[self.index_x];
        let ry = state[self.index_y];
        let rz = state[self.index_z];
        let dt = state[self.index_dt];

        let _zwd = self.index_zwd.map(|idx| state[idx]).unwrap_or(0.0);

        let dx = self.sat_pos.x - rx;
        let dy = self.sat_pos.y - ry;
        let dz = self.sat_pos.z - rz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        let expected_pr = dist + dt - self.sat_clock_bias + self.tropo_dry_delay; // Ignore ZWD for now

        DVector::from_vec(vec![self.measured_pr - expected_pr])
    }

    fn jacobian(&self, state: &DVector<f64>) -> DMatrix<f64> {
        let rx = state[self.index_x];
        let ry = state[self.index_y];
        let rz = state[self.index_z];

        let dx = self.sat_pos.x - rx;
        let dy = self.sat_pos.y - ry;
        let dz = self.sat_pos.z - rz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        // jacobian w.r.t [x, y, z, ..., dt, ...]
        let mut jac = DMatrix::zeros(1, state.len());

        if dist > 1e-6 {
            jac[(0, self.index_x)] = dx / dist; // negative of derivative of expected_pr
            jac[(0, self.index_y)] = dy / dist;
            jac[(0, self.index_z)] = dz / dist;
            jac[(0, self.index_dt)] = -1.0;
            if let Some(idx) = self.index_zwd {
                jac[(0, idx)] = 0.0; // Disable ZWD estimation
            }
        }

        jac
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_element(1, 1, 1.0 / self.variance.max(1e-9))
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(self.robust_threshold)
    }

    fn is_cauchy_rejectable(&self) -> bool {
        true
    }
}

/// Factor for a Carrier Phase measurement.
pub struct CarrierPhaseFactor {
    pub sat_pos: Vector3<f64>,
    pub measured_cp: f64,
    pub variance: f64,
    pub sat_clock_bias: f64,
    pub tropo_dry_delay: f64,
    pub map_wet: f64,
    pub wavelength: f64,
    pub index_x: usize,
    pub index_y: usize,
    pub index_z: usize,
    pub index_dt: usize,
    pub index_zwd: Option<usize>,
    pub index_amb: usize, // index of the ambiguity state (in cycles)
    pub robust_threshold: f64,
}

impl Factor for CarrierPhaseFactor {
    fn residual(&self, state: &DVector<f64>) -> DVector<f64> {
        let rx = state[self.index_x];
        let ry = state[self.index_y];
        let rz = state[self.index_z];
        let dt = state[self.index_dt];
        let amb = state[self.index_amb];
        let is_iono_free = false;
        let _zwd = self.index_zwd.map(|idx| state[idx]).unwrap_or(0.0);
        let _var_cp = if is_iono_free { 0.001 } else { 0.01 }; // 10cm stddev for phase
        let _huber_cp = 3.0; // 3-sigma (30cm) to reject cycle slips.

        let dx = self.sat_pos.x - rx;
        let dy = self.sat_pos.y - ry;
        let dz = self.sat_pos.z - rz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        let expected_cp =
            dist + dt - self.sat_clock_bias + self.tropo_dry_delay + amb * self.wavelength; // Ignore ZWD

        DVector::from_vec(vec![self.measured_cp - expected_cp])
    }

    fn jacobian(&self, state: &DVector<f64>) -> DMatrix<f64> {
        let rx = state[self.index_x];
        let ry = state[self.index_y];
        let rz = state[self.index_z];

        let dx = self.sat_pos.x - rx;
        let dy = self.sat_pos.y - ry;
        let dz = self.sat_pos.z - rz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        let mut jac = DMatrix::zeros(1, state.len());

        if dist > 1e-6 {
            jac[(0, self.index_x)] = dx / dist;
            jac[(0, self.index_y)] = dy / dist;
            jac[(0, self.index_z)] = dz / dist;
            jac[(0, self.index_dt)] = -1.0;
            if let Some(idx) = self.index_zwd {
                jac[(0, idx)] = 0.0;
            }
            jac[(0, self.index_amb)] = -self.wavelength;
        }

        jac
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_element(1, 1, 1.0 / self.variance.max(1e-9))
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(self.robust_threshold)
    }

    fn is_cauchy_rejectable(&self) -> bool {
        true
    }
}

/// Error-State Factor for a Pseudorange measurement.
pub struct ErrorStatePseudorangeFactor {
    pub sat_pos: Vector3<f64>,
    pub measured_pr: f64,
    pub variance: f64,
    pub sat_clock_bias: f64,
    pub tropo_dry_delay: f64,
    pub map_wet: f64,

    pub nominal_rx: f64,
    pub nominal_ry: f64,
    pub nominal_rz: f64,
    pub nominal_dt: f64,
    pub nominal_dt_gal: f64,
    pub nominal_dt_bds: f64,
    pub nominal_dt_glo: f64,
    pub nominal_zwd: f64,

    pub index_x: usize,
    pub index_y: usize,
    pub index_z: usize,
    pub index_dt: usize,
    pub index_dt_gal: Option<usize>,
    pub index_dt_bds: Option<usize>,
    pub index_dt_glo: Option<usize>,
    pub index_zwd: Option<usize>,
    pub sat_id: gneiss_core::sat::SatelliteId,
    pub robust_threshold: f64,
}

impl Factor for ErrorStatePseudorangeFactor {
    fn residual(&self, delta: &DVector<f64>) -> DVector<f64> {
        let (rx, ry, rz) = (
            self.nominal_rx + delta[self.index_x],
            self.nominal_ry + delta[self.index_y],
            self.nominal_rz + delta[self.index_z],
        );
        let dt = compute_clock_delta!(self, delta);
        let zwd = self.nominal_zwd + self.index_zwd.map(|i| delta[i]).unwrap_or(0.0);
        let dist = ((self.sat_pos.x - rx).powi(2)
            + (self.sat_pos.y - ry).powi(2)
            + (self.sat_pos.z - rz).powi(2))
        .sqrt();
        let expected_pr =
            dist + dt - self.sat_clock_bias + self.tropo_dry_delay + zwd * self.map_wet;
        DVector::from_vec(vec![self.measured_pr - expected_pr])
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(self.robust_threshold)
    }
    fn is_cauchy_rejectable(&self) -> bool {
        true
    }

    fn jacobian(&self, delta: &DVector<f64>) -> DMatrix<f64> {
        let (rx, ry, rz) = (
            self.nominal_rx + delta[self.index_x],
            self.nominal_ry + delta[self.index_y],
            self.nominal_rz + delta[self.index_z],
        );
        let (dx, dy, dz) = (
            self.sat_pos.x - rx,
            self.sat_pos.y - ry,
            self.sat_pos.z - rz,
        );
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        let mut jac = DMatrix::zeros(1, delta.len());
        if dist > 1e-6 {
            jac[(0, self.index_x)] = dx / dist;
            jac[(0, self.index_y)] = dy / dist;
            jac[(0, self.index_z)] = dz / dist;
            apply_clock_jacobian!(self, jac);
            if let Some(idx) = self.index_zwd {
                jac[(0, idx)] = -self.map_wet;
            }
        }
        jac
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_element(1, 1, 1.0 / self.variance.max(1e-9))
    }
}

/// Error-State Factor for a Carrier Phase measurement.
pub struct ErrorStateCarrierPhaseFactor {
    pub sat_pos: Vector3<f64>,
    pub measured_cp: f64,
    pub variance: f64,
    pub sat_clock_bias: f64,
    pub tropo_dry_delay: f64,
    pub map_wet: f64,
    pub wavelength: f64,

    pub nominal_rx: f64,
    pub nominal_ry: f64,
    pub nominal_rz: f64,
    pub nominal_dt: f64,
    pub nominal_dt_gal: f64,
    pub nominal_dt_bds: f64,
    pub nominal_dt_glo: f64,
    pub nominal_zwd: f64,
    pub nominal_amb: f64,

    pub index_x: usize,
    pub index_y: usize,
    pub index_z: usize,
    pub index_dt: usize,
    pub index_dt_gal: Option<usize>,
    pub index_dt_bds: Option<usize>,
    pub index_dt_glo: Option<usize>,
    pub index_zwd: Option<usize>,
    pub index_amb: usize,
    pub sat_id: gneiss_core::sat::SatelliteId,
    pub robust_threshold: f64,
}

impl Factor for ErrorStateCarrierPhaseFactor {
    fn residual(&self, delta: &DVector<f64>) -> DVector<f64> {
        let (rx, ry, rz) = (
            self.nominal_rx + delta[self.index_x],
            self.nominal_ry + delta[self.index_y],
            self.nominal_rz + delta[self.index_z],
        );
        let dt = compute_clock_delta!(self, delta);
        let amb = self.nominal_amb + delta[self.index_amb];
        let zwd = self.nominal_zwd + self.index_zwd.map(|idx| delta[idx]).unwrap_or(0.0);
        let dist = ((self.sat_pos.x - rx).powi(2)
            + (self.sat_pos.y - ry).powi(2)
            + (self.sat_pos.z - rz).powi(2))
        .sqrt();
        let expected_cp = dist + dt - self.sat_clock_bias
            + self.tropo_dry_delay
            + zwd * self.map_wet
            + amb * self.wavelength;
        DVector::from_vec(vec![self.measured_cp - expected_cp])
    }

    fn jacobian(&self, delta: &DVector<f64>) -> DMatrix<f64> {
        let (rx, ry, rz) = (
            self.nominal_rx + delta[self.index_x],
            self.nominal_ry + delta[self.index_y],
            self.nominal_rz + delta[self.index_z],
        );
        let (dx, dy, dz) = (
            self.sat_pos.x - rx,
            self.sat_pos.y - ry,
            self.sat_pos.z - rz,
        );
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        let mut jac = DMatrix::zeros(1, delta.len());
        if dist > 1e-6 {
            jac[(0, self.index_x)] = dx / dist;
            jac[(0, self.index_y)] = dy / dist;
            jac[(0, self.index_z)] = dz / dist;
            apply_clock_jacobian!(self, jac);
            if let Some(idx) = self.index_zwd {
                jac[(0, idx)] = -self.map_wet;
            }
            jac[(0, self.index_amb)] = -self.wavelength;
        }
        jac
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_element(1, 1, 1.0 / self.variance.max(1e-9))
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(self.robust_threshold) // threshold to reject cycle slips.
    }

    fn is_cauchy_rejectable(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper to calculate numerical jacobian via central differencing
    fn numerical_jacobian<F>(factor: &F, state: &DVector<f64>) -> DMatrix<f64>
    where
        F: Factor,
    {
        let n = state.len();
        let mut jac = DMatrix::zeros(1, n);
        let eps = 1e-4;

        for i in 0..n {
            let mut state_plus = state.clone();
            state_plus[i] += eps;
            let res_plus = factor.residual(&state_plus);

            let mut state_minus = state.clone();
            state_minus[i] -= eps;
            let res_minus = factor.residual(&state_minus);

            jac[(0, i)] = (res_plus[0] - res_minus[0]) / (2.0 * eps);
        }

        jac
    }

    #[test]
    fn test_error_state_pseudorange_jacobian() {
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(20000000.0, 10000000.0, 5000000.0),
            measured_pr: 22000000.0,
            variance: 1.0,
            sat_clock_bias: 0.0001,
            tropo_dry_delay: 2.3,
            map_wet: 3.1,
            nominal_rx: 1000.0,
            nominal_ry: 2000.0,
            nominal_rz: 3000.0,
            nominal_dt: 0.0002,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };

        let state = DVector::from_vec(vec![5.0, -3.0, 2.0, 0.1, 0.05]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);

        println!("Anal PR: {}\nNum PR: {}", anal_jac, num_jac);
        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "Analytical PR jacobian diverges from numerical"
        );
    }

    #[test]
    fn test_error_state_carrier_phase_jacobian() {
        let factor = ErrorStateCarrierPhaseFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(20000000.0, 10000000.0, 5000000.0),
            measured_cp: 120000000.0,
            variance: 1.0,
            sat_clock_bias: 0.0001,
            tropo_dry_delay: 2.3,
            map_wet: 3.1,
            wavelength: 0.19,
            nominal_rx: 1000.0,
            nominal_ry: 2000.0,
            nominal_rz: 3000.0,
            nominal_dt: 0.0002,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            nominal_amb: 50.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            index_amb: 5,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };

        let state = DVector::from_vec(vec![5.0, -3.0, 2.0, 0.1, 0.05, -2.0]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);

        println!("Anal CP: {}\nNum CP: {}", anal_jac, num_jac);
        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "Analytical CP jacobian diverges from numerical"
        );
    }

    // ==================== PseudorangeFactor tests ====================

    #[test]
    fn test_pseudorange_residual_zero() {
        let factor = PseudorangeFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.0,
            variance: 2.5,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
        let res = factor.residual(&state);
        assert!(
            res[0].abs() < 1e-10,
            "PseudorangeFactor residual should be 0 for correct state. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_pseudorange_residual_with_bias() {
        // satellite at [2,3,6], receiver at [0,0,0] -> dist=7.0
        // dt=0.1, clock_bias=0.5, tropo=0.3 -> expected=7.0+0.1-0.5+0.3=6.9
        // measured_pr=6.9 -> residual=0
        let factor = PseudorangeFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 6.9,
            variance: 1.0,
            sat_clock_bias: 0.5,
            tropo_dry_delay: 0.3,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.1]);
        let res = factor.residual(&state);
        assert!(
            res[0].abs() < 1e-10,
            "PseudorangeFactor residual with bias should be 0. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_pseudorange_jacobian_numerical() {
        let factor = PseudorangeFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.0,
            variance: 2.5,
            sat_clock_bias: 0.5,
            tropo_dry_delay: 0.3,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.1]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);
        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "PseudorangeFactor jacobian diverges from numerical"
        );
    }

    #[test]
    fn test_pseudorange_jacobian_nonzero_offset() {
        // Test jacobian when receiver is not at origin
        // sat_pos=[2,3,6], rx=[1,1,1] -> dist=sqrt(1+4+25)=sqrt(30)=5.477...
        // dx/dist=1/5.477=0.1826, dy/dist=2/5.477=0.3651, dz/dist=5/5.477=0.9129
        let factor = PseudorangeFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 20.0,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![1.0, 1.0, 1.0, 0.0]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);
        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "PseudorangeFactor jacobian diverges from numerical at nonzero offset"
        );
    }

    #[test]
    fn test_pseudorange_information() {
        let factor = PseudorangeFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.0,
            variance: 4.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            robust_threshold: 3.0,
        };
        let info = factor.information();
        assert!(
            (info[(0, 0)] - 0.25).abs() < 1e-15,
            "PseudorangeFactor information should be 1/variance=0.25"
        );
    }

    #[test]
    fn test_pseudorange_information_tiny_variance() {
        let factor = PseudorangeFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.0,
            variance: 0.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            robust_threshold: 3.0,
        };
        let info = factor.information();
        // variance.max(1e-9) = 1e-9, so info = 1.0 / 1e-9 = 1e9
        assert!(
            (info[(0, 0)] - 1e9).abs() < 1.0,
            "PseudorangeFactor information should clamp near-zero variance to 1e9"
        );
    }

    #[test]
    fn test_pseudorange_robust_methods() {
        let factor = PseudorangeFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.0,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            robust_threshold: 3.0,
        };
        assert_eq!(factor.robust_threshold(), Some(3.0));
        assert!(factor.is_cauchy_rejectable());
    }

    #[test]
    fn test_pseudorange_jacobian_dist_near_zero() {
        // When sat_pos ~ receiver_pos, dist < 1e-6, jacobian should be zeros for position
        let factor = PseudorangeFactor {
            sat_pos: Vector3::new(0.0, 0.0, 1e-7),
            measured_pr: 0.0,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
        let jac = factor.jacobian(&state);
        // All entries should be 0 since dist < 1e-6 (entire if-block is skipped)
        assert_eq!(
            jac[(0, 0)], 0.0,
            "x jacobian should be 0 when dist near 0"
        );
        assert_eq!(
            jac[(0, 1)], 0.0,
            "y jacobian should be 0 when dist near 0"
        );
        assert_eq!(
            jac[(0, 2)], 0.0,
            "z jacobian should be 0 when dist near 0"
        );
        assert_eq!(
            jac[(0, 3)], 0.0,
            "dt jacobian should be 0 when dist near 0 (entire if-block skipped)"
        );
    }

    #[test]
    fn test_pseudorange_jacobian_with_zwd() {
        // With index_zwd = Some(4), the jacobian should have a 0.0 entry at that index
        let factor = PseudorangeFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.0,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 0.0]);
        let jac = factor.jacobian(&state);
        assert_eq!(jac[(0, 4)], 0.0, "ZWD jacobian should be 0.0 (disabled)");
    }

    // ==================== CarrierPhaseFactor tests ====================

    #[test]
    fn test_carrier_phase_residual_zero() {
        let factor = CarrierPhaseFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.0,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_amb: 4,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 0.0]);
        let res = factor.residual(&state);
        assert!(
            res[0].abs() < 1e-10,
            "CarrierPhaseFactor residual should be 0 for correct state. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_carrier_phase_residual_with_ambiguity() {
        // dist=7.0, dt=0.1, amb=2.0, wavelength=0.19 -> expected = 7.0+0.1+2.0*0.19 = 7.48
        // measured_cp=7.48 -> residual=0
        let factor = CarrierPhaseFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.48,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_amb: 4,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.1, 2.0]);
        let res = factor.residual(&state);
        assert!(
            res[0].abs() < 1e-10,
            "CarrierPhaseFactor residual with ambiguity should be 0. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_carrier_phase_jacobian_numerical() {
        let factor = CarrierPhaseFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 100.0,
            variance: 0.01,
            sat_clock_bias: 0.5,
            tropo_dry_delay: 0.3,
            map_wet: 0.0,
            wavelength: 0.19,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_amb: 4,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.1, 1.0]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);
        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "CarrierPhaseFactor jacobian diverges from numerical"
        );
    }

    #[test]
    fn test_carrier_phase_jacobian_ambiguity_term() {
        // The ambiguity Jacobian should be -wavelength
        let factor = CarrierPhaseFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.0,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_amb: 4,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 1.0]);
        let jac = factor.jacobian(&state);
        assert!(
            (jac[(0, 4)] - (-0.19)).abs() < 1e-15,
            "Ambiguity jacobian should be -wavelength = -0.19"
        );
    }

    #[test]
    fn test_carrier_phase_information() {
        let factor = CarrierPhaseFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.0,
            variance: 4.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_amb: 4,
            robust_threshold: 3.0,
        };
        let info = factor.information();
        assert!(
            (info[(0, 0)] - 0.25).abs() < 1e-15,
            "CarrierPhaseFactor information should be 1/variance=0.25"
        );
    }

    #[test]
    fn test_carrier_phase_robust_methods() {
        let factor = CarrierPhaseFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.0,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_amb: 4,
            robust_threshold: 3.0,
        };
        assert_eq!(factor.robust_threshold(), Some(3.0));
        assert!(factor.is_cauchy_rejectable());
    }

    #[test]
    fn test_carrier_phase_jacobian_with_zwd() {
        // CarrierPhaseFactor jacobian with index_zwd set should have 0.0 at the ZWD index
        let factor = CarrierPhaseFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.0,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            index_amb: 5,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0]);
        let jac = factor.jacobian(&state);
        // ZWD jacobian should be 0.0 (disabled in carrier phase)
        assert!(
            (jac[(0, 4)]).abs() < 1e-15,
            "ZWD jacobian should be 0.0, got {}",
            jac[(0, 4)]
        );
        // Ambiguity jacobian should still be -wavelength
        assert!(
            (jac[(0, 5)] - (-0.19)).abs() < 1e-15,
            "Ambiguity jacobian should be -wavelength = -0.19, got {}",
            jac[(0, 5)]
        );
    }

    #[test]
    fn test_carrier_phase_jacobian_numerical_with_zwd() {
        // Full numerical jacobian check for CarrierPhaseFactor with ZWD enabled
        let factor = CarrierPhaseFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.0,
            variance: 0.01,
            sat_clock_bias: 0.5,
            tropo_dry_delay: 0.3,
            map_wet: 0.0,
            wavelength: 0.19,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            index_amb: 5,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![1.0, 2.0, 3.0, 0.1, 0.05, 2.0]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);
        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "CarrierPhaseFactor jacobian with ZWD diverges from numerical"
        );
    }

    // ==================== ErrorStatePseudorangeFactor extended tests ====================

    #[test]
    fn test_error_state_pseudorange_residual_zero() {
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.0,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
        let res = factor.residual(&delta);
        assert!(
            res[0].abs() < 1e-10,
            "ErrorStatePseudorangeFactor residual should be 0. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_error_state_pseudorange_residual_zwd() {
        // ZWD: nominal_zwd=0.5, delta[index_zwd]=0 -> zwd=0.5
        // dist=7.0, dt=0, zwd*map_wet=0.5*1.0=0.5 -> expected=7.5
        // measured_pr=7.0 -> residual=-0.5
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.0,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 1.0,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.5,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 0.5]);
        // zwd = 0.5 + 0.5 = 1.0, zwd*map_wet = 1.0
        // expected = 7.0 + 0 + 0 + 0 + 1.0 = 8.0
        // residual = 7.0 - 8.0 = -1.0
        let res = factor.residual(&delta);
        assert!(
            (res[0] - (-1.0)).abs() < 1e-10,
            "ErrorStatePseudorangeFactor ZWD residual incorrect. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_error_state_pseudorange_galileo_clock() {
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.18,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.1,
            nominal_dt_gal: 0.05,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: Some(4),
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Galileo,
                prn: 2,
            },
        };
        // delta: state=[0,0,0,0.02,0.03]
        // dt = nominal_dt(0.1) + delta[3](0.02) + nominal_dt_gal(0.05) + delta[4](0.03) = 0.20
        // dist = 7.0, expected = 7.0 + 0.20 = 7.20
        // measured_pr = 7.18, residual = 7.18 - 7.20 = -0.02
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.02, 0.03]);
        let res = factor.residual(&delta);
        assert!(
            (res[0] - (-0.02)).abs() < 1e-10,
            "ErrorStatePseudorangeFactor Galileo clock residual incorrect. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_error_state_pseudorange_beidou_clock() {
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.18,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.1,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.05,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: None,
            index_dt_bds: Some(4),
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Beidou,
                prn: 3,
            },
        };
        // dt = 0.1 + 0.02 + 0.05 + 0.03 = 0.20
        // dist = 7.0, expected = 7.20, residual = 7.18 - 7.20 = -0.02
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.02, 0.03]);
        let res = factor.residual(&delta);
        assert!(
            (res[0] - (-0.02)).abs() < 1e-10,
            "ErrorStatePseudorangeFactor Beidou clock residual incorrect. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_error_state_pseudorange_glonass_clock() {
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.18,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.1,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.05,
            nominal_zwd: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: Some(4),
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Glonass,
                prn: 4,
            },
        };
        // dt = 0.1 + 0.02 + 0.05 + 0.03 = 0.20
        // dist = 7.0, expected = 7.20, residual = 7.18 - 7.20 = -0.02
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.02, 0.03]);
        let res = factor.residual(&delta);
        assert!(
            (res[0] - (-0.02)).abs() < 1e-10,
            "ErrorStatePseudorangeFactor Glonass clock residual incorrect. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_error_state_pseudorange_jacobian_zwd() {
        // With index_zwd, the jacobian should have -map_wet at that index
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.0,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 3.1,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        let state = DVector::from_vec(vec![5.0, -3.0, 2.0, 0.1, 0.05]);
        let anal_jac = factor.jacobian(&state);
        assert!(
            (anal_jac[(0, 4)] - (-3.1)).abs() < 1e-10,
            "ZWD jacobian should be -map_wet = -3.1"
        );
    }

    #[test]
    fn test_error_state_pseudorange_information() {
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.0,
            variance: 4.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        let info = factor.information();
        assert!(
            (info[(0, 0)] - 0.25).abs() < 1e-15,
            "ErrorStatePseudorangeFactor information should be 1/variance=0.25"
        );
    }

    #[test]
    fn test_error_state_pseudorange_robust_methods() {
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 5.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.0,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        assert_eq!(factor.robust_threshold(), Some(5.0));
        assert!(factor.is_cauchy_rejectable());
    }

    // ==================== ErrorStateCarrierPhaseFactor extended tests ====================

    #[test]
    fn test_error_state_carrier_phase_residual_zero() {
        let factor = ErrorStateCarrierPhaseFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.0,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            nominal_amb: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            index_amb: 4,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 0.0]);
        let res = factor.residual(&delta);
        assert!(
            res[0].abs() < 1e-10,
            "ErrorStateCarrierPhaseFactor residual should be 0. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_error_state_carrier_phase_residual_with_amb() {
        // nominal_amb = 1.0, delta[index_amb] = 2.0 -> amb = 3.0
        // dist = 7.0, dt=0, expected = 7.0 + 3.0 * 0.19 = 7.57
        // measured_cp = 7.57 -> residual = 0
        let factor = ErrorStateCarrierPhaseFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.57,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            nominal_amb: 1.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            index_amb: 4,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 2.0]);
        let res = factor.residual(&delta);
        assert!(
            res[0].abs() < 1e-10,
            "ErrorStateCarrierPhaseFactor amb residual should be 0. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_error_state_carrier_phase_galileo_clock() {
        let factor = ErrorStateCarrierPhaseFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.2,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.1,
            nominal_dt_gal: 0.05,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            nominal_amb: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: Some(4),
            index_dt_bds: None,
            index_dt_glo: None,
            index_amb: 5,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Galileo,
                prn: 2,
            },
        };
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.02, 0.03, 0.0]);
        // dt = 0.1 + 0.02 + 0.05 + 0.03 = 0.20
        // dist = 7.0, expected = 7.0 + 0.20 = 7.20
        // measured_cp = 7.20 -> residual = 0
        let res = factor.residual(&delta);
        assert!(
            res[0].abs() < 1e-10,
            "ErrorStateCarrierPhaseFactor Galileo residual should be 0. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_error_state_carrier_phase_jacobian_matches_numerical() {
        let factor = ErrorStateCarrierPhaseFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(20000000.0, 10000000.0, 5000000.0),
            measured_cp: 120000000.0,
            variance: 1.0,
            sat_clock_bias: 0.0001,
            tropo_dry_delay: 2.3,
            map_wet: 3.1,
            wavelength: 0.19,
            nominal_rx: 1000.0,
            nominal_ry: 2000.0,
            nominal_rz: 3000.0,
            nominal_dt: 0.0002,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            nominal_amb: 50.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            index_amb: 5,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        let state = DVector::from_vec(vec![5.0, -3.0, 2.0, 0.1, 0.05, -2.0]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);
        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "ErrorStateCarrierPhaseFactor jacobian diverges from numerical"
        );
    }

    #[test]
    fn test_error_state_carrier_phase_information() {
        let factor = ErrorStateCarrierPhaseFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.0,
            variance: 10.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            nominal_amb: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            index_amb: 4,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        let info = factor.information();
        assert!(
            (info[(0, 0)] - 0.1).abs() < 1e-15,
            "ErrorStateCarrierPhaseFactor information should be 1/variance=0.1"
        );
    }

    #[test]
    fn test_error_state_carrier_phase_robust_methods() {
        let factor = ErrorStateCarrierPhaseFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.0,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            nominal_amb: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            index_amb: 4,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        assert_eq!(factor.robust_threshold(), Some(3.0));
        assert!(factor.is_cauchy_rejectable());
    }

    // ==================== ErrorStateDopplerFactor tests ====================

    #[test]
    fn test_error_state_doppler_jacobian_numerical() {
        let los = Vector3::new(0.6, 0.0, 0.8);
        let factor = ErrorStateDopplerFactor {
            robust_threshold: 3.0,
            los,
            sat_vel: Vector3::new(10.0, 0.0, 10.0),
            measured_doppler_hz: 0.0,
            variance: 2.5,
            wavelength: 0.19,
            sat_clock_drift: 0.0,
            nominal_vx: 0.0,
            nominal_vy: 0.0,
            nominal_vz: 0.0,
            nominal_cdt: 0.0,
            index_vx: 0,
            index_vy: 1,
            index_vz: 2,
            index_cdt: 3,
        };
        let delta = DVector::from_vec(vec![1.0, 2.0, 3.0, 0.1]);
        let anal_jac = factor.jacobian(&delta);
        let num_jac = numerical_jacobian(&factor, &delta);
        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "Doppler jacobian diverges from numerical"
        );
    }

    #[test]
    fn test_error_state_doppler_jacobian_components() {
        let los = Vector3::new(0.6, 0.8, 0.0);
        let factor = ErrorStateDopplerFactor {
            robust_threshold: 3.0,
            los,
            sat_vel: Vector3::new(10.0, 0.0, 5.0),
            measured_doppler_hz: -42.10526315789474, // = -(range_rate)/0.19
            variance: 1.0,
            wavelength: 0.19,
            sat_clock_drift: 0.0,
            nominal_vx: 0.0,
            nominal_vy: 0.0,
            nominal_vz: 0.0,
            nominal_cdt: 0.0,
            index_vx: 0,
            index_vy: 1,
            index_vz: 2,
            index_cdt: 3,
        };
        let delta = DVector::from_vec(vec![2.0, 0.0, 1.0, 0.0]);
        let jac = factor.jacobian(&delta);
        // Jacobian should be [los.x, los.y, los.z, -1.0] = [0.6, 0.8, 0.0, -1.0]
        assert!((jac[(0, 0)] - 0.6).abs() < 1e-15);
        assert!((jac[(0, 1)] - 0.8).abs() < 1e-15);
        assert!((jac[(0, 2)]).abs() < 1e-15);
        assert!((jac[(0, 3)] - (-1.0)).abs() < 1e-15);
    }

    #[test]
    fn test_error_state_doppler_information() {
        let factor = ErrorStateDopplerFactor {
            robust_threshold: 3.0,
            los: Vector3::new(1.0, 0.0, 0.0),
            sat_vel: Vector3::new(10.0, 0.0, 0.0),
            measured_doppler_hz: 0.0,
            variance: 4.0,
            wavelength: 0.19,
            sat_clock_drift: 0.0,
            nominal_vx: 0.0,
            nominal_vy: 0.0,
            nominal_vz: 0.0,
            nominal_cdt: 0.0,
            index_vx: 0,
            index_vy: 1,
            index_vz: 2,
            index_cdt: 3,
        };
        let info = factor.information();
        assert!(
            (info[(0, 0)] - 0.25).abs() < 1e-15,
            "Doppler information should be 1/variance=0.25"
        );
    }

    #[test]
    fn test_error_state_pseudorange_jacobian_galileo_clock() {
        // Verify that the Galileo-specific clock jacobian entry is -1.0 via numerical check
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(20000000.0, 10000000.0, 5000000.0),
            measured_pr: 22000000.0,
            variance: 1.0,
            sat_clock_bias: 0.0001,
            tropo_dry_delay: 2.3,
            map_wet: 3.1,
            nominal_rx: 1000.0,
            nominal_ry: 2000.0,
            nominal_rz: 3000.0,
            nominal_dt: 0.0002,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            index_dt_gal: Some(5),
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Galileo,
                prn: 2,
            },
        };

        let state = DVector::from_vec(vec![5.0, -3.0, 2.0, 0.1, 0.05, 0.01]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);

        // Galileo clock column should have -1.0
        assert!(
            (anal_jac[(0, 5)] - (-1.0)).abs() < 1e-10,
            "Galileo clock jacobian should be -1.0, got {}",
            anal_jac[(0, 5)]
        );

        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "Analytical PR jacobian with Galileo clock diverges from numerical"
        );
    }

    #[test]
    fn test_error_state_pseudorange_jacobian_glonass_clock() {
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(20000000.0, 10000000.0, 5000000.0),
            measured_pr: 22000000.0,
            variance: 1.0,
            sat_clock_bias: 0.0001,
            tropo_dry_delay: 2.3,
            map_wet: 3.1,
            nominal_rx: 1000.0,
            nominal_ry: 2000.0,
            nominal_rz: 3000.0,
            nominal_dt: 0.0002,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(5),
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: Some(4),
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Glonass,
                prn: 5,
            },
        };
        let state = DVector::from_vec(vec![5.0, -3.0, 2.0, 0.1, 0.01, 0.05]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);
        assert!(
            (anal_jac[(0, 4)] - (-1.0)).abs() < 1e-10,
            "Glonass clock jacobian should be -1.0 at index 4, got {}",
            anal_jac[(0, 4)]
        );
        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "PR jacobian with Glonass clock diverges from numerical"
        );
    }

    #[test]
    fn test_error_state_carrier_phase_jacobian_beidou_clock() {
        let factor = ErrorStateCarrierPhaseFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(20000000.0, 10000000.0, 5000000.0),
            measured_cp: 120000000.0,
            variance: 1.0,
            sat_clock_bias: 0.0001,
            tropo_dry_delay: 2.3,
            map_wet: 3.1,
            wavelength: 0.19,
            nominal_rx: 1000.0,
            nominal_ry: 2000.0,
            nominal_rz: 3000.0,
            nominal_dt: 0.0002,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            nominal_amb: 50.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(5),
            index_dt_gal: None,
            index_dt_bds: Some(4),
            index_dt_glo: None,
            index_amb: 6,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Beidou,
                prn: 3,
            },
        };
        let state = DVector::from_vec(vec![5.0, -3.0, 2.0, 0.1, 0.01, 0.05, -2.0]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);
        assert!(
            (anal_jac[(0, 4)] - (-1.0)).abs() < 1e-10,
            "Beidou clock jacobian should be -1.0, got {}",
            anal_jac[(0, 4)]
        );
        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "CP jacobian with Beidou clock diverges from numerical"
        );
    }

    #[test]
    fn test_error_state_carrier_phase_jacobian_galileo_clock() {
        let factor = ErrorStateCarrierPhaseFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(20000000.0, 10000000.0, 5000000.0),
            measured_cp: 120000000.0,
            variance: 1.0,
            sat_clock_bias: 0.0001,
            tropo_dry_delay: 2.3,
            map_wet: 3.1,
            wavelength: 0.19,
            nominal_rx: 1000.0,
            nominal_ry: 2000.0,
            nominal_rz: 3000.0,
            nominal_dt: 0.0002,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            nominal_amb: 50.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(6),
            index_dt_gal: Some(4),
            index_dt_bds: None,
            index_dt_glo: Some(5), // unused but present in state to check they don't interfere
            index_amb: 7,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Galileo,
                prn: 2,
            },
        };
        let state = DVector::from_vec(vec![5.0, -3.0, 2.0, 0.1, 0.01, -0.02, 0.05, -2.0]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);
        // Galileo clock column should be -1.0
        assert!(
            (anal_jac[(0, 4)] - (-1.0)).abs() < 1e-10,
            "Galileo clock jacobian should be -1.0, got {}",
            anal_jac[(0, 4)]
        );
        // Glonass clock column (index 5) should be 0 since sat is Galileo
        assert!(
            anal_jac[(0, 5)].abs() < 1e-10,
            "Glonass clock jacobian should be 0 for Galileo sat, got {}",
            anal_jac[(0, 5)]
        );
        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "CP jacobian with Galileo clock diverges from numerical"
        );
    }

    // ==================== ErrorStateDopplerFactor additional tests ====================

    #[test]
    fn test_error_state_doppler_residual_zero() {
        let los = Vector3::new(1.0, 0.0, 0.0);
        let sat_vel = Vector3::new(10.0, 0.0, 0.0);
        // True rx_vel = 2 m/s toward satellite
        // pred_rr = los.dot(sat_vel - rx_vel) + cdt = 1*(10-2) + 0 = 8 m/s
        // doppler_hz = -rr / wavelength = -8/0.19
        let wavelength = 0.19;
        let doppler_hz = -8.0 / wavelength;
        let factor = ErrorStateDopplerFactor {
            los,
            sat_vel,
            measured_doppler_hz: doppler_hz,
            variance: 1.0,
            wavelength,
            sat_clock_drift: 0.0,
            nominal_vx: 0.0,
            nominal_vy: 0.0,
            nominal_vz: 0.0,
            nominal_cdt: 0.0,
            index_vx: 0,
            index_vy: 1,
            index_vz: 2,
            index_cdt: 3,
            robust_threshold: 3.0,
        };
        let delta = DVector::from_vec(vec![2.0, 0.0, 0.0, 0.0]);
        let res = factor.residual(&delta);
        assert!(
            res[0].abs() < 1e-6,
            "Doppler residual should be 0 for correct velocity. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_error_state_doppler_residual_with_cdt() {
        let los = Vector3::new(1.0, 0.0, 0.0);
        let sat_vel = Vector3::new(10.0, 0.0, 0.0);
        // pred_rr = los.dot(sat_vel - rx_vel) + cdt = 1*(10-2) + 0.5 = 8.5
        // doppler = -8.5 / 0.19
        let wavelength = 0.19;
        let doppler_hz = -8.5 / wavelength;
        let factor = ErrorStateDopplerFactor {
            los,
            sat_vel,
            measured_doppler_hz: doppler_hz,
            variance: 1.0,
            wavelength,
            sat_clock_drift: 0.0,
            nominal_vx: 0.0,
            nominal_vy: 0.0,
            nominal_vz: 0.0,
            nominal_cdt: 0.0,
            index_vx: 0,
            index_vy: 1,
            index_vz: 2,
            index_cdt: 3,
            robust_threshold: 3.0,
        };
        let delta = DVector::from_vec(vec![2.0, 0.0, 0.0, 0.5]);
        let res = factor.residual(&delta);
        assert!(
            res[0].abs() < 1e-6,
            "Doppler residual with cdt should be 0. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_error_state_doppler_residual_mismatch() {
        let los = Vector3::new(1.0, 0.0, 0.0);
        let sat_vel = Vector3::new(10.0, 0.0, 0.0);
        // True rr = 10-2 = 8, doppler = -8/0.19 ≈ -42.105
        // But state gives rx_vel=0, so pred_rr = 1*(10-0) + 0 = 10
        // Residual = observed_rr - predicted_rr = (-8) - 10 = -18
        let wavelength = 0.19;
        let doppler_hz = -8.0 / wavelength;
        let factor = ErrorStateDopplerFactor {
            los,
            sat_vel,
            measured_doppler_hz: doppler_hz,
            variance: 1.0,
            wavelength,
            sat_clock_drift: 0.0,
            nominal_vx: 0.0,
            nominal_vy: 0.0,
            nominal_vz: 0.0,
            nominal_cdt: 0.0,
            index_vx: 0,
            index_vy: 1,
            index_vz: 2,
            index_cdt: 3,
            robust_threshold: 3.0,
        };
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0]); // Wrong velocity
        let res = factor.residual(&delta);
        // observed_rr = -doppler * wavelength = 8.0
        // predicted_rr = 1*(10-0) + 0 = 10.0
        // residual = 8.0 - 10.0 = -2.0
        assert!(
            (res[0] - (-2.0)).abs() < 1e-6,
            "Doppler residual mismatch should be -2.0. Got {}",
            res[0]
        );
    }

    #[test]
    fn test_error_state_doppler_robust_threshold() {
        let factor = ErrorStateDopplerFactor {
            robust_threshold: 3.0,
            los: Vector3::new(1.0, 0.0, 0.0),
            sat_vel: Vector3::new(10.0, 0.0, 0.0),
            measured_doppler_hz: 0.0,
            variance: 1.0,
            wavelength: 0.19,
            sat_clock_drift: 0.0,
            nominal_vx: 0.0,
            nominal_vy: 0.0,
            nominal_vz: 0.0,
            nominal_cdt: 0.0,
            index_vx: 0,
            index_vy: 1,
            index_vz: 2,
            index_cdt: 3,
        };
        assert_eq!(factor.robust_threshold(), Some(3.0));
        assert!(factor.is_cauchy_rejectable());
    }

    // ==================== Additional PseudorangeFactor coverage ====================

    #[test]
    fn test_pseudorange_residual_with_zwd_some() {
        // Cover the Some(index) branch of index_zwd in residual()
        let factor = PseudorangeFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.0,
            variance: 2.5,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 99.0]);
        let res = factor.residual(&state);
        // zwd is read but ignored in residual, so residual should be 0
        assert!(
            res[0].abs() < 1e-10,
            "Pseudorange residual with ZWD Some should be 0, got {}",
            res[0]
        );
    }

    // ==================== Additional CarrierPhaseFactor coverage ====================

    #[test]
    fn test_carrier_phase_jacobian_dist_near_zero() {
        // Cover dist < 1e-6: entire if-block in jacobian skipped, all entries 0
        let factor = CarrierPhaseFactor {
            sat_pos: Vector3::new(0.0, 0.0, 1e-7),
            measured_cp: 0.0,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_amb: 4,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 1.0]);
        let jac = factor.jacobian(&state);
        for i in 0..5 {
            assert_eq!(
                jac[(0, i)], 0.0,
                "CarrierPhase jacobian[0,{}] should be 0 when dist near 0",
                i
            );
        }
    }

    #[test]
    fn test_carrier_phase_residual_with_zwd_some() {
        // Cover the Some(index) branch of index_zwd in CP residual
        let factor = CarrierPhaseFactor {
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_cp: 7.0,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            index_amb: 5,
            robust_threshold: 3.0,
        };
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 42.0, 0.0]);
        let res = factor.residual(&state);
        // zwd is read but ignored, amb=0, dt=0, dist=7, expected=7, residual=0
        assert!(
            res[0].abs() < 1e-10,
            "CarrierPhase residual with ZWD Some should be 0, got {}",
            res[0]
        );
    }

    // ==================== Additional ErrorStatePseudorangeFactor coverage ====================

    #[test]
    fn test_error_state_pseudorange_residual_gps_clock_offset() {
        // Test GPS clock delta: dt = nominal_dt + delta[index_dt] (no ISB for GPS)
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(2.0, 3.0, 6.0),
            measured_pr: 7.5,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.1,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        // delta[3] = 0.4, so dt = 0.1 + 0.4 = 0.5
        // dist = 7.0, expected = 7.0 + 0.5 = 7.5, residual = 7.5 - 7.5 = 0
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.4]);
        let res = factor.residual(&delta);
        assert!(
            res[0].abs() < 1e-10,
            "GPS clock offset residual should be 0, got {}",
            res[0]
        );
    }

    #[test]
    fn test_error_state_pseudorange_jacobian_dist_near_zero() {
        // Cover dist < 1e-6: entire if-block in jacobian skipped
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(0.0, 0.0, 1e-7),
            measured_pr: 0.0,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
        let jac = factor.jacobian(&delta);
        for i in 0..4 {
            assert_eq!(
                jac[(0, i)], 0.0,
                "ErrorState PR jacobian[0,{}] should be 0 when dist near 0",
                i
            );
        }
    }

    #[test]
    fn test_error_state_pseudorange_jacobian_dist_near_zero_with_zwd() {
        // dist < 1e-6 with index_zwd set: entire if-block skipped
        let factor = ErrorStatePseudorangeFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(0.0, 0.0, 1e-7),
            measured_pr: 0.0,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 2.5,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 0.0]);
        let jac = factor.jacobian(&delta);
        for i in 0..5 {
            assert_eq!(
                jac[(0, i)], 0.0,
                "ErrorState PR jacobian[0,{}] should be 0 when dist near 0 with ZWD",
                i
            );
        }
    }

    // ==================== Additional ErrorStateCarrierPhaseFactor coverage ====================

    #[test]
    fn test_error_state_carrier_phase_jacobian_dist_near_zero() {
        // Cover dist < 1e-6: entire if-block in CP jacobian skipped
        let factor = ErrorStateCarrierPhaseFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(0.0, 0.0, 1e-7),
            measured_cp: 0.0,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength: 0.19,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            nominal_amb: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            index_amb: 4,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 0.0]);
        let jac = factor.jacobian(&delta);
        for i in 0..5 {
            assert_eq!(
                jac[(0, i)], 0.0,
                "ErrorState CP jacobian[0,{}] should be 0 when dist near 0",
                i
            );
        }
    }

    #[test]
    fn test_error_state_carrier_phase_jacobian_glonass_clock() {
        // Test Glonass clock ISB jacobian entry for ErrorStateCarrierPhaseFactor
        let factor = ErrorStateCarrierPhaseFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(20000000.0, 10000000.0, 5000000.0),
            measured_cp: 120000000.0,
            variance: 1.0,
            sat_clock_bias: 0.0001,
            tropo_dry_delay: 2.3,
            map_wet: 3.1,
            wavelength: 0.19,
            nominal_rx: 1000.0,
            nominal_ry: 2000.0,
            nominal_rz: 3000.0,
            nominal_dt: 0.0002,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            nominal_amb: 50.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(5),
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: Some(4),
            index_amb: 6,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Glonass,
                prn: 5,
            },
        };
        let state = DVector::from_vec(vec![5.0, -3.0, 2.0, 0.1, 0.01, 0.05, -2.0]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);
        assert!(
            (anal_jac[(0, 4)] - (-1.0)).abs() < 1e-10,
            "Glonass clock jacobian should be -1.0, got {}",
            anal_jac[(0, 4)]
        );
        assert!(
            (anal_jac - num_jac).norm() < 2e-3,
            "CP jacobian with Glonass clock diverges from numerical"
        );
    }

    #[test]
    fn test_error_state_carrier_phase_jacobian_dist_near_zero_with_zwd() {
        // dist < 1e-6 with index_zwd and index_amb set: entire if-block skipped
        let factor = ErrorStateCarrierPhaseFactor {
            robust_threshold: 3.0,
            sat_pos: Vector3::new(0.0, 0.0, 1e-7),
            measured_cp: 0.0,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 2.5,
            wavelength: 0.19,
            nominal_rx: 0.0,
            nominal_ry: 0.0,
            nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0,
            nominal_dt_bds: 0.0,
            nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            nominal_amb: 1.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: Some(4),
            index_dt_gal: None,
            index_dt_bds: None,
            index_dt_glo: None,
            index_amb: 5,
            sat_id: gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
        };
        let delta = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let jac = factor.jacobian(&delta);
        for i in 0..6 {
            assert_eq!(
                jac[(0, i)], 0.0,
                "ErrorState CP jacobian[0,{}] should be 0 when dist near 0 with ZWD",
                i
            );
        }
    }
}

/// Error-State Factor for a Doppler measurement.
pub struct ErrorStateDopplerFactor {
    pub los: Vector3<f64>, // Line of sight vector from receiver to satellite
    pub sat_vel: Vector3<f64>,
    pub measured_doppler_hz: f64,
    pub variance: f64,
    pub wavelength: f64,
    pub sat_clock_drift: f64,

    pub nominal_vx: f64,
    pub nominal_vy: f64,
    pub nominal_vz: f64,
    pub nominal_cdt: f64, // Receiver clock drift

    pub index_vx: usize,
    pub index_vy: usize,
    pub index_vz: usize,
    pub index_cdt: usize,
    pub robust_threshold: f64,
}

impl Factor for ErrorStateDopplerFactor {
    fn residual(&self, delta: &DVector<f64>) -> DVector<f64> {
        let vx = self.nominal_vx + delta[self.index_vx];
        let vy = self.nominal_vy + delta[self.index_vy];
        let vz = self.nominal_vz + delta[self.index_vz];
        let cdt = self.nominal_cdt + delta[self.index_cdt];

        let v_rx = Vector3::new(vx, vy, vz);
        let predicted_rr = self.los.dot(&(self.sat_vel - v_rx)) + cdt - self.sat_clock_drift;
        let observed_rr = -self.measured_doppler_hz * self.wavelength;

        let res = observed_rr - predicted_rr;

        DVector::from_vec(vec![res])
    }

    fn jacobian(&self, delta: &DVector<f64>) -> DMatrix<f64> {
        let mut jac = DMatrix::zeros(1, delta.len());
        // predicted_rr = los.dot(sat_vel - rx_vel) + cdt
        // predicted_rr = los.dot(sat_vel) - los.dot(rx_vel) + cdt
        // res = observed_rr - predicted_rr
        // res = observed_rr - los.dot(sat_vel) + los.dot(rx_vel) - cdt
        // d(res)/d(rx_vel) = +los
        // d(res)/d(cdt) = -1.0
        jac[(0, self.index_vx)] = self.los.x;
        jac[(0, self.index_vy)] = self.los.y;
        jac[(0, self.index_vz)] = self.los.z;
        jac[(0, self.index_cdt)] = -1.0;
        jac
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_element(1, 1, 1.0 / self.variance.max(1e-9))
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(self.robust_threshold)
    }

    fn is_cauchy_rejectable(&self) -> bool {
        true
    }
}

// ---------------------------------------------------------------------------
// Double-Differenced (DD) RTK measurement factors
// ---------------------------------------------------------------------------
// These use pre-computed DD innovations and LOS Jacobians from the EKF's
// H matrix. They are LINEAR factors: residual = H·δx - z, and the Jacobian
// is the constant H row. This avoids recomputing satellite geometry during
// factor graph optimization.

/// Linear DD pseudorange factor using pre-computed EKF innovation and LOS Jacobian.
///
/// `h_pos` is the DD line-of-sight vector `e_ref - e_sat` (unitless).
/// The residual is `h_pos · δx_pos - z` where `z` is the pre-computed
/// DD innovation (observed minus predicted at the EKF linearization point).
pub struct ErrorStateDdPseudorangeFactor {
    pub sat_id: gneiss_core::sat::SatelliteId,
    pub ref_sat_id: gneiss_core::sat::SatelliteId,
    /// Pre-computed DD innovation (meters).
    pub z: f64,
    /// DD LOS Jacobian for position: [e_ref - e_sat] (3 elements).
    pub h_pos: [f64; 3],
    /// Measurement variance (from EKF R diagonal).
    pub variance: f64,
    /// Factor graph state indices for position (px, py, pz) for this epoch.
    pub index_px: usize,
    pub index_py: usize,
    pub index_pz: usize,
    /// Total error-state dimension (for zero-padding Jacobian).
    pub total_dim: usize,
}

impl Factor for ErrorStateDdPseudorangeFactor {
    fn residual(&self, delta: &DVector<f64>) -> DVector<f64> {
        let pred = self.h_pos[0] * delta[self.index_px]
            + self.h_pos[1] * delta[self.index_py]
            + self.h_pos[2] * delta[self.index_pz];
        DVector::from_vec(vec![pred - self.z])
    }

    fn jacobian(&self, _delta: &DVector<f64>) -> DMatrix<f64> {
        let mut jac = DMatrix::zeros(1, self.total_dim);
        jac[(0, self.index_px)] = self.h_pos[0];
        jac[(0, self.index_py)] = self.h_pos[1];
        jac[(0, self.index_pz)] = self.h_pos[2];
        jac
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_element(1, 1, 1.0 / self.variance.max(1e-9))
    }
}

/// Linear DD carrier phase factor using pre-computed EKF innovation and LOS Jacobian.
///
/// The residual is `h_pos · δx_pos - z + δamb[sat] - δamb[ref]`,
/// where `δamb` is the error-state correction to the float ambiguity estimate
/// (in meters, since ambiguities are stored in meters in the EKF state).
pub struct ErrorStateDdCarrierPhaseFactor {
    pub sat_id: gneiss_core::sat::SatelliteId,
    pub ref_sat_id: gneiss_core::sat::SatelliteId,
    /// Pre-computed DD innovation (meters).
    pub z: f64,
    /// DD LOS Jacobian for position: [e_ref - e_sat] (3 elements).
    pub h_pos: [f64; 3],
    /// Measurement variance (from EKF R diagonal).
    pub variance: f64,
    /// Factor graph state indices for position (px, py, pz) for this epoch.
    pub index_px: usize,
    pub index_py: usize,
    pub index_pz: usize,
    /// Factor graph ambiguity index for the satellite.
    pub index_amb_sat: usize,
    /// Factor graph ambiguity index for the reference satellite.
    pub index_amb_ref: usize,
    /// Ionosphere state indices and scale in the factor graph ambiguity block.
    /// (sat_iono_fg_idx, ref_iono_fg_idx, iono_scale).
    pub iono_pair: Option<(usize, usize, f64)>,
    /// Total error-state dimension (for zero-padding Jacobian).
    pub total_dim: usize,
}

impl Factor for ErrorStateDdCarrierPhaseFactor {
    fn residual(&self, delta: &DVector<f64>) -> DVector<f64> {
        let pred = self.h_pos[0] * delta[self.index_px]
            + self.h_pos[1] * delta[self.index_py]
            + self.h_pos[2] * delta[self.index_pz]
            + delta[self.index_amb_sat]
            - delta[self.index_amb_ref];
        if let Some((si, ri, scale)) = self.iono_pair {
            DVector::from_vec(vec![pred - self.z + scale * (delta[si] - delta[ri])])
        } else {
            DVector::from_vec(vec![pred - self.z])
        }
    }

    fn jacobian(&self, _delta: &DVector<f64>) -> DMatrix<f64> {
        let mut jac = DMatrix::zeros(1, self.total_dim);
        jac[(0, self.index_px)] = self.h_pos[0];
        jac[(0, self.index_py)] = self.h_pos[1];
        jac[(0, self.index_pz)] = self.h_pos[2];
        jac[(0, self.index_amb_sat)] = 1.0;
        jac[(0, self.index_amb_ref)] = -1.0;
        if let Some((si, ri, scale)) = self.iono_pair {
            jac[(0, si)] = scale;
            jac[(0, ri)] = -scale;
        }
        jac
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_element(1, 1, 1.0 / self.variance.max(1e-9))
    }
}

#[cfg(test)]
mod doppler_tests {
    use super::*;
    use nalgebra::{DVector, Vector3};

    #[test]
    fn test_error_state_doppler_signs() {
        // Assume satellite at [1000, 0, 0], moving away from receiver at [0, 0, 0]
        let los = Vector3::new(1.0, 0.0, 0.0);
        let sat_vel = Vector3::new(10.0, 0.0, 0.0); // moving away at 10 m/s
        let _rx_vel = Vector3::new(2.0, 0.0, 0.0); // moving towards sat at 2 m/s

        // True relative velocity = sat_vel - rx_vel = [8.0, 0, 0]
        // Range rate = los.dot(sat_vel - rx_vel) = 8.0 m/s

        // doppler_hz = - range_rate / wavelength
        let wavelength = 0.19;
        let doppler_hz = -8.0 / wavelength;

        let factor = ErrorStateDopplerFactor {
            robust_threshold: 3.0,
            los,
            sat_vel,
            measured_doppler_hz: doppler_hz,
            variance: 1.0,
            wavelength,
            sat_clock_drift: 0.0,
            nominal_vx: 0.0,
            nominal_vy: 0.0,
            nominal_vz: 0.0,
            nominal_cdt: 0.0,
            index_vx: 0,
            index_vy: 1,
            index_vz: 2,
            index_cdt: 3,
        };

        // If the state is exactly correct (rx_vel = 2.0), residual should be 0.
        let state = DVector::from_vec(vec![2.0, 0.0, 0.0, 0.0]);
        let res = factor.residual(&state);

        assert!(
            res[0].abs() < 1e-6,
            "Doppler factor residual should be 0 for correct velocity. Got {}",
            res[0]
        );
    }
}
