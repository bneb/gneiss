use nalgebra::{DMatrix, DVector, Vector3};
use super::Factor;

/// Factor for a Pseudorange measurement.
pub struct PseudorangeFactor {
    pub sat_pos: Vector3<f64>,
    pub measured_pr: f64,
    pub variance: f64,
    pub sat_clock_bias: f64,
    pub tropo_dry_delay: f64,
    pub map_wet: f64,
    pub index_x: usize, // index of x in state
    pub index_y: usize, // index of y
    pub index_z: usize, // index of z
    pub index_dt: usize, // index of receiver clock bias
    pub index_zwd: Option<usize>, // index of zenith wet delay
}

impl Factor for PseudorangeFactor {
    fn residual(&self, state: &DVector<f64>) -> DVector<f64> {
        let rx = state[self.index_x];
        let ry = state[self.index_y];
        let rz = state[self.index_z];
        let dt = state[self.index_dt];
        
        let zwd = self.index_zwd.map(|idx| state[idx]).unwrap_or(0.0);
        
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
            jac[(0, self.index_x)] = dx / dist;  // negative of derivative of expected_pr
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
        Some(3.0) // 3-sigma threshold
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
}

impl Factor for CarrierPhaseFactor {
    fn residual(&self, state: &DVector<f64>) -> DVector<f64> {
        let rx = state[self.index_x];
        let ry = state[self.index_y];
        let rz = state[self.index_z];
        let dt = state[self.index_dt];
        let amb = state[self.index_amb];
        let is_iono_free = false;
        let zwd = self.index_zwd.map(|idx| state[idx]).unwrap_or(0.0);
        let var_cp = if is_iono_free { 0.001 } else { 0.01 }; // 10cm stddev for phase
        let huber_cp = 3.0; // 3-sigma (30cm) to reject cycle slips.

        let dx = self.sat_pos.x - rx;
        let dy = self.sat_pos.y - ry;
        let dz = self.sat_pos.z - rz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        
        let expected_cp = dist + dt - self.sat_clock_bias + self.tropo_dry_delay + amb * self.wavelength; // Ignore ZWD
        
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
        Some(3.0)
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
}

impl Factor for ErrorStatePseudorangeFactor {
    fn residual(&self, delta: &DVector<f64>) -> DVector<f64> {
        let rx = self.nominal_rx + delta[self.index_x];
        let ry = self.nominal_ry + delta[self.index_y];
        let rz = self.nominal_rz + delta[self.index_z];
        
        let mut dt = self.nominal_dt + delta[self.index_dt];
        match self.sat_id.constellation {
            gneiss_core::sat::Constellation::Galileo => {
                if let Some(idx) = self.index_dt_gal { dt += self.nominal_dt_gal + delta[idx]; }
            },
            gneiss_core::sat::Constellation::Beidou => {
                if let Some(idx) = self.index_dt_bds { dt += self.nominal_dt_bds + delta[idx]; }
            },
            gneiss_core::sat::Constellation::Glonass => {
                if let Some(idx) = self.index_dt_glo { dt += self.nominal_dt_glo + delta[idx]; }
            },
            _ => {},
        }
        
        let zwd = self.nominal_zwd + self.index_zwd.map(|i| delta[i]).unwrap_or(0.0);
        
        let dx = self.sat_pos.x - rx;
        let dy = self.sat_pos.y - ry;
        let dz = self.sat_pos.z - rz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        
        let expected_pr = dist + dt - self.sat_clock_bias + self.tropo_dry_delay + zwd * self.map_wet;
        
        DVector::from_vec(vec![self.measured_pr - expected_pr])
    }
    
    fn robust_threshold(&self) -> Option<f64> {
        Some(3.0) // 3-sigma threshold
    }
    
    fn is_cauchy_rejectable(&self) -> bool {
        true
    }
    
    fn jacobian(&self, delta: &DVector<f64>) -> DMatrix<f64> {
        let rx = self.nominal_rx + delta[self.index_x];
        let ry = self.nominal_ry + delta[self.index_y];
        let rz = self.nominal_rz + delta[self.index_z];
        
        let dx = self.sat_pos.x - rx;
        let dy = self.sat_pos.y - ry;
        let dz = self.sat_pos.z - rz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        
        let mut jac = DMatrix::zeros(1, delta.len());
        if dist > 1e-6 {
            jac[(0, self.index_x)] = dx / dist;
            jac[(0, self.index_y)] = dy / dist;
            jac[(0, self.index_z)] = dz / dist;
            jac[(0, self.index_dt)] = -1.0;
            
            match self.sat_id.constellation {
                gneiss_core::sat::Constellation::Galileo => {
                    if let Some(idx) = self.index_dt_gal { jac[(0, idx)] = -1.0; }
                },
                gneiss_core::sat::Constellation::Beidou => {
                    if let Some(idx) = self.index_dt_bds { jac[(0, idx)] = -1.0; }
                },
                gneiss_core::sat::Constellation::Glonass => {
                    if let Some(idx) = self.index_dt_glo { jac[(0, idx)] = -1.0; }
                },
                _ => {},
            }
            
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
}

impl Factor for ErrorStateCarrierPhaseFactor {
    fn residual(&self, delta: &DVector<f64>) -> DVector<f64> {
        let rx = self.nominal_rx + delta[self.index_x];
        let ry = self.nominal_ry + delta[self.index_y];
        let rz = self.nominal_rz + delta[self.index_z];
        
        let mut dt = self.nominal_dt + delta[self.index_dt];
        match self.sat_id.constellation {
            gneiss_core::sat::Constellation::Galileo => {
                if let Some(idx) = self.index_dt_gal { dt += self.nominal_dt_gal + delta[idx]; }
            },
            gneiss_core::sat::Constellation::Beidou => {
                if let Some(idx) = self.index_dt_bds { dt += self.nominal_dt_bds + delta[idx]; }
            },
            gneiss_core::sat::Constellation::Glonass => {
                if let Some(idx) = self.index_dt_glo { dt += self.nominal_dt_glo + delta[idx]; }
            },
            _ => {},
        }
        
        let amb = self.nominal_amb + delta[self.index_amb];
        let zwd = self.nominal_zwd + self.index_zwd.map(|idx| delta[idx]).unwrap_or(0.0);
        
        let dx = self.sat_pos.x - rx;
        let dy = self.sat_pos.y - ry;
        let dz = self.sat_pos.z - rz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        
        let expected_cp = dist + dt - self.sat_clock_bias + self.tropo_dry_delay + zwd * self.map_wet + amb * self.wavelength;
        DVector::from_vec(vec![self.measured_cp - expected_cp])
    }
    

    fn jacobian(&self, delta: &DVector<f64>) -> DMatrix<f64> {
        let rx = self.nominal_rx + delta[self.index_x];
        let ry = self.nominal_ry + delta[self.index_y];
        let rz = self.nominal_rz + delta[self.index_z];
        
        let dx = self.sat_pos.x - rx;
        let dy = self.sat_pos.y - ry;
        let dz = self.sat_pos.z - rz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        
        let mut jac = DMatrix::zeros(1, delta.len());
        if dist > 1e-6 {
            jac[(0, self.index_x)] = dx / dist;
            jac[(0, self.index_y)] = dy / dist;
            jac[(0, self.index_z)] = dz / dist;
            jac[(0, self.index_dt)] = -1.0;
            
            match self.sat_id.constellation {
                gneiss_core::sat::Constellation::Galileo => {
                    if let Some(idx) = self.index_dt_gal { jac[(0, idx)] = -1.0; }
                },
                gneiss_core::sat::Constellation::Beidou => {
                    if let Some(idx) = self.index_dt_bds { jac[(0, idx)] = -1.0; }
                },
                gneiss_core::sat::Constellation::Glonass => {
                    if let Some(idx) = self.index_dt_glo { jac[(0, idx)] = -1.0; }
                },
                _ => {},
            }
            
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
        Some(3.0) // 3-sigma threshold to reject cycle slips.
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
    where F: Factor {
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
            nominal_dt_gal: 0.0, nominal_dt_bds: 0.0, nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            index_x: 0, index_y: 1, index_z: 2, index_dt: 3, index_zwd: Some(4),
            index_dt_gal: None, index_dt_bds: None, index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId { constellation: gneiss_core::sat::Constellation::Gps, prn: 1 },
        };

        let state = DVector::from_vec(vec![5.0, -3.0, 2.0, 0.1, 0.05]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);

        println!("Anal PR: {}\nNum PR: {}", anal_jac, num_jac);        assert!((anal_jac - num_jac).norm() < 2e-3, "Analytical PR jacobian diverges from numerical");
    }

    #[test]
    fn test_error_state_carrier_phase_jacobian() {
        let factor = ErrorStateCarrierPhaseFactor {
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
            nominal_dt_gal: 0.0, nominal_dt_bds: 0.0, nominal_dt_glo: 0.0,
            nominal_zwd: 0.1,
            nominal_amb: 50.0,
            index_x: 0, index_y: 1, index_z: 2, index_dt: 3, index_zwd: Some(4), index_amb: 5,
            index_dt_gal: None, index_dt_bds: None, index_dt_glo: None,
            sat_id: gneiss_core::sat::SatelliteId { constellation: gneiss_core::sat::Constellation::Gps, prn: 1 },
        };

        let state = DVector::from_vec(vec![5.0, -3.0, 2.0, 0.1, 0.05, -2.0]);
        let anal_jac = factor.jacobian(&state);
        let num_jac = numerical_jacobian(&factor, &state);

        println!("Anal CP: {}\nNum CP: {}", anal_jac, num_jac);        assert!((anal_jac - num_jac).norm() < 2e-3, "Analytical CP jacobian diverges from numerical");
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
        Some(3.0)
    }
    
    fn is_cauchy_rejectable(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod doppler_tests {
    use super::*;
    use nalgebra::{Vector3, DVector};

    #[test]
    fn test_error_state_doppler_signs() {
        // Assume satellite at [1000, 0, 0], moving away from receiver at [0, 0, 0]
        let los = Vector3::new(1.0, 0.0, 0.0);
        let sat_vel = Vector3::new(10.0, 0.0, 0.0); // moving away at 10 m/s
        let rx_vel = Vector3::new(2.0, 0.0, 0.0);   // moving towards sat at 2 m/s
        
        // True relative velocity = sat_vel - rx_vel = [8.0, 0, 0]
        // Range rate = los.dot(sat_vel - rx_vel) = 8.0 m/s
        
        // doppler_hz = - range_rate / wavelength
        let wavelength = 0.19;
        let doppler_hz = -8.0 / wavelength;
        
        let factor = ErrorStateDopplerFactor {
            los, sat_vel, measured_doppler_hz: doppler_hz, variance: 1.0, wavelength,
            sat_clock_drift: 0.0, nominal_vx: 0.0, nominal_vy: 0.0, nominal_vz: 0.0, nominal_cdt: 0.0,
            index_vx: 0, index_vy: 1, index_vz: 2, index_cdt: 3,
        };
        
        // If the state is exactly correct (rx_vel = 2.0), residual should be 0.
        let state = DVector::from_vec(vec![2.0, 0.0, 0.0, 0.0]);
        let res = factor.residual(&state);
        
        assert!(res[0].abs() < 1e-6, "Doppler factor residual should be 0 for correct velocity. Got {}", res[0]);
    }
}
