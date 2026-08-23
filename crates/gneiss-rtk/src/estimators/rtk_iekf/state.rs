//! RTK State Vector & Covariance Management for DD-IEKF.

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use gneiss_core::time::GpsTime;

/// Key uniquely identifying a double-difference satellite pair on a frequency band.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct DoubleDiffKey {
    pub constellation_id: u8,
    pub sat: u16,
    pub ref_sat: u16,
    pub freq_band: u8,
}

/// Dynamic state vector for Double-Difference IEKF.
///
/// State vector layout:
///   x[0..3] — Position ECEF (m)
///   x[3..6] — Velocity ECEF (m/s)
///   x[6..]  — Double-difference carrier phase ambiguities (cycles)
///   x[last] — Rover ZWD residual (m), only when enabled (long-baseline
///     mode): non-dispersive wet-delay mis-modelling, random walk.
#[derive(Debug, Clone)]
pub struct RtkState {
    pub time: GpsTime,
    pub pos_ecef: Vector3<f64>,
    pub vel_ecef: Vector3<f64>,
    pub ambiguities: Vec<(DoubleDiffKey, f64)>,
    pub zwd_enabled: bool,
    pub zwd_m: f64,
    pub cov: DMatrix<f64>,
}

impl RtkState {
    /// Create a new RTK state initialized at a seed position.
    pub fn new(initial_pos: Vector3<f64>, initial_time: GpsTime) -> Self {
        let dim = 6;
        let mut cov = DMatrix::zeros(dim, dim);
        for i in 0..3 {
            cov[(i, i)] = 100.0 * 100.0; // 100m initial position uncertainty
        }
        for i in 3..6 {
            cov[(i, i)] = 10.0 * 10.0; // 10 m/s initial velocity uncertainty
        }
        Self {
            time: initial_time,
            pos_ecef: initial_pos,
            vel_ecef: Vector3::zeros(),
            ambiguities: Vec::new(),
            zwd_enabled: false,
            zwd_m: 0.0,
            cov,
        }
    }

    /// Total dimension of the state vector.
    pub fn dim(&self) -> usize {
        6 + self.ambiguities.len() + self.zwd_enabled as usize
    }

    /// Column offset applied to ambiguity indices (1 when ZWD present).
    /// Layout with ZWD: [pos, vel, zwd@6, ambs@7..]; without: legacy layout.
    pub fn amb_offset(&self) -> usize {
        6 + self.zwd_enabled as usize
    }

    /// Column index of the ZWD state, when present.
    pub fn zwd_idx(&self) -> Option<usize> {
        self.zwd_enabled.then_some(6)
    }

    /// Enable the rover ZWD residual state. Must be called before any
    /// ambiguities exist (i.e., at engine construction): the column is
    /// reserved at index 6 and every ambiguity shifts behind it.
    pub fn enable_zwd(&mut self, init_var: f64) {
        if self.zwd_enabled || !self.ambiguities.is_empty() {
            return;
        }
        self.zwd_enabled = true;
        let old_dim = self.cov.nrows();
        let new_dim = old_dim + 1;
        let mut new_cov = DMatrix::zeros(new_dim, new_dim);
        new_cov.view_range_mut(0..old_dim, 0..old_dim).copy_from(&self.cov);
        new_cov[(old_dim, old_dim)] = init_var.max(1e-6);
        self.cov = new_cov;
    }

    /// Find index of ambiguity in the state vector.
    pub fn get_amb_idx(&self, key: &DoubleDiffKey) -> Option<usize> {
        self.ambiguities.iter().position(|(k, _)| k == key)
            .map(|idx| self.amb_offset() + idx)
    }

    /// Pack current state into a flat DVector.
    pub fn to_dvector(&self) -> DVector<f64> {
        let mut vec = DVector::zeros(self.dim());
        vec[0] = self.pos_ecef.x;
        vec[1] = self.pos_ecef.y;
        vec[2] = self.pos_ecef.z;
        vec[3] = self.vel_ecef.x;
        vec[4] = self.vel_ecef.y;
        vec[5] = self.vel_ecef.z;
        let off = self.amb_offset();
        for (i, (_, val)) in self.ambiguities.iter().enumerate() {
            vec[off + i] = *val;
        }
        if let Some(i) = self.zwd_idx() {
            vec[i] = self.zwd_m;
        }
        vec
    }

    /// Unpack a flat DVector into the structured state fields.
    pub fn update_from_dvector(&mut self, vec: &DVector<f64>) {
        self.pos_ecef = Vector3::new(vec[0], vec[1], vec[2]);
        self.vel_ecef = Vector3::new(vec[3], vec[4], vec[5]);
        let off = self.amb_offset();
        for (i, (_, val)) in self.ambiguities.iter_mut().enumerate() {
            *val = vec[off + i];
        }
        if let Some(i) = self.zwd_idx() {
            self.zwd_m = vec[i];
        }
    }

    /// Ensure an ambiguity state exists, adding it with given seed variance if new.
    pub fn ensure_ambiguity(&mut self, key: DoubleDiffKey, initial_val: f64, initial_var: f64) {
        if self.get_amb_idx(&key).is_some() {
            return;
        }
        self.ambiguities.push((key, initial_val));
        let old_dim = self.cov.nrows();
        let new_dim = old_dim + 1;
        let mut new_cov = DMatrix::zeros(new_dim, new_dim);
        new_cov.view_range_mut(0..old_dim, 0..old_dim).copy_from(&self.cov);
        new_cov[(old_dim, old_dim)] = initial_var.max(1.0);
        self.cov = new_cov;
    }

    /// Reset an ambiguity variance (e.g. after detected cycle slip).
    pub fn reset_ambiguity(&mut self, key: &DoubleDiffKey, initial_val: f64, initial_var: f64) {
        if let Some(idx) = self.get_amb_idx(key) {
            let rel = idx - self.amb_offset();
            self.ambiguities[rel].1 = initial_val;
            // Clear cross-covariances for this ambiguity
            for r in 0..self.cov.nrows() {
                self.cov[(r, idx)] = 0.0;
                self.cov[(idx, r)] = 0.0;
            }
            self.cov[(idx, idx)] = initial_var.max(1.0);
        }
    }

    /// Remove stale ambiguities not in the active set.
    pub fn retain_active_ambiguities(&mut self, active_keys: &[DoubleDiffKey]) {
        let mut keep_indices = Vec::with_capacity(6 + active_keys.len());
        for i in 0..self.amb_offset() {
            keep_indices.push(i);
        }

        let mut new_ambs = Vec::new();
        for (i, (key, val)) in self.ambiguities.iter().enumerate() {
            if active_keys.contains(key) {
                keep_indices.push(self.amb_offset() + i);
                new_ambs.push((*key, *val));
            }
        }

        if keep_indices.len() == self.cov.nrows() {
            return;
        }

        let new_dim = keep_indices.len();
        let mut new_cov = DMatrix::zeros(new_dim, new_dim);
        for (new_r, &old_r) in keep_indices.iter().enumerate() {
            for (new_c, &old_c) in keep_indices.iter().enumerate() {
                new_cov[(new_r, new_c)] = self.cov[(old_r, old_c)];
            }
        }
        self.ambiguities = new_ambs;
        self.cov = new_cov;
    }

    /// Extract 3x3 position covariance matrix.
    pub fn extract_pos_cov(&self) -> Matrix3<f64> {
        let mut p = Matrix3::zeros();
        for r in 0..3 {
            for c in 0..3 {
                p[(r, c)] = self.cov[(r, c)];
            }
        }
        p
    }

    /// Extract float ambiguity vector and its covariance submatrix for LAMBDA.
    pub fn extract_amb_block(&self) -> (DVector<f64>, DMatrix<f64>) {
        let off = self.amb_offset();
        let n_amb = self.ambiguities.len();
        let mut a = DVector::zeros(n_amb);
        let mut q = DMatrix::zeros(n_amb, n_amb);
        for i in 0..n_amb {
            a[i] = self.ambiguities[i].1;
            for j in 0..n_amb {
                q[(i, j)] = self.cov[(off + i, off + j)];
            }
        }
        (a, q)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rtk_state_lifecycle() {
        let mut state = RtkState::new(Vector3::new(100.0, 200.0, 300.0), GpsTime::new(2000, 100.0));
        assert_eq!(state.dim(), 6);

        let k1 = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        let k2 = DoubleDiffKey { constellation_id: 0, sat: 3, ref_sat: 1, freq_band: 1 };

        state.ensure_ambiguity(k1, 5.0, 100.0);
        state.ensure_ambiguity(k2, 10.0, 100.0);
        assert_eq!(state.dim(), 8);
        assert_eq!(state.get_amb_idx(&k1), Some(6));
        assert_eq!(state.get_amb_idx(&k2), Some(7));

        let vec = state.to_dvector();
        assert_eq!(vec[6], 5.0);
        assert_eq!(vec[7], 10.0);

        state.retain_active_ambiguities(&[k2]);
        assert_eq!(state.dim(), 7);
        assert_eq!(state.get_amb_idx(&k2), Some(6));
    }
}
