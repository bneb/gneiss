//! Time-Differenced Carrier Phase (TDCP) delta-position solver.
//!
//! TDCP differences DD carrier phase between consecutive epochs to cancel
//! integer ambiguities, yielding mm-level precise position changes that are
//! immune to code multipath. Combined with long-arc PR for absolute anchor,
//! this breaks the RTK code-multipath deadlock.
//!
//! # Theory
//!
//! At epoch k, the DD carrier phase observation for satellite pair (s, r) is:
//!
//! ```text
//!   DD_CP_k = geom_DD(r_k) + λ·N + ε_k
//! ```
//!
//! Time-differencing between epochs k and k-1 cancels the ambiguity N:
//!
//! ```text
//!   ΔDD_CP = DD_CP_k - DD_CP_{k-1}
//!          = geom_DD(r_k) - geom_DD(r_{k-1}) + Δε
//! ```
//!
//! Linearizing about the predicted position r_k^pred and using the previous
//! estimate r_{k-1}^est gives a measurement of the position correction δr_k:
//!
//! ```text
//!   z_s ≈ (e_{s,k} - e_{r,k})ᵀ · δr_k
//! ```
//!
//! Stacking all valid satellite pairs yields a linear system solved via WLS.

use nalgebra::{DMatrix, DVector, Vector3};
use std::collections::VecDeque;

use crate::filter::DdObservation;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

// ---------------------------------------------------------------------------
// Stored previous-epoch carrier phase data
// ---------------------------------------------------------------------------

/// DD carrier phase data for one satellite pair from the previous epoch.
#[derive(Clone, Debug)]
struct PrevCpDd {
    /// Non-reference satellite ID.
    sat: SatelliteId,
    /// Reference satellite ID used for DD at the previous epoch.
    ref_sat: SatelliteId,
    /// DD CP L1 value in meters.
    dd_cp_l1_m: f64,
    /// Geometric DD from rover at previous epoch (meters).
    geom_dd_rov: f64,
    /// Geometric DD from base at previous epoch (meters).
    geom_dd_base: f64,
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Compute satellite position with approximate range for signal travel time.
fn compute_sat_position(
    eph: &Ephemeris,
    rx_pos: Vector3<f64>,
    time: GpsTime,
) -> (Vector3<f64>, Vector3<f64>) {
    // Use the approximate range from rover to satellite for signal travel time.
    // First call with pr=0 gives rough position; compute range; call again.
    let (sat_approx, _) =
        crate::engine::measurement_math::get_sat_state(eph, 0.0, 0.0, time, rx_pos);
    let range_approx = (sat_approx - rx_pos).norm();
    crate::engine::measurement_math::get_sat_state(eph, range_approx, 0.0, time, rx_pos)
}

/// Get L1 wavelength for a satellite (handles GLONASS FDMA via freq_num).
fn l1_wavelength(sat: SatelliteId, freq_num: i8) -> f64 {
    let (f1, _f2) = gneiss_core::signal::satellite_frequencies(sat, freq_num);
    gneiss_core::constants::SPEED_OF_LIGHT_M_S / f1
}

// ---------------------------------------------------------------------------
// TDCP Delta-Position Solver
// ---------------------------------------------------------------------------

/// TDCP delta-position solver.
///
/// Computes precise position change between consecutive epochs using
/// time-differenced DD carrier phase. Ambiguities cancel in the difference,
/// yielding mm-level delta positions without requiring AR.
pub struct TdcpSolver {
    /// Previous-epoch DD CP data.
    prev_data: Option<Vec<PrevCpDd>>,
    /// Previous-epoch rover position (ECEF).
    prev_rover_pos: Vector3<f64>,
    /// Whether the solver has valid previous-epoch data.
    has_prev: bool,
    /// Minimum number of TDCP measurements required for a solution.
    min_measurements: usize,
}

impl TdcpSolver {
    /// Create a new TDCP solver.
    pub fn new() -> Self {
        Self {
            prev_data: None,
            prev_rover_pos: Vector3::zeros(),
            has_prev: false,
            min_measurements: 4,
        }
    }

    /// Store current-epoch data for time-differencing at the next epoch.
    pub fn store_epoch(
        &mut self,
        matched_obs: &[(DdObservation, DdObservation)],
        ref_sats: &[(Constellation, SatelliteId)],
        rover_pos: Vector3<f64>,
        base_pos: Vector3<f64>,
        ephemerides: &[Ephemeris],
        time: GpsTime,
    ) {
        let mut data = Vec::with_capacity(matched_obs.len());

        for (rov_obs, base_obs) in matched_obs {
            let constellation = rov_obs.sat.constellation;
            let ref_sat = match ref_sats.iter().find(|(c, _)| *c == constellation) {
                Some((_, rs)) => *rs,
                None => continue,
            };
            if rov_obs.sat == ref_sat {
                continue;
            }

            // Find reference satellite observations
            let ref_rov = match matched_obs.iter().find(|(r, _)| r.sat == ref_sat) {
                Some((r, _)) => r,
                None => continue,
            };
            let ref_base = match matched_obs.iter().find(|(_, b)| b.sat == ref_sat) {
                Some((_, b)) => b,
                None => continue,
            };

            // Require valid CP L1 on all four observations
            let (Some(rov_cp), Some(ref_rov_cp)) = (rov_obs.cp_l1, ref_rov.cp_l1) else {
                continue;
            };
            let (Some(base_cp), Some(ref_base_cp)) = (base_obs.cp_l1, ref_base.cp_l1) else {
                continue;
            };

            // Find ephemerides
            let eph_sat = match crate::engine::measurement::geometry::find_ephemeris(
                ephemerides,
                rov_obs.sat,
                time.tow,
            ) {
                Some(e) => e,
                None => continue,
            };
            let eph_ref = match crate::engine::measurement::geometry::find_ephemeris(
                ephemerides,
                ref_sat,
                time.tow,
            ) {
                Some(e) => e,
                None => continue,
            };

            // Wavelengths (handles GLONASS FDMA)
            let lam_sat = l1_wavelength(rov_obs.sat, eph_sat.freq_num());
            let lam_ref = l1_wavelength(ref_sat, eph_ref.freq_num());

            // DD CP L1 in meters
            let dd_cp_l1_m = (rov_cp * lam_sat - ref_rov_cp * lam_ref)
                - (base_cp * lam_sat - ref_base_cp * lam_ref);

            // Compute satellite positions from common reference (base pos)
            // so compute_delta sees consistent satellite geometry.
            let (sat_pos, _) = compute_sat_position(eph_sat, base_pos, time);
            let (ref_sat_pos, _) = compute_sat_position(eph_ref, base_pos, time);

            // Geometric DD at rover position using common satellite positions
            let geom_dd_rov =
                (rover_pos - sat_pos).norm() - (rover_pos - ref_sat_pos).norm();

            // Geometric DD at base position using same satellite positions
            let geom_dd_base =
                (base_pos - sat_pos).norm() - (base_pos - ref_sat_pos).norm();

            data.push(PrevCpDd {
                sat: rov_obs.sat,
                ref_sat,
                dd_cp_l1_m,
                geom_dd_rov,
                geom_dd_base,
            });
        }

        self.prev_data = Some(data);
        self.prev_rover_pos = rover_pos;
        self.has_prev = true;
    }

    /// Compute TDCP delta position between the previous and current epoch.
    ///
    /// Returns `(delta_pos_ecef, covariance_3x3)` where `delta_pos` is the
    /// estimated position change r_k - r_{k-1} in ECEF meters, or `None` if
    /// insufficient valid TDCP measurements are available.
    pub fn compute_delta(
        &self,
        matched_obs_curr: &[(DdObservation, DdObservation)],
        ref_sats_curr: &[(Constellation, SatelliteId)],
        rover_pos_curr: Vector3<f64>,
        base_pos: Vector3<f64>,
        ephemerides: &[Ephemeris],
        time_curr: GpsTime,
    ) -> Option<(Vector3<f64>, DMatrix<f64>)> {
        let prev_data = self.prev_data.as_ref()?;
        if !self.has_prev {
            return None;
        }

        // CP measurement noise (meters): σ ≈ 5mm for short-baseline L1.
        // Time-differencing doubles the variance.
        let sigma_cp = 0.005_f64;
        let var_tdcp = 2.0 * sigma_cp * sigma_cp;

        let mut h_rows: Vec<Vector3<f64>> = Vec::with_capacity(prev_data.len());
        let mut z_vals: Vec<f64> = Vec::with_capacity(prev_data.len());

        for prev in prev_data {
            // Find current-epoch observations for the same satellite
            let rov_curr = match matched_obs_curr.iter().find(|(r, _)| r.sat == prev.sat) {
                Some((r, _)) => r,
                None => continue,
            };
            let base_curr = match matched_obs_curr.iter().find(|(_, b)| b.sat == prev.sat) {
                Some((_, b)) => b,
                None => continue,
            };

            // Must use the same reference satellite as previous epoch
            let ref_sat_curr = match ref_sats_curr
                .iter()
                .find(|(c, _)| *c == prev.sat.constellation)
            {
                Some((_, rs)) => *rs,
                None => continue,
            };
            if ref_sat_curr != prev.ref_sat {
                continue;
            }

            let ref_rov_curr =
                match matched_obs_curr.iter().find(|(r, _)| r.sat == ref_sat_curr) {
                    Some((r, _)) => r,
                    None => continue,
                };
            let ref_base_curr =
                match matched_obs_curr.iter().find(|(_, b)| b.sat == ref_sat_curr) {
                    Some((_, b)) => b,
                    None => continue,
                };

            // Require valid CP L1 on all four current observations
            let (Some(cp_sat), Some(cp_ref)) = (rov_curr.cp_l1, ref_rov_curr.cp_l1) else {
                continue;
            };
            let (Some(cp_base_sat), Some(cp_base_ref)) = (base_curr.cp_l1, ref_base_curr.cp_l1)
            else {
                continue;
            };

            // Find ephemerides for current epoch
            let eph_sat = match crate::engine::measurement::geometry::find_ephemeris(
                ephemerides,
                prev.sat,
                time_curr.tow,
            ) {
                Some(e) => e,
                None => continue,
            };
            let eph_ref = match crate::engine::measurement::geometry::find_ephemeris(
                ephemerides,
                prev.ref_sat,
                time_curr.tow,
            ) {
                Some(e) => e,
                None => continue,
            };

            // Wavelengths
            let lam_sat = l1_wavelength(prev.sat, eph_sat.freq_num());
            let lam_ref = l1_wavelength(prev.ref_sat, eph_ref.freq_num());

            // Current-epoch DD CP L1 in meters
            let dd_cp_curr = (cp_sat * lam_sat - cp_ref * lam_ref)
                - (cp_base_sat * lam_sat - cp_base_ref * lam_ref);

            // Time-differenced DD CP
            let delta_dd_cp = dd_cp_curr - prev.dd_cp_l1_m;

            // Compute satellite positions from a COMMON reference point
            // (the base station) so that satellite motion cancels exactly
            // between rover and base geometric DDs.  Using different ref
            // points for rover vs base introduces a ~5cm satellite position
            // difference (from signal travel time), which, when differenced
            // across epochs with ~160m of satellite motion, produces ~10m
            // RMS z-variation.
            //
            // The rover range is still computed from the rover position;
            // the ~5cm satellite position error from using the base ref
            // point is constant across epochs and cancels in the time
            // difference.
            let (sat_pos_curr, _) = compute_sat_position(eph_sat, base_pos, time_curr);
            let (ref_sat_pos_curr, _) = compute_sat_position(eph_ref, base_pos, time_curr);

            // Geometric DD at current predicted position using common sat pos
            let geom_dd_rov_curr = (rover_pos_curr - sat_pos_curr).norm()
                - (rover_pos_curr - ref_sat_pos_curr).norm();

            // Base geometric DD at current epoch using same sat pos
            let geom_dd_base_curr = (base_pos - sat_pos_curr).norm()
                - (base_pos - ref_sat_pos_curr).norm();

            // TDCP innovation with common satellite positions:
            //
            // ΔDD_CP has two components: satellite motion + rover motion.
            // Δgeom_base captures the satellite motion (base is stationary).
            // Δgeom_rov captures satellite motion + predicted rover motion.
            //
            // With common sat pos (both computed from base_pos), the satellite
            // motion terms are identical in all three deltas, so:
            //   (ΔDD_CP - Δgeom_base)     = CP-measured rover motion change
            //   (Δgeom_rov - Δgeom_base)  = predicted rover motion change
            //   z = (ΔDD_CP - Δgeom_base) - (Δgeom_rov - Δgeom_base)
            //     = ΔDD_CP - Δgeom_rov
            //
            // This cancels satellite motion (~160m/epoch) leaving only the
            // rover motion residual (~0.1-1m).
            let delta_geom_rov = geom_dd_rov_curr - prev.geom_dd_rov;
            let delta_geom_base = geom_dd_base_curr - prev.geom_dd_base;

            // Satellite-motion-free innovation
            let z = delta_dd_cp - delta_geom_rov;

            // Debug: verify satellite motion cancellation

            // Reject cycle-slipped pairs: genuine rover motion produces
            // |z| < 1.0m per epoch (0.2s × 30m/s × 0.1 LOS ≈ 0.6m max).
            // Larger values indicate a cycle slip or reference satellite
            // change that broke the CP continuity assumption.
            if z.abs() > 5.0 {
                continue;
            }

            // Jacobian: H = (e_sat - e_ref)ᵀ at current epoch
            let los_sat_curr = (sat_pos_curr - rover_pos_curr).normalize();
            let los_ref_curr = (ref_sat_pos_curr - rover_pos_curr).normalize();
            let h_row = los_sat_curr - los_ref_curr;

            h_rows.push(h_row);
            z_vals.push(z);
        }

        if h_rows.len() < self.min_measurements {
            return None;
        }

        let n = h_rows.len();
        let mut h_mat = DMatrix::zeros(n, 3);
        let z_vec = DVector::from_vec(z_vals);

        for i in 0..n {
            h_mat[(i, 0)] = h_rows[i].x;
            h_mat[(i, 1)] = h_rows[i].y;
            h_mat[(i, 2)] = h_rows[i].z;
        }

        // Solve: δr = (HᵀH)⁻¹ Hᵀz
        let ht_h = &h_mat.transpose() * &h_mat;
        let ht_z = &h_mat.transpose() * &z_vec;

        // Check conditioning
        let ht_h_eig = ht_h.symmetric_eigenvalues();
        if ht_h_eig[0] < 1e-6 {
            return None;
        }

        let ht_h_inv = ht_h.try_inverse()?;
        let delta_dv = &ht_h_inv * ht_z;
        let delta_pos = Vector3::new(delta_dv[0], delta_dv[1], delta_dv[2]);

        // Covariance: (HᵀH)⁻¹ · var_tdcp
        let cov = ht_h_inv * var_tdcp;

        // Sanity check: with the corrected formulation (common sat ref
        // + outlier rejection), deltas should be ~0.1-2m per epoch.
        // Anything above 5m indicates uncorrected cycle slips.
        let delta_norm = delta_pos.norm();
        if delta_norm > 5.0 {
            let z_rms = (z_vec.norm() / (n as f64).sqrt()).sqrt();
            tracing::debug!(
                "TDCP: rejected delta {:.1}m ({} pairs, z_rms={:.1}m)",
                delta_norm, n, z_rms
            );
            return None;
        }

        Some((delta_pos, cov))
    }

    /// Returns true if the solver has valid previous-epoch data.
    pub fn has_prev_epoch(&self) -> bool {
        self.has_prev
    }
}

impl Default for TdcpSolver {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// TDCP Trajectory Accumulator
// ---------------------------------------------------------------------------

/// Sliding-window TDCP trajectory accumulator.
///
/// Integrates TDCP delta positions over a configurable window to produce a
/// smooth relative trajectory. Combined with a long-arc PR absolute anchor,
/// this provides precise absolute positioning independent of instantaneous
/// code multipath.
pub struct TdcpTrajectory {
    deltas: VecDeque<(Vector3<f64>, f64)>,
    cumulative: Vector3<f64>,
    window_size: usize,
}

impl TdcpTrajectory {
    /// Create a new trajectory accumulator with the given window size.
    pub fn new(window_size: usize) -> Self {
        Self {
            deltas: VecDeque::with_capacity(window_size.max(1)),
            cumulative: Vector3::zeros(),
            window_size: window_size.max(1),
        }
    }

    /// Push a new TDCP delta position.
    ///
    /// If the window is full, the oldest delta is removed and subtracted from
    /// the cumulative sum.
    pub fn push(&mut self, delta: Vector3<f64>, cov: &DMatrix<f64>) {
        let trace = cov[(0, 0)] + cov[(1, 1)] + cov[(2, 2)];

        if self.deltas.len() >= self.window_size {
            if let Some((old_delta, _)) = self.deltas.pop_front() {
                self.cumulative -= old_delta;
            }
        }

        self.cumulative += delta;
        self.deltas.push_back((delta, trace));
    }

    /// Current cumulative position change from the anchor epoch.
    pub fn cumulative_delta(&self) -> Vector3<f64> {
        self.cumulative
    }

    /// Smoothed delta position and estimated 3D accuracy over the window.
    ///
    /// Returns `(smoothed_delta, sigma_3d)` or `None` if the window is empty.
    pub fn smoothed_delta(&self) -> Option<(Vector3<f64>, f64)> {
        if self.deltas.is_empty() {
            return None;
        }
        let total_trace: f64 = self.deltas.iter().map(|(_, t)| t).sum();
        let sigma_3d = total_trace.sqrt();
        Some((self.cumulative, sigma_3d))
    }

    /// Number of epochs currently in the window.
    pub fn len(&self) -> usize {
        self.deltas.len()
    }

    /// Whether the window is empty.
    pub fn is_empty(&self) -> bool {
        self.deltas.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solver_new_has_no_prev() {
        let solver = TdcpSolver::new();
        assert!(!solver.has_prev_epoch());
    }

    #[test]
    fn test_trajectory_empty() {
        let traj = TdcpTrajectory::new(60);
        assert!(traj.is_empty());
        assert_eq!(traj.len(), 0);
        assert!(traj.smoothed_delta().is_none());
        assert_eq!(traj.cumulative_delta(), Vector3::zeros());
    }

    #[test]
    fn test_trajectory_single_delta() {
        let mut traj = TdcpTrajectory::new(60);
        let delta = Vector3::new(1.0, 0.5, -0.2);
        let cov = DMatrix::from_diagonal(&DVector::from_vec(vec![0.01, 0.01, 0.01]));
        traj.push(delta, &cov);

        assert_eq!(traj.len(), 1);
        let cum = traj.cumulative_delta();
        assert!((cum.x - 1.0).abs() < 1e-10);
        assert!((cum.y - 0.5).abs() < 1e-10);
        assert!((cum.z + 0.2).abs() < 1e-10);

        let (smoothed, sigma) = traj.smoothed_delta().unwrap();
        assert_eq!(smoothed.x, 1.0);
        assert!((sigma - 0.03_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_trajectory_window_rolls() {
        let mut traj = TdcpTrajectory::new(3);
        let cov = DMatrix::from_diagonal(&DVector::from_vec(vec![0.01, 0.01, 0.01]));

        traj.push(Vector3::new(1.0, 0.0, 0.0), &cov);
        traj.push(Vector3::new(2.0, 0.0, 0.0), &cov);
        traj.push(Vector3::new(3.0, 0.0, 0.0), &cov);
        assert_eq!(traj.len(), 3);
        assert!((traj.cumulative_delta().x - 6.0).abs() < 1e-10);

        traj.push(Vector3::new(4.0, 0.0, 0.0), &cov);
        assert_eq!(traj.len(), 3);
        assert!((traj.cumulative_delta().x - 9.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_delta_no_prev_returns_none() {
        let solver = TdcpSolver::new();
        let result = solver.compute_delta(
            &[],
            &[],
            Vector3::zeros(),
            Vector3::zeros(),
            &[],
            GpsTime::new(0, 0.0),
        );
        assert!(result.is_none());
    }

    #[test]
    fn test_solver_default() {
        let solver = TdcpSolver::default();
        assert!(!solver.has_prev_epoch());
    }

    #[test]
    fn test_trajectory_window_size_one() {
        let mut traj = TdcpTrajectory::new(1);
        let cov = DMatrix::from_diagonal(&DVector::from_vec(vec![0.01, 0.01, 0.01]));

        traj.push(Vector3::new(1.0, 2.0, 3.0), &cov);
        assert_eq!(traj.len(), 1);

        traj.push(Vector3::new(4.0, 5.0, 6.0), &cov);
        assert_eq!(traj.len(), 1);
        assert!((traj.cumulative_delta().x - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_l1_wavelength_gps() {
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let lam = l1_wavelength(sat, 0);
        // GPS L1: 1575.42 MHz → λ ≈ 0.1903 m
        assert!((lam - 0.1903).abs() < 0.001);
    }
}
