//! Virtual Reference Station (VRS) Network Processing Engine.
//!
//! Synthesizes a virtual reference station at the rover's approximate
//! coordinate by interpolating atmospheric delay surfaces (ZWD and slant
//! ionosphere) across a multi-station CORS network cluster.

use nalgebra::{DMatrix, DVector, Vector3};
use std::collections::HashMap;

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::sat::SatelliteId;
use gneiss_core::time::GpsTime;

/// Station metadata for a network reference base.
#[derive(Debug, Clone)]
pub struct NetworkStation {
    pub id: String,
    pub position_ecef: Vector3<f64>,
    pub epochs: Vec<EpochObs>,
}

/// Atmospheric surface interpolation model across network baselines.
#[derive(Debug, Clone)]
pub struct NetworkAtmosphereSurface {
    pub master_id: String,
    pub master_pos: Vector3<f64>,
    /// Troposphere zenith wet delay gradient [dZWD/dEast, dZWD/dNorth]
    pub tropo_gradient: Vector3<f64>,
    /// Per-satellite ionosphere slant delay gradients [dI/dEast, dI/dNorth]
    pub iono_gradients: HashMap<SatelliteId, Vector3<f64>>,
}

impl NetworkAtmosphereSurface {
    /// Create a new atmosphere surface from master base coordinate.
    pub fn new(master_id: &str, master_pos: Vector3<f64>) -> Self {
        Self {
            master_id: master_id.to_string(),
            master_pos,
            tropo_gradient: Vector3::zeros(),
            iono_gradients: HashMap::new(),
        }
    }

    /// Estimate spatial gradients from secondary station residuals.
    pub fn estimate_gradients(
        &mut self,
        station_offsets_enu: &[Vector3<f64>],
        tropo_residuals: &[f64],
        iono_residuals: &HashMap<SatelliteId, Vec<f64>>,
    ) {
        if station_offsets_enu.len() >= 2 && station_offsets_enu.len() == tropo_residuals.len() {
            self.tropo_gradient = fit_plane_gradient(station_offsets_enu, tropo_residuals);
        }

        for (&sat, residuals) in iono_residuals {
            if station_offsets_enu.len() >= 2 && station_offsets_enu.len() == residuals.len() {
                let grad = fit_plane_gradient(station_offsets_enu, residuals);
                self.iono_gradients.insert(sat, grad);
            }
        }
    }

    /// Interpolate atmospheric corrections at a target coordinate.
    pub fn interpolate_delays(
        &self,
        target_enu_offset: Vector3<f64>,
        sat: SatelliteId,
    ) -> (f64, f64) {
        let delta_tropo = self.tropo_gradient.dot(&target_enu_offset);
        let delta_iono = self.iono_gradients
            .get(&sat)
            .map_or(0.0, |g| g.dot(&target_enu_offset));
        (delta_tropo, delta_iono)
    }
}

/// Fit a first-order linear gradient [d/dEast, d/dNorth, 0] to regional residuals.
fn fit_plane_gradient(offsets: &[Vector3<f64>], residuals: &[f64]) -> Vector3<f64> {
    let n = offsets.len();
    if n < 2 {
        return Vector3::zeros();
    }
    let mut a = DMatrix::<f64>::zeros(n, 2);
    let mut b = DVector::<f64>::zeros(n);
    for i in 0..n {
        a[(i, 0)] = offsets[i].x; // East
        a[(i, 1)] = offsets[i].y; // North
        b[i] = residuals[i];
    }
    match (a.transpose() * &a).try_inverse() {
        Some(ata_inv) => {
            let grad = ata_inv * (a.transpose() * b);
            Vector3::new(grad[0], grad[1], 0.0)
        }
        None => Vector3::zeros(),
    }
}

/// Synthesize a Virtual Reference Station (VRS) epoch at `vrs_position`.
pub fn synthesize_vrs_epoch(
    master_epoch: &EpochObs,
    master_pos: Vector3<f64>,
    vrs_pos: Vector3<f64>,
    surface: &NetworkAtmosphereSurface,
    ephemerides: &[Ephemeris],
) -> EpochObs {
    let delta_ecef = vrs_pos - master_pos;
    let master_llh = gneiss_core::coords::ecef_to_llh(master_pos);
    let r_enu = gneiss_core::coords::ecef_to_ned_matrix(master_llh);
    let delta_enu = r_enu * delta_ecef;

    let mut vrs_sats = Vec::with_capacity(master_epoch.satellites.len());

    for sat_obs in &master_epoch.satellites {
        if let Some(syn_sat) = synthesize_sat_obs(
            sat_obs, master_pos, vrs_pos, delta_enu, surface, ephemerides, master_epoch.time,
        ) {
            vrs_sats.push(syn_sat);
        }
    }

    EpochObs {
        time: master_epoch.time,
        satellites: vrs_sats,
    }
}

fn synthesize_sat_obs(
    master_obs: &SatObs,
    master_pos: Vector3<f64>,
    vrs_pos: Vector3<f64>,
    delta_enu: Vector3<f64>,
    surface: &NetworkAtmosphereSurface,
    ephemerides: &[Ephemeris],
    time: GpsTime,
) -> Option<SatObs> {
    let eph = ephemerides.iter().find(|e| e.sat() == master_obs.sat)?;
    let (sat_pos, _, _, _) = eph.position(time);

    let rho_master = (sat_pos - master_pos).norm();
    let rho_vrs = (sat_pos - vrs_pos).norm();
    let delta_geom = rho_vrs - rho_master;

    let (delta_tropo, delta_iono) = surface.interpolate_delays(delta_enu, master_obs.sat);

    let mut syn = master_obs.clone();
    shift_observables(&mut syn, delta_geom, delta_tropo, delta_iono, master_obs.sat);
    Some(syn)
}

fn shift_observables(
    sat_obs: &mut SatObs,
    delta_geom: f64,
    delta_tropo: f64,
    delta_iono: f64,
    sat: SatelliteId,
) {
    let (f1, f2) = gneiss_core::signal::satellite_frequencies(sat, 0);
    let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    let lambda1 = c / f1;
    let lambda2 = c / f2;
    let gamma = (f1 / f2).powi(2);

    for obs in &mut sat_obs.observations {
        let band = obs.code.signal.freq_band;
        match obs.code.obs_type {
            gneiss_core::obs::ObsType::Pseudorange => {
                let iono = if band == 2 { delta_iono * gamma } else { delta_iono };
                obs.value += delta_geom + delta_tropo + iono;
            }
            gneiss_core::obs::ObsType::CarrierPhase => {
                let (lambda, iono) = if band == 2 { (lambda2, delta_iono * gamma) } else { (lambda1, delta_iono) };
                obs.value += (delta_geom + delta_tropo - iono) / lambda;
            }
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::sat::Constellation;

    #[test]
    fn test_plane_gradient_fit() {
        let offsets = vec![
            Vector3::new(10_000.0, 0.0, 0.0),
            Vector3::new(0.0, 10_000.0, 0.0),
            Vector3::new(10_000.0, 10_000.0, 0.0),
        ];
        let residuals = vec![0.10, 0.05, 0.15];
        let grad = fit_plane_gradient(&offsets, &residuals);
        assert!((grad.x - 1e-5).abs() < 1e-7);
        assert!((grad.y - 5e-6).abs() < 1e-7);
    }

    #[test]
    fn test_vrs_surface_interpolation() {
        let mut surface = NetworkAtmosphereSurface::new(
            "P181",
            Vector3::new(-2694123.0, -4298123.0, 3854123.0),
        );
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let offsets = vec![
            Vector3::new(10_000.0, 0.0, 0.0),
            Vector3::new(0.0, 10_000.0, 0.0),
        ];
        let tropo = vec![0.02, 0.01];
        let mut iono = HashMap::new();
        iono.insert(sat, vec![0.04, 0.02]);
        surface.estimate_gradients(&offsets, &tropo, &iono);

        let target_enu = Vector3::new(5_000.0, 5_000.0, 0.0);
        let (dt, di) = surface.interpolate_delays(target_enu, sat);
        assert!((dt - 0.015).abs() < 1e-4);
        assert!((di - 0.030).abs() < 1e-4);
    }
}
