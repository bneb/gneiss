//! Virtual Reference Station (VRS) Network Atmospheric Engine.
//!
//! Synthesizes localized Virtual Reference Station observables at the rover's
//! approximate coordinate by interpolating multi-station CORS atmospheric surfaces
//! (tropospheric ZWD and satellite-specific ionospheric pierce points) via 2D
//! Delaunay triangulation with barycentric interpolation.

use std::collections::HashMap;
use nalgebra::{DMatrix, DVector, Vector2, Vector3};

use gneiss_core::coords::{az_el, ecef_delta_to_enu, ecef_to_llh};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::sat::SatelliteId;
use gneiss_core::time::GpsTime;

use crate::post_process::network_adj::{CorsStation, NetworkAdjuster};
use crate::spatial::delaunay::{Delaunay2D, EngineError, Point2D};

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
    pub tropo_gradient: Vector3<f64>,
    pub iono_gradients: HashMap<SatelliteId, Vector3<f64>>,
    pub delaunay_model: Option<DelaunayAtmosphereModel>,
}

/// Spatial 2D Delaunay atmospheric model for CORS regional networks.
#[derive(Debug, Clone)]
pub struct DelaunayAtmosphereModel {
    pub station_mesh: Delaunay2D,
    pub station_zwd: Vec<f64>,
    pub sat_meshes: HashMap<SatelliteId, Delaunay2D>,
    pub sat_iono_delays: HashMap<SatelliteId, Vec<f64>>,
}

impl NetworkAtmosphereSurface {
    /// Create a new atmosphere surface from master base coordinate.
    pub fn new(master_id: &str, master_pos: Vector3<f64>) -> Self {
        Self {
            master_id: master_id.to_string(),
            master_pos,
            tropo_gradient: Vector3::zeros(),
            iono_gradients: HashMap::new(),
            delaunay_model: None,
        }
    }

    /// Estimate planar spatial gradients from secondary station residuals.
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

    /// Attach a 2D Delaunay triangulation model to the surface.
    pub fn set_delaunay_model(&mut self, model: DelaunayAtmosphereModel) {
        self.delaunay_model = Some(model);
    }

    /// Interpolate atmospheric corrections at target coordinate for given satellite.
    pub fn interpolate_delays(
        &self,
        target_enu: Vector3<f64>,
        sat: SatelliteId,
        sat_pos: Option<Vector3<f64>>,
        target_ecef: Option<Vector3<f64>>,
    ) -> (f64, f64) {
        if let Some(ref model) = self.delaunay_model {
            return model.interpolate(target_enu, sat, sat_pos, target_ecef);
        }
        let delta_tropo = self.tropo_gradient.dot(&target_enu);
        let delta_iono = self.iono_gradients
            .get(&sat)
            .map_or(0.0, |g| g.dot(&target_enu));
        (delta_tropo, delta_iono)
    }
}

impl DelaunayAtmosphereModel {
    /// Build Delaunay atmosphere model from station positions, ZWDs, and IPPs.
    pub fn build(
        station_enus: &[Vector2<f64>],
        station_zwds: &[f64],
        sat_ipps: &HashMap<SatelliteId, Vec<(Vector2<f64>, f64)>>,
    ) -> Result<Self, EngineError> {
        let station_mesh = Delaunay2D::new(station_enus)?;
        let mut sat_meshes = HashMap::new();
        let mut sat_iono_delays = HashMap::new();

        for (&sat, ipp_list) in sat_ipps {
            if ipp_list.len() >= 3 {
                let pts: Vec<Vector2<f64>> = ipp_list.iter().map(|(p, _)| *p).collect();
                let vals: Vec<f64> = ipp_list.iter().map(|(_, v)| *v).collect();
                if let Ok(mesh) = Delaunay2D::new(&pts) {
                    sat_meshes.insert(sat, mesh);
                    sat_iono_delays.insert(sat, vals);
                }
            }
        }

        Ok(Self {
            station_mesh,
            station_zwd: station_zwds.to_vec(),
            sat_meshes,
            sat_iono_delays,
        })
    }

    /// Interpolate ZWD and slant iono using Delaunay meshes.
    pub fn interpolate(
        &self,
        target_enu: Vector3<f64>,
        sat: SatelliteId,
        sat_pos: Option<Vector3<f64>>,
        target_ecef: Option<Vector3<f64>>,
    ) -> (f64, f64) {
        let pt2d = Point2D::new(target_enu.x, target_enu.y);
        let tropo_total = self.station_mesh.interpolate(&self.station_zwd, pt2d);
        let master_zwd = self.station_zwd.first().copied().unwrap_or(0.15);
        let delta_zwd = tropo_total - master_zwd;

        let delta_tropo = if let (Some(s_pos), Some(t_pos)) = (sat_pos, target_ecef) {
            let t_llh = ecef_to_llh(t_pos);
            let (_, el) = az_el(t_llh, t_pos, s_pos);
            delta_zwd / el.max(0.1).sin()
        } else {
            delta_zwd
        };

        let iono = self.interpolate_iono(sat, sat_pos, target_ecef);
        (delta_tropo, iono)
    }

    fn interpolate_iono(
        &self,
        sat: SatelliteId,
        sat_pos: Option<Vector3<f64>>,
        target_ecef: Option<Vector3<f64>>,
    ) -> f64 {
        let (Some(s_pos), Some(t_pos)) = (sat_pos, target_ecef) else { return 0.0 };
        let (Some(mesh), Some(delays)) = (self.sat_meshes.get(&sat), self.sat_iono_delays.get(&sat)) else {
            return 0.0;
        };
        let t_llh = ecef_to_llh(t_pos);
        let (az, el) = az_el(t_llh, t_pos, s_pos);
        let (lat_ipp, lon_ipp) = compute_ipp(t_llh.x, t_llh.y, az, el);
        mesh.interpolate(delays, Point2D::new(lon_ipp, lat_ipp))
    }
}

/// Compute Ionospheric Pierce Point (IPP) for single-layer shell model at H = 350 km.
pub fn compute_ipp(rec_lat_deg: f64, rec_lon_deg: f64, az_rad: f64, el_rad: f64) -> (f64, f64) {
    let re = 6371.0;
    let h = 350.0;
    let el = el_rad.max(0.05);
    let ratio = (re / (re + h) * libm::cos(el)).min(1.0);
    let psi = (libm::asin(ratio) - (std::f64::consts::PI * 0.5 - el)).max(0.0);
    let lat_r = rec_lat_deg.to_radians();
    let sin_ipp = lat_r.sin() * libm::cos(psi) + lat_r.cos() * libm::sin(psi) * libm::cos(az_rad);
    let ipp_lat = libm::asin(sin_ipp.clamp(-1.0, 1.0)).to_degrees();
    let cos_ipp = ipp_lat.to_radians().cos().max(1e-6);
    let sin_dlon = libm::sin(psi) * libm::sin(az_rad) / cos_ipp;
    let ipp_lon = rec_lon_deg + libm::asin(sin_dlon.clamp(-1.0, 1.0)).to_degrees();
    (ipp_lat, ipp_lon)
}

/// VRS Synthesizer for localized virtual reference stations.
pub struct VrsSynthesizer {
    pub master_id: String,
    pub master_pos: Vector3<f64>,
}

impl VrsSynthesizer {
    pub fn new(master_id: &str, master_pos: Vector3<f64>) -> Self {
        Self { master_id: master_id.to_string(), master_pos }
    }

    /// Synthesize localized VRS observation stream across multiple epochs.
    pub fn synthesize_vrs_stream(
        &self,
        stations: &[CorsStation],
        rover_approx_ecef: Vector3<f64>,
        ephemerides: &[Ephemeris],
    ) -> Result<Vec<EpochObs>, EngineError> {
        let master = stations.iter().find(|s| s.id == self.master_id)
            .ok_or_else(|| EngineError::DegenerateMesh("Master not found".into()))?;
        let adjuster = NetworkAdjuster::new(stations, &self.master_id)?;
        let num_epochs = master.epochs.len();
        let mut vrs_epochs = Vec::with_capacity(num_epochs);

        for epoch_idx in 0..num_epochs {
            let adj_res = adjuster.adjust_epoch(stations, epoch_idx, ephemerides).ok();
            let surface = build_epoch_surface(&adj_res, stations, &self.master_id, self.master_pos, epoch_idx, ephemerides);
            let m_epoch = &master.epochs[epoch_idx];
            let vrs_epoch = synthesize_vrs_epoch(m_epoch, self.master_pos, rover_approx_ecef, &surface, ephemerides);
            vrs_epochs.push(vrs_epoch);
        }
        Ok(vrs_epochs)
    }
}

/// Construct atmosphere surface with Delaunay model from adjustment result.
fn build_epoch_surface(
    adj: &Option<crate::post_process::network_adj::NetworkAdjustmentResult>,
    stations: &[CorsStation],
    master_id: &str,
    master_pos: Vector3<f64>,
    epoch_idx: usize,
    eph: &[Ephemeris],
) -> NetworkAtmosphereSurface {
    let mut surface = NetworkAtmosphereSurface::new(master_id, master_pos);
    let master_llh = ecef_to_llh(master_pos);

    let mut st_enus = Vec::new();
    let mut st_zwds = Vec::new();
    let mut sat_ipps: HashMap<SatelliteId, Vec<(Vector2<f64>, f64)>> = HashMap::new();

    for s in stations {
        let enu = ecef_delta_to_enu(s.pos_ecef, master_pos, master_llh);
        st_enus.push(Vector2::new(enu.x, enu.y));
        let zwd = adj.as_ref().and_then(|a| a.station_atmospheres.get(&s.id)).map_or(0.15, |at| at.zwd_m);
        st_zwds.push(zwd);

        if let Some(epoch) = s.epochs.get(epoch_idx) {
            collect_station_ipps(&mut sat_ipps, s, epoch, eph, adj);
        }
    }

    if let Ok(model) = DelaunayAtmosphereModel::build(&st_enus, &st_zwds, &sat_ipps) {
        surface.set_delaunay_model(model);
    }
    surface
}

/// Collect IPP coordinates and iono delays for all satellites tracked by station.
fn collect_station_ipps(
    sat_ipps: &mut HashMap<SatelliteId, Vec<(Vector2<f64>, f64)>>,
    st: &CorsStation,
    epoch: &EpochObs,
    eph: &[Ephemeris],
    adj: &Option<crate::post_process::network_adj::NetworkAdjustmentResult>,
) {
    let s_llh = ecef_to_llh(st.pos_ecef);
    for sat_obs in &epoch.satellites {
        if let Some(e) = eph.iter().find(|e| e.sat() == sat_obs.sat) {
            let (sat_p, _, _, _) = e.position(epoch.time);
            let (az, el) = az_el(s_llh, st.pos_ecef, sat_p);
            let (lat_ipp, lon_ipp) = compute_ipp(s_llh.x, s_llh.y, az, el);
            let delay = adj.as_ref()
                .and_then(|a| a.station_atmospheres.get(&st.id))
                .and_then(|at| at.iono_slant_m.get(&sat_obs.sat).copied())
                .unwrap_or(0.0);
            sat_ipps.entry(sat_obs.sat).or_default().push((Vector2::new(lon_ipp, lat_ipp), delay));
        }
    }
}

/// Fit a first-order linear gradient [d/dEast, d/dNorth, 0] to regional residuals.
fn fit_plane_gradient(offsets: &[Vector3<f64>], residuals: &[f64]) -> Vector3<f64> {
    let n = offsets.len();
    if n < 2 { return Vector3::zeros(); }
    let mut a = DMatrix::<f64>::zeros(n, 2);
    let mut b = DVector::<f64>::zeros(n);
    for i in 0..n {
        a[(i, 0)] = offsets[i].x;
        a[(i, 1)] = offsets[i].y;
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
    let master_llh = ecef_to_llh(master_pos);
    let delta_enu = ecef_delta_to_enu(vrs_pos, master_pos, master_llh);
    let mut vrs_sats = Vec::with_capacity(master_epoch.satellites.len());

    for sat_obs in &master_epoch.satellites {
        if let Some(syn_sat) = synthesize_sat_obs(
            sat_obs, master_pos, vrs_pos, delta_enu, surface, ephemerides, master_epoch.time,
        ) {
            vrs_sats.push(syn_sat);
        }
    }
    EpochObs { time: master_epoch.time, satellites: vrs_sats }
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
    let sat_pos_m = compute_sat_pos_sagnac(eph, time, master_pos);
    let sat_pos_v = compute_sat_pos_sagnac(eph, time, vrs_pos);
    let rho_master = (sat_pos_m - master_pos).norm();
    let rho_vrs = (sat_pos_v - vrs_pos).norm();
    let delta_geom = rho_vrs - rho_master;

    let (delta_tropo, delta_iono) = surface.interpolate_delays(
        delta_enu, master_obs.sat, Some(sat_pos_v), Some(vrs_pos),
    );

    let mut syn = master_obs.clone();
    shift_observables(&mut syn, delta_geom, delta_tropo, delta_iono, master_obs.sat);
    Some(syn)
}

fn compute_sat_pos_sagnac(eph: &Ephemeris, t_rx: GpsTime, rx_pos: Vector3<f64>) -> Vector3<f64> {
    use gneiss_core::constants::{EARTH_ROTATION_RATE_RAD_S, SPEED_OF_LIGHT_M_S};
    let (p0, _, clk0, _) = eph.position(t_rx);
    let tau0 = (rx_pos - p0).norm() / SPEED_OF_LIGHT_M_S;
    let t_tx = t_rx - tau0 - clk0;
    let (ptx, _, _, _) = eph.position(t_tx);
    let wt = EARTH_ROTATION_RATE_RAD_S * tau0;
    let (sw, cw) = libm::sincos(wt);
    Vector3::new(
        ptx.x * cw + ptx.y * sw,
        -ptx.x * sw + ptx.y * cw,
        ptx.z,
    )
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
    fn test_compute_ipp() {
        let (lat, lon) = compute_ipp(37.5, -122.2, 0.5, 0.7);
        assert!((lat - 37.5).abs() < 5.0);
        assert!((lon - -122.2).abs() < 5.0);
    }

    #[test]
    fn test_vrs_delaunay_surface_interpolation() {
        let enus = vec![
            Vector2::new(0.0, 0.0),
            Vector2::new(20_000.0, 0.0),
            Vector2::new(0.0, 20_000.0),
        ];
        let zwds = vec![0.12, 0.16, 0.14];
        let ipps = HashMap::new();
        let model = DelaunayAtmosphereModel::build(&enus, &zwds, &ipps).expect("build model");
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let (tropo, _) = model.interpolate(Vector3::new(5_000.0, 5_000.0, 0.0), sat, None, None);
        assert!((tropo - 0.015).abs() < 1e-3);
    }
}
