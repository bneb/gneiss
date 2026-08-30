//! Atmospheric delay and mapping models (Troposphere and Ionosphere).

pub mod ionosphere;
pub mod mapping;
pub mod troposphere;

#[cfg(test)]
mod tests;

use crate::time::GpsTime;
use alloc::vec::Vec;
use nalgebra::Vector3;

/// Ionospheric delay model parameters (Klobuchar).
#[derive(Debug, Clone, Copy)]
pub struct KlobucharParams {
    pub alpha: [f64; 4],
    pub beta: [f64; 4],
}

impl Default for KlobucharParams {
    fn default() -> Self {
        Self {
            alpha: [0.1118E-07, -0.7451E-08, -0.5960E-07, 0.1192E-06],
            beta: [0.1167E+06, -0.2294E+06, -0.1311E+06, 0.1049E+07],
        }
    }
}

/// Tropospheric delay model parameters.
#[derive(Debug, Clone, Copy)]
pub struct TropoParams {
    pub temp_k: f64,
    pub press_hpa: f64,
    pub hum_rel: f64,
}

impl Default for TropoParams {
    fn default() -> Self {
        Self {
            temp_k: 288.15, // 15 C
            press_hpa: 1013.25,
            hum_rel: 0.5,
        }
    }
}

/// Models atmospheric delays for a specific satellite-receiver geometry.
pub struct AtmosphereModel;

impl AtmosphereModel {
    pub fn iono_klobuchar(
        params: &KlobucharParams,
        pos_llh: Vector3<f64>,
        az: f64,
        el: f64,
        time: GpsTime,
    ) -> f64 {
        ionosphere::iono_klobuchar(params, pos_llh, az, el, time)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn iono_ionex(
        tec_maps: &[(GpsTime, &Vec<Vec<f64>>)],
        grid_lat1: f64,
        grid_lat2: f64,
        grid_dlat: f64,
        grid_lon1: f64,
        grid_lon2: f64,
        grid_dlon: f64,
        height_km: f64,
        pos_llh: Vector3<f64>,
        az: f64,
        el: f64,
        time: GpsTime,
    ) -> f64 {
        ionosphere::iono_ionex(
            tec_maps, grid_lat1, grid_lat2, grid_dlat, grid_lon1, grid_lon2, grid_dlon,
            height_km, pos_llh, az, el, time,
        )
    }

    pub fn tropo_saastamoinen(params: &TropoParams, el: f64, height: f64) -> f64 {
        troposphere::tropo_saastamoinen(params, el, height)
    }

    pub fn tropo_rtklib_saastamoinen(params: &TropoParams, pos_llh: Vector3<f64>, el: f64) -> f64 {
        troposphere::tropo_rtklib_saastamoinen(params, pos_llh, el)
    }

    pub fn nmf_mapping_functions(pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> (f64, f64) {
        mapping::nmf_impl(pos_llh, el, time)
    }

    pub fn tropo_nmf(params: &TropoParams, pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> f64 {
        troposphere::tropo_nmf(params, pos_llh, el, time)
    }
}
