//! Tropospheric delay models (Saastamoinen, RTKLIB, NMF).

use crate::time::GpsTime;
use nalgebra::Vector3;
use super::mapping::nmf_impl;
use super::TropoParams;

/// Computes Tropospheric delay in meters using the Saastamoinen model.
pub fn tropo_saastamoinen(params: &TropoParams, el: f64, height: f64) -> f64 {
    let z = core::f64::consts::FRAC_PI_2 - el;
    let p = params.press_hpa * libm::pow(1.0 - 0.000022557 * height, 5.2568);
    let t = params.temp_k - 0.0065 * height;
    let e = 6.108 * libm::exp((17.15 * t - 4684.0) / (t - 38.45)) * params.hum_rel;

    0.002277 / libm::cos(z) * (p + (1255.0 / t + 0.05) * e - libm::tan(z) * libm::tan(z))
}

/// Computes Tropospheric delay in meters using the EXACT RTKLIB Saastamoinen model.
pub fn tropo_rtklib_saastamoinen(params: &TropoParams, pos_llh: Vector3<f64>, el: f64) -> f64 {
    if pos_llh.z < -100.0 || pos_llh.z > 10000.0 || el <= 0.0 {
        return 0.0;
    }
    let z = core::f64::consts::FRAC_PI_2 - el;
    let height = pos_llh.z;
    let p = params.press_hpa * libm::pow(1.0 - 0.000022557 * height, 5.2568);
    let t = params.temp_k - 0.0065 * height;
    let e = 6.108 * libm::exp((17.15 * t - 4684.0) / (t - 38.45)) * params.hum_rel;

    let trph = 0.0022768 * p
        / (1.0 - 0.00266 * libm::cos(2.0 * pos_llh.x) - 0.00028 * height / 1000.0);
    let trpw = 0.002277 * (1255.0 / t + 0.05) * e;

    (trph + trpw) / libm::cos(z)
}

/// Computes Tropospheric delay in meters using Saastamoinen zenith delay mapped with NMF.
pub fn tropo_nmf(params: &TropoParams, pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> f64 {
    if el <= 0.0 {
        return 0.0;
    }

    let hgt = pos_llh.z;
    let (m_h, m_w) = nmf_impl(pos_llh, el, time);

    let z_dry = 0.0022768 * params.press_hpa
        / (1.0 - 0.00266 * libm::cos(2.0 * pos_llh.x) - 0.00028 * hgt / 1000.0);

    let e = 6.108
        * libm::exp((17.15 * params.temp_k - 4684.0) / (params.temp_k - 38.45))
        * params.hum_rel;
    let z_wet = 0.002277 * (1255.0 / params.temp_k + 0.05) * e;

    z_dry * m_h + z_wet * m_w
}
