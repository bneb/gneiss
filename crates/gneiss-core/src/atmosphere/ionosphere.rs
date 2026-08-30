//! Ionospheric delay models (Klobuchar, IONEX grid interpolation).

use crate::time::GpsTime;
use alloc::vec::Vec;
use nalgebra::Vector3;
use super::KlobucharParams;

/// Computes Ionospheric delay in meters using the Klobuchar model.
pub fn iono_klobuchar(
    params: &KlobucharParams,
    pos_llh: Vector3<f64>,
    az: f64,
    el: f64,
    time: GpsTime,
) -> f64 {
    let pi = core::f64::consts::PI;
    let f = 1.0 + 16.0 * libm::pow(0.53 - el / pi, 3.0);
    let psi = 0.0137 / (el / pi + 0.11) - 0.022;
    let phi_u = pos_llh.x / pi;
    let phi_i = (phi_u + psi * libm::cos(az)).clamp(-0.416, 0.416);
    let lambda_u = pos_llh.y / pi;
    let lambda_i = lambda_u + psi * libm::sin(az) / libm::cos(phi_i * pi);
    let phi_m = phi_i + 0.064 * libm::cos(lambda_i * pi - 1.617);

    let mut t = 43200.0 * lambda_i + time.tow;
    t %= 86400.0;
    if t < 0.0 {
        t += 86400.0;
    }

    let mut a = params.alpha[0]
        + params.alpha[1] * phi_m
        + params.alpha[2] * phi_m * phi_m
        + params.alpha[3] * phi_m * phi_m * phi_m;
    if a < 0.0 {
        a = 0.0;
    }

    let mut p = params.beta[0]
        + params.beta[1] * phi_m
        + params.beta[2] * phi_m * phi_m
        + params.beta[3] * phi_m * phi_m * phi_m;
    if p < 72000.0 {
        p = 72000.0;
    }

    let x = 2.0 * pi * (t - 50400.0) / p;

    let delay = if libm::fabs(x) < 1.57 {
        5e-9 + a * (1.0 - x * x / 2.0 + x * x * x * x / 24.0)
    } else {
        5e-9
    };

    delay * f * crate::constants::SPEED_OF_LIGHT_M_S
}

fn bilinear_interpolate_tec(
    map: &[Vec<f64>],
    lat: f64,
    lon: f64,
    grid_lat1: f64,
    grid_dlat: f64,
    grid_lon1: f64,
    grid_dlon: f64,
) -> f64 {
    let nlat = map.len();
    if nlat == 0 {
        return 0.0;
    }
    let nlon = map[0].len();
    if nlon == 0 {
        return 0.0;
    }

    let lat_frac = (lat - grid_lat1) / grid_dlat;
    let i0 = (libm::floor(lat_frac) as isize).clamp(0, nlat as isize - 1) as usize;
    let i1 = (i0 + 1).min(nlat - 1);
    let lat_w1 = lat_frac - i0 as f64;
    let lat_w0 = 1.0 - lat_w1;

    let mut lon_deg = lon % 360.0;
    if lon_deg < grid_lon1 {
        lon_deg += 360.0;
    }
    let lon_frac = (lon_deg - grid_lon1) / grid_dlon;
    let j0 = (libm::floor(lon_frac) as isize).clamp(0, nlon as isize - 1) as usize;
    let j1 = if j0 + 1 < nlon { j0 + 1 } else { 0 };
    let lon_w1 = lon_frac - j0 as f64;
    let lon_w0 = 1.0 - lon_w1;

    let v00 = map[i0][j0];
    let v01 = map[i0][j1];
    let v10 = map[i1][j0];
    let v11 = map[i1][j1];

    lat_w0 * (lon_w0 * v00 + lon_w1 * v01) + lat_w1 * (lon_w0 * v10 + lon_w1 * v11)
}

/// Computes ionospheric delay in meters from an IONEX TEC grid.
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
    if tec_maps.is_empty() || el <= 0.0 || tec_maps[0].1.is_empty() || tec_maps[0].1[0].is_empty() {
        return 0.0;
    }

    let re = 6371.0;
    let h = height_km;
    let psi = (libm::asin(re / (re + h) * libm::cos(el)) - (core::f64::consts::PI / 2.0 - el)).max(0.0);
    let lat_r = pos_llh.x.to_degrees();
    let lon_r = pos_llh.y.to_degrees();
    let ipp_lat = libm::asin(
        libm::sin(lat_r.to_radians()) * libm::cos(psi)
            + libm::cos(lat_r.to_radians()) * libm::sin(psi) * libm::cos(az),
    ).to_degrees().clamp(grid_lat1.min(grid_lat2), grid_lat1.max(grid_lat2));
    let ipp_lon = (lon_r + libm::asin(libm::sin(psi) * libm::sin(az) / libm::cos(ipp_lat.to_radians())).to_degrees())
        .clamp(grid_lon1.min(grid_lon2), grid_lon1.max(grid_lon2));

    let n = tec_maps.len();
    let (t1_idx, t2_idx, frac) = if time.tow <= tec_maps[0].0.tow {
        (0, 0, 0.0)
    } else if time.tow >= tec_maps[n - 1].0.tow {
        (n - 1, n - 1, 0.0)
    } else {
        match tec_maps.binary_search_by(|m| m.0.tow.partial_cmp(&time.tow).expect("TOW is never NaN")) {
            Ok(idx) => (idx, idx, 0.0),
            Err(idx) if idx >= n => (n - 1, n - 1, 0.0),
            Err(0) => (0, 0, 0.0),
            Err(idx) => {
                let dt = tec_maps[idx].0.tow - tec_maps[idx - 1].0.tow;
                let fr = if dt > 0.0 { (time.tow - tec_maps[idx - 1].0.tow) / dt } else { 0.0 };
                (idx - 1, idx, fr)
            }
        }
    };

    let vtec1 = bilinear_interpolate_tec(tec_maps[t1_idx].1, ipp_lat, ipp_lon, grid_lat1, grid_dlat, grid_lon1, grid_dlon);
    let vtec2 = bilinear_interpolate_tec(tec_maps[t2_idx].1, ipp_lat, ipp_lon, grid_lat1, grid_dlat, grid_lon1, grid_dlon);
    let vtec = vtec1 * (1.0 - frac) + vtec2 * frac;

    let chi = libm::asin(re / (re + h) * libm::cos(el));
    let mf = 1.0 / libm::cos(chi);
    (0.1624 * vtec * mf).max(0.0)
}
