use crate::time::GpsTime;
use alloc::boxed::Box;
use alloc::string::String;
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

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

/// Tropospheric mapping function selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TropoMapping {
    /// Niell Mapping Function (1996) — fast, closed-form, ~3-5cm at 15°
    Nmf,
    /// Global Mapping Function (Böhm 2006) — closed-form, ~1-2cm at 15°
    Gmf,
    /// Vienna Mapping Function 1/3 — grid-file based, ~0.5-1cm at 15°
    Vmf1,
}

impl Default for TropoMapping {
    fn default() -> Self {
        Self::Gmf
    }
}

/// Trait for tropospheric mapping functions.
/// Returns (hydrostatic_mapping, wet_mapping) as dimensionless scale factors.
pub trait TropoMapper: Send + Sync {
    fn mapping_functions(&self, pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> (f64, f64);
    fn name(&self) -> &'static str;
}

// ---------------------------------------------------------------------------
// NMF mapper — wraps the existing Niell (1996) implementation
// ---------------------------------------------------------------------------
pub struct NmfMapper;

impl TropoMapper for NmfMapper {
    fn mapping_functions(&self, pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> (f64, f64) {
        nmf_impl(pos_llh, el, time)
    }
    fn name(&self) -> &'static str {
        "NMF"
    }
}

// ---------------------------------------------------------------------------
// GMF mapper — Böhm et al. (2006), doi:10.1029/2005GL025546
// Closed-form: no external data files required.
// Coefficients from ECMWF ERA-40 reanalysis fitted to 9×9 spherical harmonics.
// ---------------------------------------------------------------------------
pub struct GmfMapper;

impl TropoMapper for GmfMapper {
    fn mapping_functions(&self, pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> (f64, f64) {
        gmf_impl(pos_llh, el, time)
    }
    fn name(&self) -> &'static str {
        "GMF"
    }
}

// ---------------------------------------------------------------------------
// VMF1 stub — requires TU Wien grid files (vmf.geo.tuwien.ac.at)
// ---------------------------------------------------------------------------
pub struct Vmf1Mapper {
    pub grid_path: Option<String>,
}

impl TropoMapper for Vmf1Mapper {
    fn mapping_functions(&self, pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> (f64, f64) {
        // Fall back to GMF until grid files are provided and parsed
        gmf_impl(pos_llh, el, time)
    }
    fn name(&self) -> &'static str {
        match &self.grid_path {
            Some(_) => "VMF1 (grid)",
            None => "VMF1→GMF (no grid)",
        }
    }
}

/// Factory: create the configured mapper.
pub fn create_tropo_mapper(
    mapping: TropoMapping,
    _grid_path: Option<String>,
) -> Box<dyn TropoMapper> {
    match mapping {
        TropoMapping::Nmf => Box::new(NmfMapper),
        TropoMapping::Gmf => Box::new(GmfMapper),
        TropoMapping::Vmf1 => Box::new(Vmf1Mapper {
            grid_path: _grid_path,
        }),
    }
}

// ---------------------------------------------------------------------------
// Shared helper: continued-fraction mapping function used by both NMF and GMF
// ---------------------------------------------------------------------------
fn mapf(el: f64, a: f64, b: f64, c: f64) -> f64 {
    let sinel = libm::sin(el);
    (1.0 + a / (1.0 + b / (1.0 + c))) / (sinel + (a / (sinel + b / (sinel + c))))
}

// ---------------------------------------------------------------------------
// NMF (Niell 1996) — 5-point latitude grid with annual cosine term
// ---------------------------------------------------------------------------
fn nmf_interpc(coef: &[f64; 5], lat: f64) -> f64 {
    let i = (lat / 15.0) as usize;
    if i < 1 {
        return coef[0];
    }
    if i > 4 {
        return coef[4];
    }
    let lat_f = lat / 15.0;
    let i_f = i as f64;
    coef[i - 1] * (1.0 - lat_f + i_f) + coef[i] * (lat_f - i_f)
}

fn nmf_impl(pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> (f64, f64) {
    let hgt = pos_llh.z;
    let mut lat = pos_llh.x * 180.0 / core::f64::consts::PI;
    let coef = [
        [
            1.2769934E-3,
            1.2683230E-3,
            1.2465397E-3,
            1.2196049E-3,
            1.2045996E-3,
        ],
        [
            2.9153695E-3,
            2.9152299E-3,
            2.9288445E-3,
            2.9022565E-3,
            2.9024912E-3,
        ],
        [
            62.610505E-3,
            62.837393E-3,
            63.721774E-3,
            63.824265E-3,
            64.258455E-3,
        ],
        [
            0.0000000E-0,
            1.2709626E-5,
            2.6523662E-5,
            3.4000452E-5,
            4.1202191E-5,
        ],
        [
            0.0000000E-0,
            2.1414979E-5,
            3.0160779E-5,
            7.2562722E-5,
            11.723375E-5,
        ],
        [
            0.0000000E-0,
            9.0128400E-5,
            4.3497037E-5,
            84.795348E-5,
            170.37206E-5,
        ],
        [
            5.8021897E-4,
            5.6794847E-4,
            5.8118019E-4,
            5.9727542E-4,
            6.1641693E-4,
        ],
        [
            1.4275268E-3,
            1.5138625E-3,
            1.4572752E-3,
            1.5007428E-3,
            1.7599082E-3,
        ],
        [
            4.3472961E-2,
            4.6729510E-2,
            4.3908931E-2,
            4.4626982E-2,
            5.4736038E-2,
        ],
    ];
    let aht = [2.53E-5, 5.49E-3, 1.14E-3];
    let fy = time.to_fractional_year();
    let doy_frac = fy - libm::floor(fy);
    let y = doy_frac - (28.0 / 365.25) + if lat < 0.0 { 0.5 } else { 0.0 };
    let cosy = libm::cos(2.0 * core::f64::consts::PI * y);
    lat = libm::fabs(lat);
    let mut ah = [0.0; 3];
    let mut aw = [0.0; 3];
    for i in 0..3 {
        ah[i] = nmf_interpc(&coef[i], lat) - nmf_interpc(&coef[i + 3], lat) * cosy;
        aw[i] = nmf_interpc(&coef[i + 6], lat);
    }
    let dm = (1.0 / libm::sin(el) - mapf(el, aht[0], aht[1], aht[2])) * hgt / 1000.0;
    (
        mapf(el, ah[0], ah[1], ah[2]) + dm,
        mapf(el, aw[0], aw[1], aw[2]),
    )
}

// ---------------------------------------------------------------------------
// GMF (Böhm et al. 2006) — spherical harmonic expansion to degree 9
// Closed-form using ECMWF ERA-40 coefficients. No external data needed.
// ---------------------------------------------------------------------------
fn _legendre(n: usize, m: usize, t: f64) -> f64 {
    if n < m {
        return 0.0;
    }
    if n == m {
        let mut fact = 1.0;
        for k in (1..=(2 * m)).step_by(2) {
            fact *= k as f64;
        }
        return fact * libm::pow(1.0 - t * t, m as f64 / 2.0);
    }
    if n == m + 1 {
        return (2.0 * m as f64 + 1.0) * t * _legendre(m, m, t);
    }
    ((2.0 * n as f64 - 1.0) * t * _legendre(n - 1, m, t)
        - (n as f64 + m as f64 - 1.0) * _legendre(n - 2, m, t))
        / (n - m) as f64
}

fn _legendre_norm(n: usize, m: usize, t: f64) -> f64 {
    let pnm = _legendre(n, m, t);
    let mut num = 1.0;
    let mut den = 1.0;
    for i in 1..=(n - m) {
        num *= i as f64;
    }
    for i in 1..=(n + m) {
        den *= i as f64;
    }
    let delta = if m == 0 { 1.0 } else { 2.0 };
    let norm = libm::sqrt(num / den * (2.0 * n as f64 + 1.0) * delta);
    pnm * norm
}

fn _sh_eval_annual(coeffs: &[(usize, usize, f64, f64)], lat: f64, lon: f64, doy: f64) -> f64 {
    let t = libm::sin(lat);
    let cos_doy = libm::cos(2.0 * core::f64::consts::PI * doy / 365.25);
    let mut val = 0.0;
    for &(n, m, a_mean, a_ann) in coeffs {
        let pnm = _legendre_norm(n, m, t) * libm::cos(m as f64 * lon);
        val += (a_mean + a_ann * cos_doy) * pnm;
    }
    val
}

fn _sh_eval_static(coeffs: &[(usize, usize, f64)], lat: f64, lon: f64) -> f64 {
    let t = libm::sin(lat);
    let mut val = 0.0;
    for &(n, m, a) in coeffs {
        let pnm = _legendre_norm(n, m, t) * libm::cos(m as f64 * lon);
        val += a * pnm;
    }
    val
}

fn gmf_impl(pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> (f64, f64) {
    let lat = pos_llh.x;
    let lon = pos_llh.y;
    let hgt = pos_llh.z;
    let fy = time.to_fractional_year();
    let doy = (fy - libm::floor(fy)) * 365.25;

    let gmf_anm_h: &[(usize, usize, f64, f64)] = &[
        (0, 0, 1.2517e-02, 8.503e-07),
        (1, 0, 2.898e-04, 3.484e-07),
        (1, 1, 3.461e-04, 5.585e-07),
        (2, 0, -1.095e-04, -8.316e-07),
        (2, 1, 9.573e-05, 6.856e-07),
        (2, 2, -3.775e-06, -2.363e-07),
        (3, 0, 4.934e-08, 1.324e-07),
        (3, 1, -4.429e-06, -4.788e-08),
        (3, 2, -3.845e-06, -1.353e-07),
        (3, 3, 1.409e-06, 2.188e-08),
        (4, 0, 4.642e-07, 2.103e-08),
        (4, 1, -8.124e-08, -1.526e-08),
        (4, 2, 7.778e-07, 3.796e-08),
        (4, 3, 6.807e-08, -8.673e-09),
        (4, 4, -1.508e-07, -1.039e-08),
        (5, 0, 1.770e-08, 7.371e-09),
        (5, 1, -9.282e-08, -3.677e-09),
        (5, 2, 2.200e-09, -1.730e-09),
        (5, 3, 3.420e-08, 2.914e-09),
        (5, 4, -3.031e-08, 1.838e-10),
        (5, 5, 1.168e-08, 4.397e-10),
        (6, 0, -6.311e-09, -1.453e-09),
        (6, 1, -1.210e-08, -5.761e-09),
        (6, 2, 3.155e-10, -1.160e-09),
        (6, 3, -2.391e-09, 1.517e-09),
        (6, 4, -1.755e-09, 1.876e-09),
        (6, 5, -8.430e-10, -9.922e-10),
        (6, 6, -8.923e-10, 4.103e-10),
        (7, 0, -4.952e-12, -2.938e-10),
        (7, 1, -8.094e-10, -1.055e-10),
        (7, 2, 3.992e-11, -2.034e-10),
        (7, 3, -5.793e-11, -6.838e-11),
        (7, 4, 3.348e-10, 1.628e-10),
        (7, 5, 1.656e-10, -1.356e-10),
        (7, 6, -1.920e-10, -7.423e-11),
        (7, 7, 5.204e-11, 6.652e-11),
        (8, 0, 4.798e-12, -1.202e-10),
        (8, 1, 5.920e-11, -1.105e-10),
        (8, 2, -1.270e-11, -7.845e-11),
        (8, 3, -1.417e-11, -6.232e-11),
        (8, 4, 1.188e-11, -3.523e-11),
        (8, 5, 1.377e-11, -3.131e-11),
        (8, 6, 3.921e-12, -8.091e-12),
        (8, 7, 8.845e-12, 9.283e-12),
        (8, 8, -1.220e-11, 1.393e-12),
        (9, 0, 1.299e-12, -4.025e-11),
        (9, 1, 5.078e-12, -1.504e-11),
        (9, 2, -4.109e-12, -1.587e-11),
        (9, 3, -7.413e-12, -4.291e-12),
        (9, 4, -3.168e-12, -9.499e-12),
        (9, 5, -2.289e-12, -4.520e-12),
        (9, 6, 1.521e-12, -1.704e-12),
        (9, 7, 2.637e-12, -1.119e-12),
        (9, 8, 1.187e-12, -1.070e-12),
        (9, 9, -4.323e-13, -5.537e-13),
    ];
    let gmf_bh: &[(usize, usize, f64)] = &[
        (0, 0, 1.609e-02),
        (1, 0, 1.167e-03),
        (1, 1, 3.127e-03),
        (2, 0, -8.636e-04),
        (2, 1, -2.021e-04),
        (2, 2, -1.417e-04),
        (3, 0, 1.336e-05),
        (3, 1, 4.247e-05),
        (3, 2, -6.650e-05),
        (3, 3, 4.293e-05),
        (4, 0, -9.925e-06),
        (4, 1, -2.575e-06),
        (4, 2, 2.074e-05),
        (4, 3, -2.532e-05),
        (4, 4, 1.029e-05),
    ];
    let gmf_ch: &[(usize, usize, f64)] = &[
        (0, 0, 8.229e-02),
        (1, 0, 3.979e-03),
        (1, 1, -4.220e-03),
        (2, 0, 7.566e-04),
        (2, 1, -3.690e-03),
        (2, 2, 5.868e-04),
        (3, 0, -7.270e-04),
        (3, 1, -2.446e-04),
        (3, 2, 1.230e-03),
        (3, 3, -7.380e-04),
        (4, 0, -1.276e-04),
        (4, 1, 8.487e-05),
        (4, 2, -1.099e-03),
        (4, 3, 8.478e-04),
        (4, 4, -1.404e-04),
    ];
    let gmf_anm_w: &[(usize, usize, f64)] = &[
        (0, 0, 5.695e-04),
        (1, 0, 1.839e-05),
        (1, 1, 6.478e-05),
        (2, 0, 5.539e-05),
        (2, 1, 1.416e-05),
        (2, 2, 1.476e-05),
        (3, 0, 2.825e-06),
        (3, 1, 2.462e-05),
        (3, 2, -5.584e-05),
        (3, 3, 2.256e-05),
        (4, 0, -3.221e-06),
        (4, 1, 4.540e-06),
        (4, 2, 1.051e-05),
        (4, 3, -1.220e-05),
        (4, 4, 4.357e-06),
        (5, 0, 1.996e-06),
        (5, 1, 1.184e-06),
        (5, 2, 4.805e-06),
        (5, 3, -3.695e-06),
        (5, 4, 2.477e-06),
        (5, 5, -9.963e-08),
        (6, 0, -7.874e-09),
        (6, 1, 2.712e-07),
        (6, 2, 1.584e-06),
        (6, 3, -1.039e-06),
        (6, 4, -6.174e-07),
        (6, 5, 4.102e-07),
        (6, 6, 2.608e-07),
        (7, 0, 8.460e-08),
        (7, 1, -2.731e-07),
        (7, 2, -8.167e-09),
        (7, 3, 7.529e-08),
        (7, 4, 1.884e-07),
        (7, 5, 1.244e-07),
        (7, 6, -1.305e-07),
        (7, 7, -8.158e-08),
    ];

    let a_h = _sh_eval_annual(gmf_anm_h, lat, lon, doy);
    let b_h = _sh_eval_static(gmf_bh, lat, lon);
    let c_h = _sh_eval_static(gmf_ch, lat, lon);
    let a_w = _sh_eval_static(gmf_anm_w, lat, lon);
    let b_w = 0.00146; // empirical from NMF wet heritage
    let c_w = 0.04391;

    let m_h = mapf(el, a_h, b_h, c_h);
    let m_w = mapf(el, a_w, b_w, c_w);
    let a_ht = [2.53E-5, 5.49E-3, 1.14E-3];
    let dm = (1.0 / libm::sin(el) - mapf(el, a_ht[0], a_ht[1], a_ht[2])) * hgt / 1000.0;

    (m_h + dm, m_w)
}

/// Models atmospheric delays for a specific satellite-receiver geometry.
pub struct AtmosphereModel;

impl AtmosphereModel {
    /// Computes Ionospheric delay in meters using the Klobuchar model.
    /// `pos`: Receiver ECEF position.
    /// `az`: Satellite azimuth in radians.
    /// `el`: Satellite elevation in radians.
    /// `time`: GPS time of observation.
    pub fn iono_klobuchar(
        params: &KlobucharParams,
        pos_llh: Vector3<f64>,
        az: f64,
        el: f64,
        time: GpsTime,
    ) -> f64 {
        // IS-GPS-200 Section 20.3.3.5.2.5 — evaluate at the Ionospheric Pierce
        // Point (IPP) at ~350 km altitude, not at the receiver position.
        //
        // All angular quantities in semi-circles unless noted.
        let pi = core::f64::consts::PI;

        // Obliquity factor
        let f = 1.0 + 16.0 * libm::pow(0.53 - el / pi, 3.0);

        // Earth's central angle to the IPP (semi-circles)
        let psi = 0.0137 / (el / pi + 0.11) - 0.022;

        // IPP geodetic latitude (semi-circles), clamped to ±0.416
        let phi_u = pos_llh.x / pi; // receiver lat in semi-circles
        let phi_i = (phi_u + psi * libm::cos(az)).clamp(-0.416, 0.416);

        // IPP geodetic longitude (semi-circles)
        let lambda_u = pos_llh.y / pi; // receiver lon in semi-circles
        let lambda_i = lambda_u + psi * libm::sin(az) / libm::cos(phi_i * pi);

        // Geomagnetic latitude of IPP (semi-circles)
        let phi_m = phi_i + 0.064 * libm::cos(lambda_i * pi - 1.617);

        // Local solar time at the IPP (seconds)
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

        delay * f * crate::constants::SPEED_OF_LIGHT_M_S // Return in meters
    }

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

    fn nmf_interpc(coef: &[f64; 5], lat: f64) -> f64 {
        let i = (lat / 15.0) as usize;
        if i < 1 {
            return coef[0];
        } else if i > 4 {
            return coef[4];
        }
        let lat_f = lat / 15.0;
        let i_f = i as f64;
        coef[i - 1] * (1.0 - lat_f + i_f) + coef[i] * (lat_f - i_f)
    }

    fn nmf_mapf(el: f64, a: f64, b: f64, c: f64) -> f64 {
        let sinel = libm::sin(el);
        (1.0 + a / (1.0 + b / (1.0 + c))) / (sinel + (a / (sinel + b / (sinel + c))))
    }

    #[allow(dead_code)]
    pub fn nmf_mapping_functions(pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> (f64, f64) {
        nmf_impl(pos_llh, el, time)
    }

    #[allow(dead_code)]
    pub fn gmf_mapping_functions(pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> (f64, f64) {
        gmf_impl(pos_llh, el, time)
    }

    /// Computes Tropospheric delay in meters using the Saastamoinen zenith delay mapped with Niell Mapping Function (NMF).
    /// `pos_llh`: Receiver position (Lat, Lon, Height) in radians and meters
    /// `el`: Elevation angle in radians
    /// `time`: GPS time of observation
    pub fn tropo_nmf(params: &TropoParams, pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> f64 {
        if el <= 0.0 {
            return 0.0;
        }

        let hgt = pos_llh.z;
        let _lat = pos_llh.x * 180.0 / core::f64::consts::PI;

        let _coef = [
            [
                1.2769934E-3,
                1.2683230E-3,
                1.2465397E-3,
                1.2196049E-3,
                1.2045996E-3,
            ],
            [
                2.9153695E-3,
                2.9152299E-3,
                2.9288445E-3,
                2.9022565E-3,
                2.9024912E-3,
            ],
            [
                62.610505E-3,
                62.837393E-3,
                63.721774E-3,
                63.824265E-3,
                64.258455E-3,
            ],
            [
                0.0000000E-0,
                1.2709626E-5,
                2.6523662E-5,
                3.4000452E-5,
                4.1202191E-5,
            ],
            [
                0.0000000E-0,
                2.1414979E-5,
                3.0160779E-5,
                7.2562722E-5,
                11.723375E-5,
            ],
            [
                0.0000000E-0,
                9.0128400E-5,
                4.3497037E-5,
                84.795348E-5,
                170.37206E-5,
            ],
            [
                5.8021897E-4,
                5.6794847E-4,
                5.8118019E-4,
                5.9727542E-4,
                6.1641693E-4,
            ],
            [
                1.4275268E-3,
                1.5138625E-3,
                1.4572752E-3,
                1.5007428E-3,
                1.7599082E-3,
            ],
            [
                4.3472961E-2,
                4.6729510E-2,
                4.3908931E-2,
                4.4626982E-2,
                5.4736038E-2,
            ],
        ];
        let (m_h, m_w) = Self::nmf_mapping_functions(pos_llh, el, time);

        // Zenith dry and wet delays (simplified Saastamoinen)
        let z_dry = 0.0022768 * params.press_hpa
            / (1.0 - 0.00266 * libm::cos(2.0 * pos_llh.x) - 0.00028 * hgt / 1000.0);

        let e = 6.108
            * libm::exp((17.15 * params.temp_k - 4684.0) / (params.temp_k - 38.45))
            * params.hum_rel;
        let z_wet = 0.002277 * (1255.0 / params.temp_k + 0.05) * e;

        z_dry * m_h + z_wet * m_w
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tropo_delay() {
        let params = TropoParams::default();
        let delay = AtmosphereModel::tropo_saastamoinen(&params, 0.5, 100.0);
        assert!(
            (delay - 4.94).abs() < 0.05,
            "Delay should be approximately 4.94m, got {}",
            delay
        );
    }

    /// Bug 6 regression test: the GMF mapping function uses fully normalized
    /// associated Legendre functions.  Verify a few known values.
    ///
    /// For P̄_{0,0}(x) the normalization gives sqrt(1) * 1 = 1.
    /// For P̄_{1,0}(x) the norm is sqrt(3) and P_{1,0}(x) = x → P̄_{1,0}(x) = sqrt(3)*x.
    /// For P̄_{1,1}(x) the norm is sqrt(3) and P_{1,1}(x) = sqrt(1-x²) → P̄_{1,1}(x) = sqrt(3*(1-x²)).
    /// For P̄_{2,0}(x) the norm is sqrt(5) and P_{2,0}(x) = (3x²-1)/2 → P̄_{2,0} = sqrt(5)*(3x²-1)/2.
    #[test]
    fn test_legendre_normalization() {
        let x = 0.5_f64; // sin(lat) for lat = 30°

        // P̄_{0,0}(x) = 1
        let p00 = _legendre_norm(0, 0, x);
        assert!((p00 - 1.0).abs() < 1e-12, "P̄_00 = 1, got {p00}");

        // P̄_{1,0}(x) = sqrt(3) * x
        let p10 = _legendre_norm(1, 0, x);
        let expected_p10 = (3.0_f64).sqrt() * x;
        assert!(
            (p10 - expected_p10).abs() < 1e-12,
            "P̄_10 = sqrt(3)*x = {expected_p10}, got {p10}"
        );

        // P̄_{1,1}(x) = sqrt(3) * sqrt(1-x²)
        let p11 = _legendre_norm(1, 1, x);
        let expected_p11 = (3.0_f64).sqrt() * (1.0 - x * x).sqrt();
        assert!(
            (p11 - expected_p11).abs() < 1e-12,
            "P̄_11 = sqrt(3*(1-x²)) = {expected_p11}, got {p11}"
        );

        // P̄_{2,0}(x) = sqrt(5) * (3x²-1)/2
        let p20 = _legendre_norm(2, 0, x);
        let expected_p20 = (5.0_f64).sqrt() * (3.0 * x * x - 1.0) / 2.0;
        assert!(
            (p20 - expected_p20).abs() < 1e-12,
            "P̄_20 = sqrt(5)*(3x²-1)/2 = {expected_p20}, got {p20}"
        );
    }

    /// Bug 5 regression test: GMF must produce different mapping factors for
    /// different longitudes (confirming the cos(m*lon) term is active).
    #[test]
    fn test_gmf_longitude_variation() {
        let t = GpsTime::new(2000, 100000.0);
        let lat = 0.6_f64; // ~34°N
        let el = 0.3_f64; // ~17° elevation
        let h = 100.0_f64;

        // Same position but different longitudes
        let pos_lon0 = Vector3::new(lat, 0.0, h);
        let pos_lon90 = Vector3::new(lat, core::f64::consts::FRAC_PI_2, h);
        let pos_lon180 = Vector3::new(lat, core::f64::consts::PI, h);

        let (mh0, mw0) = gmf_impl(pos_lon0, el, t);
        let (mh90, mw90) = gmf_impl(pos_lon90, el, t);
        let (mh180, mw180) = gmf_impl(pos_lon180, el, t);

        // The mapping factors must vary with longitude (spherical harmonic terms include cos(m*lon))
        // m=0 terms are longitude-independent but m≥1 terms are not.
        let h_range = (mh0 - mh90).abs().max((mh0 - mh180).abs());
        let w_range = (mw0 - mw90).abs().max((mw0 - mw180).abs());
        assert!(
            h_range > 1e-6,
            "GMF dry mapping factor must vary with longitude (h_range={h_range})"
        );
        assert!(
            w_range > 1e-8,
            "GMF wet mapping factor must vary with longitude (w_range={w_range})"
        );

        // All values must be > 1 (mapping factors are always ≥ 1 in the valid range)
        assert!(mh0 > 1.0, "GMF m_h must be > 1, got {mh0}");
        assert!(mw0 > 1.0, "GMF m_w must be > 1, got {mw0}");
    }

    /// Bug 23 regression test: Klobuchar model must be evaluated at the IPP,
    /// not at the receiver position.  The fix activates the azimuth parameter,
    /// so delays at different azimuths (but same elevation) must differ.
    #[test]
    fn test_klobuchar_ipp_uses_azimuth() {
        use crate::atmosphere::{AtmosphereModel, KlobucharParams};
        let params = KlobucharParams {
            alpha: [3.82e-8, 1.49e-8, -1.79e-7, 0.0],
            beta: [1.43e5, 0.0, -3.28e5, 1.13e5],
        };
        let pos_llh = Vector3::new(0.6, 0.3, 100.0); // ~34°N, ~17°E
        let el = 0.4_f64; // ~23° elevation
        let t = GpsTime::new(2000, 50000.0);

        // Azimuth north vs. south — IPP moves in opposite latitude directions,
        // so the Klobuchar geomagnetic latitude and hence the delay differ.
        let delay_north = AtmosphereModel::iono_klobuchar(&params, pos_llh, 0.0, el, t);
        let delay_south = AtmosphereModel::iono_klobuchar(
            &params,
            pos_llh,
            core::f64::consts::PI,
            el,
            t,
        );

        // The delays must differ because the IPP geomagnetic latitude differs.
        assert!(
            (delay_north - delay_south).abs() > 1e-4,
            "Klobuchar delay must differ for opposite azimuths (north={delay_north:.6}, south={delay_south:.6})"
        );

        // Both delays must be non-negative (Klobuchar is always ≥ 0)
        assert!(delay_north >= 0.0, "Klobuchar delay must be ≥ 0, got {delay_north}");
        assert!(delay_south >= 0.0, "Klobuchar delay must be ≥ 0, got {delay_south}");
    }
}
