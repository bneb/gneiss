//! Tropospheric mapping functions (NMF, GMF).

use crate::time::GpsTime;
use nalgebra::Vector3;

/// Continued-fraction mapping function used by both NMF and GMF.
pub(crate) fn mapf(el: f64, a: f64, b: f64, c: f64) -> f64 {
    let sinel = libm::sin(el);
    (1.0 + a / (1.0 + b / (1.0 + c))) / (sinel + (a / (sinel + b / (sinel + c))))
}

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

pub fn nmf_impl(pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> (f64, f64) {
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

#[cfg(test)]
pub(crate) fn _legendre(n: usize, m: usize, t: f64) -> f64 {
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

#[cfg(test)]
pub(crate) fn _legendre_norm(n: usize, m: usize, t: f64) -> f64 {
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

#[cfg(test)]
pub(crate) fn _sh_eval_annual(
    coeffs_cos: &[(usize, usize, f64, f64)],
    coeffs_sin: Option<&[(usize, usize, f64, f64)]>,
    lat: f64,
    lon: f64,
    doy: f64,
) -> f64 {
    let t = libm::sin(lat);
    let cos_doy = libm::cos(2.0 * core::f64::consts::PI * doy / 365.25);
    let mut val = 0.0;
    for &(n, m, a_mean, a_ann) in coeffs_cos {
        let pnm = _legendre_norm(n, m, t);
        val += (a_mean + a_ann * cos_doy) * pnm * libm::cos(m as f64 * lon);
    }
    if let Some(sin_coeffs) = coeffs_sin {
        for &(n, m, b_mean, b_ann) in sin_coeffs {
            if m == 0 {
                continue;
            }
            let pnm = _legendre_norm(n, m, t);
            val += (b_mean + b_ann * cos_doy) * pnm * libm::sin(m as f64 * lon);
        }
    }
    val
}

#[cfg(test)]
pub(crate) fn _sh_eval_static(
    coeffs_cos: &[(usize, usize, f64)],
    coeffs_sin: Option<&[(usize, usize, f64)]>,
    lat: f64,
    lon: f64,
) -> f64 {
    let t = libm::sin(lat);
    let mut val = 0.0;
    for &(n, m, a) in coeffs_cos {
        let pnm = _legendre_norm(n, m, t);
        val += a * pnm * libm::cos(m as f64 * lon);
    }
    if let Some(sin_coeffs) = coeffs_sin {
        for &(n, m, b) in sin_coeffs {
            if m == 0 {
                continue;
            }
            let pnm = _legendre_norm(n, m, t);
            val += b * pnm * libm::sin(m as f64 * lon);
        }
    }
    val
}

#[cfg(test)]
pub fn gmf_impl(pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> (f64, f64) {
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

    let a_h = _sh_eval_annual(gmf_anm_h, None, lat, lon, doy);
    let b_h = _sh_eval_static(gmf_bh, None, lat, lon);
    let c_h = _sh_eval_static(gmf_ch, None, lat, lon);
    let a_w = _sh_eval_static(gmf_anm_w, None, lat, lon);
    let b_w = 0.00146;
    let c_w = 0.04391;

    let m_h = mapf(el, a_h, b_h, c_h);
    let m_w = mapf(el, a_w, b_w, c_w);
    let a_ht = [2.53E-5, 5.49E-3, 1.14E-3];
    let dm = (1.0 / libm::sin(el) - mapf(el, a_ht[0], a_ht[1], a_ht[2])) * hgt / 1000.0;

    (m_h + dm, m_w)
}
