//! Shared Keplerian orbit propagation, used by every broadcast-ephemeris
//! constellation except GLONASS (which integrates numerically instead --
//! see `ephemeris.rs`'s `GlonassEphemeris`). Split out of `ephemeris.rs`
//! (CLAUDE.md's 500-line file standard): GPS, Galileo, BeiDou, and QZSS
//! each call [`calc_keplerian`] from their own `position()` method with
//! their own broadcast parameters -- this is genuinely shared orbital
//! mechanics, not a per-constellation concern. See
//! `docs/PROJECT_STATUS.md` Sprint 13.

use nalgebra::Vector3;

use crate::ephemeris::F;
use crate::time::GpsTime;

fn solve_eccentric_anomaly(mk: f64, e: f64) -> f64 {
    let mut ek = mk;
    for _ in 0..10 {
        ek = mk + e * libm::sin(ek);
    }
    ek
}

#[allow(clippy::too_many_arguments)]
fn compute_orbit_plane(
    ek: f64,
    e: f64,
    omega: f64,
    cus: f64,
    cuc: f64,
    crs: f64,
    crc: f64,
    cis: f64,
    cic: f64,
    a: f64,
    i0: f64,
    idot: f64,
    tk: f64,
) -> (f64, f64, f64) {
    let cos_ek = libm::cos(ek);
    let sin_ek = libm::sin(ek);
    let vk = libm::atan2(libm::sqrt(1.0 - e * e) * sin_ek, cos_ek - e);
    let uk = omega + vk;
    let sin_2uk = libm::sin(2.0 * uk);
    let cos_2uk = libm::cos(2.0 * uk);
    let u = uk + cus * sin_2uk + cuc * cos_2uk;
    let r = a * (1.0 - e * cos_ek) + crs * sin_2uk + crc * cos_2uk;
    let i = i0 + cis * sin_2uk + cic * cos_2uk + idot * tk;
    (u, r, i)
}

fn compute_orbit_velocities(a: f64, e: f64, ek: f64, n: f64, u: f64, r: f64) -> (f64, f64) {
    let cos_ek = libm::cos(ek);
    let rk_dot = a * e * libm::sin(ek) * n / (1.0 - e * cos_ek);
    let uk_dot = (libm::sqrt(1.0 - e * e) / (1.0 - e * cos_ek)) * n;
    let xk_prime_dot = rk_dot * libm::cos(u) - r * libm::sin(u) * uk_dot;
    let yk_prime_dot = rk_dot * libm::sin(u) + r * libm::cos(u) * uk_dot;
    (xk_prime_dot, yk_prime_dot)
}

#[allow(clippy::too_many_arguments)]
fn rotate_to_ecef(
    xk_p: f64,
    yk_p: f64,
    vx_p: f64,
    vy_p: f64,
    omegak: f64,
    om_dot: f64,
    i: f64,
    idot: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    let cos_ok = libm::cos(omegak);
    let sin_ok = libm::sin(omegak);
    let cos_i = libm::cos(i);
    let sin_i = libm::sin(i);
    let x = xk_p * cos_ok - yk_p * cos_i * sin_ok;
    let y = xk_p * sin_ok + yk_p * cos_i * cos_ok;
    let z = yk_p * sin_i;
    let vx = vx_p * cos_ok - vy_p * cos_i * sin_ok - y * om_dot - yk_p * sin_i * idot * sin_ok;
    let vy = vx_p * sin_ok + vy_p * cos_i * cos_ok + x * om_dot + yk_p * sin_i * idot * cos_ok;
    let vz = vy_p * sin_i + yk_p * cos_i * idot;
    (x, y, z, vx, vy, vz)
}

#[allow(clippy::too_many_arguments)]
fn apply_bds_geo_rotation(
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    tk: f64,
    omega_e: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    let sin_5 = libm::sin(-5.0f64.to_radians());
    let cos_5 = libm::cos(-5.0f64.to_radians());
    let sin_oet = libm::sin(omega_e * tk);
    let cos_oet = libm::cos(omega_e * tk);
    let xg = x * cos_oet + y * sin_oet * cos_5 + z * sin_oet * sin_5;
    let yg = -x * sin_oet + y * cos_oet * cos_5 + z * cos_oet * sin_5;
    let zg = -y * sin_5 + z * cos_5;
    let vxg = vx * cos_oet - x * omega_e * sin_oet
        + vy * sin_oet * cos_5
        + y * omega_e * cos_oet * cos_5
        + vz * sin_oet * sin_5
        + z * omega_e * cos_oet * sin_5;
    let vyg = -vx * sin_oet - x * omega_e * cos_oet + vy * cos_oet * cos_5
        - y * omega_e * sin_oet * cos_5
        + vz * cos_oet * sin_5
        - z * omega_e * sin_oet * sin_5;
    let vzg = -vy * sin_5 + vz * cos_5;
    (xg, yg, zg, vxg, vyg, vzg)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_keplerian(
    t: GpsTime,
    toe: GpsTime,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    crs: f64,
    crc: f64,
    cuc: f64,
    cus: f64,
    cic: f64,
    cis: f64,
    m0: f64,
    e: f64,
    sqrt_a: f64,
    delta_n: f64,
    omega0: f64,
    omega_dot: f64,
    i0: f64,
    idot: f64,
    omega: f64,
    tgd: f64,
    mu: f64,
    omega_e: f64,
    is_bds_geo: bool,
) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
    let tk = t - toe;
    let a = sqrt_a * sqrt_a;
    let n = libm::sqrt(mu / (a * a * a)) + delta_n;
    let ek = solve_eccentric_anomaly(m0 + n * tk, e);

    let (u, r, i) =
        compute_orbit_plane(ek, e, omega, cus, cuc, crs, crc, cis, cic, a, i0, idot, tk);
    let (vx_p, vy_p) = compute_orbit_velocities(a, e, ek, n, u, r);

    let xk_p = r * libm::cos(u);
    let yk_p = r * libm::sin(u);

    let omegak = if is_bds_geo {
        omega0 + omega_dot * tk - omega_e * toe.tow
    } else {
        omega0 + (omega_dot - omega_e) * tk - omega_e * toe.tow
    };
    let om_dot = if is_bds_geo {
        omega_dot
    } else {
        omega_dot - omega_e
    };

    let (mut x, mut y, mut z, mut vx, mut vy, mut vz) =
        rotate_to_ecef(xk_p, yk_p, vx_p, vy_p, omegak, om_dot, i, idot);

    if is_bds_geo {
        let rot = apply_bds_geo_rotation(x, y, z, vx, vy, vz, tk, omega_e);
        x = rot.0;
        y = rot.1;
        z = rot.2;
        vx = rot.3;
        vy = rot.4;
        vz = rot.5;
    }

    let tc = t - toc;
    let clk_err = af0 + af1 * tc + af2 * tc * tc + F * e * sqrt_a * libm::sin(ek) - tgd;
    let clk_drift = af1 + 2.0 * af2 * tc;

    (
        Vector3::new(x, y, z),
        Vector3::new(vx, vy, vz),
        clk_err,
        clk_drift,
    )
}
