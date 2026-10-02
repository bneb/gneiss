//! Synthetic end-to-end Doppler fixtures for the velocity estimator.
//!
//! Sign convention under test (the classic trap in Doppler velocimetry).
//! The receiver measures carrier frequency OFFSET from the carrier, so
//!     f_D = -(1/lambda) * d(rho)/dt ,   hence   d(rho)/dt = -lambda * f_D .
//! A receiver receding from the satellite at speed s therefore reports a
//! POSITIVE range rate of +s, and one closing on it reports -s. The estimator
//! builds its measurement as `rho_meas = -lambda * f_D` and its predicted
//! satellite-side rate as `rho_sat = e_los . v_sat - c * d(sat clock)/dt`,
//! giving the innovation
//!     y = d(rho)/dt_measured - d(rho)/dt_predicted = -(e_los . v_rx) + c*dtdt_rx
//! whose design-matrix row is exactly `[-e_los, +1]` -- see `solve_wls`.

use super::*;
use gneiss_core::ephemeris::GpsEphemeris;
use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::SatelliteId;

pub const C: f64 = SPEED_OF_LIGHT_M_S;
pub const OMEGA: f64 = 7.292_115_146_7e-5;

/// Receiver at the equator on the prime meridian: ECEF Up = +x,
/// North = +z, East = +y.
pub fn receiver() -> Vector3<f64> {
    Vector3::new(6_378_137.0, 0.0, 0.0)
}

/// A valid GPS broadcast ephemeris with distinct orbital elements per PRN.
pub fn ephemeris(prn: u8, t: GpsTime) -> Ephemeris {
    Ephemeris::Gps(GpsEphemeris {
        sat: SatelliteId { constellation: Constellation::Gps, prn },
        toe: t,
        toc: t,
        af0: 0.0,
        af1: 0.0,
        af2: 0.0,
        crs: 0.0,
        crc: 0.0,
        cuc: 0.0,
        cus: 0.0,
        cic: 0.0,
        cis: 0.0,
        m0: f64::from(prn) * core::f64::consts::PI / 3.0,
        e: 0.01,
        sqrt_a: 5153.6,
        delta_n: 0.0,
        omega0: f64::from(prn) * core::f64::consts::PI / 2.0,
        omega_dot: 0.0,
        i0: 0.95,
        idot: 0.0,
        omega: 0.0,
        tgd: 0.0,
        iode: 1,
        iodc: 1,
    })
}

fn sagnac(v: Vector3<f64>, tau: f64) -> Vector3<f64> {
    let th = OMEGA * tau;
    let (c, s) = (th.cos(), th.sin());
    Vector3::new(v.x * c + v.y * s, -v.x * s + v.y * c, v.z)
}

/// Builds one satellite whose Doppler encodes a receiver velocity `v_rx`.
///
/// Everything here is the textbook forward model:
///
/// 1. converge the light travel time tau = |r_sat(t_rx - tau) - r_rx| / c;
/// 2. rotate the satellite state from transmit ECEF to receive ECEF,
///    accounting for Earth rotation during flight;
/// 3. rho_dot = e_los . (v_sat - v_rx) - c * sat_clock_drift;
/// 4. f_D = -rho_dot / lambda.
///
/// The satellite velocity is a central difference with a 1 s step (the
/// estimator uses 0.05 s), so this is an independent numerical path.
pub fn sat_obs(eph: &Ephemeris, v_rx: Vector3<f64>) -> SatObs {
    let t_rx = eph.toe();
    let r_rcv = receiver();
    let lambda = get_wavelength(eph.sat(), 1, eph.freq_num());
    let mut tau = 0.0;
    for _ in 0..8 {
        let p = eph.position(GpsTime::new(t_rx.week, t_rx.tow - tau)).0;
        tau = (p - r_rcv).norm() / C;
    }
    let t_tx = GpsTime::new(t_rx.week, t_rx.tow - tau);
    let p = eph.position(t_tx).0;
    let (_, _, _, drift) = eph.position(t_tx);
    // 1 s central difference (the estimator uses 0.05 s).
    let v_sat = (eph.position(GpsTime::new(t_tx.week, t_tx.tow + 1.0)).0
        - eph.position(GpsTime::new(t_tx.week, t_tx.tow - 1.0)).0)
        / 2.0;
    let p = sagnac(p, tau);
    let v_sat = sagnac(v_sat, tau);
    let e_los = (p - r_rcv).normalize();
    let rho_dot = e_los.dot(&(v_sat - v_rx)) - C * drift;
    SatObs {
        sat: eph.sat(),
        observations: vec![Observation {
            code: ObsCode {
                obs_type: ObsType::Doppler,
                signal: SignalCode { freq_band: 1, attribute: 'D' },
            },
            value: -rho_dot / lambda,
            lock_time: None,
            lli: None,
        }],
    }
}

pub fn epoch(sats: Vec<SatObs>, t: GpsTime) -> EpochObs {
    EpochObs { time: t, satellites: sats }
}