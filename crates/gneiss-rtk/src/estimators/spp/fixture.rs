#![allow(clippy::unwrap_used)]

//! Shared synthetic SPP fixtures.
//!
//! Published carrier frequencies used by the tests below (MHz):
//!   GPS L1 = 1575.42, GPS L2 = 1227.60, Galileo E1 = 1575.42,
//!   Galileo E5b = 1207.14, BeiDou B1I = 1561.098, BeiDou B2I = 1207.14.

use super::*;
use gneiss_core::ephemeris::GpsEphemeris;
use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
use nalgebra::Vector3;

pub const L1_HZ: f64 = 1575.42e6;
pub const L2_HZ: f64 = 1227.60e6;
pub const E5B_HZ: f64 = 1207.14e6;
pub const B1I_HZ: f64 = 1561.098e6;
pub const B2I_HZ: f64 = 1207.14e6;

/// Nominal epoch used by every fixture.
pub fn epoch_time() -> GpsTime {
    GpsTime::new(2000, 100_000.0)
}

/// A valid GPS broadcast ephemeris with distinct orbital elements per PRN.
pub fn gps_ephemeris(prn: u8, t: GpsTime) -> Ephemeris {
    Ephemeris::Gps(GpsEphemeris {
        sat: gneiss_core::sat::SatelliteId { constellation: gneiss_core::sat::Constellation::Gps, prn },
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

/// One pseudorange observation on the given RINEX band.
pub fn pseudorange(band: u8, attr: char, value: f64) -> Observation {
    Observation {
        code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: band, attribute: attr } },
        value,
        lock_time: None,
        lli: None,
    }
}

pub fn sat_obs(constellation: gneiss_core::sat::Constellation, prn: u8, obs: Vec<Observation>) -> SatObs {
    SatObs { sat: gneiss_core::sat::SatelliteId { constellation, prn }, observations: obs }
}

/// Config that lets every satellite through: no elevation masking and a very
/// loose geometry gate, so the test exercises the estimator rather than gating.
pub fn open_config() -> SppConfig {
    SppConfig {
        enable_sagnac: false,
        enable_tropo: false,
        enable_iono: false,
        geometry_variance_threshold: 1e9,
        elevation_mask_rad: -core::f64::consts::PI,
        ..Default::default()
    }
}

/// Builds `n` GPS satellites plus a matching ephemeris list. Pseudoranges are
/// generated from the true geometry so the least-squares fix converges on
/// `true_pos`.
pub fn four_gps_scene(true_pos: Vector3<f64>, true_cdt: f64, t: GpsTime) -> (EpochObs, Vec<Ephemeris>) {
    let ephems: Vec<Ephemeris> = (1..=6u8).map(|p| gps_ephemeris(p, t)).collect();
    let satellites = ephems
        .iter()
        .map(|eph| {
            let raw_pr = solve_pseudorange(eph, true_pos, true_cdt, t);
            sat_obs(eph.sat().constellation, eph.sat().prn, vec![pseudorange(1, 'C', raw_pr)])
        })
        .collect();
    (EpochObs { time: t, satellites }, ephems)
}

/// Iterate the light-time / clock-bias loop to the geometric range the receiver
/// would measure, plus the receiver clock offset.
pub fn solve_pseudorange(eph: &Ephemeris, true_pos: Vector3<f64>, true_cdt: f64, t: GpsTime) -> f64 {
    let c = LIGHT_SPEED;
    let mut raw_pr = 20_000_000.0 + true_cdt;
    for _ in 0..6 {
        let t_tx = GpsTime::new(t.week, t.tow - raw_pr / c);
        let (_, _, clk_rough, _) = eph.position(t_tx);
        let t_tx = GpsTime::new(t.week, t.tow - raw_pr / c - clk_rough);
        let (sat_pos, _, clk, _) = eph.position(t_tx);
        let d = true_pos - sat_pos;
        raw_pr = d.norm() + true_cdt - clk * c;
    }
    raw_pr
}
/// A valid Galileo broadcast ephemeris (needed because `build_single_measurement`
/// matches observations to ephemerides by satellite id).
pub fn galileo_ephemeris(prn: u8, t: GpsTime) -> Ephemeris {
    Ephemeris::Galileo(gneiss_core::ephemeris::GalileoEphemeris {
        sat: gneiss_core::sat::SatelliteId { constellation: gneiss_core::sat::Constellation::Galileo, prn },
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
        m0: f64::from(prn) * core::f64::consts::PI / 4.0,
        e: 0.01,
        sqrt_a: 5153.6,
        delta_n: 0.0,
        omega0: f64::from(prn) * core::f64::consts::PI / 3.0,
        omega_dot: 0.0,
        i0: 0.96,
        idot: 0.0,
        omega: 0.0,
        bgd_e1_e5a: 0.0,
        bgd_e1_e5b: 0.0,
        iod_nav: 1,
    })
}

/// A valid BeiDou broadcast ephemeris (BDT-2 shape).
pub fn beidou_ephemeris(prn: u8, t: GpsTime) -> Ephemeris {
    Ephemeris::Beidou(gneiss_core::ephemeris::BeidouEphemeris {
        sat: gneiss_core::sat::SatelliteId { constellation: gneiss_core::sat::Constellation::Beidou, prn },
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
        m0: f64::from(prn) * core::f64::consts::PI / 5.0,
        e: 0.01,
        sqrt_a: 5282.0,
        delta_n: 0.0,
        omega0: f64::from(prn) * core::f64::consts::PI / 2.5,
        omega_dot: 0.0,
        i0: 0.96,
        idot: 0.0,
        omega: 0.0,
        tgd1: 0.0,
        tgd2: 0.0,
        aode: 1,
        aodc: 1,
    })
}
