//! Ephemeris builder functions mapping raw RINEX broadcast values into typed Ephemeris structs.

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

pub(crate) fn build_ephemeris(
    constellation: Constellation,
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<Ephemeris> {
    match constellation {
        Constellation::Glonass => build_glonass_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Gps => build_gps_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Galileo => build_galileo_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Beidou => build_beidou_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Qzss => build_qzss_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Navic => build_gps_ephemeris(sat, toc, af0, af1, af2, vals),
        _ => None,
    }
}

fn build_glonass_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<Ephemeris> {
    Some(Ephemeris::Glonass(
        gneiss_core::ephemeris::GlonassEphemeris {
            sat,
            toe: toc,
            freq_num: vals[7] as i8,
            tau_n: af0,
            gamma_n: af1,
            delta_tau_n: af2,
            x: vals[0] * 1000.0,
            y: vals[4] * 1000.0,
            z: vals[8] * 1000.0,
            vx: vals[1] * 1000.0,
            vy: vals[5] * 1000.0,
            vz: vals[9] * 1000.0,
            ax: vals[2] * 1000.0,
            ay: vals[6] * 1000.0,
            az: vals[10] * 1000.0,
        },
    ))
}

fn build_gps_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<Ephemeris> {
    Some(Ephemeris::Gps(
        gneiss_core::ephemeris::GpsEphemeris {
            sat,
            toc,
            toe: GpsTime::new(toc.week, vals[8]),
            af0,
            af1,
            af2,
            iode: vals[0] as u32,
            crs: vals[1],
            delta_n: vals[2],
            m0: vals[3],
            cuc: vals[4],
            e: vals[5],
            cus: vals[6],
            sqrt_a: vals[7],
            cic: vals[9],
            omega0: vals[10],
            cis: vals[11],
            i0: vals[12],
            crc: vals[13],
            omega: vals[14],
            omega_dot: vals[15],
            idot: vals[16],
            tgd: vals[22],
            iodc: vals[23] as u32,
        },
    ))
}

fn build_galileo_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<Ephemeris> {
    Some(Ephemeris::Galileo(
        gneiss_core::ephemeris::GalileoEphemeris {
            sat,
            toc,
            toe: GpsTime::new(toc.week, vals[8]),
            af0,
            af1,
            af2,
            iod_nav: vals[0] as u32,
            crs: vals[1],
            delta_n: vals[2],
            m0: vals[3],
            cuc: vals[4],
            e: vals[5],
            cus: vals[6],
            sqrt_a: vals[7],
            cic: vals[9],
            omega0: vals[10],
            cis: vals[11],
            i0: vals[12],
            crc: vals[13],
            omega: vals[14],
            omega_dot: vals[15],
            idot: vals[16],
            bgd_e1_e5a: vals[22],
            bgd_e1_e5b: vals[23],
        },
    ))
}

fn build_beidou_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<Ephemeris> {
    Some(Ephemeris::Beidou(
        gneiss_core::ephemeris::BeidouEphemeris {
            sat,
            toc,
            toe: GpsTime::new(toc.week, vals[8]),
            af0,
            af1,
            af2,
            aode: vals[0] as u32,
            crs: vals[1],
            delta_n: vals[2],
            m0: vals[3],
            cuc: vals[4],
            e: vals[5],
            cus: vals[6],
            sqrt_a: vals[7],
            cic: vals[9],
            omega0: vals[10],
            cis: vals[11],
            i0: vals[12],
            crc: vals[13],
            omega: vals[14],
            omega_dot: vals[15],
            idot: vals[16],
            tgd1: vals[22],
            tgd2: vals[23],
            aodc: vals[25] as u32,
        },
    ))
}

fn build_qzss_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<Ephemeris> {
    Some(Ephemeris::Qzss(
        gneiss_core::ephemeris::QzssEphemeris {
            sat,
            toc,
            toe: GpsTime::new(toc.week, vals[8]),
            af0,
            af1,
            af2,
            iode: vals[0] as u32,
            crs: vals[1],
            delta_n: vals[2],
            m0: vals[3],
            cuc: vals[4],
            e: vals[5],
            cus: vals[6],
            sqrt_a: vals[7],
            cic: vals[9],
            omega0: vals[10],
            cis: vals[11],
            i0: vals[12],
            crc: vals[13],
            omega: vals[14],
            omega_dot: vals[15],
            idot: vals[16],
            tgd: vals[22],
            iodc: vals[23] as u32,
        },
    ))
}
