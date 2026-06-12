use super::*;
use gneiss_core::time::GpsTime;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::obs::{EpochObs, SatObs, Observation, ObsCode, SignalCode, ObsType};
use nalgebra::Vector3;

fn make_ephemeris(constellation: Constellation, prn: u16, t: GpsTime, pos: (f64, f64, f64)) -> Ephemeris {
    let sat_id = SatelliteId { constellation, prn };
    if constellation == Constellation::Gps {
        Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: sat_id, toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.00, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0,
            omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        })
    } else {
        // Mock others if needed, for now just GPS
        panic!("Unsupported test constellation")
    }
}
//...
