//! GNSS broadcast ephemeris orbital computation for GPS, GLONASS, Galileo, BeiDou, and QZSS.

pub mod keplerian;
pub mod glonass;
#[cfg(test)]
mod tests;

pub use keplerian::{BeidouEphemeris, GalileoEphemeris, GpsEphemeris, QzssEphemeris};
pub use glonass::GlonassEphemeris;

use crate::sat::SatelliteId;
use crate::time::GpsTime;
use nalgebra::Vector3;

pub(crate) const MU_GPS: f64 = 3.986005e14;
pub(crate) const MU_GAL: f64 = 3.986004418e14;
pub(crate) const MU_BDS: f64 = 3.986004418e14;
pub(crate) const MU_GLO: f64 = 3.9860044e14;
pub(crate) const OMEGA_E_GPS: f64 = crate::constants::EARTH_ROTATION_RATE_RAD_S;
pub(crate) const OMEGA_E_GAL: f64 = crate::constants::EARTH_ROTATION_RATE_RAD_S;
pub(crate) const OMEGA_E_BDS: f64 = 7.292115e-5;
pub(crate) const OMEGA_E_GLO: f64 = 7.292115e-5;
pub(crate) const J2_GLO: f64 = 1.0826257e-3;
pub(crate) const RADIUS_GLO: f64 = 6378136.0;
pub(crate) const F: f64 = -4.442807633e-10;

/// Broadcast ephemeris variant supporting all major constellations.
#[derive(Debug, Clone, PartialEq)]
pub enum Ephemeris {
    Gps(GpsEphemeris),
    Galileo(GalileoEphemeris),
    Beidou(BeidouEphemeris),
    Qzss(QzssEphemeris),
    Glonass(GlonassEphemeris),
}

impl Ephemeris {
    pub fn sat(&self) -> SatelliteId {
        match self {
            Ephemeris::Gps(e) => e.sat,
            Ephemeris::Galileo(e) => e.sat,
            Ephemeris::Beidou(e) => e.sat,
            Ephemeris::Qzss(e) => e.sat,
            Ephemeris::Glonass(e) => e.sat,
        }
    }

    pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        match self {
            Ephemeris::Gps(e) => e.position(t),
            Ephemeris::Galileo(e) => e.position(t),
            Ephemeris::Beidou(e) => e.position(t),
            Ephemeris::Qzss(e) => e.position(t),
            Ephemeris::Glonass(e) => e.position(t),
        }
    }

    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        match self {
            Ephemeris::Gps(e) => e.position_iono_free(t),
            Ephemeris::Galileo(e) => e.position_iono_free(t),
            Ephemeris::Beidou(e) => e.position_iono_free(t),
            Ephemeris::Qzss(e) => e.position_iono_free(t),
            Ephemeris::Glonass(e) => e.position(t),
        }
    }

    pub fn position_e5b(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        match self {
            Ephemeris::Galileo(e) => e.position_e5b(t),
            other => other.position(t),
        }
    }

    pub fn toe(&self) -> GpsTime {
        match self {
            Ephemeris::Gps(e) => e.toe,
            Ephemeris::Galileo(e) => e.toe,
            Ephemeris::Beidou(e) => e.toe,
            Ephemeris::Qzss(e) => e.toe,
            Ephemeris::Glonass(e) => e.toe,
        }
    }

    pub fn freq_num(&self) -> i8 {
        match self {
            Ephemeris::Glonass(e) => e.freq_num,
            _ => 0,
        }
    }

    pub fn tgd(&self) -> f64 {
        match self {
            Ephemeris::Gps(e) => e.tgd,
            Ephemeris::Galileo(e) => e.bgd_e1_e5a,
            Ephemeris::Beidou(e) => e.tgd1,
            Ephemeris::Qzss(e) => e.tgd,
            Ephemeris::Glonass(_) => 0.0,
        }
    }

    pub fn bgd_e5b(&self) -> f64 {
        match self {
            Ephemeris::Galileo(e) => e.bgd_e1_e5b,
            other => other.tgd(),
        }
    }
}
