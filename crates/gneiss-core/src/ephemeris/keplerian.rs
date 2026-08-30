//! Keplerian orbit models for GPS, Galileo, BeiDou, and QZSS broadcast ephemeris.

use crate::keplerian::calc_keplerian;
use crate::sat::SatelliteId;
use crate::time::GpsTime;
use nalgebra::Vector3;
use super::{MU_BDS, MU_GAL, MU_GPS, OMEGA_E_BDS, OMEGA_E_GAL, OMEGA_E_GPS};

#[derive(Debug, Clone, PartialEq)]
pub struct GpsEphemeris {
    pub sat: SatelliteId,
    pub toe: GpsTime,
    pub toc: GpsTime,
    pub af0: f64,
    pub af1: f64,
    pub af2: f64,
    pub crs: f64,
    pub crc: f64,
    pub cuc: f64,
    pub cus: f64,
    pub cic: f64,
    pub cis: f64,
    pub m0: f64,
    pub e: f64,
    pub sqrt_a: f64,
    pub delta_n: f64,
    pub omega0: f64,
    pub omega_dot: f64,
    pub i0: f64,
    pub idot: f64,
    pub omega: f64,
    pub tgd: f64,
    pub iode: u32,
    pub iodc: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GalileoEphemeris {
    pub sat: SatelliteId,
    pub toe: GpsTime,
    pub toc: GpsTime,
    pub af0: f64,
    pub af1: f64,
    pub af2: f64,
    pub crs: f64,
    pub crc: f64,
    pub cuc: f64,
    pub cus: f64,
    pub cic: f64,
    pub cis: f64,
    pub m0: f64,
    pub e: f64,
    pub sqrt_a: f64,
    pub delta_n: f64,
    pub omega0: f64,
    pub omega_dot: f64,
    pub i0: f64,
    pub idot: f64,
    pub omega: f64,
    pub bgd_e1_e5a: f64,
    pub bgd_e1_e5b: f64,
    pub iod_nav: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BeidouEphemeris {
    pub sat: SatelliteId,
    pub toe: GpsTime,
    pub toc: GpsTime,
    pub af0: f64,
    pub af1: f64,
    pub af2: f64,
    pub crs: f64,
    pub crc: f64,
    pub cuc: f64,
    pub cus: f64,
    pub cic: f64,
    pub cis: f64,
    pub m0: f64,
    pub e: f64,
    pub sqrt_a: f64,
    pub delta_n: f64,
    pub omega0: f64,
    pub omega_dot: f64,
    pub i0: f64,
    pub idot: f64,
    pub omega: f64,
    pub tgd1: f64,
    pub tgd2: f64,
    pub aode: u32,
    pub aodc: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct QzssEphemeris {
    pub sat: SatelliteId,
    pub toe: GpsTime,
    pub toc: GpsTime,
    pub af0: f64,
    pub af1: f64,
    pub af2: f64,
    pub crs: f64,
    pub crc: f64,
    pub cuc: f64,
    pub cus: f64,
    pub cic: f64,
    pub cis: f64,
    pub m0: f64,
    pub e: f64,
    pub sqrt_a: f64,
    pub delta_n: f64,
    pub omega0: f64,
    pub omega_dot: f64,
    pub i0: f64,
    pub idot: f64,
    pub omega: f64,
    pub tgd: f64,
    pub iode: u32,
    pub iodc: u32,
}

impl GpsEphemeris {
    pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(
            t, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n,
            self.omega0, self.omega_dot, self.i0, self.idot, self.omega,
            self.tgd, MU_GPS, OMEGA_E_GPS, false,
        )
    }

    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(
            t, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n,
            self.omega0, self.omega_dot, self.i0, self.idot, self.omega,
            0.0, MU_GPS, OMEGA_E_GPS, false,
        )
    }
}

impl GalileoEphemeris {
    pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(
            t, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n,
            self.omega0, self.omega_dot, self.i0, self.idot, self.omega,
            self.bgd_e1_e5a, MU_GAL, OMEGA_E_GAL, false,
        )
    }

    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(
            t, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n,
            self.omega0, self.omega_dot, self.i0, self.idot, self.omega,
            0.0, MU_GAL, OMEGA_E_GAL, false,
        )
    }

    pub fn position_e5b(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(
            t, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n,
            self.omega0, self.omega_dot, self.i0, self.idot, self.omega,
            self.bgd_e1_e5b, MU_GAL, OMEGA_E_GAL, false,
        )
    }
}

impl BeidouEphemeris {
    pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
        let is_bds_geo = self.sat.prn <= 5 || self.sat.prn >= 59;
        calc_keplerian(
            t_bdt, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n,
            self.omega0, self.omega_dot, self.i0, self.idot, self.omega,
            self.tgd1, MU_BDS, OMEGA_E_BDS, is_bds_geo,
        )
    }

    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
        let is_bds_geo = self.sat.prn <= 5 || self.sat.prn >= 59;
        calc_keplerian(
            t_bdt, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n,
            self.omega0, self.omega_dot, self.i0, self.idot, self.omega,
            0.0, MU_BDS, OMEGA_E_BDS, is_bds_geo,
        )
    }
}

impl QzssEphemeris {
    pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(
            t, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n,
            self.omega0, self.omega_dot, self.i0, self.idot, self.omega,
            self.tgd, MU_GPS, OMEGA_E_GPS, false,
        )
    }

    pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(
            t, self.toe, self.toc, self.af0, self.af1, self.af2,
            self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
            self.m0, self.e, self.sqrt_a, self.delta_n,
            self.omega0, self.omega_dot, self.i0, self.idot, self.omega,
            0.0, MU_GPS, OMEGA_E_GPS, false,
        )
    }
}
