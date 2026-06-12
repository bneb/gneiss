use crate::sat::SatelliteId;
use crate::time::GpsTime;
use nalgebra::Vector3;

const MU_GPS: f64 = 3.986005e14;
const MU_GAL: f64 = 3.986004418e14;
const MU_BDS: f64 = 3.986004418e14;
const MU_GLO: f64 = 3.9860044e14;
const OMEGA_E_GPS: f64 = crate::constants::EARTH_ROTATION_RATE_RAD_S;
const OMEGA_E_GAL: f64 = crate::constants::EARTH_ROTATION_RATE_RAD_S; // Same as GPS
const OMEGA_E_BDS: f64 = 7.292115e-5;
const OMEGA_E_GLO: f64 = 7.292115e-5;
const J2_GLO: f64 = 1.0826257e-3;
const RADIUS_GLO: f64 = 6378136.0;
const F: f64 = -4.442807633e-10; // Relativistic constant

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
}

#[derive(Debug, Clone, PartialEq)]
pub struct GpsEphemeris {
    pub sat: SatelliteId, pub toe: GpsTime, pub toc: GpsTime,
    pub af0: f64, pub af1: f64, pub af2: f64,
    pub crs: f64, pub crc: f64, pub cuc: f64, pub cus: f64, pub cic: f64, pub cis: f64,
    pub m0: f64, pub e: f64, pub sqrt_a: f64, pub delta_n: f64,
    pub omega0: f64, pub omega_dot: f64, pub i0: f64, pub idot: f64, pub omega: f64, pub tgd: f64,
    pub iode: u32, pub iodc: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GalileoEphemeris {
    pub sat: SatelliteId, pub toe: GpsTime, pub toc: GpsTime,
    pub af0: f64, pub af1: f64, pub af2: f64,
    pub crs: f64, pub crc: f64, pub cuc: f64, pub cus: f64, pub cic: f64, pub cis: f64,
    pub m0: f64, pub e: f64, pub sqrt_a: f64, pub delta_n: f64,
    pub omega0: f64, pub omega_dot: f64, pub i0: f64, pub idot: f64, pub omega: f64, pub bgd_e1_e5a: f64,
    pub iod_nav: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BeidouEphemeris {
    pub sat: SatelliteId, pub toe: GpsTime, pub toc: GpsTime,
    pub af0: f64, pub af1: f64, pub af2: f64,
    pub crs: f64, pub crc: f64, pub cuc: f64, pub cus: f64, pub cic: f64, pub cis: f64,
    pub m0: f64, pub e: f64, pub sqrt_a: f64, pub delta_n: f64,
    pub omega0: f64, pub omega_dot: f64, pub i0: f64, pub idot: f64, pub omega: f64, pub tgd1: f64,
    pub aode: u32, pub aodc: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct QzssEphemeris {
    pub sat: SatelliteId, pub toe: GpsTime, pub toc: GpsTime,
    pub af0: f64, pub af1: f64, pub af2: f64,
    pub crs: f64, pub crc: f64, pub cuc: f64, pub cus: f64, pub cic: f64, pub cis: f64,
    pub m0: f64, pub e: f64, pub sqrt_a: f64, pub delta_n: f64,
    pub omega0: f64, pub omega_dot: f64, pub i0: f64, pub idot: f64, pub omega: f64, pub tgd: f64,
    pub iode: u32, pub iodc: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GlonassEphemeris {
    pub sat: SatelliteId, pub toe: GpsTime,
    pub freq_num: i8,
    pub tau_n: f64, pub gamma_n: f64, pub delta_tau_n: f64,
    pub x: f64, pub y: f64, pub z: f64,
    pub vx: f64, pub vy: f64, pub vz: f64,
    pub ax: f64, pub ay: f64, pub az: f64,
}

fn glonass_derivatives(state: &[f64; 6], acc: &[f64; 3]) -> [f64; 6] {
    let r2 = state[0] * state[0] + state[1] * state[1] + state[2] * state[2];
    let r = libm::sqrt(r2);
    let r3 = r2 * r;
    let ae2 = RADIUS_GLO * RADIUS_GLO;
    let factor = 1.5 * J2_GLO * MU_GLO * ae2 / (r2 * r3);
    let z2_r2 = state[2] * state[2] / r2;

    let ax = -MU_GLO * state[0] / r3 - factor * state[0] * (1.0 - 5.0 * z2_r2) + OMEGA_E_GLO * OMEGA_E_GLO * state[0] + 2.0 * OMEGA_E_GLO * state[4] + acc[0];
    let ay = -MU_GLO * state[1] / r3 - factor * state[1] * (1.0 - 5.0 * z2_r2) + OMEGA_E_GLO * OMEGA_E_GLO * state[1] - 2.0 * OMEGA_E_GLO * state[3] + acc[1];
    let az = -MU_GLO * state[2] / r3 - factor * state[2] * (3.0 - 5.0 * z2_r2) + acc[2];

    [state[3], state[4], state[5], ax, ay, az]
}

fn rk4_step(state: &[f64; 6], acc: &[f64; 3], h: f64) -> [f64; 6] {
    let k1 = glonass_derivatives(state, acc);
    
    let mut s2 = [0.0; 6];
    for i in 0..6 { s2[i] = state[i] + 0.5 * h * k1[i]; }
    let k2 = glonass_derivatives(&s2, acc);
    
    let mut s3 = [0.0; 6];
    for i in 0..6 { s3[i] = state[i] + 0.5 * h * k2[i]; }
    let k3 = glonass_derivatives(&s3, acc);
    
    let mut s4 = [0.0; 6];
    for i in 0..6 { s4[i] = state[i] + h * k3[i]; }
    let k4 = glonass_derivatives(&s4, acc);

    let mut next_state = [0.0; 6];
    for i in 0..6 {
        next_state[i] = state[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    next_state
}

impl GlonassEphemeris {
    pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        let dt = t - self.toe;
        let mut state = [self.x, self.y, self.z, self.vx, self.vy, self.vz];
        let acc = [self.ax, self.ay, self.az];
        
        let step = if dt < 0.0 { -30.0 } else { 30.0 };
        let mut t_rem = dt;

        while libm::fabs(t_rem) > 1e-14 {
            let h = if libm::fabs(t_rem) < libm::fabs(step) { t_rem } else { step };
            state = rk4_step(&state, &acc, h);
            t_rem -= h;
        }

        // Relativistic effect is already absorbed into GLONASS tau_n usually, but we compute clock error:
        let clk_err = -self.tau_n + self.gamma_n * dt;
        let clk_drift = self.gamma_n;

        (
            Vector3::new(state[0], state[1], state[2]),
            Vector3::new(state[3], state[4], state[5]),
            clk_err,
            clk_drift
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn calc_keplerian(
    t: GpsTime, toe: GpsTime, toc: GpsTime,
    af0: f64, af1: f64, af2: f64, crs: f64, crc: f64, cuc: f64, cus: f64, cic: f64, cis: f64,
    m0: f64, e: f64, sqrt_a: f64, delta_n: f64, omega0: f64, omega_dot: f64, i0: f64, idot: f64, omega: f64, tgd: f64,
    mu: f64, omega_e: f64, is_bds_geo: bool
) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
    let tk = t - toe;
    let a = sqrt_a * sqrt_a;
    let n0 = libm::sqrt(mu / (a * a * a));
    let n = n0 + delta_n;
    let mk = m0 + n * tk;

    let mut ek = mk;
    for _ in 0..10 {
        ek = mk + e * libm::sin(ek);
    }

    let cos_ek = libm::cos(ek);
    let sin_ek = libm::sin(ek);

    let vk = libm::atan2(libm::sqrt(1.0 - e * e) * sin_ek, cos_ek - e);
    let uk = omega + vk;

    let sin_2uk = libm::sin(2.0 * uk);
    let cos_2uk = libm::cos(2.0 * uk);

    let duk = cus * sin_2uk + cuc * cos_2uk;
    let drk = crs * sin_2uk + crc * cos_2uk;
    let dik = cis * sin_2uk + cic * cos_2uk;

    let u = uk + duk;
    let r = a * (1.0 - e * cos_ek) + drk;
    let i = i0 + dik + idot * tk;

    let xk_prime = r * libm::cos(u);
    let yk_prime = r * libm::sin(u);

    let omegak = if is_bds_geo {
        omega0 + omega_dot * tk - omega_e * toe.tow
    } else {
        omega0 + (omega_dot - omega_e) * tk - omega_e * toe.tow
    };

    let cos_omegak = libm::cos(omegak);
    let sin_omegak = libm::sin(omegak);
    let sin_ik = libm::sin(i);
    let cos_ik = libm::cos(i);

    let x_orb = xk_prime * cos_omegak - yk_prime * cos_ik * sin_omegak;
    let y_orb = xk_prime * sin_omegak + yk_prime * cos_ik * cos_omegak;
    let z_orb = yk_prime * sin_ik;

    let rk_dot = sqrt_a * sqrt_a * e * sin_ek * n / (1.0 - e * cos_ek);
    let uk_dot = (libm::sqrt(1.0 - e * e) / (1.0 - e * cos_ek)) * n;
    
    let xk_prime_dot = rk_dot * libm::cos(u) - r * libm::sin(u) * uk_dot;
    let yk_prime_dot = rk_dot * libm::sin(u) + r * libm::cos(u) * uk_dot;
    
    let omegak_dot = if is_bds_geo { omega_dot } else { omega_dot - omega_e };
    
    let x_orb_dot = xk_prime_dot * cos_omegak - yk_prime_dot * cos_ik * sin_omegak - y_orb * omegak_dot - yk_prime * sin_ik * idot * sin_omegak;
    let y_orb_dot = xk_prime_dot * sin_omegak + yk_prime_dot * cos_ik * cos_omegak + x_orb * omegak_dot + yk_prime * sin_ik * idot * cos_omegak;
    let z_orb_dot = yk_prime_dot * sin_ik + yk_prime * cos_ik * idot;

    let (x, y, z, vx, vy, vz) = if is_bds_geo {
        let sin_5 = libm::sin(-5.0f64.to_radians());
        let cos_5 = libm::cos(-5.0f64.to_radians());
        let sin_oet = libm::sin(omega_e * tk);
        let cos_oet = libm::cos(omega_e * tk);
        
        let xg = x_orb * cos_oet + y_orb * sin_oet * cos_5 + z_orb * sin_oet * sin_5;
        let yg = -x_orb * sin_oet + y_orb * cos_oet * cos_5 + z_orb * cos_oet * sin_5;
        let zg = -y_orb * sin_5 + z_orb * cos_5;

        let vxg = x_orb_dot * cos_oet - x_orb * omega_e * sin_oet 
                  + y_orb_dot * sin_oet * cos_5 + y_orb * omega_e * cos_oet * cos_5 
                  + z_orb_dot * sin_oet * sin_5 + z_orb * omega_e * cos_oet * sin_5;
        
        let vyg = -x_orb_dot * sin_oet - x_orb * omega_e * cos_oet 
                  + y_orb_dot * cos_oet * cos_5 - y_orb * omega_e * sin_oet * cos_5 
                  + z_orb_dot * cos_oet * sin_5 - z_orb * omega_e * sin_oet * sin_5;
        
        let vzg = -y_orb_dot * sin_5 + z_orb_dot * cos_5;
        
        (xg, yg, zg, vxg, vyg, vzg)
    } else {
        (x_orb, y_orb, z_orb, x_orb_dot, y_orb_dot, z_orb_dot)
    };

    let tc = t - toc;
    
    let dt_rel = F * e * sqrt_a * sin_ek;
    let clk_err = af0 + af1 * tc + af2 * tc * tc + dt_rel - tgd;
    let clk_drift = af1 + 2.0 * af2 * tc;

    (Vector3::new(x, y, z), Vector3::new(vx, vy, vz), clk_err, clk_drift)
}

impl GpsEphemeris {
    pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(t, self.toe, self.toc, self.af0, self.af1, self.af2, self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis, self.m0, self.e, self.sqrt_a, self.delta_n, self.omega0, self.omega_dot, self.i0, self.idot, self.omega, self.tgd, MU_GPS, OMEGA_E_GPS, false)
    }
}

impl GalileoEphemeris {
    pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(t, self.toe, self.toc, self.af0, self.af1, self.af2, self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis, self.m0, self.e, self.sqrt_a, self.delta_n, self.omega0, self.omega_dot, self.i0, self.idot, self.omega, self.bgd_e1_e5a, MU_GAL, OMEGA_E_GAL, false)
    }
}

impl BeidouEphemeris {
    pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        // BDT is 14 seconds behind GPS Time (GPST = BDT + 14s).
        // Since t is passed in GPST, we convert it to BDT for keplerian projection.
        let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
        
        let is_bds_geo = self.sat.prn <= 5 || self.sat.prn >= 59;
        
        calc_keplerian(t_bdt, self.toe, self.toc, self.af0, self.af1, self.af2, self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis, self.m0, self.e, self.sqrt_a, self.delta_n, self.omega0, self.omega_dot, self.i0, self.idot, self.omega, self.tgd1, MU_BDS, OMEGA_E_BDS, is_bds_geo)
    }
}

impl QzssEphemeris {
    pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        calc_keplerian(t, self.toe, self.toc, self.af0, self.af1, self.af2, self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis, self.m0, self.e, self.sqrt_a, self.delta_n, self.omega0, self.omega_dot, self.i0, self.idot, self.omega, self.tgd, MU_GPS, OMEGA_E_GPS, false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sat::Constellation;

    #[test]
    fn test_gps_position_calculation() {
        let eph = GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100010.0),
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: -2.0e-9, i0: 0.95, idot: 0.0,
            omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        };

        let (pos, vel, _clk_err, _clk_drift) = eph.position(GpsTime::new(2000, 100000.0));
        assert!(pos.norm() > 20_000_000.0);
        assert!(pos.norm() < 30_000_000.0);
        assert!(vel.norm() > 1000.0); // GPS satellites move at ~3.9 km/s
    }

    #[test]
    fn test_glonass_rk4_physics() {
        // Test GLONASS Cartesian RK4 numerical integrator
        let eph = GlonassEphemeris {
            sat: SatelliteId { constellation: Constellation::Glonass, prn: 1 },
            toe: GpsTime::new(2000, 100000.0),
            freq_num: 1,
            tau_n: 1e-5, gamma_n: 1e-9, delta_tau_n: 0.0,
            x: 10_000_000.0, y: 15_000_000.0, z: 20_000_000.0,
            vx: -2000.0, vy: 1500.0, vz: 1000.0,
            ax: 0.0, ay: 0.0, az: 0.0, // Solar/lunar accels
        };
        
        let (pos, vel, clk_err, _clk_drift) = eph.position(GpsTime::new(2000, 100060.0)); // Propagate 60 seconds
        
        // 60s at roughly 2.5km/s gives about 150km change
        let dist_moved = (pos - Vector3::new(10_000_000.0, 15_000_000.0, 20_000_000.0)).norm();
        assert!(dist_moved > 100_000.0 && dist_moved < 200_000.0);
        assert!(vel.norm() > 1000.0);
        assert!((clk_err - (-1e-5 + 1e-9 * 60.0)).abs() < 1e-12);
    }

    #[test]
    fn test_ephemeris_enum_dispatch() {
        let gal_eph = GalileoEphemeris {
            sat: SatelliteId { constellation: Constellation::Galileo, prn: 2 },
            toe: GpsTime::new(2000, 100000.0), toc: GpsTime::new(2000, 100010.0),
            af0: 0.0, af1: 0.0, af2: 0.0, crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.01, sqrt_a: 5440.6, delta_n: 0.0, omega0: 0.0, omega_dot: -2.0e-9,
            i0: 0.95, idot: 0.0, omega: 0.0, bgd_e1_e5a: 0.0, iod_nav: 1,
        };
        
        let enum_eph = Ephemeris::Galileo(gal_eph.clone());
        assert_eq!(enum_eph.sat().constellation, Constellation::Galileo);
        
        let (pos1, _, _, _) = gal_eph.position(GpsTime::new(2000, 100000.0));
        let (pos2, _, _, _) = enum_eph.position(GpsTime::new(2000, 100000.0));
        assert_eq!(pos1, pos2);
    }

    #[test]
    fn test_ephemeris_exact() {
        let gps_eph = GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            toe: GpsTime::new(2000, 100000.0), toc: GpsTime::new(2000, 100010.0),
            af0: 1.0, af1: 2.0, af2: 3.0, crs: 4.0, crc: 5.0, cuc: 6.0, cus: 7.0, cic: 8.0, cis: 9.0,
            m0: 0.1, e: 0.01, sqrt_a: 5153.6, delta_n: 0.001, omega0: 0.2, omega_dot: -2.0e-9,
            i0: 0.95, idot: 0.01, omega: 0.3, tgd: 0.02, iode: 1, iodc: 1,
        };
        let gal_eph = GalileoEphemeris {
            sat: SatelliteId { constellation: Constellation::Galileo, prn: 2 },
            toe: GpsTime::new(2000, 100000.0), toc: GpsTime::new(2000, 100010.0),
            af0: 1.0, af1: 2.0, af2: 3.0, crs: 4.0, crc: 5.0, cuc: 6.0, cus: 7.0, cic: 8.0, cis: 9.0,
            m0: 0.1, e: 0.01, sqrt_a: 5153.6, delta_n: 0.001, omega0: 0.2, omega_dot: -2.0e-9,
            i0: 0.95, idot: 0.01, omega: 0.3, bgd_e1_e5a: 0.02, iod_nav: 1,
        };
        let bds_geo_eph = BeidouEphemeris {
            sat: SatelliteId { constellation: Constellation::Beidou, prn: 1 }, // Geo
            toe: GpsTime::new(2000, 100000.0), toc: GpsTime::new(2000, 100010.0),
            af0: 1.0, af1: 2.0, af2: 3.0, crs: 4.0, crc: 5.0, cuc: 6.0, cus: 7.0, cic: 8.0, cis: 9.0,
            m0: 0.1, e: 0.01, sqrt_a: 5153.6, delta_n: 0.001, omega0: 0.2, omega_dot: -2.0e-9,
            i0: 0.95, idot: 0.01, omega: 0.3, tgd1: 0.02, aode: 1, aodc: 1,
        };
        let bds_igso_eph = BeidouEphemeris {
            sat: SatelliteId { constellation: Constellation::Beidou, prn: 6 }, // Not Geo
            toe: GpsTime::new(2000, 100000.0), toc: GpsTime::new(2000, 100010.0),
            af0: 1.0, af1: 2.0, af2: 3.0, crs: 4.0, crc: 5.0, cuc: 6.0, cus: 7.0, cic: 8.0, cis: 9.0,
            m0: 0.1, e: 0.01, sqrt_a: 5153.6, delta_n: 0.001, omega0: 0.2, omega_dot: -2.0e-9,
            i0: 0.95, idot: 0.01, omega: 0.3, tgd1: 0.02, aode: 1, aodc: 1,
        };
        let qzss_eph = QzssEphemeris {
            sat: SatelliteId { constellation: Constellation::Qzss, prn: 4 },
            toe: GpsTime::new(2000, 100000.0), toc: GpsTime::new(2000, 100010.0),
            af0: 1.0, af1: 2.0, af2: 3.0, crs: 4.0, crc: 5.0, cuc: 6.0, cus: 7.0, cic: 8.0, cis: 9.0,
            m0: 0.1, e: 0.01, sqrt_a: 5153.6, delta_n: 0.001, omega0: 0.2, omega_dot: -2.0e-9,
            i0: 0.95, idot: 0.01, omega: 0.3, tgd: 0.02, iode: 1, iodc: 1,
        };
        let glo_eph_fwd = GlonassEphemeris {
            sat: SatelliteId { constellation: Constellation::Glonass, prn: 5 },
            toe: GpsTime::new(2000, 100000.0), freq_num: 7,
            tau_n: 1e-5, gamma_n: 1e-9, delta_tau_n: 0.0,
            x: 10_000_000.0, y: 15_000_000.0, z: 20_000_000.0,
            vx: -2000.0, vy: 1500.0, vz: 1000.0, ax: 0.1, ay: 0.2, az: 0.3,
        };
        let glo_eph_bwd = glo_eph_fwd.clone();

        let t = GpsTime::new(2000, 100060.0);
        let e_gps = Ephemeris::Gps(gps_eph.clone());
        let e_gal = Ephemeris::Galileo(gal_eph.clone());
        let e_bds = Ephemeris::Beidou(bds_igso_eph.clone());
        let e_qzss = Ephemeris::Qzss(qzss_eph.clone());
        let e_glo = Ephemeris::Glonass(glo_eph_fwd.clone());
        
        assert_eq!(e_gps.sat(), SatelliteId { constellation: Constellation::Gps, prn: 1 });
        assert_eq!(e_gps.toe(), GpsTime::new(2000, 100000.0));
        assert_eq!(e_gps.freq_num(), 0);

        assert_eq!(e_gal.sat(), SatelliteId { constellation: Constellation::Galileo, prn: 2 });
        assert_eq!(e_gal.toe(), GpsTime::new(2000, 100000.0));
        assert_eq!(e_gal.freq_num(), 0);

        assert_eq!(e_bds.sat(), SatelliteId { constellation: Constellation::Beidou, prn: 6 });
        assert_eq!(e_bds.toe(), GpsTime::new(2000, 100000.0));
        assert_eq!(e_bds.freq_num(), 0);

        assert_eq!(e_qzss.sat(), SatelliteId { constellation: Constellation::Qzss, prn: 4 });
        assert_eq!(e_qzss.toe(), GpsTime::new(2000, 100000.0));
        assert_eq!(e_qzss.freq_num(), 0);

        assert_eq!(e_glo.sat(), SatelliteId { constellation: Constellation::Glonass, prn: 5 });
        assert_eq!(e_glo.toe(), GpsTime::new(2000, 100000.0));
        assert_eq!(e_glo.freq_num(), 7);

        // Calculate positions
        let p_gps = gps_eph.position(t);
        let p_gal = gal_eph.position(t);
        let p_bds_geo = bds_geo_eph.position(t);
        let p_bds_igso = bds_igso_eph.position(t);
        let p_qzss = qzss_eph.position(t);
        let p_glo_fwd = glo_eph_fwd.position(t);
        let p_glo_bwd = glo_eph_bwd.position(GpsTime::new(2000, 99940.0));
        
        // Assert exact values to kill arithmetic mutants
        use nalgebra::Vector3;
        
        let p_gps = gps_eph.position(t);
        let p_gal = gal_eph.position(t);
        let p_bds_geo = bds_geo_eph.position(t);
        let p_bds_igso = bds_igso_eph.position(t);
        let p_qzss = qzss_eph.position(t);
        let p_glo_fwd = glo_eph_fwd.position(t);
        let p_glo_bwd = glo_eph_bwd.position(GpsTime::new(2000, 99940.0));
        
        // macro to assert approx equal for positions to avoid precision issue on diff architectures, but strict enough to catch mutation
        macro_rules! assert_vec_eq {
            ($a:expr, $b:expr) => {
                assert!(($a - $b).norm() < 1e-10, "Vectors not equal: {:?} and {:?}", $a, $b);
            };
        }
        
        assert_vec_eq!(p_gps.0, Vector3::new(-20111686.91105315, 16173621.585368361, -5050844.794698619));
        assert_vec_eq!(p_gps.1, Vector3::new(-42915.9866782879, -49950.87879784308, -59327.05664687349));
        extern crate std;
        std::println!("p_gps.2 = {:.15}", p_gps.2);
        assert!((p_gps.2 - 7600.979999996116).abs() < 1e-12);
        assert!((p_gps.3 - 302.0).abs() < 1e-12);

        assert_vec_eq!(p_gal.0, Vector3::new(-20111686.90739729, 16173621.585960612, -5050844.807207882));
        assert_vec_eq!(p_gal.1, Vector3::new(-42915.98666471951, -49950.87870889442, -59327.056414597435));
        assert!((p_gal.2 - 7600.979999996116).abs() < 1e-12);
        assert!((p_gal.3 - 302.0).abs() < 1e-12);

        assert_vec_eq!(p_bds_geo.0, Vector3::new(-20530859.400528494, 16091202.889367301, -3331282.3000197206));
        assert_vec_eq!(p_bds_geo.1, Vector3::new(-42143.353117060804, -43745.76008367346, -67353.04008649733));
        assert!((p_bds_geo.2 - 3960.9799999964825).abs() < 1e-12);
        assert!((p_bds_geo.3 - 218.0).abs() < 1e-12);

        assert_vec_eq!(p_bds_igso.0, Vector3::new(-20532037.822592802, 15739895.005709056, -4715036.370416497));
        assert_vec_eq!(p_bds_igso.1, Vector3::new(-42188.08236021536, -49443.13151989598, -63140.92184141916));
        assert!((p_bds_igso.2 - 3960.9799999964825).abs() < 1e-12);
        assert!((p_bds_igso.3 - 218.0).abs() < 1e-12);

        assert_vec_eq!(p_qzss.0, Vector3::new(-20111686.91105315, 16173621.585368361, -5050844.794698619));
        assert_vec_eq!(p_qzss.1, Vector3::new(-42915.9866782879, -49950.87879784308, -59327.05664687349));
        assert!((p_qzss.2 - 7600.979999996116).abs() < 1e-12);
        assert!((p_qzss.3 - 302.0).abs() < 1e-12);

        assert_vec_eq!(p_glo_fwd.0, Vector3::new(9880305.169245299, 15090476.713825395, 20059805.549791936));
        assert_vec_eq!(p_glo_fwd.1, Vector3::new(-1989.7750855530473, 1515.879262912734, 993.5292177774875));
        assert!((p_glo_fwd.2 - -9.940000000000001e-6).abs() < 1e-12);
        assert!((p_glo_fwd.3 - 1e-9).abs() < 1e-12);

        assert_vec_eq!(p_glo_bwd.0, Vector3::new(10120298.841760032, 14910478.053419847, 19939804.282132257));
        assert_vec_eq!(p_glo_bwd.1, Vector3::new(-2009.9085434547917, 1484.0537579365496, 1006.5341630852408));
        assert!((p_glo_bwd.2 - -1.006e-5).abs() < 1e-12);
        assert!((p_glo_bwd.3 - 1e-9).abs() < 1e-12);
    }


    #[test]
    fn test_bds_boundaries() {
        let mut eph = BeidouEphemeris {
            sat: SatelliteId { constellation: Constellation::Beidou, prn: 5 },
            toe: GpsTime::new(2000, 100000.0), toc: GpsTime::new(2000, 100010.0),
            af0: 1.0, af1: 2.0, af2: 3.0, crs: 4.0, crc: 5.0, cuc: 6.0, cus: 7.0, cic: 8.0, cis: 9.0,
            m0: 0.1, e: 0.01, sqrt_a: 5153.6, delta_n: 0.001, omega0: 0.2, omega_dot: -2.0e-9,
            i0: 0.95, idot: 0.01, omega: 0.3, tgd1: 0.02, aode: 1, aodc: 1,
        };
        let t = GpsTime::new(2000, 100060.0);
        let p5 = eph.position(t); // geo
        
        eph.sat.prn = 59;
        let p59 = eph.position(t); // geo
        
        eph.sat.prn = 6;
        let p6 = eph.position(t); // igso
        
        // Ensure 5 and 59 behave like geo, 6 like igso
        assert_eq!(p5.0, p59.0); // Wait, position will be identical because prn is only used for geo check
        assert_ne!(p5.0, p6.0);
    }


    #[test]
    fn test_glonass_partial_step() {
        let glo_eph = GlonassEphemeris {
            sat: SatelliteId { constellation: Constellation::Glonass, prn: 5 },
            toe: GpsTime::new(2000, 100000.0), freq_num: 7,
            tau_n: 1e-5, gamma_n: 1e-9, delta_tau_n: 0.0,
            x: 10_000_000.0, y: 15_000_000.0, z: 20_000_000.0,
            vx: -2000.0, vy: 1500.0, vz: 1000.0, ax: 0.1, ay: 0.2, az: 0.3,
        };
        let t = GpsTime::new(2000, 100050.0);
        let p = glo_eph.position(t);
        // We just assert on exactly what it produces so it locks it in.
        
        let expected_x = 9900211.557692498;
        let expected_y = 15075331.129011473;
        let expected_z = 20049864.88968714;
        
        // Use a tiny delta to lock it in
        assert!((p.0.x - expected_x).abs() < 1e-4);
        assert!((p.0.y - expected_y).abs() < 1e-4);
        assert!((p.0.z - expected_z).abs() < 1e-4);
    }

}
