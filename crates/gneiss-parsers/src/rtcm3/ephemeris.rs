use bitvec::prelude::*;
use gneiss_core::ephemeris::GpsEphemeris;
use gneiss_core::sat::{SatelliteId, Constellation};
use gneiss_core::time::GpsTime;
use super::RtcmParseError;
use super::{sign_extend_i16, sign_extend_i32};

pub fn parse_1019(payload: &[u8]) -> Result<GpsEphemeris, RtcmParseError> {
    let bits = payload.view_bits::<Msb0>();
    if bits.len() < 488 {
        return Err(RtcmParseError::Incomplete);
    }

    let mut offset = 0;
    
    macro_rules! next_u8 {
        ($len:expr) => {{
            let v = bits[offset..offset + $len].load_be::<u8>();
            offset += $len;
            v
        }}
    }
    macro_rules! next_u16 {
        ($len:expr) => {{
            let v = bits[offset..offset + $len].load_be::<u16>();
            offset += $len;
            v
        }}
    }
    macro_rules! next_u32 {
        ($len:expr) => {{
            let v = bits[offset..offset + $len].load_be::<u32>();
            offset += $len;
            v
        }}
    }

    let msg_num = next_u16!(12);
    if msg_num != 1019 {
        return Err(RtcmParseError::UnsupportedMsmType);
    }

    let prn = next_u8!(6);
    let _week = next_u16!(10);
    let _ura = next_u8!(4);
    let _code_l2 = next_u8!(2);

    let idot_raw = next_u16!(14);
    let idot = sign_extend_i16(idot_raw, 14) as f64 * f64::powf(2.0, -43.0) * core::f64::consts::PI;

    let iode = next_u8!(8) as u32;

    let toc_raw = next_u16!(16);
    let toc = GpsTime::new(0, toc_raw as f64 * 16.0); // We ignore full week resolution here for simplicity

    let af2_raw = next_u8!(8);
    let af2 = sign_extend_i32(af2_raw as u32, 8) as f64 * f64::powf(2.0, -55.0);

    let af1_raw = next_u16!(16);
    let af1 = sign_extend_i16(af1_raw, 16) as f64 * f64::powf(2.0, -43.0);

    let af0_raw = next_u32!(22);
    let af0 = sign_extend_i32(af0_raw, 22) as f64 * f64::powf(2.0, -31.0);

    let iodc = next_u16!(10) as u32;

    let crs_raw = next_u16!(16);
    let crs = sign_extend_i16(crs_raw, 16) as f64 * f64::powf(2.0, -5.0);

    let delta_n_raw = next_u16!(16);
    let delta_n = sign_extend_i16(delta_n_raw, 16) as f64 * f64::powf(2.0, -43.0) * core::f64::consts::PI;

    let m0_raw = next_u32!(32);
    let m0 = sign_extend_i32(m0_raw, 32) as f64 * f64::powf(2.0, -31.0) * core::f64::consts::PI;

    let cuc_raw = next_u16!(16);
    let cuc = sign_extend_i16(cuc_raw, 16) as f64 * f64::powf(2.0, -29.0);

    let e_raw = next_u32!(32);
    let e = e_raw as f64 * f64::powf(2.0, -33.0);

    let cus_raw = next_u16!(16);
    let cus = sign_extend_i16(cus_raw, 16) as f64 * f64::powf(2.0, -29.0);

    let sqrt_a_raw = next_u32!(32);
    let sqrt_a = sqrt_a_raw as f64 * f64::powf(2.0, -19.0);

    let toe_raw = next_u16!(16);
    let toe = GpsTime::new(0, toe_raw as f64 * 16.0);

    let cic_raw = next_u16!(16);
    let cic = sign_extend_i16(cic_raw, 16) as f64 * f64::powf(2.0, -29.0);

    let omega0_raw = next_u32!(32);
    let omega0 = sign_extend_i32(omega0_raw, 32) as f64 * f64::powf(2.0, -31.0) * core::f64::consts::PI;

    let cis_raw = next_u16!(16);
    let cis = sign_extend_i16(cis_raw, 16) as f64 * f64::powf(2.0, -29.0);

    let i0_raw = next_u32!(32);
    let i0 = sign_extend_i32(i0_raw, 32) as f64 * f64::powf(2.0, -31.0) * core::f64::consts::PI;

    let crc_raw = next_u16!(16);
    let crc = sign_extend_i16(crc_raw, 16) as f64 * f64::powf(2.0, -5.0);

    let omega_raw = next_u32!(32);
    let omega = sign_extend_i32(omega_raw, 32) as f64 * f64::powf(2.0, -31.0) * core::f64::consts::PI;

    let omega_dot_raw = next_u32!(24);
    let omega_dot = sign_extend_i32(omega_dot_raw, 24) as f64 * f64::powf(2.0, -43.0) * core::f64::consts::PI;

    let tgd_raw = next_u8!(8);
    let tgd = sign_extend_i32(tgd_raw as u32, 8) as f64 * f64::powf(2.0, -31.0);

    Ok(GpsEphemeris {
        sat: SatelliteId {
            constellation: Constellation::Gps,
            prn,
        },
        toe,
        toc,
        af0,
        af1,
        af2,
        crs,
        crc,
        cuc,
        cus,
        cic,
        cis,
        m0,
        e,
        sqrt_a,
        delta_n,
        omega0,
        omega_dot,
        i0,
        idot,
        omega,
        tgd,
        iode,
        iodc,
    })
}
