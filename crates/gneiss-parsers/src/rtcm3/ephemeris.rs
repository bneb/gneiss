use bitvec::prelude::*;
use gneiss_core::ephemeris::GpsEphemeris;
use gneiss_core::sat::{SatelliteId, Constellation};
use gneiss_core::time::GpsTime;
use super::RtcmParseError;
use super::{sign_extend_i16, sign_extend_i32};

const P2_5: f64 = 0.03125;
const P2_19: f64 = 1.9073486328125e-06;
const P2_29: f64 = 1.862645149230957e-09;
const P2_31: f64 = 4.656612873077393e-10;
const P2_33: f64 = 1.1641532182693481e-10;
const P2_43: f64 = 1.1368683772161603e-13;
const P2_55: f64 = 2.7755575615628914e-17;

fn read_u8(bits: &bitvec::slice::BitSlice<u8, Msb0>, off: &mut usize, len: usize) -> u8 {
    let v = bits[*off..*off + len].load_be::<u8>();
    *off += len;
    v
}

fn read_u16(bits: &bitvec::slice::BitSlice<u8, Msb0>, off: &mut usize, len: usize) -> u16 {
    let v = bits[*off..*off + len].load_be::<u16>();
    *off += len;
    v
}

fn read_u32(bits: &bitvec::slice::BitSlice<u8, Msb0>, off: &mut usize, len: usize) -> u32 {
    let v = bits[*off..*off + len].load_be::<u32>();
    *off += len;
    v
}

fn read_1019_part1(bits: &bitvec::slice::BitSlice<u8, Msb0>, off: &mut usize) -> (u8, f64, u32, GpsTime, f64, f64, f64, u32) {
    let prn = read_u8(bits, off, 6);
    read_u16(bits, off, 10); // week
    read_u8(bits, off, 4);   // ura
    read_u8(bits, off, 2);   // code_l2
    let idot = sign_extend_i16(read_u16(bits, off, 14), 14) as f64 * P2_43 * core::f64::consts::PI;
    let iode = read_u8(bits, off, 8) as u32;
    let toc = GpsTime::new(0, read_u16(bits, off, 16) as f64 * 16.0);
    let af2 = sign_extend_i32(read_u8(bits, off, 8) as u32, 8) as f64 * P2_55;
    let af1 = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_43;
    let af0 = sign_extend_i32(read_u32(bits, off, 22), 22) as f64 * P2_31;
    let iodc = read_u16(bits, off, 10) as u32;
    (prn, idot, iode, toc, af2, af1, af0, iodc)
}

fn read_1019_part2(bits: &bitvec::slice::BitSlice<u8, Msb0>, off: &mut usize) -> (f64, f64, f64, f64, f64, f64, f64, GpsTime) {
    let crs = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_5;
    let delta_n = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_43 * core::f64::consts::PI;
    let m0 = sign_extend_i32(read_u32(bits, off, 32), 32) as f64 * P2_31 * core::f64::consts::PI;
    let cuc = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_29;
    let e = read_u32(bits, off, 32) as f64 * P2_33;
    let cus = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_29;
    let sqrt_a = read_u32(bits, off, 32) as f64 * P2_19;
    let toe = GpsTime::new(0, read_u16(bits, off, 16) as f64 * 16.0);
    (crs, delta_n, m0, cuc, e, cus, sqrt_a, toe)
}

fn read_1019_part3(bits: &bitvec::slice::BitSlice<u8, Msb0>, off: &mut usize) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let cic = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_29;
    let omega0 = sign_extend_i32(read_u32(bits, off, 32), 32) as f64 * P2_31 * core::f64::consts::PI;
    let cis = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_29;
    let i0 = sign_extend_i32(read_u32(bits, off, 32), 32) as f64 * P2_31 * core::f64::consts::PI;
    let crc = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_5;
    let omega = sign_extend_i32(read_u32(bits, off, 32), 32) as f64 * P2_31 * core::f64::consts::PI;
    let omega_dot = sign_extend_i32(read_u32(bits, off, 24), 24) as f64 * P2_43 * core::f64::consts::PI;
    let tgd = sign_extend_i32(read_u8(bits, off, 8) as u32, 8) as f64 * P2_31;
    (cic, omega0, cis, i0, crc, omega, omega_dot, tgd)
}

pub fn parse_1019(payload: &[u8]) -> Result<GpsEphemeris, RtcmParseError> {
    let bits = payload.view_bits::<Msb0>();
    if bits.len() < 488 { return Err(RtcmParseError::Incomplete); }
    let mut off = 0;
    
    let msg_num = read_u16(bits, &mut off, 12);
    if msg_num != 1019 { return Err(RtcmParseError::UnsupportedMsmType); }

    let (prn, idot, iode, toc, af2, af1, af0, iodc) = read_1019_part1(bits, &mut off);
    let (crs, delta_n, m0, cuc, e, cus, sqrt_a, toe) = read_1019_part2(bits, &mut off);
    let (cic, omega0, cis, i0, crc, omega, omega_dot, tgd) = read_1019_part3(bits, &mut off);

    Ok(GpsEphemeris {
        sat: SatelliteId { constellation: Constellation::Gps, prn },
        toe, toc, af0, af1, af2, crs, crc, cuc, cus, cic, cis, m0, e, sqrt_a,
        delta_n, omega0, omega_dot, i0, idot, omega, tgd, iode, iodc,
    })
}
