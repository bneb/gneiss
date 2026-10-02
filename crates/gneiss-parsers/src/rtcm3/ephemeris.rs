//! RTCM 3 message type 1019 — GPS Ephemeris Data.
//!
//! Field offsets below are payload-relative (bit 0 = first payload bit, i.e. the
//! most significant bit of the first payload byte) and were derived by hand from
//! the reference decoder in RTKLIB `rtcm3.c` `decode_type1019`, which indexes the
//! same fields from the start of the frame at bit 24 (= payload bit 0):
//!
//! ```text
//! offset  width  field            offset  width  field
//!      0     12  message number     192     16  C_uc
//!     12      6  PRN                208     32  e
//!     18     10  GPS week           240     16  C_us
//!     28      4  SV accuracy        256     32  sqrt(A)
//!     32      2  L2 code flag       288     16  T_oe
//!     34     14  IDOT               304     16  C_ic
//!     48      8  IODE               320     32  Omega_0
//!     56     16  T_oc               352     16  C_is
//!     72      8  a_f2               368     32  i_0
//!     80     16  a_f1               400     16  C_rc
//!     96     22  a_f0               416     32  omega
//!    118     10  IODC               448     24  Omega_dot
//!    128     16  C_rs               472      8  T_gd
//!    144     16  delta_n            480      8  SV health (6), L2 flag (1), fit (1)
//!    160     32  M_0
//! ```

use super::RtcmParseError;
use super::{sign_extend_i16, sign_extend_i32};
use bitvec::prelude::*;
use gneiss_core::ephemeris::GpsEphemeris;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

const P2_5: f64 = 0.03125;
const P2_19: f64 = 1.9073486328125e-06;
const P2_29: f64 = 1.862645149230957e-09;
const P2_31: f64 = 4.656612873077393e-10;
const P2_33: f64 = 1.1641532182693481e-10;
const P2_43: f64 = 1.1368683772161603e-13;
const P2_55: f64 = 2.7755575615628914e-17;

/// Payload size of message 1019: 12 bits of message number + 476 data bits.
/// RTKLIB requires 512 frame bits, i.e. 512 - 24 header bits = 488 payload bits.
pub const TYPE_1019_PAYLOAD_BITS: usize = 488;

/// GPS week counter modulus. The broadcast week field is 10 bits wide, so it
/// counts modulo 1024 weeks and must be resolved against a reference epoch.
const GPS_WEEK_MODULUS: i64 = 1024;

/// Resolves the truncated 10-bit broadcast week to a full GPS week number.
///
/// `raw_week` is a 10-bit field; the returned week is the unique value `w` with
/// `w mod 1024 == raw_week` that lies closest to `reference_week` (ties break
/// forward). This is the same discontinuity handling RTKLIB performs in
/// `adjgpsweek()`, but stateless: the caller supplies the reference instead of
/// relying on a compile-time constant.
pub fn adjust_gps_week(raw_week: u32, reference_week: u32) -> u32 {
    let raw = (raw_week as i64) & 0x3FF;
    let delta = reference_week as i64 - raw;
    let cycles = (delta + GPS_WEEK_MODULUS / 2).div_euclid(GPS_WEEK_MODULUS);
    (raw + cycles * GPS_WEEK_MODULUS).max(0) as u32
}

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

/// Parses a 1019 payload, resolving the 10-bit week against a reference epoch.
pub fn parse_1019_at(payload: &[u8], reference: GpsTime) -> Result<GpsEphemeris, RtcmParseError> {
    let bits = payload.view_bits::<Msb0>();
    if bits.len() < TYPE_1019_PAYLOAD_BITS {
        return Err(RtcmParseError::Incomplete);
    }
    let mut off = 0;

    let msg_num = read_u16(bits, &mut off, 12);
    if msg_num != 1019 {
        return Err(RtcmParseError::UnsupportedMsmType);
    }

    let p1 = read_1019_part1(bits, &mut off);
    let p2 = read_1019_part2(bits, &mut off);
    let p3 = read_1019_part3(bits, &mut off);

    let week = adjust_gps_week(p1.week, reference.week);

    Ok(GpsEphemeris {
        sat: SatelliteId {
            constellation: Constellation::Gps,
            prn: p1.prn,
        },
        toe: GpsTime::new(week, p2.toe_s),
        toc: GpsTime::new(week, p1.toc_s),
        af0: p1.af0,
        af1: p1.af1,
        af2: p1.af2,
        crs: p2.crs,
        crc: p3.crc,
        cuc: p2.cuc,
        cus: p2.cus,
        cic: p3.cic,
        cis: p3.cis,
        m0: p2.m0,
        e: p2.e,
        sqrt_a: p2.sqrt_a,
        delta_n: p2.delta_n,
        omega0: p3.omega0,
        omega_dot: p3.omega_dot,
        i0: p3.i0,
        idot: p1.idot,
        omega: p3.omega,
        tgd: p3.tgd,
        iode: p1.iode,
        iodc: p1.iodc,
    })
}

/// Parses a 1019 payload using GPS week 2048 (2019-04-06) as the rollover
/// reference, i.e. the week of the first LNAV rollover *after* the 2019 handover.
///
/// Callers that know the observation time should prefer [`parse_1019_at`], whose
/// reference is exact instead of assumed.
pub fn parse_1019(payload: &[u8]) -> Result<GpsEphemeris, RtcmParseError> {
    parse_1019_at(payload, GpsTime::new(2048, 0.0))
}

struct Part1 {
    prn: u8,
    week: u32,
    idot: f64,
    iode: u32,
    toc_s: f64,
    af2: f64,
    af1: f64,
    af0: f64,
    iodc: u32,
}

fn read_1019_part1(bits: &bitvec::slice::BitSlice<u8, Msb0>, off: &mut usize) -> Part1 {
    let prn = read_u8(bits, off, 6);
    let week = read_u16(bits, off, 10) as u32;
    read_u8(bits, off, 4); // SV accuracy
    read_u8(bits, off, 2); // L2 code flag
    let idot = sign_extend_i16(read_u16(bits, off, 14), 14) as f64 * P2_43 * core::f64::consts::PI;
    let iode = read_u8(bits, off, 8) as u32;
    let toc_s = read_u16(bits, off, 16) as f64 * 16.0;
    let af2 = sign_extend_i32(read_u8(bits, off, 8) as u32, 8) as f64 * P2_55;
    let af1 = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_43;
    let af0 = sign_extend_i32(read_u32(bits, off, 22), 22) as f64 * P2_31;
    let iodc = read_u16(bits, off, 10) as u32;
    Part1 { prn, week, idot, iode, toc_s, af2, af1, af0, iodc }
}

struct Part2 {
    crs: f64,
    delta_n: f64,
    m0: f64,
    cuc: f64,
    e: f64,
    cus: f64,
    sqrt_a: f64,
    toe_s: f64,
}

fn read_1019_part2(bits: &bitvec::slice::BitSlice<u8, Msb0>, off: &mut usize) -> Part2 {
    let crs = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_5;
    let delta_n =
        sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_43 * core::f64::consts::PI;
    let m0 = sign_extend_i32(read_u32(bits, off, 32), 32) as f64 * P2_31 * core::f64::consts::PI;
    let cuc = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_29;
    let e = read_u32(bits, off, 32) as f64 * P2_33;
    let cus = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_29;
    let sqrt_a = read_u32(bits, off, 32) as f64 * P2_19;
    let toe_s = read_u16(bits, off, 16) as f64 * 16.0;
    Part2 { crs, delta_n, m0, cuc, e, cus, sqrt_a, toe_s }
}

struct Part3 {
    cic: f64,
    omega0: f64,
    cis: f64,
    i0: f64,
    crc: f64,
    omega: f64,
    omega_dot: f64,
    tgd: f64,
}

fn read_1019_part3(bits: &bitvec::slice::BitSlice<u8, Msb0>, off: &mut usize) -> Part3 {
    let cic = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_29;
    let omega0 =
        sign_extend_i32(read_u32(bits, off, 32), 32) as f64 * P2_31 * core::f64::consts::PI;
    let cis = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_29;
    let i0 = sign_extend_i32(read_u32(bits, off, 32), 32) as f64 * P2_31 * core::f64::consts::PI;
    let crc = sign_extend_i16(read_u16(bits, off, 16), 16) as f64 * P2_5;
    let omega = sign_extend_i32(read_u32(bits, off, 32), 32) as f64 * P2_31 * core::f64::consts::PI;
    let omega_dot =
        sign_extend_i32(read_u32(bits, off, 24), 24) as f64 * P2_43 * core::f64::consts::PI;
    let tgd = sign_extend_i32(read_u8(bits, off, 8) as u32, 8) as f64 * P2_31;
    Part3 { cic, omega0, cis, i0, crc, omega, omega_dot, tgd }
}

#[cfg(test)]
mod tests;
