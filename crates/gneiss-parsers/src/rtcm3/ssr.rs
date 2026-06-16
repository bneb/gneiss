use bitvec::prelude::*;
use super::{RtcmParseError, sign_extend_i32};

const RES_ORBIT_RADIAL: f64 = 0.0001;
const RES_ORBIT_TRACK: f64 = 0.0004;
const RES_ORBIT_DOT_RADIAL: f64 = 0.000001;
const RES_ORBIT_DOT_TRACK: f64 = 0.000004;

const RES_CLOCK_C0: f64 = 0.0001;
const RES_CLOCK_C1: f64 = 0.000001;
const RES_CLOCK_C2: f64 = 0.00000002;

const RES_CODE_BIAS: f64 = 0.01;

#[derive(Debug, Clone, PartialEq)]
pub struct SsrHeader {
    pub message_number: u16,
    pub epoch_time: u32,
    pub update_interval: u8,
    pub multiple_message_indicator: bool,
    pub satellite_reference_datum: bool,
    pub iod_ssr: u8,
    pub provider_id: u16,
    pub solution_id: u8,
    pub num_satellites: u8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SsrOrbitSat {
    pub sat_id: u8,
    pub iode: u8,
    pub delta_radial: f64,
    pub delta_along_track: f64,
    pub delta_cross_track: f64,
    pub dot_delta_radial: f64,
    pub dot_delta_along_track: f64,
    pub dot_delta_cross_track: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SsrClockSat {
    pub sat_id: u8,
    pub delta_clock_c0: f64,
    pub delta_clock_c1: f64,
    pub delta_clock_c2: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SsrCodeBiasSat {
    pub sat_id: u8,
    pub num_biases: u8,
    pub biases: Vec<SsrSignalBias>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SsrSignalBias {
    pub signal_and_tracking_mode: u8,
    pub bias: f64,
}

pub fn parse_ssr_header(bits: &BitSlice<u8, Msb0>) -> Result<(&BitSlice<u8, Msb0>, SsrHeader), RtcmParseError> {
    if bits.len() < 68 {
        return Err(RtcmParseError::Incomplete);
    }

    let header = SsrHeader {
        message_number: bits[0..12].load_be::<u16>(),
        epoch_time: bits[12..32].load_be::<u32>(),
        update_interval: bits[32..36].load_be::<u8>(),
        multiple_message_indicator: bits[36],
        satellite_reference_datum: bits[37],
        iod_ssr: bits[38..42].load_be::<u8>(),
        provider_id: bits[42..58].load_be::<u16>(),
        solution_id: bits[58..62].load_be::<u8>(),
        num_satellites: bits[62..68].load_be::<u8>(),
    };

    Ok((&bits[68..], header))
}

pub fn parse_ssr_orbit(payload: &[u8]) -> Result<(SsrHeader, Vec<SsrOrbitSat>), RtcmParseError> {
    let bits = payload.view_bits::<Msb0>();
    let (mut bits, header) = parse_ssr_header(bits)?;

    let mut sats = Vec::with_capacity(header.num_satellites as usize);
    for _ in 0..header.num_satellites {
        let (next_bits, sat) = parse_orbit_sat(bits)?;
        sats.push(sat);
        bits = next_bits;
    }

    Ok((header, sats))
}

fn parse_orbit_sat(bits: &BitSlice<u8, Msb0>) -> Result<(&BitSlice<u8, Msb0>, SsrOrbitSat), RtcmParseError> {
    if bits.len() < 135 {
        return Err(RtcmParseError::Incomplete);
    }

    let sat = SsrOrbitSat {
        sat_id: bits[0..6].load_be::<u8>(),
        iode: bits[6..14].load_be::<u8>(),
        delta_radial: sign_extend_i32(bits[14..36].load_be::<u32>(), 22) as f64 * RES_ORBIT_RADIAL,
        delta_along_track: sign_extend_i32(bits[36..56].load_be::<u32>(), 20) as f64 * RES_ORBIT_TRACK,
        delta_cross_track: sign_extend_i32(bits[56..76].load_be::<u32>(), 20) as f64 * RES_ORBIT_TRACK,
        dot_delta_radial: sign_extend_i32(bits[76..97].load_be::<u32>(), 21) as f64 * RES_ORBIT_DOT_RADIAL,
        dot_delta_along_track: sign_extend_i32(bits[97..116].load_be::<u32>(), 19) as f64 * RES_ORBIT_DOT_TRACK,
        dot_delta_cross_track: sign_extend_i32(bits[116..135].load_be::<u32>(), 19) as f64 * RES_ORBIT_DOT_TRACK,
    };

    Ok((&bits[135..], sat))
}

pub fn parse_ssr_clock(payload: &[u8]) -> Result<(SsrHeader, Vec<SsrClockSat>), RtcmParseError> {
    let bits = payload.view_bits::<Msb0>();
    let (mut bits, header) = parse_ssr_header(bits)?;

    let mut sats = Vec::with_capacity(header.num_satellites as usize);
    for _ in 0..header.num_satellites {
        let (next_bits, sat) = parse_clock_sat(bits)?;
        sats.push(sat);
        bits = next_bits;
    }

    Ok((header, sats))
}

fn parse_clock_sat(bits: &BitSlice<u8, Msb0>) -> Result<(&BitSlice<u8, Msb0>, SsrClockSat), RtcmParseError> {
    if bits.len() < 76 {
        return Err(RtcmParseError::Incomplete);
    }

    let sat = SsrClockSat {
        sat_id: bits[0..6].load_be::<u8>(),
        delta_clock_c0: sign_extend_i32(bits[6..28].load_be::<u32>(), 22) as f64 * RES_CLOCK_C0,
        delta_clock_c1: sign_extend_i32(bits[28..49].load_be::<u32>(), 21) as f64 * RES_CLOCK_C1,
        delta_clock_c2: sign_extend_i32(bits[49..76].load_be::<u32>(), 27) as f64 * RES_CLOCK_C2,
    };

    Ok((&bits[76..], sat))
}

pub fn parse_ssr_code_bias(payload: &[u8]) -> Result<(SsrHeader, Vec<SsrCodeBiasSat>), RtcmParseError> {
    let bits = payload.view_bits::<Msb0>();
    let (mut bits, header) = parse_ssr_header(bits)?;

    let mut sats = Vec::with_capacity(header.num_satellites as usize);
    for _ in 0..header.num_satellites {
        let (next_bits, sat) = parse_code_bias_sat(bits)?;
        sats.push(sat);
        bits = next_bits;
    }

    Ok((header, sats))
}

fn parse_code_bias_sat(bits: &BitSlice<u8, Msb0>) -> Result<(&BitSlice<u8, Msb0>, SsrCodeBiasSat), RtcmParseError> {
    if bits.len() < 11 {
        return Err(RtcmParseError::Incomplete);
    }

    let sat_id = bits[0..6].load_be::<u8>();
    let num_biases = bits[6..11].load_be::<u8>();
    let mut current_bits = &bits[11..];

    let mut biases = Vec::with_capacity(num_biases as usize);
    for _ in 0..num_biases {
        if current_bits.len() < 19 {
            return Err(RtcmParseError::Incomplete);
        }
        biases.push(SsrSignalBias {
            signal_and_tracking_mode: current_bits[0..5].load_be::<u8>(),
            bias: sign_extend_i32(current_bits[5..19].load_be::<u32>(), 14) as f64 * RES_CODE_BIAS,
        });
        current_bits = &current_bits[19..];
    }

    Ok((current_bits, SsrCodeBiasSat { sat_id, num_biases, biases }))
}
