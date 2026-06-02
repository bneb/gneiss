
/// A raw RTCM3 frame containing the verified payload (message type and data).
#[derive(Debug, Clone, PartialEq)]
pub struct RtcmFrame<'a> {
    pub payload: &'a [u8],
}

#[derive(Debug, Clone, PartialEq)]
pub enum RtcmParseError {
    Incomplete,
    InvalidPreamble,
    CrcMismatch,
    UnsupportedMsmType,
}

const CRC24_POLY: u32 = 0x864CFB;

/// Computes the CRC-24Q checksum for RTCM3 validation.
pub fn crc24q(data: &[u8]) -> u32 {
    let mut crc: u32 = 0;
    for &byte in data {
        crc ^= (byte as u32) << 16;
        for _ in 0..8 {
            if (crc & 0x800000) != 0 {
                crc = (crc << 1) ^ CRC24_POLY;
            } else {
                crc <<= 1;
            }
        }
    }
    crc & 0xFFFFFF
}

/// Parses a single RTCM3 frame from the input stream.
/// Returns the remaining bytes and the valid frame.
pub fn parse_rtcm3_frame(input: &[u8]) -> Result<(&[u8], RtcmFrame<'_>), RtcmParseError> {
    if input.is_empty() {
        return Err(RtcmParseError::Incomplete);
    }

    if input[0] != 0xD3 {
        return Err(RtcmParseError::InvalidPreamble);
    }

    if input.len() < 3 {
        return Err(RtcmParseError::Incomplete);
    }

    let len_bits = u16::from_be_bytes([input[1], input[2]]);
    let payload_len = (len_bits & 0x03FF) as usize;

    let frame_len = 3 + payload_len + 3;

    if input.len() < frame_len {
        return Err(RtcmParseError::Incomplete);
    }

    let expected_crc = crc24q(&input[0..(3 + payload_len)]);
    
    let crc_bytes = &input[(3 + payload_len)..frame_len];
    let actual_crc = ((crc_bytes[0] as u32) << 16) | ((crc_bytes[1] as u32) << 8) | (crc_bytes[2] as u32);

    if expected_crc != actual_crc {
        return Err(RtcmParseError::CrcMismatch);
    }

    let payload = &input[3..(3 + payload_len)];
    let remaining = &input[frame_len..];

    Ok((remaining, RtcmFrame { payload }))
}

pub fn sign_extend_i16(value: u16, bits: u32) -> i16 {
    let shift = 16 - bits;
    (value << shift) as i16 >> shift
}

pub fn sign_extend_i32(value: u32, bits: u32) -> i32 {
    let shift = 32 - bits;
    (value << shift) as i32 >> shift
}

pub mod ephemeris;
pub mod msm;

pub use ephemeris::*;
pub use msm::*;