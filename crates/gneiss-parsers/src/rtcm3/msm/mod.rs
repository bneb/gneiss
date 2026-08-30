//! RTCM 3 Multiple Signal Message (MSM4, MSM5, MSM6, MSM7) decoding.

pub mod decoder;
pub mod signals;
#[cfg(test)]
mod tests;

use super::{sign_extend_i16, sign_extend_i32, RtcmParseError};
use bitvec::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MsmType {
    Msm4,
    Msm5,
    Msm6,
    Msm7,
}

impl MsmType {
    /// Infers the MSM type from the message number (e.g., 1074 -> Msm4)
    pub fn from_message_number(num: u16) -> Option<Self> {
        match num % 10 {
            4 => Some(MsmType::Msm4),
            5 => Some(MsmType::Msm5),
            6 => Some(MsmType::Msm6),
            7 => Some(MsmType::Msm7),
            _ => None,
        }
    }
}

/// A complete, decoded Multiple Signal Message (MSM).
#[derive(Debug, Clone, PartialEq)]
pub struct MsmMessage {
    pub msm_type: MsmType,
    pub header: MsmHeader,
    pub masks: MsmMasks,
    pub satellite_data: MsmSatelliteData,
    pub signal_data: MsmSignalData,
}

/// Common header for all MSM messages (MSM1 - MSM7)
#[derive(Debug, Clone, PartialEq)]
pub struct MsmHeader {
    pub message_number: u16,
    pub station_id: u16,
    pub epoch_time: u32,
    pub multiple_message: bool,
    pub iods: u8,
    pub clock_steering: u8,
    pub external_clock: u8,
    pub smoothing_indicator: bool,
    pub smoothing_interval: u8,
}

/// Masks determining which satellites, signals, and specific cells are present in the MSM message.
#[derive(Debug, Clone, PartialEq)]
pub struct MsmMasks {
    pub satellite_mask: u64,
    pub signal_mask: u32,
    pub cell_mask: Vec<bool>,
}

/// Data provided for each active satellite.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MsmSatelliteData {
    pub rough_range_int_ms: Vec<u8>,
    pub extended_sat_info: Vec<u8>,
    pub rough_ranges: Vec<u16>,
    pub rough_phase_range_rates: Vec<i16>,
}

/// Data provided for each active signal cell.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MsmSignalData {
    pub fine_pseudoranges: Vec<i32>,
    pub fine_phase_ranges: Vec<i32>,
    pub lock_time_indicators: Vec<u16>,
    pub half_cycle_ambiguities: Vec<bool>,
    pub cnrs: Vec<u16>,
    pub fine_phase_range_rates: Vec<i16>,
}

/// Parses a complete MSM message from the raw RTCM3 payload bytes.
pub fn parse_msm_message(payload: &[u8]) -> Result<MsmMessage, RtcmParseError> {
    let bits = payload.view_bits::<Msb0>();
    let (bits, header) = parse_msm_header(bits)?;

    let msm_type = MsmType::from_message_number(header.message_number)
        .ok_or(RtcmParseError::UnsupportedMsmType)?;

    let (bits, masks) = parse_msm_masks(bits)?;

    let n_sat = masks.satellite_mask.count_ones() as usize;
    let (bits, satellite_data) = parse_satellite_data(bits, n_sat, msm_type)?;

    let n_cell = masks.cell_mask.iter().filter(|&&b| b).count();
    let (_bits, signal_data) = parse_signal_data(bits, n_cell, msm_type)?;

    Ok(MsmMessage {
        msm_type,
        header,
        masks,
        satellite_data,
        signal_data,
    })
}

/// Parses the 73-bit common MSM header from the bit stream.
pub fn parse_msm_header(
    bits: &BitSlice<u8, Msb0>,
) -> Result<(&BitSlice<u8, Msb0>, MsmHeader), RtcmParseError> {
    if bits.len() < 73 {
        return Err(RtcmParseError::Incomplete);
    }

    let message_number = bits[0..12].load_be::<u16>();
    let station_id = bits[12..24].load_be::<u16>();
    let epoch_time = bits[24..54].load_be::<u32>();
    let multiple_message = bits[54];
    let iods = bits[55..58].load_be::<u8>();
    let clock_steering = bits[65..67].load_be::<u8>();
    let external_clock = bits[67..69].load_be::<u8>();
    let smoothing_indicator = bits[69];
    let smoothing_interval = bits[70..73].load_be::<u8>();

    let header = MsmHeader {
        message_number,
        station_id,
        epoch_time,
        multiple_message,
        iods,
        clock_steering,
        external_clock,
        smoothing_indicator,
        smoothing_interval,
    };

    Ok((&bits[73..], header))
}

/// Parses the Satellite, Signal, and Cell masks.
pub fn parse_msm_masks(
    bits: &BitSlice<u8, Msb0>,
) -> Result<(&BitSlice<u8, Msb0>, MsmMasks), RtcmParseError> {
    if bits.len() < 96 {
        return Err(RtcmParseError::Incomplete);
    }

    let satellite_mask = bits[0..64].load_be::<u64>();
    let signal_mask = bits[64..96].load_be::<u32>();

    let n_sat = satellite_mask.count_ones() as usize;
    let n_sig = signal_mask.count_ones() as usize;
    let num_cells = n_sat * n_sig;

    let end_of_masks = 96 + num_cells;
    if bits.len() < end_of_masks {
        return Err(RtcmParseError::Incomplete);
    }

    let mut cell_mask = Vec::with_capacity(num_cells);
    for i in 0..num_cells {
        cell_mask.push(bits[96 + i]);
    }

    let masks = MsmMasks {
        satellite_mask,
        signal_mask,
        cell_mask,
    };

    Ok((&bits[end_of_masks..], masks))
}

/// Parses the Satellite Data Section.
pub fn parse_satellite_data(
    bits: &BitSlice<u8, Msb0>,
    n_sat: usize,
    msm_type: MsmType,
) -> Result<(&BitSlice<u8, Msb0>, MsmSatelliteData), RtcmParseError> {
    let mut offset = 0;
    let mut data = MsmSatelliteData::default();

    let req_bits = n_sat * 8;
    if bits.len() < offset + req_bits {
        return Err(RtcmParseError::Incomplete);
    }
    for _ in 0..n_sat {
        data.rough_range_int_ms
            .push(bits[offset..offset + 8].load_be::<u8>());
        offset += 8;
    }

    if matches!(msm_type, MsmType::Msm5 | MsmType::Msm7) {
        let req_bits = n_sat * 4;
        if bits.len() < offset + req_bits {
            return Err(RtcmParseError::Incomplete);
        }
        for _ in 0..n_sat {
            data.extended_sat_info
                .push(bits[offset..offset + 4].load_be::<u8>());
            offset += 4;
        }
    }

    let req_bits = n_sat * 10;
    if bits.len() < offset + req_bits {
        return Err(RtcmParseError::Incomplete);
    }
    for _ in 0..n_sat {
        data.rough_ranges
            .push(bits[offset..offset + 10].load_be::<u16>());
        offset += 10;
    }

    if matches!(msm_type, MsmType::Msm5 | MsmType::Msm7) {
        let req_bits = n_sat * 14;
        if bits.len() < offset + req_bits {
            return Err(RtcmParseError::Incomplete);
        }
        for _ in 0..n_sat {
            let val = bits[offset..offset + 14].load_be::<u16>();
            data.rough_phase_range_rates.push(sign_extend_i16(val, 14));
            offset += 14;
        }
    }

    Ok((&bits[offset..], data))
}

/// Parses the Signal Data Section.
pub fn parse_signal_data(
    bits: &BitSlice<u8, Msb0>,
    n_cell: usize,
    msm_type: MsmType,
) -> Result<(&BitSlice<u8, Msb0>, MsmSignalData), RtcmParseError> {
    let mut offset = 0;
    let mut data = MsmSignalData::default();

    let (pr_bits, ph_bits, lock_bits, cnr_bits) = match msm_type {
        MsmType::Msm4 | MsmType::Msm5 => (15, 22, 4, 6),
        MsmType::Msm6 | MsmType::Msm7 => (20, 24, 10, 10),
    };

    let req_bits = n_cell * pr_bits;
    if bits.len() < offset + req_bits {
        return Err(RtcmParseError::Incomplete);
    }
    for _ in 0..n_cell {
        let val = bits[offset..offset + pr_bits].load_be::<u32>();
        data.fine_pseudoranges
            .push(sign_extend_i32(val, pr_bits as u32));
        offset += pr_bits;
    }

    let req_bits = n_cell * ph_bits;
    if bits.len() < offset + req_bits {
        return Err(RtcmParseError::Incomplete);
    }
    for _ in 0..n_cell {
        let val = bits[offset..offset + ph_bits].load_be::<u32>();
        data.fine_phase_ranges
            .push(sign_extend_i32(val, ph_bits as u32));
        offset += ph_bits;
    }

    let req_bits = n_cell * lock_bits;
    if bits.len() < offset + req_bits {
        return Err(RtcmParseError::Incomplete);
    }
    for _ in 0..n_cell {
        data.lock_time_indicators
            .push(bits[offset..offset + lock_bits].load_be::<u16>());
        offset += lock_bits;
    }

    let req_bits = n_cell;
    if bits.len() < offset + req_bits {
        return Err(RtcmParseError::Incomplete);
    }
    for _ in 0..n_cell {
        data.half_cycle_ambiguities.push(bits[offset]);
        offset += 1;
    }

    let req_bits = n_cell * cnr_bits;
    if bits.len() < offset + req_bits {
        return Err(RtcmParseError::Incomplete);
    }
    for _ in 0..n_cell {
        data.cnrs
            .push(bits[offset..offset + cnr_bits].load_be::<u16>());
        offset += cnr_bits;
    }

    if matches!(msm_type, MsmType::Msm5 | MsmType::Msm7) {
        let req_bits = n_cell * 15;
        if bits.len() < offset + req_bits {
            return Err(RtcmParseError::Incomplete);
        }
        for _ in 0..n_cell {
            let val = bits[offset..offset + 15].load_be::<u16>();
            data.fine_phase_range_rates.push(sign_extend_i16(val, 15));
            offset += 15;
        }
    }

    Ok((&bits[offset..], data))
}
