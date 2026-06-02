use bitvec::prelude::*;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::sat::{SatelliteId, Constellation};
use gneiss_core::time::GpsTime;
use super::{RtcmParseError, sign_extend_i16, sign_extend_i32};

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
            _ => None, // MSM1, 2, 3 are considered legacy/unsupported for modern RTK
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

impl MsmMessage {
    pub fn into_epoch_obs(&self) -> EpochObs {
        let time = GpsTime::new(0, self.header.epoch_time as f64 / 1000.0); // Simple time mapping
        
        let constellation = match self.header.message_number / 10 {
            107 => Constellation::Gps,
            108 => Constellation::Glonass,
            109 => Constellation::Galileo,
            110 => Constellation::Sbas,
            111 => Constellation::Qzss,
            112 => Constellation::Beidou,
            _ => Constellation::Gps, // Default fallback
        };

        let mut satellites = Vec::new();
        let mut sat_idx = 0;
        let _cell_idx = 0;

        for i in 0..64 {
            // Check if satellite 'i' is present in the mask
            if (self.masks.satellite_mask & (1 << (63 - i))) != 0 {
                let prn = (i + 1) as u8;
                let sat = SatelliteId {
                    constellation,
                    prn,
                };

                let observations = Vec::new();
                let _rough_range = self.satellite_data.rough_ranges[sat_idx] as f64; // Stub: need proper scaling

                for j in 0..32 {
                    // Check if signal 'j' is present in the signal mask
                    if (self.masks.signal_mask & (1 << (31 - j))) != 0 {
                        // Check if this specific cell (sat i, signal j) is active
                        let _cell_offset = sat_idx * self.masks.signal_mask.count_ones() as usize + j;
                        // wait, cell_mask is sparsely populated only for true signal bits.
                        // Actually, cell_mask has size = N_sat * N_sig.
                        // We need to track the active signal index
                        
                        // We will simplify this loop structure to just linearly read cells.
                    }
                }

                satellites.push(SatObs { sat, observations });
                sat_idx += 1;
            }
        }

        EpochObs { time, satellites }
    }
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

/// Parses the 73-bit common MSM header from the bit stream.
pub fn parse_msm_header(bits: &BitSlice<u8, Msb0>) -> Result<(&BitSlice<u8, Msb0>, MsmHeader), RtcmParseError> {
    if bits.len() < 73 {
        return Err(RtcmParseError::Incomplete);
    }
    
    let message_number = bits[0..12].load_be::<u16>();
    let station_id = bits[12..24].load_be::<u16>();
    let epoch_time = bits[24..54].load_be::<u32>();
    let multiple_message = bits[54];
    let iods = bits[55..58].load_be::<u8>();
    // bits 58..65 are reserved (7 bits)
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

/// Masks determining which satellites, signals, and specific cells are present in the MSM message.
#[derive(Debug, Clone, PartialEq)]
pub struct MsmMasks {
    pub satellite_mask: u64,
    pub signal_mask: u32,
    /// Boolean mask indicating which cells (satellite + signal combination) are present.
    /// Length is equal to `n_sat * n_sig`.
    pub cell_mask: Vec<bool>,
}

/// Parses the Satellite, Signal, and Cell masks.
pub fn parse_msm_masks(bits: &BitSlice<u8, Msb0>) -> Result<(&BitSlice<u8, Msb0>, MsmMasks), RtcmParseError> {
    if bits.len() < 96 { // 64 (sat) + 32 (sig)
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

/// Data provided for each active satellite.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MsmSatelliteData {
    pub rough_ranges: Vec<u16>,             // 10 bits
    pub extended_sat_info: Vec<u8>,         // 4 bits (MSM4+)
    pub rough_phase_range_rates: Vec<i16>,  // 14 bits (MSM5, 7 only)
}

/// Parses the Satellite Data Section.
pub fn parse_satellite_data(
    bits: &BitSlice<u8, Msb0>,
    n_sat: usize,
    msm_type: MsmType,
) -> Result<(&BitSlice<u8, Msb0>, MsmSatelliteData), RtcmParseError> {
    let mut offset = 0;
    let mut data = MsmSatelliteData::default();

    // 1. Rough Ranges (10 bits)
    let req_bits = n_sat * 10;
    if bits.len() < offset + req_bits { return Err(RtcmParseError::Incomplete); }
    for _ in 0..n_sat {
        data.rough_ranges.push(bits[offset..offset+10].load_be::<u16>());
        offset += 10;
    }

    // 2. Extended Sat Info (4 bits) - MSM4, 5, 6, 7
    if matches!(msm_type, MsmType::Msm4 | MsmType::Msm5 | MsmType::Msm6 | MsmType::Msm7) {
        let req_bits = n_sat * 4;
        if bits.len() < offset + req_bits { return Err(RtcmParseError::Incomplete); }
        for _ in 0..n_sat {
            data.extended_sat_info.push(bits[offset..offset+4].load_be::<u8>());
            offset += 4;
        }
    }

    // 3. Rough PhaseRange Rates (14 bits) - MSM5, 7
    if matches!(msm_type, MsmType::Msm5 | MsmType::Msm7) {
        let req_bits = n_sat * 14;
        if bits.len() < offset + req_bits { return Err(RtcmParseError::Incomplete); }
        for _ in 0..n_sat {
            let val = bits[offset..offset+14].load_be::<u16>();
            data.rough_phase_range_rates.push(sign_extend_i16(val, 14));
            offset += 14;
        }
    }

    Ok((&bits[offset..], data))
}

/// Data provided for each active signal cell.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MsmSignalData {
    pub fine_pseudoranges: Vec<i32>,        // 15b/20b
    pub fine_phase_ranges: Vec<i32>,        // 22b/24b
    pub lock_time_indicators: Vec<u16>,     // 4b/10b
    pub half_cycle_ambiguities: Vec<bool>,  // 1b
    pub cnrs: Vec<u16>,                     // 6b/10b
    pub fine_phase_range_rates: Vec<i16>,   // 15b (MSM5, 7 only)
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

    // Fine Pseudoranges
    let req_bits = n_cell * pr_bits;
    if bits.len() < offset + req_bits { return Err(RtcmParseError::Incomplete); }
    for _ in 0..n_cell {
        let val = bits[offset..offset+pr_bits].load_be::<u32>();
        data.fine_pseudoranges.push(sign_extend_i32(val, pr_bits as u32));
        offset += pr_bits;
    }

    // Fine PhaseRanges
    let req_bits = n_cell * ph_bits;
    if bits.len() < offset + req_bits { return Err(RtcmParseError::Incomplete); }
    for _ in 0..n_cell {
        let val = bits[offset..offset+ph_bits].load_be::<u32>();
        data.fine_phase_ranges.push(sign_extend_i32(val, ph_bits as u32));
        offset += ph_bits;
    }

    // Lock Time Indicators
    let req_bits = n_cell * lock_bits;
    if bits.len() < offset + req_bits { return Err(RtcmParseError::Incomplete); }
    for _ in 0..n_cell {
        data.lock_time_indicators.push(bits[offset..offset+lock_bits].load_be::<u16>());
        offset += lock_bits;
    }

    // Half-cycle Ambiguities (1 bit)
    let req_bits = n_cell;
    if bits.len() < offset + req_bits { return Err(RtcmParseError::Incomplete); }
    for _ in 0..n_cell {
        data.half_cycle_ambiguities.push(bits[offset]);
        offset += 1;
    }

    // CNRs
    let req_bits = n_cell * cnr_bits;
    if bits.len() < offset + req_bits { return Err(RtcmParseError::Incomplete); }
    for _ in 0..n_cell {
        data.cnrs.push(bits[offset..offset+cnr_bits].load_be::<u16>());
        offset += cnr_bits;
    }

    // Fine PhaseRangeRates (15 bits, signed) - MSM5/7
    if matches!(msm_type, MsmType::Msm5 | MsmType::Msm7) {
        let req_bits = n_cell * 15;
        if bits.len() < offset + req_bits { return Err(RtcmParseError::Incomplete); }
        for _ in 0..n_cell {
            let val = bits[offset..offset+15].load_be::<u16>();
            data.fine_phase_range_rates.push(sign_extend_i16(val, 15));
            offset += 15;
        }
    }

    Ok((&bits[offset..], data))
}

