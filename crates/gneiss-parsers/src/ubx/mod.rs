use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

/// A raw UBX frame.
#[derive(Debug, Clone, PartialEq)]
pub struct UbxFrame<'a> {
    pub class: u8,
    pub id: u8,
    pub payload: &'a [u8],
}

#[derive(Debug, Clone, PartialEq)]
pub enum UbxParseError {
    Incomplete,
    InvalidSync,
    ChecksumMismatch,
    InvalidLength,
}

/// Computes the UBX 8-bit Fletcher checksum (RFC1145).
pub fn ubx_checksum(data: &[u8]) -> (u8, u8) {
    let mut a: u8 = 0;
    let mut b: u8 = 0;
    for &byte in data {
        a = a.wrapping_add(byte);
        b = b.wrapping_add(a);
    }
    (a, b)
}

/// Parses a single UBX frame.
pub fn parse_ubx_frame(input: &[u8]) -> Result<(&[u8], UbxFrame<'_>), UbxParseError> {
    if input.len() < 8 {
        return Err(UbxParseError::Incomplete);
    }

    if input[0] != 0xB5 || input[1] != 0x62 {
        return Err(UbxParseError::InvalidSync);
    }

    let class = input[2];
    let id = input[3];
    let payload_len = u16::from_le_bytes([input[4], input[5]]) as usize;

    if input.len() < 6 + payload_len + 2 {
        return Err(UbxParseError::Incomplete);
    }

    let (ck_a, ck_b) = ubx_checksum(&input[2..(6 + payload_len)]);
    let actual_a = input[6 + payload_len];
    let actual_b = input[7 + payload_len];

    if ck_a != actual_a || ck_b != actual_b {
        return Err(UbxParseError::ChecksumMismatch);
    }

    let payload = &input[6..(6 + payload_len)];
    let remaining = &input[(6 + payload_len + 2)..];

    Ok((remaining, UbxFrame { class, id, payload }))
}

/// Individual measurement block in UBX-RXM-RAWX
#[derive(Debug, Clone, PartialEq)]
pub struct RxmRawxMeas {
    pub pr_mes: f64,
    pub cp_mes: f64,
    pub do_mes: f32,
    pub gnss_id: u8,
    pub sv_id: u8,
    pub sig_id: u8,
    pub freq_id: u8,
    pub locktime: u16,
    pub cno: u8,
    pub pr_stdev: f64,
    pub cp_stdev: f64,
    pub do_stdev: f64,
    pub pr_valid: bool,
    pub cp_valid: bool,
    pub half_cycle_valid: bool,
    pub sub_half_cycle: bool,
}

/// Decoded UBX-RXM-SFRBX message (Broadcast Navigation Data Subframe)
#[derive(Debug, Clone, PartialEq)]
pub struct UbxRxmSfrbx {
    pub gnss_id: u8,
    pub sv_id: u8,
    pub sig_id: u8,
    pub freq_id: u8,
    pub num_words: u8,
    pub chn: u8,
    pub version: u8,
    pub words: Vec<u32>,
}

/// Parses the payload of a UBX-RXM-SFRBX message.
pub fn parse_rxm_sfrbx(payload: &[u8]) -> Result<UbxRxmSfrbx, UbxParseError> {
    if payload.len() < 8 {
        return Err(UbxParseError::InvalidLength);
    }

    let gnss_id = payload[0];
    let sv_id = payload[1];
    let sig_id = payload[2];
    let freq_id = payload[3];
    let num_words = payload[4];
    let chn = payload[5];
    let version = payload[6];
    // payload[7] is reserved

    let expected_len = 8 + (num_words as usize) * 4;
    if payload.len() != expected_len {
        return Err(UbxParseError::InvalidLength);
    }

    let mut words = Vec::with_capacity(num_words as usize);
    for i in 0..num_words as usize {
        let offset = 8 + i * 4;
        let word = u32::from_le_bytes(
            payload[offset..offset + 4]
                .try_into()
                .map_err(|_| UbxParseError::InvalidLength)?,
        );
        words.push(word);
    }

    Ok(UbxRxmSfrbx {
        gnss_id,
        sv_id,
        sig_id,
        freq_id,
        num_words,
        chn,
        version,
        words,
    })
}
#[derive(Debug, Clone, PartialEq)]
pub struct UbxRxmRawx {
    pub rcv_tow: f64,
    pub week: u16,
    pub leap_s: i8,
    pub num_meas: u8,
    pub rec_stat: u8,
    pub version: u8,
    pub measurements: Vec<RxmRawxMeas>,
}

impl UbxRxmRawx {
    pub fn into_epoch_obs(&self) -> EpochObs {
        let time = GpsTime::new(self.week as u32, self.rcv_tow);
        use std::collections::HashMap;
        let mut sat_map: HashMap<SatelliteId, Vec<Observation>> = HashMap::new();

        for meas in &self.measurements {
            let constellation = match meas.gnss_id {
                0 => Constellation::Gps,
                1 => Constellation::Sbas,
                2 => Constellation::Galileo,
                3 => Constellation::Beidou,
                5 => Constellation::Qzss,
                6 => Constellation::Glonass,
                _ => continue, // Unknown constellation, skip
            };

            let sat = SatelliteId {
                constellation,
                prn: meas.sv_id,
            };

            let observations = sat_map.entry(sat).or_default();

            let freq_band = if meas.sig_id == 0 { 1 } else { 2 };
            let attribute = 'C';

            if meas.pr_valid {
                observations.push(Observation {
                    code: ObsCode {
                        obs_type: ObsType::Pseudorange,
                        signal: SignalCode {
                            freq_band,
                            attribute,
                        },
                    },
                    value: meas.pr_mes,
                    lock_time: None,
                    lli: None,
                });
            }

            if meas.cp_valid && meas.half_cycle_valid {
                // Emit Carrier Phase in pure cycles. The engine will scale to meters using the correct satellite frequency.
                observations.push(Observation {
                    code: ObsCode {
                        obs_type: ObsType::CarrierPhase,
                        signal: SignalCode {
                            freq_band,
                            attribute,
                        },
                    },
                    value: meas.cp_mes,
                    lock_time: Some(meas.locktime),
                    lli: None,
                });
            }

            observations.push(Observation {
                code: ObsCode {
                    obs_type: ObsType::Doppler,
                    signal: SignalCode {
                        freq_band,
                        attribute,
                    },
                },
                value: meas.do_mes as f64,
                lock_time: None,
                lli: None,
            });

            observations.push(Observation {
                code: ObsCode {
                    obs_type: ObsType::Snr,
                    signal: SignalCode {
                        freq_band,
                        attribute,
                    },
                },
                value: meas.cno as f64,
                lock_time: None,
                lli: None,
            });
        }

        let mut satellites = Vec::new();
        for (sat, observations) in sat_map {
            satellites.push(SatObs { sat, observations });
        }

        EpochObs { time, satellites }
    }
}

/// Parses the payload of a UBX-RXM-RAWX message.
pub fn parse_rxm_rawx(payload: &[u8]) -> Result<UbxRxmRawx, UbxParseError> {
    if payload.len() < 16 {
        return Err(UbxParseError::InvalidLength);
    }

    let num_meas = payload[11];
    let expected_len = 16 + (num_meas as usize) * 32;
    if payload.len() != expected_len {
        return Err(UbxParseError::InvalidLength);
    }

    let mut measurements = Vec::with_capacity(num_meas as usize);
    for i in 0..num_meas as usize {
        measurements.push(parse_rxm_rawx_block(
            &payload[16 + i * 32..16 + i * 32 + 32],
        )?);
    }

    Ok(UbxRxmRawx {
        rcv_tow: f64::from_le_bytes(
            payload[0..8]
                .try_into()
                .map_err(|_| UbxParseError::InvalidLength)?,
        ),
        week: u16::from_le_bytes(
            payload[8..10]
                .try_into()
                .map_err(|_| UbxParseError::InvalidLength)?,
        ),
        leap_s: payload[10] as i8,
        num_meas,
        rec_stat: payload[12],
        version: payload[13],
        measurements,
    })
}

fn parse_rxm_rawx_block(block: &[u8]) -> Result<RxmRawxMeas, UbxParseError> {
    let trk_stat = block[30];
    Ok(RxmRawxMeas {
        pr_mes: f64::from_le_bytes(
            block[0..8]
                .try_into()
                .map_err(|_| UbxParseError::InvalidLength)?,
        ),
        cp_mes: f64::from_le_bytes(
            block[8..16]
                .try_into()
                .map_err(|_| UbxParseError::InvalidLength)?,
        ),
        do_mes: f32::from_le_bytes(
            block[16..20]
                .try_into()
                .map_err(|_| UbxParseError::InvalidLength)?,
        ),
        gnss_id: block[20],
        sv_id: block[21],
        sig_id: block[22],
        freq_id: block[23],
        locktime: u16::from_le_bytes(
            block[24..26]
                .try_into()
                .map_err(|_| UbxParseError::InvalidLength)?,
        ),
        cno: block[26],
        pr_stdev: 0.01 * f64::powi(2.0, block[27] as i32),
        cp_stdev: 0.004 * f64::powi(2.0, block[28] as i32),
        do_stdev: 0.002 * f64::powi(2.0, block[29] as i32),
        pr_valid: (trk_stat & 0x01) != 0,
        cp_valid: (trk_stat & 0x02) != 0,
        half_cycle_valid: (trk_stat & 0x04) != 0,
        sub_half_cycle: (trk_stat & 0x08) != 0,
    })
}

/// Struct representing a single sensor measurement inside UBX-ESF-MEAS.
#[derive(Debug, Clone, PartialEq)]
pub struct EsfMeasData {
    /// Sensor data (24-bit signed, sign-extended to 32-bit)
    pub data: i32,
    /// Type of the data (4 bits, e.g. 5 = z-axis gyro, 14 = x-axis accel)
    pub data_type: u8,
}

impl EsfMeasData {
    /// Returns the correctly scaled sensor value in SI units (m/s^2 for accel, rad/s for gyro).
    /// Returns the raw value for unhandled types.
    pub fn scaled_value(&self) -> f64 {
        let val = self.data as f64;
        match self.data_type {
            5 | 13 | 14 => {
                // Gyroscope: 2^-12 deg/s -> convert to rad/s
                val * 2.0f64.powi(-12) * (std::f64::consts::PI / 180.0)
            }
            16..=18 => {
                // Accelerometer: 2^-10 m/s^2
                val * 2.0f64.powi(-10)
            }
            _ => val,
        }
    }
}

/// Raw UBX-ESF-MEAS (External Sensor Fusion Measurements) message.
#[derive(Debug, Clone, PartialEq)]
pub struct UbxEsfMeas {
    /// Time tag of the measurement (typically ms or ticks).
    pub time_tag: u32,
    /// Flags indicating time mark availability and calibration status.
    pub flags: u16,
    /// Identification number of the data provider.
    pub id: u16,
    /// List of measurements in this packet.
    pub measurements: Vec<EsfMeasData>,
}

/// Parses UBX-ESF-MEAS payload (Class 0x10, ID 0x02)
pub fn parse_esf_meas(payload: &[u8]) -> Result<UbxEsfMeas, UbxParseError> {
    if payload.len() < 8 {
        return Err(UbxParseError::InvalidLength);
    }

    let time_tag = u32::from_le_bytes(
        payload[0..4]
            .try_into()
            .map_err(|_| UbxParseError::InvalidLength)?,
    );
    let flags = u16::from_le_bytes(
        payload[4..6]
            .try_into()
            .map_err(|_| UbxParseError::InvalidLength)?,
    );
    let id = u16::from_le_bytes(
        payload[6..8]
            .try_into()
            .map_err(|_| UbxParseError::InvalidLength)?,
    );

    let num_bytes = payload.len() - 8;
    if !num_bytes.is_multiple_of(4) {
        return Err(UbxParseError::InvalidLength);
    }
    let num_meas = num_bytes / 4;

    let mut measurements = Vec::with_capacity(num_meas);
    for i in 0..num_meas {
        let offset = 8 + i * 4;
        let raw = u32::from_le_bytes(
            payload[offset..offset + 4]
                .try_into()
                .map_err(|_| UbxParseError::InvalidLength)?,
        );

        let mut data = raw & 0x00FFFFFF;
        if (data & 0x00800000) != 0 {
            data |= 0xFF000000;
        }
        let data_i32 = data as i32;
        let data_type = (raw >> 24) as u8;

        measurements.push(EsfMeasData {
            data: data_i32,
            data_type,
        });
    }

    Ok(UbxEsfMeas {
        time_tag,
        flags,
        id,
        measurements,
    })
}

/// Raw UBX-ESF-STATUS (External Sensor Fusion Status) message.
#[derive(Debug, Clone, PartialEq)]
pub struct UbxEsfStatus {
    /// Time tag of the measurement (typically ms or ticks).
    pub time_tag: u32,
    /// Fusion mode (0=Init, 1=Fusion, etc.)
    pub fusion_mode: u8,
    pub num_sensors: u8,
}

/// Parses UBX-ESF-STATUS payload (Class 0x10, ID 0x10)
pub fn parse_esf_status(payload: &[u8]) -> Result<UbxEsfStatus, UbxParseError> {
    if payload.len() < 16 {
        return Err(UbxParseError::InvalidLength);
    }

    let time_tag = u32::from_le_bytes(
        payload[0..4]
            .try_into()
            .map_err(|_| UbxParseError::InvalidLength)?,
    );
    let fusion_mode = payload[5];
    let num_sensors = payload[11];

    Ok(UbxEsfStatus {
        time_tag,
        fusion_mode,
        num_sensors,
    })
}


#[cfg(test)]
mod tests;
