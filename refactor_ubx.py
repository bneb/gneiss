import sys

with open('crates/gneiss-parsers/src/ubx.rs', 'r') as f:
    content = f.read()

# Replace into_epoch_obs and parse_rxm_rawx
# Find UbxRxmRawx impl block start
impl_start = content.find('impl UbxRxmRawx {')
impl_end = content.find('pub fn parse_rxm_rawx')

if impl_start == -1 or impl_end == -1:
    print("Could not find impl block")
    sys.exit(1)

parse_rawx_start = content.find('pub fn parse_rxm_rawx')
parse_esf_meas_start = content.find('/// Struct representing a single sensor measurement inside UBX-ESF-MEAS.')

pre_content = content[:impl_start]
post_content = content[parse_esf_meas_start:]

new_content = """impl UbxRxmRawx {
    pub fn into_epoch_obs(&self) -> EpochObs {
        let time = GpsTime::new(self.week as u32, self.rcv_tow);
        use std::collections::HashMap;
        let mut sat_map: HashMap<SatelliteId, Vec<Observation>> = HashMap::new();

        for meas in &self.measurements {
            if let Some((sat, mut new_obs)) = Self::convert_meas(meas) {
                sat_map.entry(sat).or_default().append(&mut new_obs);
            }
        }

        let satellites = sat_map.into_iter().map(|(sat, observations)| SatObs { sat, observations }).collect();
        EpochObs { time, satellites }
    }

    fn convert_meas(meas: &RxmRawxMeas) -> Option<(SatelliteId, Vec<Observation>)> {
        let constellation = match meas.gnss_id {
            0 => Constellation::Gps,
            1 => Constellation::Sbas,
            2 => Constellation::Galileo,
            3 => Constellation::Beidou,
            5 => Constellation::Qzss,
            6 => Constellation::Glonass,
            _ => return None,
        };

        let sat = SatelliteId { constellation, prn: meas.sv_id };
        let freq_band = if meas.sig_id == 0 { 1 } else { 2 };
        let attribute = 'C';
        let mut observations = Vec::new();

        if meas.pr_valid {
            observations.push(Observation {
                code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band, attribute } },
                value: meas.pr_mes,
                lock_time: None,
            });
        }

        if meas.cp_valid {
            observations.push(Observation {
                code: ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band, attribute } },
                value: meas.cp_mes,
                lock_time: Some(meas.locktime),
            });
        }

        observations.push(Observation {
            code: ObsCode { obs_type: ObsType::Doppler, signal: SignalCode { freq_band, attribute } },
            value: meas.do_mes as f64,
            lock_time: None,
        });

        observations.push(Observation {
            code: ObsCode { obs_type: ObsType::Snr, signal: SignalCode { freq_band, attribute } },
            value: meas.cno as f64,
            lock_time: None,
        });

        Some((sat, observations))
    }
}

/// Parses the payload of a UBX-RXM-RAWX message.
pub fn parse_rxm_rawx(payload: &[u8]) -> Result<UbxRxmRawx, UbxParseError> {
    if payload.len() < 16 {
        return Err(UbxParseError::InvalidLength);
    }

    let rcv_tow = f64::from_le_bytes(payload[0..8].try_into().map_err(|_| UbxParseError::InvalidLength)?);
    let week = u16::from_le_bytes(payload[8..10].try_into().map_err(|_| UbxParseError::InvalidLength)?);
    let leap_s = payload[10] as i8;
    let num_meas = payload[11];
    let rec_stat = payload[12];
    let version = payload[13];

    let expected_len = 16 + (num_meas as usize) * 32;
    if payload.len() != expected_len {
        return Err(UbxParseError::InvalidLength);
    }

    let mut measurements = Vec::with_capacity(num_meas as usize);
    for i in 0..num_meas as usize {
        let offset = 16 + i * 32;
        let block = &payload[offset..offset + 32];
        measurements.push(parse_rxm_rawx_meas(block)?);
    }

    Ok(UbxRxmRawx { rcv_tow, week, leap_s, num_meas, rec_stat, version, measurements })
}

fn parse_rxm_rawx_meas(block: &[u8]) -> Result<RxmRawxMeas, UbxParseError> {
    let pr_mes = f64::from_le_bytes(block[0..8].try_into().map_err(|_| UbxParseError::InvalidLength)?);
    let cp_mes = f64::from_le_bytes(block[8..16].try_into().map_err(|_| UbxParseError::InvalidLength)?);
    let do_mes = f32::from_le_bytes(block[16..20].try_into().map_err(|_| UbxParseError::InvalidLength)?);
    let gnss_id = block[20];
    let sv_id = block[21];
    let sig_id = block[22];
    let freq_id = block[23];
    let locktime = u16::from_le_bytes(block[24..26].try_into().map_err(|_| UbxParseError::InvalidLength)?);
    let cno = block[26];
    
    let pr_stdev_n = block[27];
    let cp_stdev_n = block[28];
    let do_stdev_n = block[29];
    let trk_stat = block[30];

    let pr_valid = (trk_stat & 0x01) != 0;
    let cp_valid = (trk_stat & 0x02) != 0;
    let half_cycle_valid = (trk_stat & 0x04) != 0;
    let sub_half_cycle = (trk_stat & 0x08) != 0;

    let pr_stdev = 0.01 * f64::powi(2.0, pr_stdev_n as i32);
    let cp_stdev = 0.004 * f64::powi(2.0, cp_stdev_n as i32);
    let do_stdev = 0.002 * f64::powi(2.0, do_stdev_n as i32);

    Ok(RxmRawxMeas {
        pr_mes, cp_mes, do_mes, gnss_id, sv_id, sig_id, freq_id, locktime, cno,
        pr_stdev, cp_stdev, do_stdev, pr_valid, cp_valid, half_cycle_valid, sub_half_cycle,
    })
}

"""

with open('crates/gneiss-parsers/src/ubx.rs', 'w') as f:
    f.write(pre_content + new_content + post_content)

print("Successfully refactored ubx.rs")
