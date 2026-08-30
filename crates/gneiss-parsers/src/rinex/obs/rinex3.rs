//! RINEX 3 observation data parser.

use gneiss_core::obs::{EpochObs, Observation, SatObs};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use std::collections::HashMap;
use super::header::parse_rinex_3_header;
use super::{map_rinex_type, parse_rinex_obs_lli, RinexObsHeader};

pub fn parse_rinex_3_obs<I: Iterator<Item = String>>(
    first_line: String,
    lines: &mut I,
) -> Result<(Vec<EpochObs>, RinexObsHeader), String> {
    let mut epochs = Vec::new();
    let (const_obs_types, header) = parse_rinex_3_header(first_line, lines)?;

    while let Some(line) = lines.next() {
        if !line.starts_with('>') || line.len() < 35 {
            continue;
        }

        let year = line[2..6].trim().parse::<i32>().unwrap_or(0);
        let month = line[7..9].trim().parse::<i32>().unwrap_or(0);
        let day = line[10..12].trim().parse::<i32>().unwrap_or(0);
        let hour = line[13..15].trim().parse::<i32>().unwrap_or(0);
        let min = line[16..18].trim().parse::<i32>().unwrap_or(0);
        let sec = line[19..29].trim().parse::<f64>().unwrap_or(0.0);

        let num_sats = line[32..35].trim().parse::<usize>().unwrap_or(0);
        let time = GpsTime::from_calendar(year, month, day, hour, min, sec);
        let mut satellites = Vec::with_capacity(num_sats);

        for _ in 0..num_sats {
            if let Some(obs_line) = lines.next() {
                if let Some(sat_obs) = parse_rinex_3_obs_line(&obs_line, &const_obs_types) {
                    satellites.push(sat_obs);
                }
            }
        }
        epochs.push(EpochObs { time, satellites });
    }
    Ok((epochs, header))
}

pub fn parse_rinex_3_obs_line(
    line: &str,
    const_obs_types: &HashMap<Constellation, Vec<String>>,
) -> Option<SatObs> {
    if line.len() < 3 {
        return None;
    }
    let sat_str = &line[0..3];
    let constellation_char = sat_str.chars().next().unwrap_or(' ');
    let constellation = match constellation_char {
        'G' => Constellation::Gps,
        'R' => Constellation::Glonass,
        'E' => Constellation::Galileo,
        'C' => Constellation::Beidou,
        'J' => Constellation::Qzss,
        'S' => Constellation::Sbas,
        _ => return None,
    };
    let prn = sat_str[1..3].trim().parse::<u8>().unwrap_or(0);
    let sat = SatelliteId { constellation, prn };

    let obs_types = const_obs_types.get(&constellation)?;
    let mut observations = Vec::new();

    for (i, type_str) in obs_types.iter().enumerate() {
        if let Some(obs) = parse_rinex_3_obs_value(line, i, type_str) {
            observations.push(obs);
        }
    }
    Some(SatObs { sat, observations })
}

fn parse_rinex_3_obs_value(obs_line: &str, type_idx: usize, type_str: &str) -> Option<Observation> {
    let start = 3 + type_idx * 16;
    if start >= obs_line.len() {
        return None;
    }
    let end = (start + 14).min(obs_line.len());
    let val_str = obs_line[start..end].trim();
    if val_str.is_empty() {
        return None;
    }
    let val = val_str.parse::<f64>().ok()?;
    let lli = parse_rinex_obs_lli(obs_line, start);
    map_rinex_type(type_str, val, lli)
}
