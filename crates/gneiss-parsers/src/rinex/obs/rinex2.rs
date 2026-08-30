//! RINEX 2 observation data parser.

use gneiss_core::obs::{EpochObs, Observation, SatObs};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use super::header::parse_rinex_2_header;
use super::{map_rinex_type, parse_rinex_obs_lli, RinexObsHeader};

pub fn parse_rinex_2_obs<I: Iterator<Item = String>>(
    first_line: String,
    lines: &mut I,
) -> Result<(Vec<EpochObs>, RinexObsHeader), String> {
    let mut epochs = Vec::new();
    let (obs_types, header) = parse_rinex_2_header(first_line, lines)?;

    while let Some(line) = lines.next() {
        if line.trim().is_empty() || line.len() < 32 {
            continue;
        }

        let year_str = line[1..3].trim();
        let year_val = year_str.parse::<i32>().unwrap_or(0);
        let year = if year_val >= 80 {
            1900 + year_val
        } else {
            2000 + year_val
        };

        let month = line[4..6].trim().parse::<i32>().unwrap_or(0);
        let day = line[7..9].trim().parse::<i32>().unwrap_or(0);
        let hour = line[10..12].trim().parse::<i32>().unwrap_or(0);
        let min = line[13..15].trim().parse::<i32>().unwrap_or(0);
        let sec = line[16..26].trim().parse::<f64>().unwrap_or(0.0);

        let flag = line[26..29].trim().parse::<i32>().unwrap_or(0);
        if flag > 1 {
            let num_skip = line[29..32].trim().parse::<usize>().unwrap_or(0);
            for _ in 0..num_skip {
                lines.next();
            }
            continue;
        }

        let num_sats = line[29..32].trim().parse::<usize>().unwrap_or(0);
        let mut sat_list = Vec::new();
        let mut sat_str = if line.len() > 32 {
            line[32..].to_string()
        } else {
            "".to_string()
        };

        while sat_list.len() < num_sats {
            let mut offset = 0;
            while offset + 3 <= sat_str.len() && sat_list.len() < num_sats {
                sat_list.push(sat_str[offset..offset + 3].to_string());
                offset += 3;
            }
            if sat_list.len() < num_sats {
                if let Some(next_line) = lines.next() {
                    sat_str = if next_line.len() > 32 {
                        next_line[32..].to_string()
                    } else {
                        next_line.to_string()
                    };
                } else {
                    break;
                }
            }
        }

        let time = GpsTime::from_calendar(year, month, day, hour, min, sec);
        let mut satellites = Vec::new();

        for sat_id_str in sat_list {
            if let Some(sat_obs) = parse_rinex_2_obs_sat(&sat_id_str, &obs_types, lines) {
                satellites.push(sat_obs);
            }
        }
        epochs.push(EpochObs { time, satellites });
    }
    Ok((epochs, header))
}

fn parse_rinex_2_obs_sat<I: Iterator<Item = String>>(
    sat_id_str: &str,
    obs_types: &[String],
    lines: &mut I,
) -> Option<SatObs> {
    let constellation_char = sat_id_str.chars().next().unwrap_or('G');
    let prn = sat_id_str[1..3].trim().parse::<u8>().unwrap_or(0);
    let constellation = match constellation_char {
        'G' | ' ' => Constellation::Gps,
        'R' => Constellation::Glonass,
        'E' => Constellation::Galileo,
        'C' => Constellation::Beidou,
        'S' => Constellation::Sbas,
        'J' => Constellation::Qzss,
        _ => return None,
    };
    let sat = SatelliteId { constellation, prn };

    let mut observations = Vec::new();
    let num_val_lines = (obs_types.len() as f64 / 5.0).ceil() as usize;
    let mut val_idx = 0;

    for _ in 0..num_val_lines {
        let Some(obs_line) = lines.next() else { continue };
        val_idx = parse_rinex_2_obs_line(&obs_line, obs_types, val_idx, &mut observations);
    }
    Some(SatObs { sat, observations })
}

fn parse_rinex_2_obs_line(
    obs_line: &str,
    obs_types: &[String],
    mut val_idx: usize,
    observations: &mut Vec<Observation>,
) -> usize {
    for col in 0..5 {
        if val_idx >= obs_types.len() {
            break;
        }
        if let Some(obs) = parse_rinex_2_obs_value(obs_line, col, &obs_types[val_idx]) {
            observations.push(obs);
        }
        val_idx += 1;
    }
    val_idx
}

fn parse_rinex_2_obs_value(obs_line: &str, col: usize, type_str: &str) -> Option<Observation> {
    let start = col * 16;
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
