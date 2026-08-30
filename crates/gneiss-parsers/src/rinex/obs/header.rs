//! RINEX observation file header parsers for RINEX 2 and RINEX 3.

use gneiss_core::sat::Constellation;
use std::collections::HashMap;
use super::{parse_rinex_f14, RinexObsHeader};

pub fn parse_rinex_2_header<I: Iterator<Item = String>>(
    first_line: String,
    lines: &mut I,
) -> Result<(Vec<String>, RinexObsHeader), String> {
    let mut obs_types: Vec<String> = Vec::new();
    let mut num_obs = 0;
    let mut header = RinexObsHeader::default();

    let mut current_line = first_line;
    loop {
        if current_line.contains("APPROX POSITION XYZ") && current_line.len() >= 42 {
            let x = parse_rinex_f14(&current_line[0..14]);
            let y = parse_rinex_f14(&current_line[14..28]);
            let z = parse_rinex_f14(&current_line[28..42]);
            if let (Some(x), Some(y), Some(z)) = (x, y, z) {
                header.approx_position = Some([x, y, z]);
            }
        }
        if current_line.contains("ANTENNA: DELTA H/E/N") && current_line.len() >= 42 {
            let h = parse_rinex_f14(&current_line[0..14]);
            let e = parse_rinex_f14(&current_line[14..28]);
            let n = parse_rinex_f14(&current_line[28..42]);
            if let (Some(h), Some(e), Some(n)) = (h, e, n) {
                header.antenna_delta = Some([h, e, n]);
            }
        }
        if current_line.contains("MARKER NAME") {
            header.marker_name = Some(current_line[0..60].trim().to_string());
        }
        if current_line.contains("# / TYPES OF OBSERV") {
            if num_obs == 0 {
                num_obs = current_line[0..6].trim().parse::<usize>().unwrap_or(0);
            }
            if current_line.len() >= 60 {
                let types_str = &current_line[6..60];
                for chunk in types_str.as_bytes().chunks(6) {
                    let t = core::str::from_utf8(chunk).unwrap_or("").trim();
                    if !t.is_empty() {
                        obs_types.push(t.into());
                    }
                }
            }
        }
        if current_line.contains("END OF HEADER") {
            break;
        }
        if let Some(next_line) = lines.next() {
            current_line = next_line;
        } else {
            break;
        }
    }

    if obs_types.is_empty() {
        return Err("No observation types found in header".into());
    }
    Ok((obs_types, header))
}

pub fn parse_rinex_3_header<I: Iterator<Item = String>>(
    first_line: String,
    lines: &mut I,
) -> Result<(HashMap<Constellation, Vec<String>>, RinexObsHeader), String> {
    let mut const_obs_types = HashMap::new();
    let mut header = RinexObsHeader::default();

    let mut current_line = first_line;
    loop {
        if current_line.contains("SYS / # / OBS TYPES") {
            let constellation_char = current_line.chars().next().unwrap_or(' ');
            let constellation = match constellation_char {
                'G' => Constellation::Gps,
                'R' => Constellation::Glonass,
                'E' => Constellation::Galileo,
                'C' => Constellation::Beidou,
                'J' => Constellation::Qzss,
                'S' => Constellation::Sbas,
                _ => {
                    if let Some(next_line) = lines.next() {
                        current_line = next_line;
                        continue;
                    } else {
                        break;
                    }
                }
            };

            let count = current_line[3..6].trim().parse::<usize>().unwrap_or(0);
            let mut types_str = if current_line.len() >= 60 {
                current_line[7..60].to_string()
            } else {
                "".to_string()
            };

            let types =
                parse_rinex_3_obs_types_list(count, &mut types_str, lines, &mut current_line)?;
            const_obs_types.insert(constellation, types);
        }
        if current_line.contains("APPROX POSITION XYZ") && current_line.len() >= 42 {
            let x = parse_rinex_f14(&current_line[0..14]);
            let y = parse_rinex_f14(&current_line[14..28]);
            let z = parse_rinex_f14(&current_line[28..42]);
            if let (Some(x), Some(y), Some(z)) = (x, y, z) {
                header.approx_position = Some([x, y, z]);
            }
        }
        if current_line.contains("ANTENNA: DELTA H/E/N") && current_line.len() >= 42 {
            let h = parse_rinex_f14(&current_line[0..14]);
            let e = parse_rinex_f14(&current_line[14..28]);
            let n = parse_rinex_f14(&current_line[28..42]);
            if let (Some(h), Some(e), Some(n)) = (h, e, n) {
                header.antenna_delta = Some([h, e, n]);
            }
        }
        if current_line.contains("MARKER NAME") {
            header.marker_name = Some(current_line[0..60].trim().to_string());
        }
        if current_line.contains("END OF HEADER") {
            break;
        }
        if let Some(next_line) = lines.next() {
            current_line = next_line;
        } else {
            break;
        }
    }

    if const_obs_types.is_empty() {
        return Err("No observation types found in RINEX 3 header".into());
    }
    Ok((const_obs_types, header))
}

pub(crate) fn parse_rinex_3_obs_types_list<I: Iterator<Item = String>>(
    count: usize,
    types_str: &mut String,
    lines: &mut I,
    current_line: &mut String,
) -> Result<Vec<String>, String> {
    let mut types = Vec::new();
    while types.len() < count {
        for chunk in types_str.as_bytes().chunks(4) {
            let t = core::str::from_utf8(chunk).unwrap_or("").trim();
            if !t.is_empty() && types.len() < count {
                types.push(t.into());
            }
        }
        if types.len() < count {
            if let Some(next_line) = lines.next() {
                *current_line = next_line;
                if !current_line.contains("SYS / # / OBS TYPES") {
                    return Err("Expected continuation of SYS / # / OBS TYPES".into());
                }
                *types_str = if current_line.len() >= 60 {
                    current_line[7..60].to_string()
                } else {
                    "".to_string()
                };
            } else {
                break;
            }
        }
    }
    Ok(types)
}
