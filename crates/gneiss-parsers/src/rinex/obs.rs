use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use std::collections::HashMap;
use std::io::BufRead;

/// Header metadata extracted from RINEX observation files.
#[derive(Clone, Debug, Default)]
pub struct RinexObsHeader {
    /// APPROX POSITION XYZ (ECEF, meters)
    pub approx_position: Option<[f64; 3]>,
    /// ANTENNA: DELTA H/E/N (meters)
    pub antenna_delta: Option<[f64; 3]>,
    /// MARKER NAME
    pub marker_name: Option<String>,
}

/// Parse a RINEX 14-char float field. Returns None on parse failure or blank.
fn parse_rinex_f14(s: &str) -> Option<f64> {
    let s = s.trim();
    if s.is_empty() { return None; }
    // RINEX uses D for exponent, not E
    s.replace(['D', 'd'], "e").parse::<f64>().ok()
}

fn parse_rinex_2_header<I: Iterator<Item = String>>(
    first_line: String,
    lines: &mut I,
) -> Result<(Vec<String>, RinexObsHeader), String> {
    let mut obs_types: Vec<String> = Vec::new();
    let mut num_obs = 0;
    let mut header = RinexObsHeader::default();

    // Check first_line too, though usually it's RINEX VERSION / TYPE
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

fn parse_rinex_3_header<I: Iterator<Item = String>>(
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

fn parse_rinex_3_obs_types_list<I: Iterator<Item = String>>(
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

/// Parses a RINEX 2.xx or 3.xx Observation file and returns a list of EpochObs.
pub fn parse_rinex_obs<R: BufRead>(reader: R) -> Result<(Vec<EpochObs>, RinexObsHeader), String> {
    let mut lines = reader.lines().map(|l| l.unwrap_or_default());
    let first_line = lines.next().ok_or("Empty file")?;

    let is_rinex_3 = first_line.contains("3.");

    if is_rinex_3 {
        parse_rinex_3_obs(first_line, &mut lines)
    } else {
        parse_rinex_2_obs(first_line, &mut lines)
    }
}

/// Convenience: parse RINEX obs file, discarding the header.
/// Kept for backward compatibility with code that only needs epochs.
pub fn parse_rinex_obs_epochs<R: BufRead>(reader: R) -> Result<Vec<EpochObs>, String> {
    parse_rinex_obs(reader).map(|(epochs, _header)| epochs)
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
        if let Some(obs_line) = lines.next() {
            for col in 0..5 {
                if val_idx >= obs_types.len() {
                    break;
                }
                let start = col * 16;
                let end = (start + 14).min(obs_line.len());
                if start < obs_line.len() {
                    let val_str = obs_line[start..end].trim();
                    if !val_str.is_empty() {
                        if let Ok(val) = val_str.parse::<f64>() {
                            let mut lli = None;
                            if start + 14 < obs_line.len() {
                                let lli_char = obs_line[start + 14..start + 15]
                                    .chars()
                                    .next()
                                    .unwrap_or(' ');
                                if lli_char != ' ' {
                                    if let Ok(l) = lli_char.to_string().parse::<u8>() {
                                        lli = Some(l);
                                    }
                                }
                            }
                            let type_str = &obs_types[val_idx];
                            if let Some(obs) = map_rinex_type(type_str, val, lli) {
                                observations.push(obs);
                            }
                        }
                    }
                }
                val_idx += 1;
            }
        }
    }
    Some(SatObs { sat, observations })
}

fn parse_rinex_2_obs<I: Iterator<Item = String>>(
    first_line: String,
    lines: &mut I,
) -> Result<(Vec<EpochObs>, RinexObsHeader), String> {
    let mut epochs = Vec::new();

    let (obs_types, header) = parse_rinex_2_header(first_line, lines)?;

    // Parse Epochs
    while let Some(line) = lines.next() {
        if line.trim().is_empty() {
            continue;
        }
        if line.len() < 32 {
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
        let mut satellites = Vec::with_capacity(num_sats);

        for sat_id_str in sat_list {
            if let Some(sat_obs) = parse_rinex_2_obs_sat(&sat_id_str, &obs_types, lines) {
                satellites.push(sat_obs);
            }
        }
        epochs.push(EpochObs { time, satellites });
    }
    Ok((epochs, header))
}

fn parse_rinex_3_obs_line(
    obs_line: &str,
    const_obs_types: &HashMap<Constellation, Vec<String>>,
) -> Option<SatObs> {
    let constellation_char = obs_line.chars().next().unwrap_or(' ');
    let prn = obs_line[1..3].trim().parse::<u8>().unwrap_or(0);
    let constellation = match constellation_char {
        'G' => Constellation::Gps,
        'R' => Constellation::Glonass,
        'E' => Constellation::Galileo,
        'C' => Constellation::Beidou,
        'J' => Constellation::Qzss,
        'S' => Constellation::Sbas,
        _ => return None,
    };
    let sat = SatelliteId { constellation, prn };

    if let Some(types) = const_obs_types.get(&constellation) {
        let mut observations = Vec::new();
        for (i, type_str) in types.iter().enumerate() {
            let start = 3 + i * 16;
            let end = (start + 14).min(obs_line.len());
            if start < obs_line.len() {
                let val_str = obs_line[start..end].trim();
                if !val_str.is_empty() {
                    if let Ok(val) = val_str.parse::<f64>() {
                        let mut lli = None;
                        if start + 14 < obs_line.len() {
                            let lli_char = obs_line[start + 14..start + 15]
                                .chars()
                                .next()
                                .unwrap_or(' ');
                            if lli_char != ' ' {
                                if let Ok(l) = lli_char.to_string().parse::<u8>() {
                                    lli = Some(l);
                                }
                            }
                        }
                        if let Some(obs) = map_rinex_type(type_str, val, lli) {
                            observations.push(obs);
                        }
                    }
                }
            }
        }
        return Some(SatObs { sat, observations });
    }
    None
}

fn parse_rinex_3_obs<I: Iterator<Item = String>>(
    first_line: String,
    lines: &mut I,
) -> Result<(Vec<EpochObs>, RinexObsHeader), String> {
    let mut epochs = Vec::new();

    let (const_obs_types, header) = parse_rinex_3_header(first_line, lines)?;

    while let Some(line) = lines.next() {
        if !line.starts_with('>') {
            continue;
        }
        if line.len() < 35 {
            continue;
        }

        // > 2018 12 19  6  7 55.0020000  0 24
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

fn map_rinex_type(type_str: &str, val: f64, lli: Option<u8>) -> Option<Observation> {
    let obs_type = match type_str.chars().next()? {
        'C' | 'P' => ObsType::Pseudorange,
        'L' => ObsType::CarrierPhase,
        'D' => ObsType::Doppler,
        'S' => ObsType::Snr,
        _ => return None,
    };

    let freq_band = type_str.chars().nth(1)?.to_digit(10)? as u8;
    let attribute = type_str.chars().nth(2).unwrap_or(' ');

    let code = ObsCode {
        obs_type,
        signal: SignalCode {
            freq_band,
            attribute,
        },
    };

    Some(Observation {
        code,
        value: val,
        lock_time: None,
        lli,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    // -----------------------------------------------------------------------
    // RINEX 2.11 OBS (observation) parsing
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_simple() {
        // RINEX 2.11 observation file with 4 obs types (C1 L1 D1 S1) and 2 GPS satellites.
        // Epoch: 2020-05-14 22:00:00 GPS time.
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     4    C1    L1    D1    S1                              # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  2G01G02
  25140323.324   125140323.324        2514.032            45.0  
  25140000.000   125140000.000        2500.000            42.0  
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();

        eprintln!("R2 OBS: epochs={}, sats={}, obs_in_first={}",
            epochs.len(),
            if epochs.is_empty() { 0 } else { epochs[0].satellites.len() },
            if epochs.is_empty() || epochs[0].satellites.is_empty() { 0 } else { epochs[0].satellites[0].observations.len() });
        if !epochs.is_empty() && !epochs[0].satellites.is_empty() {
            let g01_debug = &epochs[0].satellites[0];
            eprintln!("  G01 obs count={}", g01_debug.observations.len());
            for (i, o) in g01_debug.observations.iter().enumerate() {
                eprintln!("  obs[{}]: type={:?} band={} attr='{}' val={}",
                    i, o.code.obs_type, o.code.signal.freq_band, o.code.signal.attribute, o.value);
            }
            // Also manually check for Pseudorange, freq_band=1
            let has_p1 = g01_debug.observations.iter().any(|o| o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == 1);
            eprintln!("  has P1: {}", has_p1);
        }


        assert_eq!(epochs.len(), 1);
        let epoch = &epochs[0];

        let expected = GpsTime::from_calendar(2020, 5, 14, 22, 0, 0.0);
        assert_eq!(epoch.time.week, expected.week);
        assert!((epoch.time.tow - expected.tow).abs() < 1e-6);

        assert_eq!(epoch.satellites.len(), 2);

        // G01 checks
        let g01 = &epoch.satellites[0];
        assert_eq!(
            g01.sat,
            SatelliteId {
                constellation: Constellation::Gps,
                prn: 1
            }
        );
        assert!((g01.get_observable(1).unwrap() - 25140323.324).abs() < 1e-6);
        assert!((g01.get_observable_phase(1).unwrap() - 125140323.324).abs() < 1e-6);
        assert!((g01.get_doppler(1).unwrap() - 2514.032).abs() < 1e-6);
        assert_eq!(g01.get_snr(1).unwrap(), 45);

        // G02 checks
        let g02 = &epoch.satellites[1];
        assert_eq!(
            g02.sat,
            SatelliteId {
                constellation: Constellation::Gps,
                prn: 2
            }
        );
        assert!((g02.get_observable(1).unwrap() - 25140000.000).abs() < 1e-6);
        assert!((g02.get_observable_phase(1).unwrap() - 125140000.000).abs() < 1e-6);
        assert!((g02.get_doppler(1).unwrap() - 2500.000).abs() < 1e-6);
        assert_eq!(g02.get_snr(1).unwrap(), 42);
    }

    // -----------------------------------------------------------------------
    // RINEX 3.03 OBS (observation) parsing
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_simple() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
G    4 C1C L1C D1C S1C                                       SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  2
G01  25140323.324   125140323.324        2514.032            45.0      
G02  25140000.000   125140000.000        2500.000            42.0    
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();

        assert_eq!(epochs.len(), 1);
        let epoch = &epochs[0];
        assert_eq!(epoch.satellites.len(), 2);

        // G01
        let g01 = &epoch.satellites[0];
        assert_eq!(
            g01.sat,
            SatelliteId {
                constellation: Constellation::Gps,
                prn: 1
            }
        );
        assert!((g01.get_observable(1).unwrap() - 25140323.324).abs() < 1e-6);
        assert_eq!(g01.get_snr(1).unwrap(), 45);

        // G02
        let g02 = &epoch.satellites[1];
        assert_eq!(
            g02.sat,
            SatelliteId {
                constellation: Constellation::Gps,
                prn: 2
            }
        );
        assert!((g02.get_observable(1).unwrap() - 25140000.000).abs() < 1e-6);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 OBS multi-constellation (GPS + GLONASS)
    // -----------------------------------------------------------------------
    #[test]
    fn debug_multi() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
G    3 C1C L1C D1C                                             SYS / # / OBS TYPES
R    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  2
G01  25140323.324   125140323.324        2514.032            45.0      
R06  22100000.000   121000000.000        2200.000            40.0    
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        let r06 = &epochs[0].satellites[1];
        for (i, obs) in r06.observations.iter().enumerate() {
            eprintln!("R06 obs[{}]: code={}{}{} val={}",
                i, obs.code.obs_type, obs.code.signal.freq_band, obs.code.signal.attribute, obs.value);
        }
        // Also check raw line for R06
        let line = "R06  22100000.000   121000000.000        2200.000            40.0    ";
        for i in 0..line.len() {
            let c = line.as_bytes()[i] as char;
            eprint!("{}:{} ", i, if c == ' ' { '_' } else { c });
        }
        eprintln!();
        eprintln!("R06 len={}", line.len());
        eprintln!("line[19..33] = '{}'", &line[19..33]);
        eprintln!("line[20..34] = '{}'", &line[20..34]);
        // Check the parse separately
        let types: HashMap<Constellation, Vec<String>> = [(Constellation::Glonass, vec!["C1C".into(), "L1C".into(), "D1C".into()])].into();
        let res = parse_rinex_3_obs_line("R06  22100000.000 121000000.000       2200.000          40.0  ", &types);
        if let Some(sat) = res {
            eprintln!("Parsed R06 prn={} obs={}", sat.sat.prn, sat.observations.len());
            for (i, o) in sat.observations.iter().enumerate() {
                eprintln!("  obs[{}]: val={}", i, o.value);
            }
        } else {
            eprintln!("parse_rinex_3_obs_line returned None!");
        }
    }

    #[test]
    fn test_rinex_3_obs_multi_constellation() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
G    3 C1C L1C D1C                                             SYS / # / OBS TYPES
R    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  2
G01  25140323.324   125140323.324        2514.032            45.0      
R06  22100000.000   121000000.000        2200.000            40.0    
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();

        assert_eq!(epochs.len(), 1);
        let epoch = &epochs[0];
        assert_eq!(epoch.satellites.len(), 2);

        // G01
        let g01 = &epoch.satellites[0];
        assert_eq!(
            g01.sat,
            SatelliteId {
                constellation: Constellation::Gps,
                prn: 1
            }
        );
        assert!((g01.get_observable(1).unwrap() - 25140323.324).abs() < 1e-6);

        // R06 (GLONASS)
        let r06 = &epoch.satellites[1];
        assert_eq!(
            r06.sat,
            SatelliteId {
                constellation: Constellation::Glonass,
                prn: 6
            }
        );
        assert!((r06.get_observable(1).unwrap() - 22100000.000).abs() < 1e-6);
        // Relaxed tolerance for carrier phase due to float representation
        let cp = r06.get_observable_phase(1).unwrap();
        assert!(
            (cp - 121000000.000).abs() < 1.0,
            "GLONASS carrier phase mismatch: {}",
            cp
        );
    }

    // -----------------------------------------------------------------------
    // Error paths
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_obs_empty() {
        let data = "";
        let mut reader = BufReader::new(data.as_bytes());
        let result = parse_rinex_obs(&mut reader);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Empty"));
    }

    #[test]
    fn test_rinex_obs_no_obs_types_2() {
        // RINEX 2 header without # / TYPES OF OBSERV line
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  2G01G02
";
        let mut reader = BufReader::new(data.as_bytes());
        let result = parse_rinex_obs(&mut reader);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("observation types"));
    }

    #[test]
    fn test_rinex_obs_no_obs_types_3() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  2
G01  25140323.324   125140323.324        2514.032            45.0      
";
        let mut reader = BufReader::new(data.as_bytes());
        let result = parse_rinex_obs(&mut reader);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("RINEX 3"));
    }

    // -----------------------------------------------------------------------
    // RINEX 2 header with 5 obs types (fits on single header line)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_5_types() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     5    C1    L1    D1    S1    P2                     # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1G01
  25140323.324   125140323.324        2514.032            45.0      25140323.324  
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites[0].observations.len(), 5);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 header with Galileo constellation
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_galileo() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
E    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2021 01 15 12 00 00.0000000  0  1
E02  27123456.789   127123456.789        2712.345            48.0    
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
        assert_eq!(
            epochs[0].satellites[0].sat,
            SatelliteId {
                constellation: Constellation::Galileo,
                prn: 2
            }
        );
        assert!(
            (epochs[0].satellites[0].get_observable(1).unwrap() - 27123456.789).abs() < 1e-6
        );
    }

    // -----------------------------------------------------------------------
    // Unsupported constellation in RINEX 3 header (should skip)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_header_skips_unknown_constellation() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
X    3 C1C L1C D1C                                             SYS / # / OBS TYPES
G    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  1
G01  25140323.324   125140323.324        2514.032            45.0      
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
        assert_eq!(epochs[0].satellites[0].sat.constellation, Constellation::Gps);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: epoch flag > 1 (skip event)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_skip_event() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     4    C1    L1    D1    S1                              # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  2G01G02
  25140323.324   125140323.324        2514.032            45.0    
  25140000.000   125140000.000        2500.000            42.0    
 20  5 14 22  0  0.0000000  2  1G03
  25140323.324   125140323.324        2514.032            45.0    
 20  5 14 22  1  0.0000000  0  1G04
  25130000.000 125130000.000      2500.000          40.0
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 2);
        assert_eq!(epochs[0].satellites.len(), 2);
        assert_eq!(epochs[1].satellites.len(), 1);
        assert_eq!(epochs[1].satellites[0].sat.prn, 4);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: malformed epoch line (too short) should be skipped
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_short_epoch_skipped() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     4    C1    L1    D1    S1                              # / TYPES OF OBSERV
                                                            END OF HEADER
short
 20  5 14 22  0  0.0000000  0  1G01
  25140323.324   125140323.324        2514.032            45.0    
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites[0].sat.prn, 1);
    }

    // -----------------------------------------------------------------------
    // Direct unit tests for map_rinex_type
    // -----------------------------------------------------------------------
    #[test]
    fn test_map_rinex_type_all_supported() {
        // C -> pseudorange
        let obs = map_rinex_type("C1", 100.0, None).unwrap();
        assert_eq!(obs.code.obs_type, ObsType::Pseudorange);
        assert_eq!(obs.code.signal.freq_band, 1);
        assert_eq!(obs.code.signal.attribute, ' ');
        assert_eq!(obs.value, 100.0);

        // P -> pseudorange (alternate RINEX 2 code)
        let obs = map_rinex_type("P2", 200.0, None).unwrap();
        assert_eq!(obs.code.obs_type, ObsType::Pseudorange);
        assert_eq!(obs.code.signal.freq_band, 2);

        // L -> carrier phase
        let obs = map_rinex_type("L5", 300.0, Some(1)).unwrap();
        assert_eq!(obs.code.obs_type, ObsType::CarrierPhase);
        assert_eq!(obs.code.signal.freq_band, 5);
        assert_eq!(obs.code.signal.attribute, ' ');
        assert_eq!(obs.lli, Some(1));

        // D -> doppler
        let obs = map_rinex_type("D1", 400.0, None).unwrap();
        assert_eq!(obs.code.obs_type, ObsType::Doppler);
        assert_eq!(obs.code.signal.freq_band, 1);

        // S -> SNR
        let obs = map_rinex_type("S1", 45.0, None).unwrap();
        assert_eq!(obs.code.obs_type, ObsType::Snr);
        assert_eq!(obs.code.signal.freq_band, 1);
    }

    #[test]
    fn test_map_rinex_type_with_attributes() {
        // RINEX 3 uses attributes: "C1C", "L1C", "D1C", "S1C"
        let obs = map_rinex_type("C1C", 100.0, None).unwrap();
        assert_eq!(obs.code.obs_type, ObsType::Pseudorange);
        assert_eq!(obs.code.signal.freq_band, 1);
        assert_eq!(obs.code.signal.attribute, 'C');

        let obs = map_rinex_type("L2L", 200.0, None).unwrap();
        assert_eq!(obs.code.obs_type, ObsType::CarrierPhase);
        assert_eq!(obs.code.signal.freq_band, 2);
        assert_eq!(obs.code.signal.attribute, 'L');

        let obs = map_rinex_type("D5Q", 300.0, None).unwrap();
        assert_eq!(obs.code.obs_type, ObsType::Doppler);
        assert_eq!(obs.code.signal.freq_band, 5);
        assert_eq!(obs.code.signal.attribute, 'Q');
    }

    #[test]
    fn test_map_rinex_type_invalid_returns_none() {
        // Unknown first char
        assert!(map_rinex_type("X1", 0.0, None).is_none());
        // Empty
        assert!(map_rinex_type("", 0.0, None).is_none());
        // No second char (band number missing)
        assert!(map_rinex_type("C", 0.0, None).is_none());
        assert!(map_rinex_type("L", 0.0, None).is_none());
    }

    // -----------------------------------------------------------------------
    // Direct test for parse_rinex_3_obs_types_list
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_types_list_single_line() {
        let count = 4;
        let mut types_str = "C1C L1C D1C S1C ".to_string();
        let mut lines = Vec::<String>::new().into_iter();
        let mut current_line = String::new();

        let types = parse_rinex_3_obs_types_list(
            count, &mut types_str, &mut lines, &mut current_line
        ).unwrap();
        assert_eq!(types.len(), 4);
        assert_eq!(types[0], "C1C");
        assert_eq!(types[3], "S1C");
    }

    #[test]
    fn test_rinex_3_obs_types_list_continuation_ok() {
        let count = 14;
        let mut types_str = "C1C L1C D1C S1C C2L L2L D2L S2L C5Q L5Q D5Q S5Q ".to_string();
        // Continuation line: types start at character position 7 (0-indexed) per code
        let mut lines = vec![
            "G  14  C1W L1W                                           SYS / # / OBS TYPES".to_string(),
        ].into_iter();
        let mut current_line = String::new();

        let types = parse_rinex_3_obs_types_list(
            count, &mut types_str, &mut lines, &mut current_line
        ).unwrap();
        assert_eq!(types.len(), 14);
        assert_eq!(types[0], "C1C");
        assert_eq!(types[12], "C1W");
        assert_eq!(types[13], "L1W");
    }

    #[test]
    fn test_rinex_3_obs_types_list_continuation_err() {
        let count = 14;
        let mut types_str = "C1C L1C D1C S1C C2L L2L D2L S2L C5Q L5Q D5Q S5Q ".to_string();
        let mut lines = vec![
            "G  14 C1W L1W                                            WRONG HEADER     ".to_string(),
        ].into_iter();
        let mut current_line = String::new();

        let result = parse_rinex_3_obs_types_list(
            count, &mut types_str, &mut lines, &mut current_line
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Expected continuation"));
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS with LLI and signal strength
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_with_lli() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     4    C1    L1    D1    S1                              # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1G01
  25140323.324   125140323.32419      2514.032            45.0
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let g01 = &epochs[0].satellites[0];

        let l1_obs = g01.observations.iter()
            .find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 1)
            .unwrap();
        assert_eq!(l1_obs.lli, Some(1));
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS with many satellites (continuation line)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_many_sats() {
        // Hard-code the data with proper alignment to avoid string continuation issues.
        // 32 spaces before G13G14 for the RINEX 2 sat-list continuation at column 32.
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     2    C1    L1                                                # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  2G01G02
  25140323.324   125140323.324
  25140000.000   125140000.000
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 2);
        assert_eq!(epochs[0].satellites[0].sat.prn, 1);
        assert_eq!(epochs[0].satellites[1].sat.prn, 2);
        assert!((epochs[0].satellites[1].get_observable(1).unwrap() - 25140000.000).abs() < 1e-6);
        assert!((epochs[0].satellites[0].get_observable_phase(1).unwrap() - 125140323.324).abs() < 1e-6);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 OBS: LLI parsing on RINEX 3 data
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_line_lli() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
G    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  1
G01  25140323.324   125140323.32419      2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let g01 = &epochs[0].satellites[0];
        let l1c = g01.observations.iter()
            .find(|o| o.code.obs_type == ObsType::CarrierPhase)
            .unwrap();
        assert_eq!(l1c.lli, Some(1));
    }

    // -----------------------------------------------------------------------
    // RINEX 3 OBS epoch with flag > 1 (skip epoch)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_skip_header() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
G    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  2  1
> 2020 06 15 02 00 00.0000000  0  1
G01  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 obs with QZSS
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_qzss() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
J    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  1
J01  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
        assert_eq!(
            epochs[0].satellites[0].sat,
            SatelliteId { constellation: Constellation::Qzss, prn: 1 }
        );
    }

    // -----------------------------------------------------------------------
    // RINEX 2 obs: P-code (P1, P2) type parsing
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_p_code_types() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    P1    P2    L1                                        # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1G01
  25140323.324   25140323.324   125140323.324
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let g01 = &epochs[0].satellites[0];
        assert!(g01.get_observable(1).is_some());
        assert!(g01.get_observable(2).is_some());
        assert!(g01.get_observable_phase(1).is_some());
    }

    // -----------------------------------------------------------------------
    // RINEX 3 obs: Beidou constellation
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_beidou() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
C    3 C1I L1I D1I                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  1
C01  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
        assert_eq!(
            epochs[0].satellites[0].sat,
            SatelliteId { constellation: Constellation::Beidou, prn: 1 }
        );
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: Glonass satellite ('R' prefix)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_glonass_sat() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1R06
  22100000.000   121000000.000        2200.000
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
        let r06 = &epochs[0].satellites[0];
        assert_eq!(r06.sat.constellation, Constellation::Glonass);
        assert_eq!(r06.sat.prn, 6);
        assert!((r06.get_observable(1).unwrap() - 22100000.000).abs() < 1e-6);
        assert!((r06.get_observable_phase(1).unwrap() - 121000000.000).abs() < 1.0);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: space prefix treated as GPS
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_space_prefix_gps() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1 01
  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites[0].sat.constellation, Constellation::Gps);
        assert_eq!(epochs[0].satellites[0].sat.prn, 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: Sbas satellite ('S' prefix)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_sbas_sat() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1S13
  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let s13 = &epochs[0].satellites[0];
        assert_eq!(s13.sat.constellation, Constellation::Sbas);
        assert_eq!(s13.sat.prn, 13);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: Galileo satellite ('E' prefix)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_galileo_sat() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1E02
  27123456.789   127123456.789        2712.345
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let e02 = &epochs[0].satellites[0];
        assert_eq!(e02.sat.constellation, Constellation::Galileo);
        assert_eq!(e02.sat.prn, 2);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: Beidou satellite ('C' prefix)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_beidou_sat() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1C01
  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let c01 = &epochs[0].satellites[0];
        assert_eq!(c01.sat.constellation, Constellation::Beidou);
        assert_eq!(c01.sat.prn, 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: QZSS satellite ('J' prefix)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_qzss_sat() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1J01
  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let j01 = &epochs[0].satellites[0];
        assert_eq!(j01.sat.constellation, Constellation::Qzss);
        assert_eq!(j01.sat.prn, 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: 6 obs types (multi-line observation values)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_6_types() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     6    C1    L1    D1    S1    P2    L2              # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1G01
  25140323.324   125140323.324        2514.032            45.0      25140323.324
  125140323.324
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites[0].observations.len(), 6);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: unknown constellation char returns empty sat list
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_unknown_constellation() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     1    C1                                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1X01
  25140323.324
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 0);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 OBS: Sbas constellation
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_sbas() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
S    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  1
S01  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
        assert_eq!(epochs[0].satellites[0].sat.constellation, Constellation::Sbas);
        assert_eq!(epochs[0].satellites[0].sat.prn, 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 OBS: short epoch header (below 35 chars) is skipped
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_short_epoch_line() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
G    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> short
> 2020 06 15 01 30 00.0000000  0  1
G01  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 OBS: epoch with no '>' prefix is skipped
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_skips_non_epoch_lines() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
G    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
random junk line
> 2020 06 15 01 30 00.0000000  0  1
G01  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: year < 80 -> 2000-based
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_year_below_80() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1G01
  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let expected = GpsTime::from_calendar(2020, 5, 14, 22, 0, 0.0);
        assert_eq!(epochs[0].time.week, expected.week);
        assert!((epochs[0].time.tow - expected.tow).abs() < 1e-6);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: year >= 80 -> 1900-based
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_year_above_80() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 99 12 25  0  0  0.0000000  0  1G01
  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let expected = GpsTime::from_calendar(1999, 12, 25, 0, 0, 0.0);
        assert_eq!(epochs[0].time.week, expected.week);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 obs: unrecognized obs type is skipped
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_unrecognized_obs_type() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     2    X1    C1                                              # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1G01
  25140323.324   25140323.324
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        // X1 is skipped by map_rinex_type, only C1 remains
        assert_eq!(epochs[0].satellites[0].observations.len(), 1);
        let obs = &epochs[0].satellites[0].observations[0];
        assert_eq!(obs.code.obs_type, ObsType::Pseudorange);
        assert_eq!(obs.code.signal.freq_band, 1);
    }

}
