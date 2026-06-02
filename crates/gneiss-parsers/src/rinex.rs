use gneiss_core::obs::{EpochObs, SatObs, Observation, ObsCode, ObsType, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use std::collections::HashMap;
use std::io::BufRead;

fn parse_rinex_2_header<I: Iterator<Item = String>>(first_line: String, lines: &mut I) -> Result<Vec<String>, String> {
    let mut obs_types: Vec<String> = Vec::new();
    let mut num_obs = 0;

    // Check first_line too, though usually it's RINEX VERSION / TYPE
    let mut current_line = first_line;
    loop {
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
    Ok(obs_types)
}

fn parse_rinex_3_header<I: Iterator<Item = String>>(first_line: String, lines: &mut I) -> Result<HashMap<Constellation, Vec<String>>, String> {
    let mut const_obs_types = HashMap::new();

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
                    if let Some(next_line) = lines.next() { current_line = next_line; continue; } else { break; }
                }
            };

            let count = current_line[3..6].trim().parse::<usize>().unwrap_or(0);
            let mut types = Vec::new();
            
            let mut types_str = if current_line.len() >= 60 { current_line[7..60].to_string() } else { "".to_string() };
            
            while types.len() < count {
                for chunk in types_str.as_bytes().chunks(4) {
                    let t = core::str::from_utf8(chunk).unwrap_or("").trim();
                    if !t.is_empty() && types.len() < count {
                        types.push(t.into());
                    }
                }
                
                if types.len() < count {
                    if let Some(next_line) = lines.next() {
                        current_line = next_line;
                        if !current_line.contains("SYS / # / OBS TYPES") {
                            return Err("Expected continuation of SYS / # / OBS TYPES".into());
                        }
                        types_str = if current_line.len() >= 60 { current_line[7..60].to_string() } else { "".to_string() };
                    } else {
                        break;
                    }
                }
            }
            const_obs_types.insert(constellation, types);
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
    Ok(const_obs_types)
}

/// Parses a RINEX 2.xx or 3.xx Observation file and returns a list of EpochObs.
pub fn parse_rinex_obs<R: BufRead>(reader: R) -> Result<Vec<EpochObs>, String> {
    let mut lines = reader.lines().map(|l| l.unwrap_or_default());
    let first_line = lines.next().ok_or("Empty file")?;
    
    let is_rinex_3 = first_line.contains("3.");
    
    if is_rinex_3 {
        parse_rinex_3_obs(first_line, &mut lines)
    } else {
        parse_rinex_2_obs(first_line, &mut lines)
    }
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
                if val_idx >= obs_types.len() { break; }
                let start = col * 16;
                let end = (start + 14).min(obs_line.len());
                if start < obs_line.len() {
                    let val_str = obs_line[start..end].trim();
                    if !val_str.is_empty() {
                        if let Ok(val) = val_str.parse::<f64>() {
                            let mut lli = None;
                            if start + 14 < obs_line.len() {
                                let lli_char = obs_line[start + 14..start + 15].chars().next().unwrap_or(' ');
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

fn parse_rinex_2_obs<I: Iterator<Item = String>>(first_line: String, lines: &mut I) -> Result<Vec<EpochObs>, String> {
    let mut epochs = Vec::new();

    let obs_types = parse_rinex_2_header(first_line, lines)?;

    // Parse Epochs
    while let Some(line) = lines.next() {
        if line.trim().is_empty() { continue; }
        if line.len() < 32 { continue; }
        
        let year_str = line[1..3].trim();
        let year_val = year_str.parse::<i32>().unwrap_or(0);
        let year = if year_val >= 80 { 1900 + year_val } else { 2000 + year_val };
        
        let month = line[4..6].trim().parse::<i32>().unwrap_or(0);
        let day = line[7..9].trim().parse::<i32>().unwrap_or(0);
        let hour = line[10..12].trim().parse::<i32>().unwrap_or(0);
        let min = line[13..15].trim().parse::<i32>().unwrap_or(0);
        let sec = line[16..26].trim().parse::<f64>().unwrap_or(0.0);
        
        let flag = line[26..29].trim().parse::<i32>().unwrap_or(0);
        if flag > 1 {
            let num_skip = line[29..32].trim().parse::<usize>().unwrap_or(0);
            for _ in 0..num_skip { lines.next(); }
            continue;
        }

        let num_sats = line[29..32].trim().parse::<usize>().unwrap_or(0);
        let mut sat_list = Vec::new();
        let mut sat_str = if line.len() > 32 { line[32..].to_string() } else { "".to_string() };
        
        while sat_list.len() < num_sats {
            let mut offset = 0;
            while offset + 3 <= sat_str.len() && sat_list.len() < num_sats {
                sat_list.push(sat_str[offset..offset+3].to_string());
                offset += 3;
            }
            if sat_list.len() < num_sats {
                if let Some(next_line) = lines.next() {
                    sat_str = if next_line.len() > 32 { next_line[32..].to_string() } else { next_line.to_string() };
                } else { break; }
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
    Ok(epochs)
}

fn parse_rinex_3_obs_line(obs_line: &str, const_obs_types: &HashMap<Constellation, Vec<String>>) -> Option<SatObs> {
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
                            let lli_char = obs_line[start + 14..start + 15].chars().next().unwrap_or(' ');
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

fn parse_rinex_3_obs<I: Iterator<Item = String>>(first_line: String, lines: &mut I) -> Result<Vec<EpochObs>, String> {
    let mut epochs = Vec::new();

    let const_obs_types = parse_rinex_3_header(first_line, lines)?;

    while let Some(line) = lines.next() {
        if !line.starts_with('>') { continue; }
        if line.len() < 35 { continue; }
        
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
    Ok(epochs)
}

fn map_rinex_type(type_str: &str, val: f64, lli: Option<u8>) -> Option<Observation> {
    let obs_char = type_str.chars().next()?;
    let freq_char = type_str.chars().nth(1)?;
    let attr_char = type_str.chars().nth(2).unwrap_or(' ');
    
    let obs_type = match obs_char {
        'C' | 'P' => ObsType::Pseudorange,
        'L' => ObsType::CarrierPhase,
        'D' => ObsType::Doppler,
        'S' => ObsType::Snr,
        _ => return None,
    };
    
    let freq_band = freq_char.to_digit(10)? as u8;
    
    let lock_time = if let Some(l) = lli {
        if l & 1 != 0 || l & 2 != 0 { Some(0) } else { None }
    } else {
        None
    };

    Some(Observation {
        code: ObsCode {
            obs_type,
            signal: SignalCode { freq_band, attribute: attr_char },
        },
        value: val,
        lock_time,
    })
}

pub fn parse_rinex_f64(s: &str) -> Result<f64, String> {
    let s = s.trim().to_string();
    if s.is_empty() { return Ok(0.0); }
    let s_clean = s.replace("D", "E").replace("d", "e");
    s_clean.parse::<f64>().map_err(|_| format!("Failed to parse RINEX f64: '{}'", s))
}

fn build_ephemeris(
    constellation: Constellation,
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32]
) -> Option<gneiss_core::ephemeris::Ephemeris> {
    use gneiss_core::ephemeris::{Ephemeris, GpsEphemeris, GalileoEphemeris, BeidouEphemeris, QzssEphemeris, GlonassEphemeris};
    match constellation {
        Constellation::Glonass => {
            Some(Ephemeris::Glonass(GlonassEphemeris {
                sat,
                toe: toc, 
                freq_num: vals[7] as i8, 
                tau_n: af0, gamma_n: af1, delta_tau_n: af2,
                x: vals[0] * 1000.0, y: vals[4] * 1000.0, z: vals[8] * 1000.0, 
                vx: vals[1] * 1000.0, vy: vals[5] * 1000.0, vz: vals[9] * 1000.0,
                ax: vals[2] * 1000.0, ay: vals[6] * 1000.0, az: vals[10] * 1000.0,
            }))
        },
        Constellation::Gps => {
            Some(Ephemeris::Gps(GpsEphemeris {
                sat, toc, toe: GpsTime::new(toc.week, vals[8]),
                af0, af1, af2,
                iode: vals[0] as u32, crs: vals[1], delta_n: vals[2], m0: vals[3],
                cuc: vals[4], e: vals[5], cus: vals[6], sqrt_a: vals[7],
                cic: vals[9], omega0: vals[10], cis: vals[11],
                i0: vals[12], crc: vals[13], omega: vals[14], omega_dot: vals[15],
                idot: vals[16], tgd: vals[22], iodc: vals[23] as u32,
            }))
        },
        Constellation::Galileo => {
            Some(Ephemeris::Galileo(GalileoEphemeris {
                sat, toc, toe: GpsTime::new(toc.week, vals[8]),
                af0, af1, af2,
                iod_nav: vals[0] as u32, crs: vals[1], delta_n: vals[2], m0: vals[3],
                cuc: vals[4], e: vals[5], cus: vals[6], sqrt_a: vals[7],
                cic: vals[9], omega0: vals[10], cis: vals[11],
                i0: vals[12], crc: vals[13], omega: vals[14], omega_dot: vals[15],
                idot: vals[16], bgd_e1_e5a: vals[22],
            }))
        },
        Constellation::Beidou => {
            Some(Ephemeris::Beidou(BeidouEphemeris {
                sat, toc, toe: GpsTime::new(toc.week, vals[8]),
                af0, af1, af2,
                aode: vals[0] as u32, crs: vals[1], delta_n: vals[2], m0: vals[3],
                cuc: vals[4], e: vals[5], cus: vals[6], sqrt_a: vals[7],
                cic: vals[9], omega0: vals[10], cis: vals[11],
                i0: vals[12], crc: vals[13], omega: vals[14], omega_dot: vals[15],
                idot: vals[16], tgd1: vals[22], aodc: vals[23] as u32,
            }))
        },
        Constellation::Qzss => {
            Some(Ephemeris::Qzss(QzssEphemeris {
                sat, toc, toe: GpsTime::new(toc.week, vals[8]),
                af0, af1, af2,
                iode: vals[0] as u32, crs: vals[1], delta_n: vals[2], m0: vals[3],
                cuc: vals[4], e: vals[5], cus: vals[6], sqrt_a: vals[7],
                cic: vals[9], omega0: vals[10], cis: vals[11],
                i0: vals[12], crc: vals[13], omega: vals[14], omega_dot: vals[15],
                idot: vals[16], tgd: vals[22], iodc: vals[23] as u32,
            }))
        },
        _ => None
    }
}

pub fn parse_rinex_nav<R: BufRead>(reader: R) -> Result<Vec<gneiss_core::ephemeris::Ephemeris>, String> {
    let mut ephemerides = Vec::new();
    let mut lines = reader.lines().map(|l| l.unwrap_or_default());

    let mut is_rinex_3 = false;
    for line in lines.by_ref() {
        if line.contains("RINEX VERSION / TYPE")
            && line.trim().starts_with('3') { is_rinex_3 = true; }
        if line.contains("END OF HEADER") { break; }
    }

    let mut current_constellation = Constellation::Gps;
    let mut current_prn = 0;
    let mut current_toc = GpsTime::new(0, 0.0);
    let mut current_af0 = 0.0;
    let mut current_af1 = 0.0;
    let mut current_af2 = 0.0;
    let mut line_idx = 0;
    let mut vals = [0.0; 32];

    for line in lines {
        if line.trim().is_empty() { continue; }
        let is_new_epoch = if is_rinex_3 {
            line.starts_with('G') || line.starts_with('R') || line.starts_with('E') || line.starts_with('C') || line.starts_with('J') || line.starts_with('S') || line.starts_with('I')
        } else {
            line_idx == 0 || (current_constellation != Constellation::Glonass && line_idx > 7) || (current_constellation == Constellation::Glonass && line_idx > 3)
        };

        if is_new_epoch {
            current_constellation = if is_rinex_3 {
                match line.chars().next().unwrap() {
                    'G' => Constellation::Gps,
                    'R' => Constellation::Glonass,
                    'E' => Constellation::Galileo,
                    'C' => Constellation::Beidou,
                    'J' => Constellation::Qzss,
                    'S' => Constellation::Sbas,
                    'I' => Constellation::Gps,
                    _ => Constellation::Gps,
                }
            } else { Constellation::Gps };

            current_prn = if is_rinex_3 {
                if line.len() >= 3 { line[1..3].trim().parse::<u8>().unwrap_or(0) } else { 0 }
            } else {
                if line.len() >= 2 { line[0..2].trim().parse::<u8>().unwrap_or(0) } else { 0 }
            };

            let year_str = if is_rinex_3 { 
                if line.len() >= 8 { &line[4..8] } else { "" }
            } else { 
                if line.len() >= 5 { &line[3..5] } else { "" }
            };
            let mut year = year_str.trim().parse::<i32>().unwrap_or(0);
            if year < 100 { year += if year > 80 { 1900 } else { 2000 }; }
            
            let month = if line.len() >= 11 { line[9..11].trim().parse::<i32>().unwrap_or(0) } else { 0 };
            let day = if line.len() >= 14 { line[12..14].trim().parse::<i32>().unwrap_or(0) } else { 0 };
            let hour = if line.len() >= 17 { line[15..17].trim().parse::<i32>().unwrap_or(0) } else { 0 };
            let min = if line.len() >= 20 { line[18..20].trim().parse::<i32>().unwrap_or(0) } else { 0 };
            let sec = if line.len() >= 22 { line[20..22].trim().parse::<f64>().unwrap_or(0.0) } else { 0.0 };

            let mut toc_gpst = GpsTime::from_calendar(year, month, day, hour, min, sec);
            match current_constellation {
                Constellation::Glonass => toc_gpst.tow += 18.0,
                Constellation::Beidou => toc_gpst.tow += 14.0,
                _ => {}
            }
            current_toc = toc_gpst;
            let (idx_af0, idx_af1, idx_af2) = if is_rinex_3 {
                (23, 42, 61)
            } else {
                (22, 41, 60)
            };
            
            current_af0 = parse_rinex_f64(if line.len() >= idx_af0 + 19 { &line[idx_af0..idx_af0+19] } else { "" })?;
            current_af1 = parse_rinex_f64(if line.len() >= idx_af1 + 19 { &line[idx_af1..idx_af1+19] } else { "" })?;
            current_af2 = parse_rinex_f64(if line.len() >= idx_af2 + 19 { &line[idx_af2..idx_af2+19] } else { "" })?;
            line_idx = 1;
            for i in 0..32 { vals[i] = 0.0; }
        } else {
            if (1..=8).contains(&line_idx) {
                let offset = (line_idx - 1) * 4;
                let (i0, i1, i2, i3) = if is_rinex_3 {
                    (4, 23, 42, 61)
                } else {
                    (3, 22, 41, 60)
                };
                
                vals[offset] = parse_rinex_f64(if line.len() >= i0 + 19 { &line[i0..i0+19] } else { "" })?;
                vals[offset+1] = parse_rinex_f64(if line.len() >= i1 + 19 { &line[i1..i1+19] } else { "" })?;
                vals[offset+2] = parse_rinex_f64(if line.len() >= i2 + 19 { &line[i2..i2+19] } else { "" })?;
                vals[offset+3] = parse_rinex_f64(if line.len() >= i3 + 19 { &line[i3..i3+19] } else { "" })?;
            }
            line_idx += 1;
            let max_lines = if current_constellation == Constellation::Glonass { 4 } else { 8 };
            if line_idx == max_lines {
                let sat = SatelliteId { constellation: current_constellation, prn: current_prn };
                if let Some(eph) = build_ephemeris(current_constellation, sat, current_toc, current_af0, current_af1, current_af2, &vals) {
                    ephemerides.push(eph);
                }
            }
        }
    }
    Ok(ephemerides)
}

#[cfg(test)]
mod test_nav_parser {
    use super::*;
    use std::fs;

    #[test]
    fn test_phone_nav() {
        let file = fs::File::open("../../datasets/rtkexplorer/rtklib-py/data/phone/Pixel4_GnssLog.nav").unwrap();
        let reader = std::io::BufReader::new(file);
        let eph = parse_rinex_nav(reader).unwrap();
        println!("Loaded {} ephemerides", eph.len());
        assert!(!eph.is_empty());
    }
}
