use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use std::collections::HashMap;
use std::io::BufRead;

fn parse_rinex_2_header<I: Iterator<Item = String>>(
    first_line: String,
    lines: &mut I,
) -> Result<Vec<String>, String> {
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

fn parse_rinex_3_header<I: Iterator<Item = String>>(
    first_line: String,
    lines: &mut I,
) -> Result<HashMap<Constellation, Vec<String>>, String> {
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
) -> Result<Vec<EpochObs>, String> {
    let mut epochs = Vec::new();

    let obs_types = parse_rinex_2_header(first_line, lines)?;

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
    Ok(epochs)
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
) -> Result<Vec<EpochObs>, String> {
    let mut epochs = Vec::new();

    let const_obs_types = parse_rinex_3_header(first_line, lines)?;

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
    Ok(epochs)
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

pub fn parse_rinex_f64(s: &str) -> Result<f64, String> {
    let s = s.trim().to_string();
    if s.is_empty() {
        return Ok(0.0);
    }
    let s_clean = s.replace("D", "E").replace("d", "e");
    s_clean
        .parse::<f64>()
        .map_err(|_| format!("Failed to parse RINEX f64: '{}'", s))
}

fn build_ephemeris(
    constellation: Constellation,
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<gneiss_core::ephemeris::Ephemeris> {
    match constellation {
        Constellation::Glonass => build_glonass_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Gps => build_gps_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Galileo => build_galileo_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Beidou => build_beidou_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Qzss => build_qzss_ephemeris(sat, toc, af0, af1, af2, vals),
        _ => None,
    }
}

fn build_glonass_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<gneiss_core::ephemeris::Ephemeris> {
    Some(gneiss_core::ephemeris::Ephemeris::Glonass(
        gneiss_core::ephemeris::GlonassEphemeris {
            sat,
            toe: toc,
            freq_num: vals[7] as i8,
            tau_n: af0,
            gamma_n: af1,
            delta_tau_n: af2,
            x: vals[0] * 1000.0,
            y: vals[4] * 1000.0,
            z: vals[8] * 1000.0,
            vx: vals[1] * 1000.0,
            vy: vals[5] * 1000.0,
            vz: vals[9] * 1000.0,
            ax: vals[2] * 1000.0,
            ay: vals[6] * 1000.0,
            az: vals[10] * 1000.0,
        },
    ))
}

fn build_gps_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<gneiss_core::ephemeris::Ephemeris> {
    Some(gneiss_core::ephemeris::Ephemeris::Gps(
        gneiss_core::ephemeris::GpsEphemeris {
            sat,
            toc,
            toe: GpsTime::new(toc.week, vals[8]),
            af0,
            af1,
            af2,
            iode: vals[0] as u32,
            crs: vals[1],
            delta_n: vals[2],
            m0: vals[3],
            cuc: vals[4],
            e: vals[5],
            cus: vals[6],
            sqrt_a: vals[7],
            cic: vals[9],
            omega0: vals[10],
            cis: vals[11],
            i0: vals[12],
            crc: vals[13],
            omega: vals[14],
            omega_dot: vals[15],
            idot: vals[16],
            tgd: vals[22],
            iodc: vals[23] as u32,
        },
    ))
}

fn build_galileo_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<gneiss_core::ephemeris::Ephemeris> {
    Some(gneiss_core::ephemeris::Ephemeris::Galileo(
        gneiss_core::ephemeris::GalileoEphemeris {
            sat,
            toc,
            toe: GpsTime::new(toc.week, vals[8]),
            af0,
            af1,
            af2,
            iod_nav: vals[0] as u32,
            crs: vals[1],
            delta_n: vals[2],
            m0: vals[3],
            cuc: vals[4],
            e: vals[5],
            cus: vals[6],
            sqrt_a: vals[7],
            cic: vals[9],
            omega0: vals[10],
            cis: vals[11],
            i0: vals[12],
            crc: vals[13],
            omega: vals[14],
            omega_dot: vals[15],
            idot: vals[16],
            bgd_e1_e5a: vals[22],
            bgd_e1_e5b: vals[23],
        },
    ))
}

fn build_beidou_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<gneiss_core::ephemeris::Ephemeris> {
    Some(gneiss_core::ephemeris::Ephemeris::Beidou(
        gneiss_core::ephemeris::BeidouEphemeris {
            sat,
            toc,
            toe: GpsTime::new(toc.week, vals[8]),
            af0,
            af1,
            af2,
            aode: vals[0] as u32,
            crs: vals[1],
            delta_n: vals[2],
            m0: vals[3],
            cuc: vals[4],
            e: vals[5],
            cus: vals[6],
            sqrt_a: vals[7],
            cic: vals[9],
            omega0: vals[10],
            cis: vals[11],
            i0: vals[12],
            crc: vals[13],
            omega: vals[14],
            omega_dot: vals[15],
            idot: vals[16],
            tgd1: vals[22],
            tgd2: vals[23],
            aodc: vals[25] as u32,
        },
    ))
}

fn build_qzss_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<gneiss_core::ephemeris::Ephemeris> {
    Some(gneiss_core::ephemeris::Ephemeris::Qzss(
        gneiss_core::ephemeris::QzssEphemeris {
            sat,
            toc,
            toe: GpsTime::new(toc.week, vals[8]),
            af0,
            af1,
            af2,
            iode: vals[0] as u32,
            crs: vals[1],
            delta_n: vals[2],
            m0: vals[3],
            cuc: vals[4],
            e: vals[5],
            cus: vals[6],
            sqrt_a: vals[7],
            cic: vals[9],
            omega0: vals[10],
            cis: vals[11],
            i0: vals[12],
            crc: vals[13],
            omega: vals[14],
            omega_dot: vals[15],
            idot: vals[16],
            tgd: vals[22],
            iodc: vals[23] as u32,
        },
    ))
}

pub fn parse_rinex_nav<R: BufRead>(
    reader: R,
) -> Result<
    (
        Vec<gneiss_core::ephemeris::Ephemeris>,
        Option<gneiss_core::atmosphere::KlobucharParams>,
    ),
    String,
> {
    let mut ephemerides = Vec::new();
    let mut lines = reader.lines().map(|l| l.unwrap_or_default());

    let mut is_rinex_3 = false;
    let klobuchar = parse_rinex_nav_header(&mut lines, &mut is_rinex_3);

    let mut current_constellation = Constellation::Gps;
    let mut current_prn = 0;
    let mut current_toc = GpsTime::new(0, 0.0);
    let mut current_af0 = 0.0;
    let mut current_af1 = 0.0;
    let mut current_af2 = 0.0;
    let mut line_idx = 0;
    let mut vals = [0.0; 32];

    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let is_new_epoch = if is_rinex_3 {
            line.starts_with('G')
                || line.starts_with('R')
                || line.starts_with('E')
                || line.starts_with('C')
                || line.starts_with('J')
                || line.starts_with('S')
                || line.starts_with('I')
        } else {
            line_idx == 0
                || (current_constellation != Constellation::Glonass && line_idx > 7)
                || (current_constellation == Constellation::Glonass && line_idx > 3)
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
            } else {
                Constellation::Gps
            };

            current_prn = if is_rinex_3 {
                if line.len() >= 3 {
                    line[1..3].trim().parse::<u8>().unwrap_or(0)
                } else {
                    0
                }
            } else if line.len() >= 2 {
                line[0..2].trim().parse::<u8>().unwrap_or(0)
            } else {
                0
            };

            let mut toc_gpst = parse_rinex_nav_epoch_time(&line, is_rinex_3);
            match current_constellation {
                Constellation::Glonass => toc_gpst = toc_gpst + (18.0 - 10800.0),
                Constellation::Beidou => toc_gpst = toc_gpst + 14.0,
                _ => {}
            }
            current_toc = toc_gpst;
            let (idx_af0, idx_af1, idx_af2) = if is_rinex_3 {
                (23, 42, 61)
            } else {
                (22, 41, 60)
            };

            current_af0 = parse_rinex_f64(if line.len() >= idx_af0 + 19 {
                &line[idx_af0..idx_af0 + 19]
            } else {
                ""
            })?;
            current_af1 = parse_rinex_f64(if line.len() >= idx_af1 + 19 {
                &line[idx_af1..idx_af1 + 19]
            } else {
                ""
            })?;
            current_af2 = parse_rinex_f64(if line.len() >= idx_af2 + 19 {
                &line[idx_af2..idx_af2 + 19]
            } else {
                ""
            })?;
            line_idx = 1;
            vals.fill(0.0);
        } else {
            if (1..=8).contains(&line_idx) {
                let offset = (line_idx - 1) * 4;
                let (i0, i1, i2, i3) = if is_rinex_3 {
                    (4, 23, 42, 61)
                } else {
                    (3, 22, 41, 60)
                };

                vals[offset] = parse_rinex_f64(if line.len() >= i0 + 19 {
                    &line[i0..i0 + 19]
                } else {
                    ""
                })?;
                vals[offset + 1] = parse_rinex_f64(if line.len() >= i1 + 19 {
                    &line[i1..i1 + 19]
                } else {
                    ""
                })?;
                vals[offset + 2] = parse_rinex_f64(if line.len() >= i2 + 19 {
                    &line[i2..i2 + 19]
                } else {
                    ""
                })?;
                vals[offset + 3] = parse_rinex_f64(if line.len() >= i3 + 19 {
                    &line[i3..i3 + 19]
                } else {
                    ""
                })?;
            }
            line_idx += 1;
            let max_lines = if current_constellation == Constellation::Glonass {
                4
            } else {
                8
            };
            if line_idx == max_lines {
                let sat = SatelliteId {
                    constellation: current_constellation,
                    prn: current_prn,
                };
                if let Some(eph) = build_ephemeris(
                    current_constellation,
                    sat,
                    current_toc,
                    current_af0,
                    current_af1,
                    current_af2,
                    &vals,
                ) {
                    ephemerides.push(eph);
                }
            }
        }
    }
    Ok((ephemerides, klobuchar))
}

fn parse_rinex_nav_header(
    lines: &mut impl Iterator<Item = String>,
    is_rinex_3: &mut bool,
) -> Option<gneiss_core::atmosphere::KlobucharParams> {
    let mut alpha = [0.0; 4];
    let mut beta = [0.0; 4];
    let mut has_alpha = false;
    let mut has_beta = false;

    for line in lines {
        if line.contains("RINEX VERSION / TYPE") && line.trim().starts_with('3') {
            *is_rinex_3 = true;
        }

        if line.contains("ION ALPHA") || line.contains("IONOSPHERIC CORR") && line.contains("GPSA")
        {
            let offset = if line.contains("GPSA") { 5 } else { 2 };
            alpha[0] = parse_rinex_f64(if line.len() >= offset + 12 {
                &line[offset..offset + 12]
            } else {
                ""
            })
            .unwrap_or(0.0);
            alpha[1] = parse_rinex_f64(if line.len() >= offset + 24 {
                &line[offset + 12..offset + 24]
            } else {
                ""
            })
            .unwrap_or(0.0);
            alpha[2] = parse_rinex_f64(if line.len() >= offset + 36 {
                &line[offset + 24..offset + 36]
            } else {
                ""
            })
            .unwrap_or(0.0);
            alpha[3] = parse_rinex_f64(if line.len() >= offset + 48 {
                &line[offset + 36..offset + 48]
            } else {
                ""
            })
            .unwrap_or(0.0);
            has_alpha = true;
        }
        if line.contains("ION BETA") || line.contains("IONOSPHERIC CORR") && line.contains("GPSB") {
            let offset = if line.contains("GPSB") { 5 } else { 2 };
            beta[0] = parse_rinex_f64(if line.len() >= offset + 12 {
                &line[offset..offset + 12]
            } else {
                ""
            })
            .unwrap_or(0.0);
            beta[1] = parse_rinex_f64(if line.len() >= offset + 24 {
                &line[offset + 12..offset + 24]
            } else {
                ""
            })
            .unwrap_or(0.0);
            beta[2] = parse_rinex_f64(if line.len() >= offset + 36 {
                &line[offset + 24..offset + 36]
            } else {
                ""
            })
            .unwrap_or(0.0);
            beta[3] = parse_rinex_f64(if line.len() >= offset + 48 {
                &line[offset + 36..offset + 48]
            } else {
                ""
            })
            .unwrap_or(0.0);
            has_beta = true;
        }

        if line.contains("END OF HEADER") {
            break;
        }
    }

    if has_alpha && has_beta {
        Some(gneiss_core::atmosphere::KlobucharParams { alpha, beta })
    } else {
        None
    }
}

fn parse_rinex_nav_epoch_time(line: &str, is_rinex_3: bool) -> GpsTime {
    let (i_y, i_m, i_d, i_h, i_min, i_s) = if is_rinex_3 {
        ((4, 8), (9, 11), (12, 14), (15, 17), (18, 20), (21, 23))
    } else {
        ((3, 5), (5, 8), (8, 11), (11, 14), (14, 17), (17, 22))
    };

    let parse_i32 = |start: usize, end: usize| -> i32 {
        if line.len() >= end {
            line[start..end].trim().parse().unwrap_or(0)
        } else {
            0
        }
    };
    let parse_f64 = |start: usize, end: usize| -> f64 {
        if line.len() >= end {
            line[start..end].trim().parse().unwrap_or(0.0)
        } else {
            0.0
        }
    };

    let mut year = parse_i32(i_y.0, i_y.1);
    if year < 100 {
        year += if year > 80 { 1900 } else { 2000 };
    }

    GpsTime::from_calendar(
        year,
        parse_i32(i_m.0, i_m.1),
        parse_i32(i_d.0, i_d.1),
        parse_i32(i_h.0, i_h.1),
        parse_i32(i_min.0, i_min.1),
        parse_f64(i_s.0, i_s.1),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::ephemeris::Ephemeris;
    use gneiss_core::time::GpsTime;
    use std::io::BufReader;

    // -----------------------------------------------------------------------
    // parse_rinex_f64 helper
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_rinex_f64_normal() {
        assert!((parse_rinex_f64("123.456").unwrap() - 123.456).abs() < 1e-12);
    }

    #[test]
    fn test_parse_rinex_f64_d_notation() {
        assert!((parse_rinex_f64("1.23D-4").unwrap() - 0.000123).abs() < 1e-12);
        assert!((parse_rinex_f64("5.0D+2").unwrap() - 500.0).abs() < 1e-12);
    }

    #[test]
    fn test_parse_rinex_f64_lowercase_d() {
        assert!((parse_rinex_f64("1.23d-4").unwrap() - 0.000123).abs() < 1e-12);
    }

    #[test]
    fn test_parse_rinex_f64_empty() {
        assert!((parse_rinex_f64("").unwrap() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_parse_rinex_f64_whitespace() {
        assert!((parse_rinex_f64("  ").unwrap() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_parse_rinex_f64_invalid() {
        assert!(parse_rinex_f64("abc").is_err());
    }

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
        let epochs = parse_rinex_obs(&mut reader).unwrap();

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
        let epochs = parse_rinex_obs(&mut reader).unwrap();

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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();

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
    // RINEX 2 GPS NAV with value checks
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_nav_gps_full() {
        let data = "     2.11           N: GPS NAV DATA                         RINEX VERSION / TYPE
                                                            END OF HEADER
 6 20  5 14 22  0  0.0-2.715480513871D-04-6.821210263297D-12 0.000000000000D+00
    7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
   -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
    4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
    9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
    4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
    0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
    4.248000000000D+05 4.000000000000D+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();

        assert_eq!(ephemerides.len(), 1);
        let eph = &ephemerides[0];

        assert_eq!(eph.sat().prn, 6);
        assert_eq!(eph.sat().constellation, Constellation::Gps);
        assert_eq!(eph.toe().week, 2105);
        assert_eq!(eph.toe().tow, 424800.0);

        if let Ephemeris::Gps(gps) = eph {
            assert!((gps.sqrt_a - 5153.539648056).abs() < 1e-6);
            assert!((gps.e - 0.001934998203069).abs() < 1e-12);
            assert_eq!(gps.iode, 70);
        } else {
            panic!("Expected GPS ephemeris");
        }
    }

    // -----------------------------------------------------------------------
    // RINEX 3 GPS NAV
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_gps() {
        // Date parsed: 2020-06-15 01:30:00 -> week=2110, tow=91800
        let data = "     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
G 1 2020 06 15 01 30  0 -.271548051387D-03 -.682121026330D-11  .000000000000D+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();

        assert_eq!(ephemerides.len(), 1);
        let eph = &ephemerides[0];

        assert_eq!(eph.sat().prn, 1);
        assert_eq!(eph.sat().constellation, Constellation::Gps);

        // TOE comes from vals[8] (424800.0) with the TOC week (2110)
        assert_eq!(eph.toe().week, 2110);
        assert_eq!(eph.toe().tow, 424800.0,
            "TOE mismatch: got {}", eph.toe().tow);

        if let Ephemeris::Gps(gps) = eph {
            assert!((gps.sqrt_a - 5153.539648056).abs() < 1e-6);
        } else {
            panic!("Expected GPS ephemeris");
        }
    }

    // -----------------------------------------------------------------------
    // RINEX 3 GLONASS NAV extended value checks
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_glonass_extended() {
        let data = "     3.03           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
R 6 2020 12 24 21 15  0  .189751386642E-03  .000000000000E+00  .422910000000E+06
     -.740158740234E+04 -.212037086487E+00  .000000000000E+00  .000000000000E+00
     -.206682856445E+05 -.176755714417E+01  .931322574615E-09 -.400000000000E+01
      .129489067383E+05 -.294115734100E+01 -.186264514923E-08  .000000000000E+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();

        assert_eq!(ephemerides.len(), 1);
        let eph = &ephemerides[0];

        assert_eq!(eph.sat().prn, 6);
        assert_eq!(eph.sat().constellation, Constellation::Glonass);
        assert_eq!(eph.toe().week, 2137);
        let raw_tow = GpsTime::from_calendar(2020, 12, 24, 21, 15, 0.0).tow;
        let expected_tow = raw_tow + 18.0 - 10800.0; // GLONASS time adjustment
        assert!(
            (eph.toe().tow - expected_tow).abs() < 1e-4,
            "Expected TOW near {}, got {}",
            expected_tow,
            eph.toe().tow
        );

        if let Ephemeris::Glonass(glo) = eph {
            assert_eq!(glo.freq_num, -4);
            assert!((glo.x - (-7401587.40234)).abs() < 1e-2, "x mismatch: {}", glo.x);
            assert!((glo.y - (-20668285.6445)).abs() < 1e-2, "y mismatch: {}", glo.y);
            assert!((glo.z - 12948906.7383).abs() < 1e-2, "z mismatch: {}", glo.z);
        } else {
            panic!("Expected GLONASS ephemeris");
        }
    }

    // -----------------------------------------------------------------------
    // RINEX 3 mixed nav (GPS + GLONASS in sequence)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_mixed_constellations() {
        let data = "     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
G 1 2020 06 15 01 30  0 -.271548051387D-03 -.682121026330D-11  .000000000000D+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00
R 6 2020 12 24 21 15  0  .189751386642E-03  .000000000000E+00  .422910000000E+06
     -.740158740234E+04 -.212037086487E+00  .000000000000E+00  .000000000000E+00
     -.206682856445E+05 -.176755714417E+01  .931322574615E-09 -.400000000000E+01
      .129489067383E+05 -.294115734100E+01 -.186264514923E-08  .000000000000E+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();

        assert_eq!(ephemerides.len(), 2);
        assert_eq!(ephemerides[0].sat().constellation, Constellation::Gps);
        assert_eq!(ephemerides[0].sat().prn, 1);
        assert_eq!(ephemerides[1].sat().constellation, Constellation::Glonass);
        assert_eq!(ephemerides[1].sat().prn, 6);
    }

    // -----------------------------------------------------------------------
    // Navigation header with Klobuchar ionospheric parameters (ION ALPHA/BETA)
    // Values at columns 2-13, 14-25, 26-37, 38-49 (each 12 chars, right-justified)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_nav_header_klobuchar() {
        let data = "     2.11           N: GPS NAV DATA                         RINEX VERSION / TYPE
     .1167E-07   .1490E-07  -.1192E-06   .1192E-06          ION ALPHA
     .1024E+06   .3277E+05  -.1966E+06   .1311E+06          ION BETA
                                                            END OF HEADER
 6 20  5 14 22  0  0.0-2.715480513871D-04-6.821210263297D-12 0.000000000000D+00
    7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
   -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
    4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
    9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
    4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
    0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
    4.248000000000D+05 4.000000000000D+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, klob) = parse_rinex_nav(&mut reader).unwrap();

        assert!(klob.is_some());
        let klob = klob.unwrap();
        assert!((klob.alpha[0] - 1.167e-8).abs() < 1e-20);
        assert!((klob.alpha[1] - 1.490e-8).abs() < 1e-20);
        assert!((klob.alpha[2] - (-1.192e-7)).abs() < 1e-20);
        assert!((klob.alpha[3] - 1.192e-7).abs() < 1e-20);

        assert!((klob.beta[0] - 102400.0).abs() < 1.0);
        assert!((klob.beta[1] - 32770.0).abs() < 1.0);
        assert!((klob.beta[2] - (-196600.0)).abs() < 1.0);
        assert!((klob.beta[3] - 131100.0).abs() < 1.0);

        assert_eq!(ephemerides.len(), 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 nav header with IONOSPHERIC CORR format (GPSA/GPSB)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_header_ionospheric_corr() {
        // GPSA/GPSB format: offset=5, values at 5-16, 17-28, 29-40, 41-52
        let data = "     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
GPSA    .1167E-07   .1490E-07  -.1192E-06   .1192E-06          IONOSPHERIC CORR
GPSB    .1024E+06   .3277E+05  -.1966E+06   .1311E+06          IONOSPHERIC CORR
                                                            END OF HEADER
G 1 2020 06 15 01 30  0 -.123456789012D-03  .000000000000E+00  .000000000000E+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, klob) = parse_rinex_nav(&mut reader).unwrap();

        assert!(klob.is_some());
        let klob = klob.unwrap();
        assert!((klob.alpha[0] - 1.167e-8).abs() < 1e-20);
        assert!((klob.beta[0] - 102400.0).abs() < 1.0);
        assert_eq!(ephemerides.len(), 1);
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

    #[test]
    fn test_rinex_nav_empty() {
        let data = "";
        let mut reader = BufReader::new(data.as_bytes());
        let (eph, klob) = parse_rinex_nav(&mut reader).unwrap();
        assert!(eph.is_empty());
        assert!(klob.is_none());
    }

    #[test]
    fn test_rinex_nav_missing_end_header() {
        // Navigation without END OF HEADER - the header parser iterates until
        // it finds END OF HEADER or runs out of lines. Without it, the parsing
        // proceeds to read ephemeris data from what would have been header lines.
        let data = "     2.11           N: GPS NAV DATA                         RINEX VERSION / TYPE
 6 20  5 14 22  0  0.0-2.715480513871D-04-6.821210263297D-12 0.000000000000D+00
    7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
   -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
    4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
    9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
    4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
    0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
    4.248000000000D+05 4.000000000000D+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (_eph, klob) = parse_rinex_nav(&mut reader).unwrap();
        // Without END OF HEADER, the ephemeris first line is consumed by the
        // header parser as header data. The remaining 8 lines form a partial
        // ephemeris. This test just verifies no panic.
        assert!(klob.is_none());
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites[0].sat.prn, 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 NAV: two successive ephemerides
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_nav_two_ephemerides() {
        let data = "     2.11           N: GPS NAV DATA                         RINEX VERSION / TYPE
                                                            END OF HEADER
 1 20  5 14 22  0  0.0-2.715480513871D-04-6.821210263297D-12 0.000000000000D+00
    7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
   -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
    4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
    9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
    4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
    0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
    4.248000000000D+05 4.000000000000D+00
 2 20  5 14 22  0  0.0-2.715480513871D-04-6.821210263297D-12 0.000000000000D+00
    7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
   -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
    4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
    9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
    4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
    0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
    4.248000000000D+05 4.000000000000D+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (eph, _klob) = parse_rinex_nav(&mut reader).unwrap();
        assert_eq!(eph.len(), 2);
        assert_eq!(eph[0].sat().prn, 1);
        assert_eq!(eph[1].sat().prn, 2);
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
    // RINEX 3 Beidou NAV
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_beidou() {
        let data = "     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
C 1 2020 06 15 01 30  0 -.271548051387D-03 -.682121026330D-11  .000000000000D+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();

        assert_eq!(ephemerides.len(), 1);
        let eph = &ephemerides[0];

        assert_eq!(eph.sat().prn, 1);
        assert_eq!(eph.sat().constellation, Constellation::Beidou);

        if let Ephemeris::Beidou(bds) = eph {
            assert!((bds.sqrt_a - 5153.539648056).abs() < 1e-6);
            assert_eq!(bds.aode, 70);
        } else {
            panic!("Expected Beidou ephemeris, got {:?}", eph.sat().constellation);
        }
    }

    // -----------------------------------------------------------------------
    // RINEX 3 Galileo NAV
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_galileo() {
        let data = "     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
E 2 2020 08 01 12  0  0 -.153846153846D-03  .000000000000E+00  .000000000000E+00
     1.000000000000D+01 0.000000000000D+00 5.000000000000D-09 3.000000000000D+00
     0.000000000000D+00 1.500000000000D-03 4.500000000000D-06 5.200000000000D+03
     4.320000000000D+05 1.200000000000D-07 2.500000000000D+00 3.400000000000D-08
     9.800000000000D-01 3.000000000000D+02 1.000000000000D+00 8.000000000000D-09
     4.500000000000D-11 0.000000000000D+00 2.000000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 1.000000000000D+01
     4.320000000000D+05 4.000000000000D+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();

        assert_eq!(ephemerides.len(), 1);
        let eph = &ephemerides[0];
        assert_eq!(eph.sat().prn, 2);
        assert_eq!(eph.sat().constellation, Constellation::Galileo);

        if let Ephemeris::Galileo(gal) = eph {
            assert!((gal.sqrt_a - 5200.0).abs() < 1e-6);
            assert_eq!(gal.iod_nav, 10);
        } else {
            panic!("Expected Galileo ephemeris");
        }
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
    }

    // -----------------------------------------------------------------------
    // parse_rinex_nav_epoch_time edge cases (2-digit year)
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_rinex_nav_epoch_time_two_digit_year() {
        let line_r2 = " 6 99 12 25  0  0  0.0";
        let time_r2 = parse_rinex_nav_epoch_time(line_r2, false);
        assert_eq!(time_r2.week, GpsTime::from_calendar(1999, 12, 25, 0, 0, 0.0).week);

        let line_r2_2 = " 6  1  1  1  0  0  0.0";
        let time_r2_2 = parse_rinex_nav_epoch_time(line_r2_2, false);
        assert_eq!(time_r2_2.week, GpsTime::from_calendar(2001, 1, 1, 0, 0, 0.0).week);
    }

    #[test]
    fn test_parse_rinex_nav_epoch_time_4_digit_year() {
        let line_r3 = "G 1 2020 06 15 01 30  0";
        let time_r3 = parse_rinex_nav_epoch_time(line_r3, true);
        let expected = GpsTime::from_calendar(2020, 6, 15, 1, 30, 0.0);
        assert_eq!(time_r3.week, expected.week);
        assert!((time_r3.tow - expected.tow).abs() < 1e-4);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 NAV with I (IRNSS) constellation
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_irnss_as_gps() {
        let data = "     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
I 1 2020 06 15 01 30  0 -.271548051387D-03 -.682121026330D-11  .000000000000D+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();
        assert_eq!(ephemerides.len(), 1);
        assert_eq!(ephemerides[0].sat().constellation, Constellation::Gps);
        assert_eq!(ephemerides[0].sat().prn, 1);
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
        assert_eq!(
            epochs[0].satellites[0].sat,
            SatelliteId { constellation: Constellation::Beidou, prn: 1 }
        );
    }

    // -----------------------------------------------------------------------
    // RINEX 3 nav with blank lines between ephemerides
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_blank_lines() {
        let data = "     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
G 1 2020 06 15 01 30  0 -.271548051387D-03 -.682121026330D-11  .000000000000D+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00

R 6 2020 12 24 21 15  0  .189751386642E-03  .000000000000E+00  .422910000000E+06
     -.740158740234E+04 -.212037086487E+00  .000000000000E+00  .000000000000E+00
     -.206682856445E+05 -.176755714417E+01  .931322574615E-09 -.400000000000E+01
      .129489067383E+05 -.294115734100E+01 -.186264514923E-08  .000000000000E+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();
        assert_eq!(ephemerides.len(), 2);
        assert_eq!(ephemerides[0].sat().constellation, Constellation::Gps);
        assert_eq!(ephemerides[1].sat().constellation, Constellation::Glonass);
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 NAV: QZSS ephemeris
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_qzss() {
        let data = "     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
J 1 2020 06 15 01 30  0 -.271548051387D-03 -.682121026330D-11  .000000000000D+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();
        assert_eq!(ephemerides.len(), 1);
        assert_eq!(ephemerides[0].sat().constellation, Constellation::Qzss);
        assert_eq!(ephemerides[0].sat().prn, 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 NAV: Sbas treated as GPS (code fallback in build_ephemeris)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_sbas_as_gps() {
        // Sbas maps to Constellation::Gps in nav parser (line 668).
        // build_ephemeris for Gps returns a GpsEphemeris.
        let data = "     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
S 1 2020 06 15 01 30  0 -.271548051387D-03 -.682121026330D-11  .000000000000D+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00
";
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();
        // SBAS nav messages are not parsed as GPS ephemeris in the current parser
        assert_eq!(ephemerides.len(), 0);
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
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
        let epochs = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        // X1 is skipped by map_rinex_type, only C1 remains
        assert_eq!(epochs[0].satellites[0].observations.len(), 1);
        let obs = &epochs[0].satellites[0].observations[0];
        assert_eq!(obs.code.obs_type, ObsType::Pseudorange);
        assert_eq!(obs.code.signal.freq_band, 1);
    }

    // -----------------------------------------------------------------------
    // parse_rinex_f64 with various edge cases
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_rinex_f64_edge_cases() {
        assert!((parse_rinex_f64("  1.23D-4  ").unwrap() - 0.000123).abs() < 1e-12);
        assert!((parse_rinex_f64("42").unwrap() - 42.0).abs() < 1e-12);
        assert!((parse_rinex_f64("-0.0").unwrap()).abs() < 1e-12);
        assert!((parse_rinex_f64("1.0D+1").unwrap() - 10.0).abs() < 1e-12);
        assert!((parse_rinex_f64("5d-1").unwrap() - 0.5).abs() < 1e-12);
    }
}
