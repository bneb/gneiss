//! RINEX navigation message parsing for GPS, GLONASS, Galileo, BeiDou, and QZSS.

pub mod builder;
#[cfg(test)]
mod tests;

use builder::build_ephemeris;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use std::io::BufRead;

pub fn parse_rinex_f64(s: &str) -> Result<f64, String> {
    let s = s.trim().to_string();
    if s.is_empty() {
        return Ok(0.0);
    }
    let s_clean = s.replace('D', "E").replace('d', "e");
    s_clean
        .parse::<f64>()
        .map_err(|_| format!("Failed to parse RINEX f64: '{}'", s))
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
                match line.chars().next().expect("RINEX line is non-empty") {
                    'G' => Constellation::Gps,
                    'R' => Constellation::Glonass,
                    'E' => Constellation::Galileo,
                    'C' => Constellation::Beidou,
                    'J' => Constellation::Qzss,
                    'S' => Constellation::Sbas,
                    'I' => Constellation::Navic,
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
                Constellation::Glonass => toc_gpst = toc_gpst + gneiss_core::gnss_time::TimeSystem::Glonass.gpst_offset(),
                Constellation::Beidou => toc_gpst = toc_gpst + gneiss_core::gnss_time::TimeSystem::Bdt.gpst_offset(),
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

        if line.contains("ION ALPHA") || line.contains("IONOSPHERIC CORR") && line.contains("GPSA") {
            let offset = if line.contains("GPSA") { 5 } else { 2 };
            alpha[0] = parse_rinex_f64(if line.len() >= offset + 12 {
                &line[offset..offset + 12]
            } else {
                ""
            }).unwrap_or(0.0);
            alpha[1] = parse_rinex_f64(if line.len() >= offset + 24 {
                &line[offset + 12..offset + 24]
            } else {
                ""
            }).unwrap_or(0.0);
            alpha[2] = parse_rinex_f64(if line.len() >= offset + 36 {
                &line[offset + 24..offset + 36]
            } else {
                ""
            }).unwrap_or(0.0);
            alpha[3] = parse_rinex_f64(if line.len() >= offset + 48 {
                &line[offset + 36..offset + 48]
            } else {
                ""
            }).unwrap_or(0.0);
            has_alpha = true;
        }
        if line.contains("ION BETA") || line.contains("IONOSPHERIC CORR") && line.contains("GPSB") {
            let offset = if line.contains("GPSB") { 5 } else { 2 };
            beta[0] = parse_rinex_f64(if line.len() >= offset + 12 {
                &line[offset..offset + 12]
            } else {
                ""
            }).unwrap_or(0.0);
            beta[1] = parse_rinex_f64(if line.len() >= offset + 24 {
                &line[offset + 12..offset + 24]
            } else {
                ""
            }).unwrap_or(0.0);
            beta[2] = parse_rinex_f64(if line.len() >= offset + 36 {
                &line[offset + 24..offset + 36]
            } else {
                ""
            }).unwrap_or(0.0);
            beta[3] = parse_rinex_f64(if line.len() >= offset + 48 {
                &line[offset + 36..offset + 48]
            } else {
                ""
            }).unwrap_or(0.0);
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
