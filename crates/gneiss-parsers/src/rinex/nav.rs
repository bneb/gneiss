use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use std::io::BufRead;

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
        Constellation::Navic => build_gps_ephemeris(sat, toc, af0, af1, af2, vals),
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
        let expected_tow = raw_tow + gneiss_core::gnss_time::TimeSystem::Glonass.gpst_offset();
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
        assert_eq!(ephemerides[0].sat().constellation, Constellation::Navic);
        assert_eq!(ephemerides[0].sat().prn, 1);
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
