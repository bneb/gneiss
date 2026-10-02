//! RINEX navigation message parsing for GPS, GLONASS, Galileo, BeiDou, and QZSS.

pub mod builder;
#[cfg(test)]
mod tests;

use builder::build_ephemeris;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use std::io::BufRead;

/// Reads one 19-byte RINEX nav value column.
///
/// RINEX 3 producers legitimately stop a record early instead of padding out
/// its trailing empty fields - the BRDC file has 9456 Galileo records whose
/// last line is 23 columns long - so a line that does not reach the column is
/// read as blank rather than treated as a truncated record.
fn field(line: &str, idx: usize) -> Result<f64, String> {
    match line.as_bytes().get(idx..idx + 19) {
        Some(bytes) => parse_rinex_f64(&String::from_utf8_lossy(bytes)),
        None => Ok(0.0),
    }
}

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
    let mut rinex_2_system = Constellation::Gps;
    let klobuchar = parse_rinex_nav_header(&mut lines, &mut is_rinex_3, &mut rinex_2_system);

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
                // RINEX 2 declares the satellite system once, in column 41 of
                // the version line; the records themselves carry only a
                // two-digit PRN. Assuming GPS here turns a GLONASS file into
                // misaligned eight-line GPS records.
                rinex_2_system
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

            current_af0 = field(&line, idx_af0)?;
            current_af1 = field(&line, idx_af1)?;
            current_af2 = field(&line, idx_af2)?;
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

                vals[offset] = field(&line, i0)?;
                vals[offset + 1] = field(&line, i1)?;
                vals[offset + 2] = field(&line, i2)?;
                vals[offset + 3] = field(&line, i3)?;
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
    rinex_2_system: &mut Constellation,
) -> Option<gneiss_core::atmosphere::KlobucharParams> {
    let mut alpha = [0.0; 4];
    let mut beta = [0.0; 4];
    let mut has_alpha = false;
    let mut has_beta = false;

    for line in lines {
        if line.contains("RINEX VERSION / TYPE") {
            if line.trim().starts_with('3') {
                *is_rinex_3 = true;
            } else if let Some(c) = line.get(40..41) {
                *rinex_2_system = match c {
                    "R" => Constellation::Glonass,
                    "E" => Constellation::Galileo,
                    "C" => Constellation::Beidou,
                    "J" => Constellation::Qzss,
                    "G" => Constellation::Gps,
                    // Blank or 'M' (mixed): GPS is the only safe default for
                    // the 2-digit-PRN records this crate already accepts.
                    _ => Constellation::Gps,
                };
            }
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
        year += if year >= 80 { 1900 } else { 2000 };
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
#[allow(clippy::unwrap_used)]
mod regression_tests {
    use super::*;
    use gneiss_core::ephemeris::Ephemeris;
    use std::io::BufReader;

    /// Two RINEX 2 GLONASS broadcast records, re-gridded from the real
    /// records of datasets/wtzr_ppp_1224/BRDC00IGS_R_20203590000_01D_MN.rnx
    /// (R01 at 00:15 and 00:45 UTC(SU)). Column 41 (byte 40) of the version
    /// line carries the satellite system: 'R' = GLONASS (RTKLIB rinex.cc:587
    /// reads `*(buff+40)`, mapping 'R' to SYS_GLO). A GLONASS record is four
    /// lines: three values on the epoch line plus 4*3 on the continuations,
    /// which reaches data index 15 (RTKLIB rinex.cc:1307 `i>=15`).
    const RINEX2_GLONASS: &str = concat!(
        "     2.11           G: GLONASS NAV DATA R                   RINEX VERSION / TYPE\n",
        "                                                            END OF HEADER\n",
        "01 20 12 24 00 15   00 7.246062159538E-05 0.000000000000E+00 3.464400000000E+05\n",
        "   -1.184028808594E+03-2.155310630798E+00 0.000000000000E+00 0.000000000000E+00\n",
        "    1.322943017578E+04-2.060539245605E+00 0.000000000000E+00 1.000000000000E+00\n",
        "    2.176970068359E+04 1.135528564453E+00-1.862645149231E-09 0.000000000000E+00\n",
        "01 20 12 24 00 45   00 7.246155291796E-05 0.000000000000E+00 3.478800000000E+05\n",
        "   -5.458924804688E+03-2.560062408447E+00 0.000000000000E+00 0.000000000000E+00\n",
        "    9.701321777344E+03-1.821957588196E+00 0.000000000000E+00 1.000000000000E+00\n",
        "    2.294521191406E+04 1.621379852295E-01-1.862645149231E-09 0.000000000000E+00\n",
    );

    /// A RINEX 2 navigation file can hold GLONASS records, and the file tells
    /// the reader so in column 41. `parse_rinex_nav` ignores that column and
    /// hard-codes `Constellation::Gps` (nav/mod.rs), then consumes eight lines
    /// per record instead of four, so the two four-line GLONASS records below
    /// are merged into one GPS ephemeris built from misaligned fields.
    #[test]
    fn rinex2_glonass_records_are_decoded_as_glonass() {
        let mut reader = BufReader::new(RINEX2_GLONASS.as_bytes());
        let (eph, _klob) = parse_rinex_nav(&mut reader)
            .expect("a valid RINEX 2 GLONASS nav file must parse");
        assert_eq!(eph.len(), 2, "four-line GLONASS records merged");
        for e in &eph {
            assert_eq!(e.sat().constellation, Constellation::Glonass);
            assert_eq!(e.sat().prn, 1, "GLONASS slot read as PRN");
        }
    }

    /// Values of the first record: PZ-90 X/Y/Z in km converted to metres, and
    /// the frequency number from field 11 of the record (RTKLIB rinex.cc:1174
    /// `geph->frq=(int)data[10]`).
    #[test]
    fn rinex2_glonass_position_and_frequency_match_the_record() {
        let mut reader = BufReader::new(RINEX2_GLONASS.as_bytes());
        let (eph, _klob) = parse_rinex_nav(&mut reader)
            .expect("a valid RINEX 2 GLONASS nav file must parse");
        let Ephemeris::Glonass(g) = &eph[0] else {
            panic!("expected GLONASS");
        };
        assert!((g.x - -1184028.808594).abs() < 1e-3, "x = {}", g.x);
        assert!((g.y - 13229430.17578).abs() < 1e-3, "y = {}", g.y);
        assert!((g.z - 21769700.68359).abs() < 1e-3, "z = {}", g.z);
        assert!((g.vx - -2155.310630798).abs() < 1e-6, "vx = {}", g.vx);
        assert_eq!(g.freq_num, 1);
        // |r| = sqrt(1184.0288^2 + 13229.4302^2 + 21769.7007^2) = 25501.76 km,
        // the GLONASS orbital radius; a misaligned read cannot land there.
        let r = (g.x * g.x + g.y * g.y + g.z * g.z).sqrt();
        assert!((r - 25501757.09).abs() < 1.0, "|r| = {r}");
    }

    /// GLONASS toc is UTC(SU). GPST = UTC(SU) - 10800 s + 18 leap seconds
    /// = UTC(SU) - 10782 s (gnss_time.rs GLONASS_TO_GPST).
    #[test]
    fn rinex2_glonass_toc_is_converted_from_utc_su() {
        let mut reader = BufReader::new(RINEX2_GLONASS.as_bytes());
        let (eph, _klob) = parse_rinex_nav(&mut reader)
            .expect("a valid RINEX 2 GLONASS nav file must parse");
        let expected = GpsTime::from_calendar(2020, 12, 24, 0, 15, 0.0) - 10_782.0;
        assert!((eph[0].toe().tow - expected.tow).abs() < 1e-6);
        assert_eq!(eph[0].toe().week, expected.week);
    }

    /// RINEX 2 carries a two-digit year. The observation parser resolves the
    /// 80..99 window to 19xx (`>= 80` in rinex2.rs); the navigation parser
    /// uses `> 80`, so the same text decodes to a century different by 100
    /// years. Both parsers must agree - this asserts the navigation side lands
    /// on the same instant as the observation side for the same file year.
    #[test]
    fn two_digit_year_80_is_the_same_century_in_nav_and_obs() {
        let nav = " 6 80 12  1  0  0  0.0-2.715480513871D-04-6.821210263297D-12 0.000000000000D+00\n";
        let t_nav = parse_rinex_nav_epoch_time(nav, false);
        let obs = " 80 12  1  0  0  0.0000000  0  0";
        let t_obs = GpsTime::from_calendar(
            {
                let y: i32 = obs[1..3].trim().parse().unwrap();
                if y >= 80 { 1900 + y } else { 2000 + y }
            },
            12,
            1,
            0,
            0,
            0.0,
        );
        assert_eq!((t_nav.week, t_nav.tow), (t_obs.week, t_obs.tow));
    }

    /// A record whose continuation lines are too short to fill the 19-byte
    /// value columns must not produce an ephemeris. `parse_rinex_nav` feeds
    /// `parse_rinex_f64("")` for every missing field, which yields 0.0, and
    /// still emits a "valid" Ephemeris.
    ///
    /// A GPS orbit has a semi-major axis of 26560 km (MEO) to 42164 km
    /// (IGSO/geo), i.e. sqrt(A) between sqrt(26560) = 162.97 km and
    /// sqrt(42164) = 205.34 km only after squaring: sqrt(A) is quoted in
    /// km^(1/2), so the broadcast value lies between 5153.6 and 6353.3. No
    /// genuine ephemeris can be outside 5100..=6400.
    #[test]
    fn truncated_record_does_not_emit_a_zero_filled_ephemeris() {
        let data = concat!(
            "     3.04           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE\n",
            "                                                            END OF HEADER\n",
            "G01 2020 12 24 00 00 00 7.914672605693E-04-5.570655048359E-12 0.000000000000E+00\n",
            "     5.100000000000E+01-2.028125000000E+01 3.843017219950E-09-1.307172517279E+00\n",
            "    -1.067295670509E-06 1.020454068203E-02\n",
            "     3.456000000000E+05 2.942979335785E-07-6.828506108746E-01 5.401670932770E-08\n",
            "     9.828138365184E-01 2.077187500000E+02 8.232702817232E-01-7.690320332710E-09\n",
            "     1.207193141583E-10 1.000000000000E+00 2.137000000000E+03 0.000000000000E+00\n",
            "     2.000000000000E+00 0.000000000000E+00 5.122274160385E-09 5.100000000000E+01\n",
            "     3.384180000000E+05 4.000000000000E+00                                           \n",
        );
        let mut reader = BufReader::new(data.as_bytes());
        let (eph, _klob) = parse_rinex_nav(&mut reader).unwrap();
        for e in &eph {
            let Ephemeris::Gps(g) = e else { panic!("expected GPS") };
            assert!(
                (5100.0..=6400.0).contains(&g.sqrt_a),
                "zero-filled ephemeris leaked out: sqrt_a = {}",
                g.sqrt_a
            );
        }
    }

    /// Differential check against
    /// datasets/wtzr_ppp_1224/BRDC00IGS_R_20203590000_01D_MN.rnx
    /// (RINEX 3.04 mixed broadcast file). Record counts come from an
    /// independent walk of the file that classifies each record by its
    /// 3-byte system+PRN prefix: GPS 418, Galileo 10427, GLONASS 1150,
    /// BeiDou 1055, QZSS 96, SBAS 8220 (SBAS has no Ephemeris variant and is
    /// dropped by design), so 13146 ephemerides must come out.
    #[test]
    fn real_rinex3_nav_record_counts_match_an_independent_walk() {
        let f = std::fs::File::open(
            "../../datasets/wtzr_ppp_1224/BRDC00IGS_R_20203590000_01D_MN.rnx",
        )
        .unwrap();
        let (eph, _klob) = parse_rinex_nav(BufReader::new(f)).unwrap();
        assert_eq!(eph.len(), 13146);
        assert_eq!(
            eph.iter().filter(|e| e.sat().constellation == Constellation::Glonass).count(),
            1150
        );
        assert_eq!(
            eph.iter().filter(|e| e.sat().constellation == Constellation::Gps).count(),
            418
        );
    }

    /// Golden vector: first record of the same file, verbatim.
    /// G01, toc 2020-12-24 00:00 GPST (week 2137, tow 345600).
    #[test]
    fn real_rinex3_first_gps_record_decodes_to_its_broadcast_fields() {
        let f = std::fs::File::open(
            "../../datasets/wtzr_ppp_1224/BRDC00IGS_R_20203590000_01D_MN.rnx",
        )
        .unwrap();
        let (eph, _klob) = parse_rinex_nav(BufReader::new(f)).unwrap();
        let g = match &eph[0] {
            Ephemeris::Gps(g) => g,
            _ => panic!("expected GPS"),
        };
        assert_eq!(g.sat.prn, 1);
        assert_eq!((g.toc.week, g.toc.tow), (2137, 345600.0));
        assert_eq!((g.toe.week, g.toe.tow), (2137, 345600.0));
        assert!((g.af0 - 7.914672605693E-04).abs() < 1e-18);
        assert!((g.af1 - -5.570655048359E-12).abs() < 1e-24);
        assert_eq!(g.iode, 51);
        assert_eq!(g.iodc, 51);
        assert!((g.sqrt_a - 5153.695047379).abs() < 1e-9);
        assert!((g.e - 1.020454068203E-02).abs() < 1e-15);
        assert!((g.m0 - -1.307172517279E+00).abs() < 1e-13);
        assert!((g.omega_dot - -7.690320332710E-09).abs() < 1e-20);
        assert!((g.idot - 1.207193141583E-10).abs() < 1e-22);
        assert!((g.tgd - 5.122274160385E-09).abs() < 1e-20);
    }
}
