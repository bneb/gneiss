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
    // The declared count, not the number of codes that happened to be printed,
    // fixes how many observation records each satellite occupies.
    let n_types = header.num_obs_declared.unwrap_or(obs_types.len()).max(obs_types.len());

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
                    sat_str = next_line.get(32..).unwrap_or("").to_string();
                } else {
                    break;
                }
            }
        }

        let time = GpsTime::from_calendar(year, month, day, hour, min, sec);
        let mut satellites = Vec::new();

        for sat_id_str in sat_list {
            if let Some(sat_obs) = parse_rinex_2_obs_sat(&sat_id_str, &obs_types, n_types, lines) {
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
    n_types: usize,
    lines: &mut I,
) -> Option<SatObs> {
    // Columns 33.. are padded with blank 3-character fields; they are padding,
    // not satellites. A blank field or an unparseable PRN must never become a
    // satellite, because each phantom entry then consumes one observation
    // record and desynchronises the rest of the file.
    if sat_id_str.trim().is_empty() {
        return None;
    }
    let constellation_char = sat_id_str.chars().next().unwrap_or('G');
    let prn = sat_id_str.get(1..3)?.trim().parse::<u8>().ok()?;
    if prn == 0 {
        return None;
    }
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
    let num_val_lines = (n_types as f64 / 5.0).ceil() as usize;
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

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::io::BufReader;

    /// Epoch record declaring 2 satellites. Only "G01" fits in columns 33-35,
    /// so the second PRN has to come from a continuation record; that
    /// continuation record is 32 bytes of blanks, i.e. it carries no PRN at
    /// all. Columns 33-68 of every satellite-list record - the epoch record
    /// and each continuation - hold the PRNs (RTKLIB rinex.cc:663-668 always
    /// restarts the column cursor at j=32, whatever the line length is).
    const BLANK_CONTINUATION: &str = concat!(
        "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE\n",
        "     2    C1    L1                                              # / TYPES OF OBSERV\n",
        "                                                            END OF HEADER\n",
        " 20  5 14 22  0  0.0000000  0  2G01\n",
        "                                \n",
        "  25140323.324   25140323.324\n",
        "                                                                \n",
        "                                                                \n",
        "                                                                \n",
    );

    /// A blank satellite-list continuation record must not invent satellites.
    /// `next_line.len() > 32` is false for a 32-byte record, so the whole
    /// record is scanned for PRNs and yields ten 3-byte "   " fields, each of
    /// which `parse_rinex_2_obs_sat` turns into GPS PRN 0 (`unwrap_or(0)` on a
    /// blank PRN). PRN 0 is not a satellite: GPS PRNs run 1..=32, so a
    /// fabricated PRN-0 entry is pure corruption that then desynchronises the
    /// rest of the file (one phantom satellite consumes one obs record).
    #[test]
    fn blank_satellite_continuation_record_invents_no_satellites() {
        let mut reader = BufReader::new(BLANK_CONTINUATION.as_bytes());
        let (epochs, _h) = super::super::parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1, "phantom PRN-0 satellites");
        assert_eq!(epochs[0].satellites[0].sat.prn, 1);
    }

    /// Every satellite identifier recovered from a real RINEX 2 file must be a
    /// real PRN: GPS 1..=32, GLONASS 1..=24, Galileo 1..=36, BeiDou 1..=63,
    /// QZSS 1..=10, SBAS 120..=141. Zero is the parser's "unparseable" default
    /// and must never reach a caller.
    #[test]
    fn real_rinex2_file_yields_no_zero_prn() {
        let f = std::fs::File::open("../../datasets/cors_baseline/rover_p123.obs").unwrap();
        let (epochs, _h) = super::super::parse_rinex_obs(BufReader::new(f)).unwrap();
        let bad = epochs
            .iter()
            .flat_map(|e| e.satellites.iter())
            .filter(|s| s.sat.prn == 0)
            .count();
        assert_eq!(bad, 0, "satellites decoded with PRN 0");
    }

    /// Differential check against datasets/cors_baseline/rover_p123.obs
    /// (RINEX 2.11, 20 obs types => 4 records per satellite, 30 s interval,
    /// 2023-04-10T00:00 to 23:59:30 UTC). Counts below come from an
    /// independent python walk of the file that only ever reads the columns
    /// named in the RINEX 2 obs record layout (epoch: year cols 2-3, month
    /// 5-6, day 8-9, hour 11-12, min 14-15, sec 17-26, flag col 29, SV count
    /// cols 30-32, PRNs from col 33 in 3-char fields; obs: 5 fields of
    /// F14.3 + LLI + strength per 80-char record).
    #[test]
    fn real_rinex2_file_matches_independent_column_walk() {
        let f = std::fs::File::open("../../datasets/cors_baseline/rover_p123.obs").unwrap();
        let (epochs, _h) = super::super::parse_rinex_obs(BufReader::new(f)).unwrap();
        let nsat: usize = epochs.iter().map(|e| e.satellites.len()).sum();
        let nobs: usize = epochs
            .iter()
            .map(|e| e.satellites.iter().map(|s| s.observations.len()).sum::<usize>())
            .sum();
        // 24 h * 3600 / 30 s = 2880 epochs.
        assert_eq!(epochs.len(), 2880);
        assert_eq!(nsat, 75464);
        assert_eq!(nobs, 825536);
    }

    /// Golden vector, first satellite record of rover_p123.obs (file line 43).
    /// Column layout of the record: each field is 16 bytes wide, made of a
    /// 14-byte F14.3 value, a one-byte LLI and a one-byte signal strength, and
    /// five fields fit in one 80-byte record.
    ///
    /// Bytes 0-13 hold " 125588294.513" (L1), byte 14 is the LLI (blank, no
    /// cycle slip) and byte 15 the strength ('6'); bytes 16-29 hold
    /// "  97861120.237" (L2).
    #[test]
    fn obs_value_and_lli_columns_match_the_file() {
        let line = " 125588294.513 6  97861120.237 6  23898625.087    23898629.285    23898624.909";
        assert_eq!(line.len(), 78);
        let obs = parse_rinex_2_obs_value(line, 0, "L1").unwrap();
        assert_eq!(obs.value, 125588294.513);
        assert_eq!(obs.lli, None, "LLI byte 14 is a blank");
        let obs = parse_rinex_2_obs_value(line, 1, "L2").unwrap();
        assert_eq!(obs.value, 97861120.237);
        assert_eq!(obs.lli, None);
        let obs = parse_rinex_2_obs_value(line, 2, "C1").unwrap();
        assert_eq!(obs.value, 23898625.087);
    }

    /// Same record, second field of rover_p123.obs line 47: bytes 30 and 31
    /// are LLI='4' (cycle slip) and strength='3'.
    #[test]
    fn lli_flag_is_read_from_byte_fourteen_of_each_field() {
        let line = " 127412365.643 6  99282462.59743  24245741.548    24245748.904    24245741.343";
        assert_eq!(line.len(), 78);
        let obs = parse_rinex_2_obs_value(line, 1, "L2").unwrap();
        assert_eq!(obs.value, 99282462.597);
        assert_eq!(obs.lli, Some(4));
    }

    /// RINEX 2 states the number of observation types in columns 1-6 of
    /// `# / TYPES OF OBSERV`; that number is what fixes how many 80-byte
    /// records each satellite occupies: ceil(N/5), five F14.3+LLI+strength
    /// fields per record. This header declares six types but only prints five
    /// codes, so every satellite still spans TWO records. `parse_rinex_2_header`
    /// reads the count into `num_obs` and then never uses it, so the parser
    /// takes ceil(5/5) = 1 record per satellite, swallows the satellite's
    /// second record as the next satellite's data, and then mistakes the
    /// leftover records for epoch records.
    #[test]
    fn declared_observation_count_fixes_records_per_satellite() {
        let data = concat!(
            "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE\n",
            "     6    C1    L1    C2    L2    S1                              # / TYPES OF OBSERV\n",
            "                                                            END OF HEADER\n",
            " 20  5 14 22  0  0.0000000  0  2G01G02\n",
            "       100.000         200.000         110.000         220.000          50.000  \n",
            "            60.000                                                        \n",
            "       101.000         201.000         111.000         221.000          51.000  \n",
            "            61.000                                                        \n",
        );
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _h) = super::super::parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1, "leftover records parsed as epochs");
        assert_eq!(epochs[0].satellites.len(), 2);
        let g2 = &epochs[0].satellites[1];
        assert_eq!(g2.sat.prn, 2);
        assert_eq!(g2.observations.len(), 5);
        assert_eq!(g2.observations[0].value, 101.0);
        assert_eq!(g2.observations[4].value, 51.0);
    }
}
