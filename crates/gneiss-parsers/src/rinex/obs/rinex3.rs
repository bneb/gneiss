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

        // Epoch flag (column 32). Flags 3, 4 and 5 mark an epoch whose
        // following records are header records, not satellites.
        let flag: usize = line[30..32].trim().parse().unwrap_or(0);
        let is_header_epoch = (3..=5).contains(&flag);

        let num_sats = line[32..35].trim().parse::<usize>().unwrap_or(0);
        let time = GpsTime::from_calendar(year, month, day, hour, min, sec);
        let mut satellites = Vec::with_capacity(num_sats);

        let n_records = if is_header_epoch { 0 } else { num_sats };
        for _ in 0..n_records {
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
        'I' => Constellation::Navic,
        _ => return None,
    };
    let prn = sat_str.get(1..3)?.trim().parse::<u8>().ok()?;
    // No constellation has PRN 0; an unparseable id is not a satellite.
    if prn == 0 {
        return None;
    }
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

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::rinex::obs::parse_rinex_obs;
    use std::io::BufReader;

    /// Minimal RINEX 3.03 file with one Navic (system 'I') observation record.
    /// Column layout: epoch record `> yyyy mm dd hh mm ss.sssssss f nnn` with
    /// the SV count at bytes 33-35 (0-based 32..35), one satellite record per
    /// SV whose 3-byte system+PRN prefix is followed by 16-byte fields
    /// (F14.3 + LLI + strength).
    const NAVIC_FILE: &str = concat!(
        "     3.03           OBSERVATION DATA    M: MIXED            RINEX VERSION / TYPE\n",
        "I    2 C5A L5A                                                        SYS / # / OBS TYPES\n",
        "                                                            END OF HEADER\n",
        "> 2020 12 24 00 00  0.0000000  0  1\n",
        "I01  23456789.123     1234567.890\n",
    );

    /// Navic/IRNSS ('I') is a first-class RINEX 3 system - the NAV parser in
    /// this crate maps it to `Constellation::Navic` (nav/mod.rs) - yet the
    /// observation side never learns about it: `header.rs` rejects the 'I'
    /// SYS / # / OBS TYPES record and `parse_rinex_3_obs_line` returns None for
    /// an 'I' satellite prefix, so the record is dropped without a word.
    #[test]
    fn navic_observation_records_are_decoded() {
        let mut reader = BufReader::new(NAVIC_FILE.as_bytes());
        let (epochs, _h) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1, "Navic record dropped");
        assert_eq!(
            epochs[0].satellites[0].sat.constellation,
            Constellation::Navic
        );
        assert_eq!(epochs[0].satellites[0].sat.prn, 1);
        assert_eq!(epochs[0].satellites[0].observations.len(), 2);
        assert_eq!(epochs[0].satellites[0].observations[0].value, 23456789.123);
    }

    /// Unit-level twin of the test above, so the failure points at the exact
    /// `match constellation_char` arm that is missing 'I'.
    #[test]
    fn navic_satellite_prefix_is_recognised() {
        let mut map = HashMap::new();
        map.insert(Constellation::Navic, vec!["C5A".to_string()]);
        let line = "I01  23456789.123     1234567.890";
        let sat = parse_rinex_3_obs_line(line, &map).expect("I-prefixed record rejected");
        assert_eq!(sat.sat.constellation, Constellation::Navic);
        assert_eq!(sat.sat.prn, 1);
    }

    /// A RINEX 3 epoch record stores its flag in byte 31 and its SV count in
    /// bytes 32..35 (RTKLIB rinex.cc:653-654). Flag 3, 4 and 5 mean "new site
    /// occupation / header information / external event follows", so the
    /// records that follow are header records, not satellite records - RTKLIB
    /// returns early for flags 3..5 (rinex.cc:659-661). The parser never reads
    /// byte 31, so it consumes the header record as a satellite: "G   16 ..."
    /// becomes GPS with a blank PRN, i.e. a fabricated GPS PRN 0.
    #[test]
    fn event_flagged_epoch_does_not_turn_header_records_into_satellites() {
        let data = concat!(
            "     3.04           OBSERVATION DATA    M: MIXED            RINEX VERSION / TYPE\n",
            "G    1 C1C                                                            SYS / # / OBS TYPES\n",
            "                                                            END OF HEADER\n",
            "> 2020 12 24 00 00  0.0000000  4  1\n",
            "G   16 C1C L1C D1C S1C                                      SYS / # / OBS TYPES\n",
            "> 2020 12 24 00 00 30.0000000  0  1\n",
            "G01  23456789.123\n",
        );
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _h) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 2);
        let zero_prn: usize = epochs
            .iter()
            .flat_map(|e| e.satellites.iter())
            .filter(|s| s.sat.prn == 0)
            .count();
        assert_eq!(zero_prn, 0, "header record decoded as a satellite");
        // The ordinary epoch that follows must still be read correctly.
        assert_eq!(epochs[1].satellites.len(), 1);
        assert_eq!(epochs[1].satellites[0].sat.prn, 1);
        assert_eq!(epochs[1].satellites[0].observations[0].value, 23456789.123);
    }

    /// Differential check against
    /// datasets/wtzr_ppp_1224/WTZR00DEU_R_20203590000_01D_30S_MO.rnx
    /// (RINEX 3.03, 30 s interval, 2020-12-24T00:00 to 23:59:30). Derived
    /// from an independent walk of the file that reads only the documented
    /// columns. Satellite-record total is the sum of the epoch records'
    /// SV counts (bytes 33-35): GPS 29668, BeiDou 28636, GLONASS 24799,
    /// Galileo 24304, Navic 8994, SBAS 6059 => 122460.
    #[test]
    fn real_rinex3_file_yields_every_declared_satellite_record() {
        let f = std::fs::File::open(
            "../../datasets/wtzr_ppp_1224/WTZR00DEU_R_20203590000_01D_30S_MO.rnx",
        )
        .unwrap();
        let (epochs, header) = parse_rinex_obs(BufReader::new(f)).unwrap();
        assert_eq!(epochs.len(), 2880);
        let nsat: usize = epochs.iter().map(|e| e.satellites.len()).sum();
        assert_eq!(nsat, 122460, "satellite records lost (Navic?)");
        assert_eq!(header.marker_name.as_deref(), Some("WTZR"));
    }

    /// Epoch times of a well-ordered 30 s RINEX 3 file must be strictly
    /// increasing and start at GPS week 2137, tow 345600
    /// (1980-01-06 + 2137*7 + 345600/86400 days = 2020-12-24T00:00:00).
    #[test]
    fn real_rinex3_epoch_times_are_ordered_and_correctly_anchored() {
        let f = std::fs::File::open(
            "../../datasets/wtzr_ppp_1224/WTZR00DEU_R_20203590000_01D_30S_MO.rnx",
        )
        .unwrap();
        let (epochs, _h) = parse_rinex_obs(BufReader::new(f)).unwrap();
        assert_eq!((epochs[0].time.week, epochs[0].time.tow), (2137, 345600.0));
        for w in epochs.windows(2) {
            assert!(w[1].time - w[0].time > 0.0, "epoch times not increasing");
        }
        assert!((epochs[2879].time - epochs[0].time - 2879.0 * 30.0).abs() < 1e-3);
    }

    /// Golden vector: first satellite record of the same file (line 39),
    /// "G08  23457408.784   123269529.49207 ...". Field 0 (C1C) occupies
    /// bytes 3..17 with a blank LLI at byte 17 and a blank strength at 18;
    /// field 1 (L1C) occupies bytes 19..33 with LLI '0' at byte 33 and
    /// strength '7' at byte 34.
    #[test]
    fn rinex3_value_and_lli_columns_match_the_file() {
        let line = "G08  23457408.784   123269529.49207      3104.737          46.300";
        assert_eq!(line.as_bytes()[17], b' ');
        assert_eq!(line.as_bytes()[33], b'0');
        let obs = parse_rinex_3_obs_value(line, 0, "C1C").unwrap();
        assert_eq!(obs.value, 23457408.784);
        assert_eq!(obs.lli, None);
        assert_eq!(obs.code.signal.freq_band, 1);
        assert_eq!(obs.code.signal.attribute, 'C');
        let obs = parse_rinex_3_obs_value(line, 1, "L1C").unwrap();
        assert_eq!(obs.value, 123269529.492);
        assert_eq!(obs.lli, Some(0));
        let obs = parse_rinex_3_obs_value(line, 2, "D1C").unwrap();
        assert_eq!(obs.value, 3104.737);
        assert_eq!(obs.lli, None);
    }
}
