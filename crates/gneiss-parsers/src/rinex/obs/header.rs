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
            // Columns 1-60 hold the name; a truncated record is still legal
            // input, so take whatever of the 60 columns is present.
            header.marker_name = Some(marker_name_field(&current_line));
        }
        if current_line.contains("# / TYPES OF OBSERV") {
            if num_obs == 0 {
                num_obs = current_line.get(0..6).unwrap_or("").trim().parse::<usize>().unwrap_or(0);
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
    if num_obs > 0 {
        header.num_obs_declared = Some(num_obs);
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
                'I' => Constellation::Navic,
                _ => {
                    if let Some(next_line) = lines.next() {
                        current_line = next_line;
                        continue;
                    } else {
                        break;
                    }
                }
            };

            let count = parse_obs_type_count(&current_line)?;
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
            // Columns 1-60 hold the name; a truncated record is still legal
            // input, so take whatever of the 60 columns is present.
            header.marker_name = Some(marker_name_field(&current_line));
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

/// Columns 1-60 of a `MARKER NAME` record hold the name. A record that does
/// not physically carry all 60 columns is still legal input, so take whatever
/// is present up to the label rather than indexing past the end.
fn marker_name_field(line: &str) -> String {
    match line.get(0..60) {
        Some(name) => name.trim().to_string(),
        None => line.split("MARKER NAME").next().unwrap_or("").trim().to_string(),
    }
}

/// Reads the mandatory "number of observation types" from columns 4-6 of a
/// `SYS / # / OBS TYPES` record. A record that does not carry a number cannot
/// be delimited, so accepting it would insert the constellation with an empty
/// type list and silently produce satellites with no observations.
fn parse_obs_type_count(line: &str) -> Result<usize, String> {
    line.get(3..6)
        .ok_or_else(|| "SYS / # / OBS TYPES record is shorter than 6 columns".to_string())?
        .trim()
        .parse::<usize>()
        .map_err(|_| format!("SYS / # / OBS TYPES record has no type count: {line:?}"))
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

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// RINEX 2 header whose MARKER NAME record is truncated to 33 bytes.
    /// Columns 1-60 carry the marker name, so a 33-byte record still names
    /// the site: everything before the trailing spaces is "P12".
    const V2_SHORT_MARKER: &str = concat!(
        "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE\n",
        "P12                              MARKER NAME\n",
        "     2    C1    L1                                              # / TYPES OF OBSERV\n",
        "                                                            END OF HEADER\n",
    );

    /// MARKER NAME shorter than 60 bytes must not panic: header.rs slices
    /// `current_line[0..60]`, which is an out-of-bounds byte-index panic for
    /// any record that does not physically carry all 60 name columns.
    #[test]
    fn rinex2_short_marker_name_record_does_not_panic() {
        let mut lines = V2_SHORT_MARKER.lines().map(str::to_string);
        let first = lines.next().unwrap();
        let (types, header) = parse_rinex_2_header(first, &mut lines).unwrap();
        assert_eq!(types, vec!["C1".to_string(), "L1".to_string()]);
        assert_eq!(header.marker_name.as_deref(), Some("P12"));
    }

    /// A SYS / # / OBS TYPES record whose mandatory columns 4-6 do not hold a
    /// number must not be accepted. `count` falls back to 0
    /// (`current_line[3..6]...unwrap_or(0)`), so the constellation is still
    /// inserted with an EMPTY observation-type list, `const_obs_types` is
    /// non-empty, and the caller never learns the header was unreadable -
    /// every satellite of that constellation then comes back with zero
    /// observations. There is no valid header behind this line, so the only
    /// acceptable outcome is an error.
    #[test]
    fn rinex3_obs_types_record_without_count_is_an_error() {
        let data = concat!(
            "     3.04           OBSERVATION DATA    M: MIXED            RINEX VERSION / TYPE\n",
            "SYS / # / OBS TYPES\n",
            "                                                            END OF HEADER\n",
        );
        let mut lines = data.lines().map(str::to_string);
        let first = lines.next().unwrap();
        assert!(parse_rinex_3_header(first, &mut lines).is_err());
    }

    /// Same defect on the RINEX 3 marker record: header.rs:122 also slices
    /// `current_line[0..60]`.
    #[test]
    fn rinex3_short_marker_name_record_does_not_panic() {
        let data = concat!(
            "     3.04           OBSERVATION DATA    M: MIXED            RINEX VERSION / TYPE\n",
            "WTZR                              MARKER NAME\n",
            "G    1 C1C                                                    SYS / # / OBS TYPES\n",
            "                                                            END OF HEADER\n",
        );
        let mut lines = data.lines().map(str::to_string);
        let first = lines.next().unwrap();
        let (types, header) = parse_rinex_3_header(first, &mut lines).unwrap();
        assert_eq!(types.get(&Constellation::Gps).map(Vec::len), Some(1));
        assert_eq!(header.marker_name.as_deref(), Some("WTZR"));
    }

    /// Golden vector for the 60-column marker-name field, taken verbatim from
    /// datasets/cors_baseline/rover_p123.obs line 21 (71 bytes: 60 name
    /// columns + "MARKER NAME" label). Columns 1-60 are "P123" + padding.
    #[test]
    fn rinex2_marker_name_takes_columns_1_to_60() {
        let line = "P123                                                        MARKER NAME";
        assert_eq!(line.len(), 71);
        let rest = vec![
            "     2    C1    L1                                              # / TYPES OF OBSERV"
                .to_string(),
            "                                                            END OF HEADER"
                .to_string(),
        ];
        let (_types, header) =
            parse_rinex_2_header(line.to_string(), &mut rest.into_iter()).unwrap();
        assert_eq!(header.marker_name.as_deref(), Some("P123"));
    }
}
