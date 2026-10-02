//! RINEX observation file parser (RINEX 2.xx and 3.xx).

pub mod header;
pub mod rinex2;
pub mod rinex3;
#[cfg(test)]
mod tests;

pub use rinex2::parse_rinex_2_obs;
pub use rinex3::parse_rinex_3_obs;
#[cfg(test)]
pub use rinex3::parse_rinex_3_obs_line;
#[cfg(test)]
pub(crate) use header::parse_rinex_3_obs_types_list;

use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SignalCode};
use std::io::BufRead;

/// Header metadata extracted from RINEX observation files.
#[derive(Clone, Debug, Default)]
pub struct RinexObsHeader {
    pub approx_position: Option<[f64; 3]>,
    pub antenna_delta: Option<[f64; 3]>,
    pub marker_name: Option<String>,
}

/// Parse a RINEX 14-char float field. Returns None on parse failure or blank.
pub(crate) fn parse_rinex_f14(s: &str) -> Option<f64> {
    let s = s.trim();
    if s.is_empty() {
        return None;
    }
    s.replace(['D', 'd'], "e").parse::<f64>().ok()
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
pub fn parse_rinex_obs_epochs<R: BufRead>(reader: R) -> Result<Vec<EpochObs>, String> {
    parse_rinex_obs(reader).map(|(epochs, _header)| epochs)
}

pub(crate) fn parse_rinex_obs_lli(obs_line: &str, start: usize) -> Option<u8> {
    if start + 14 >= obs_line.len() {
        return None;
    }
    let lli_char = obs_line[start + 14..start + 15]
        .chars()
        .next()
        .unwrap_or(' ');
    if lli_char == ' ' {
        return None;
    }
    lli_char.to_string().parse::<u8>().ok()
}

pub(crate) fn map_rinex_type(type_str: &str, val: f64, lli: Option<u8>) -> Option<Observation> {
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
#[allow(clippy::unwrap_used)]
mod helper_tests {
    use super::*;

    /// LLI lives in byte 14 of each 16-byte observation field (F14.3 + LLI +
    /// strength); byte 15 is the signal-strength flag. Golden vector from
    /// rover_p123.obs line 43, field 1: bytes 16..30 = "  97861120.237",
    /// byte 30 = ' ', byte 31 = '6'.
    #[test]
    fn lli_is_the_fifteenth_byte_of_a_sixteen_byte_field() {
        let line = " 125588294.513 6  97861120.237 6  23898625.087    23898629.285";
        // Field 0 starts at byte 0 (LLI at 14, strength at 15); field 1 starts
        // at byte 16 (LLI at 30, strength at 31).
        assert_eq!(parse_rinex_obs_lli(line, 0), None, "byte 14 is blank");
        assert_eq!(parse_rinex_obs_lli(line, 16), None, "byte 30 is blank");
        let flagged = " 125588294.5133   97861120.2377   23898625.0871 ";
        assert_eq!(flagged.len(), 48);
        assert_eq!(flagged.as_bytes()[14], b'3', "LLI of field 0");
        assert_eq!(flagged.as_bytes()[15], b' ', "strength of field 0");
        assert_eq!(flagged.as_bytes()[30], b'7', "LLI of field 1");
        assert_eq!(parse_rinex_obs_lli(flagged, 0), Some(3));
        assert_eq!(parse_rinex_obs_lli(flagged, 16), Some(7));
        assert_eq!(parse_rinex_obs_lli(flagged, 32), Some(1));
    }

    /// Byte 14 must be a character boundary. Every observation record is
    /// sliced at a fixed byte offset, so a record carrying one multi-byte
    /// character makes `obs_line[start + 14..start + 15]` (and the value slice
    /// `obs_line[start..start + 14]`) panic on a non-boundary index. RINEX
    /// records are nominally ASCII, but the parser accepts any UTF-8 line from
    /// a `BufRead` and must not abort the process on it.
    #[test]
    fn observation_field_slicing_survives_a_multibyte_character() {
        // 13 ASCII bytes then U+00E9, which occupies bytes 13 and 14: byte 14
        // lands in the middle of the character.
        let line = "abcdefghijklm\u{e9} ";
        assert_eq!(line.len(), 16);
        let _ = parse_rinex_obs_lli(line, 0);
    }

    /// `parse_rinex_f14` is the shared 14-column float reader for APPROX
    /// POSITION XYZ, ANTENNA: DELTA H/E/N and MARKER NUMBER.
    #[test]
    fn f14_reader_rejects_blank_and_accepts_fortran_exponents() {
        assert_eq!(parse_rinex_f14("   -1405299.676"), Some(-1405299.676));
        assert_eq!(parse_rinex_f14("  931853.5784"), Some(931853.5784));
        assert_eq!(parse_rinex_f14("  1.2345D+03 "), Some(1234.5));
        assert_eq!(parse_rinex_f14(""), None);
        assert_eq!(parse_rinex_f14("          "), None);
    }

    /// `map_rinex_type` splits a RINEX observation code into type, band and
    /// attribute. RINEX 3 codes are three characters, RINEX 2 codes two, and
    /// the band digit is mandatory in both.
    #[test]
    fn observation_code_splits_into_type_band_and_attribute() {
        let c = map_rinex_type("C1C", 1.0, None).unwrap();
        assert_eq!(c.code.obs_type, ObsType::Pseudorange);
        assert_eq!(c.code.signal.freq_band, 1);
        assert_eq!(c.code.signal.attribute, 'C');
        let c = map_rinex_type("L2W", 2.0, Some(5)).unwrap();
        assert_eq!(c.code.obs_type, ObsType::CarrierPhase);
        assert_eq!(c.code.signal.freq_band, 2);
        assert_eq!(c.code.signal.attribute, 'W');
        assert_eq!(c.lli, Some(5));
        // RINEX 2 two-character code: missing attribute reads as blank.
        let c = map_rinex_type("D5", 0.0, None).unwrap();
        assert_eq!(c.code.obs_type, ObsType::Doppler);
        assert_eq!(c.code.signal.freq_band, 5);
        assert_eq!(c.code.signal.attribute, ' ');
    }

    /// Codes with an unknown type letter or a missing band digit carry no
    /// decodable observation and must be dropped, not guessed at.
    #[test]
    fn unusable_observation_codes_are_dropped() {
        assert!(map_rinex_type("X1", 1.0, None).is_none());
        assert!(map_rinex_type("CZ", 1.0, None).is_none());
        assert!(map_rinex_type("", 1.0, None).is_none());
    }
}
