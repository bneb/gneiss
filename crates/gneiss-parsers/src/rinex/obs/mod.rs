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
