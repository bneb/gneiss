use core::fmt;
use core::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ObsType {
    Pseudorange, // 'C'
    CarrierPhase, // 'L'
    Doppler, // 'D'
    Snr, // 'S'
}

impl fmt::Display for ObsType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let c = match self {
            ObsType::Pseudorange => 'C',
            ObsType::CarrierPhase => 'L',
            ObsType::Doppler => 'D',
            ObsType::Snr => 'S',
        };
        write!(f, "{}", c)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SignalCode {
    pub freq_band: u8,
    pub attribute: char,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ObsCode {
    pub obs_type: ObsType,
    pub signal: SignalCode,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseObsCodeError;

impl fmt::Display for ParseObsCodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid observation code format")
    }
}

impl FromStr for ObsCode {
    type Err = ParseObsCodeError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() != 3 {
            return Err(ParseObsCodeError);
        }

        let mut chars = s.chars();
        let obs_char = chars.next().unwrap();
        let freq_char = chars.next().unwrap();
        let attr_char = chars.next().unwrap();

        let obs_type = match obs_char {
            'C' => ObsType::Pseudorange,
            'L' => ObsType::CarrierPhase,
            'D' => ObsType::Doppler,
            'S' => ObsType::Snr,
            _ => return Err(ParseObsCodeError),
        };

        let freq_band = freq_char.to_digit(10).ok_or(ParseObsCodeError)? as u8;

        Ok(ObsCode {
            obs_type,
            signal: SignalCode {
                freq_band,
                attribute: attr_char,
            },
        })
    }
}

impl fmt::Display for ObsCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}", self.obs_type, self.signal.freq_band, self.signal.attribute)
    }
}

use crate::sat::SatelliteId;
use crate::time::GpsTime;
use alloc::vec::Vec;

/// Represents a single observation (e.g. L1C carrier phase) for a satellite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observation {
    pub code: ObsCode,
    pub value: f64,
    /// Optional lock time for carrier phase measurements (in increments)
    pub lock_time: Option<u16>,
}

/// All observations for a specific satellite at a specific epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct SatObs {
    pub sat: SatelliteId,
    pub observations: Vec<Observation>,
}

impl SatObs {
    pub fn get_observable(&self, freq_band: u8) -> Option<f64> {
        self.observations.iter().find(|o| o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == freq_band).map(|o| o.value)
    }

    pub fn get_observable_phase(&self, freq_band: u8) -> Option<f64> {
        self.observations.iter().find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == freq_band).map(|o| o.value)
    }

    pub fn get_locktime(&self, freq_band: u8) -> Option<u16> {
        self.observations.iter().find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == freq_band).and_then(|o| o.lock_time)
    }

    pub fn get_snr(&self, freq_band: u8) -> Option<u8> {
        self.observations.iter().find(|o| o.code.obs_type == ObsType::Snr && o.code.signal.freq_band == freq_band).map(|o| o.value as u8)
    }
}

/// A complete epoch of observations across all tracked satellites.
#[derive(Debug, Clone, PartialEq)]
pub struct EpochObs {
    pub time: GpsTime,
    pub satellites: Vec<SatObs>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::string::ToString;

    #[test]
    fn test_obs_code_parsing() {
        let code1 = ObsCode::from_str("L1C").unwrap();
        assert_eq!(code1.obs_type, ObsType::CarrierPhase);
        assert_eq!(code1.signal.freq_band, 1);
        assert_eq!(code1.signal.attribute, 'C');

        let code2 = ObsCode::from_str("C2W").unwrap();
        assert_eq!(code2.obs_type, ObsType::Pseudorange);
        assert_eq!(code2.signal.freq_band, 2);
        assert_eq!(code2.signal.attribute, 'W');

        let code3 = ObsCode::from_str("D5Q").unwrap();
        assert_eq!(code3.obs_type, ObsType::Doppler);
        assert_eq!(code3.signal.freq_band, 5);
        assert_eq!(code3.signal.attribute, 'Q');

        let code4 = ObsCode::from_str("S8X").unwrap();
        assert_eq!(code4.obs_type, ObsType::Snr);
        assert_eq!(code4.signal.freq_band, 8);
        assert_eq!(code4.signal.attribute, 'X');

        assert!(ObsCode::from_str("X1C").is_err()); // Invalid ObsType
        assert!(ObsCode::from_str("LC").is_err()); // Too short
        assert!(ObsCode::from_str("LXC").is_err()); // Invalid freq band
    }

    #[test]
    fn test_obs_code_display() {
        let code1 = ObsCode {
            obs_type: ObsType::Pseudorange,
            signal: SignalCode { freq_band: 1, attribute: 'C' }
        };
        assert_eq!(code1.to_string(), "C1C");

        let code2 = ObsCode {
            obs_type: ObsType::CarrierPhase,
            signal: SignalCode { freq_band: 2, attribute: 'W' }
        };
        assert_eq!(code2.to_string(), "L2W");
    }
}
