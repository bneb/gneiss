use core::fmt;
use core::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ObsType {
    Pseudorange,  // 'C'
    CarrierPhase, // 'L'
    Doppler,      // 'D'
    Snr,          // 'S'
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
        let obs_char = chars.next().expect("s.len() == 3, so at least 3 chars");
        let freq_char = chars.next().expect("s.len() == 3, so at least 3 chars");
        let attr_char = chars.next().expect("s.len() == 3, so at least 3 chars");

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
        write!(
            f,
            "{}{}{}",
            self.obs_type, self.signal.freq_band, self.signal.attribute
        )
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
    /// Optional Loss of Lock Indicator (from RINEX)
    pub lli: Option<u8>,
}

/// All observations for a specific satellite at a specific epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct SatObs {
    pub sat: SatelliteId,
    pub observations: Vec<Observation>,
}

impl SatObs {
    fn matches_band(&self, o_band: u8, o_attr: char, req_band: u8) -> bool {
        if o_band == req_band {
            return true;
        }
        if self.sat.constellation == crate::sat::Constellation::Beidou && o_attr == 'I' {
            return (req_band == 1 && o_band == 2) || (req_band == 2 && o_band == 1);
        }
        false
    }

    pub fn get_observable(&self, freq_band: u8) -> Option<f64> {
        self.observations
            .iter()
            .find(|o| {
                o.code.obs_type == ObsType::Pseudorange
                    && self.matches_band(o.code.signal.freq_band, o.code.signal.attribute, freq_band)
            })
            .map(|o| o.value)
    }

    pub fn get_observable_phase(&self, freq_band: u8) -> Option<f64> {
        self.observations
            .iter()
            .find(|o| {
                o.code.obs_type == ObsType::CarrierPhase
                    && self.matches_band(o.code.signal.freq_band, o.code.signal.attribute, freq_band)
            })
            .map(|o| o.value)
    }

    pub fn get_observable_phase_lli(&self, freq_band: u8) -> Option<(f64, Option<u8>)> {
        self.observations
            .iter()
            .find(|o| {
                o.code.obs_type == ObsType::CarrierPhase
                    && self.matches_band(o.code.signal.freq_band, o.code.signal.attribute, freq_band)
            })
            .map(|o| (o.value, o.lli))
    }

    pub fn get_doppler(&self, freq_band: u8) -> Option<f64> {
        self.observations
            .iter()
            .find(|o| {
                o.code.obs_type == ObsType::Doppler
                    && self.matches_band(o.code.signal.freq_band, o.code.signal.attribute, freq_band)
            })
            .map(|o| o.value)
    }

    pub fn get_locktime(&self, freq_band: u8) -> Option<u16> {
        self.observations
            .iter()
            .find(|o| {
                o.code.obs_type == ObsType::CarrierPhase
                    && self.matches_band(o.code.signal.freq_band, o.code.signal.attribute, freq_band)
            })
            .and_then(|o| o.lock_time)
    }

    pub fn get_lli(&self, freq_band: u8) -> Option<u8> {
        self.observations
            .iter()
            .find(|o| {
                o.code.obs_type == ObsType::CarrierPhase
                    && self.matches_band(o.code.signal.freq_band, o.code.signal.attribute, freq_band)
            })
            .and_then(|o| o.lli)
    }

    pub fn get_snr(&self, freq_band: u8) -> Option<u8> {
        self.observations
            .iter()
            .find(|o| {
                o.code.obs_type == ObsType::Snr
                    && self.matches_band(o.code.signal.freq_band, o.code.signal.attribute, freq_band)
            })
            .map(|o| o.value as u8)
    }
}

/// A complete epoch of observations across all tracked satellites.
#[derive(Debug, Clone, PartialEq)]
pub struct EpochObs {
    pub time: GpsTime,
    pub satellites: Vec<SatObs>,
}

/// Retain only satellites whose constellation appears in `allowed`
/// (single-letter codes: `G`=GPS, `R`=GLONASS, `E`=Galileo, `C`=BeiDou).
/// Shared by constellation-mix experiments and CLI `--systems` filtering
/// so both apply the exact same set semantics.
pub fn filter_constellations(epochs: &mut [EpochObs], allowed: &str) {
    let keep = |c: crate::sat::Constellation| match c {
        crate::sat::Constellation::Gps => allowed.contains('G'),
        crate::sat::Constellation::Glonass => allowed.contains('R'),
        crate::sat::Constellation::Galileo => allowed.contains('E'),
        crate::sat::Constellation::Beidou => allowed.contains('C'),
        _ => false,
    };
    for e in epochs.iter_mut() {
        e.satellites.retain(|s| keep(s.sat.constellation));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sat::{Constellation, SatelliteId};
    use crate::time::GpsTime;
    use alloc::string::ToString;
    use alloc::vec;

    fn epoch_with(constellations: &[Constellation]) -> EpochObs {
        EpochObs {
            time: GpsTime::new(2000, 0.0),
            satellites: constellations
                .iter()
                .enumerate()
                .map(|(i, &constellation)| SatObs {
                    sat: SatelliteId { constellation, prn: i as u8 + 1 },
                    observations: Vec::new(),
                })
                .collect(),
        }
    }

    #[test]
    fn filter_constellations_keeps_only_allowed_letters() {
        let mut epochs = vec![epoch_with(&[
            Constellation::Gps,
            Constellation::Glonass,
            Constellation::Galileo,
            Constellation::Beidou,
        ])];
        filter_constellations(&mut epochs, "GE");
        let kept: alloc::vec::Vec<Constellation> =
            epochs[0].satellites.iter().map(|s| s.sat.constellation).collect();
        assert_eq!(kept, alloc::vec![Constellation::Gps, Constellation::Galileo]);
    }

    #[test]
    fn filter_constellations_empty_allowed_drops_everything() {
        let mut epochs = vec![epoch_with(&[Constellation::Gps, Constellation::Glonass])];
        filter_constellations(&mut epochs, "");
        assert!(epochs[0].satellites.is_empty());
    }

    #[test]
    fn filter_constellations_applies_across_all_epochs() {
        let mut epochs = vec![
            epoch_with(&[Constellation::Gps, Constellation::Glonass]),
            epoch_with(&[Constellation::Glonass, Constellation::Galileo]),
        ];
        filter_constellations(&mut epochs, "G");
        assert_eq!(epochs[0].satellites.len(), 1);
        assert_eq!(epochs[1].satellites.len(), 0);
    }

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
            signal: SignalCode {
                freq_band: 1,
                attribute: 'C',
            },
        };
        assert_eq!(code1.to_string(), "C1C");

        let code2 = ObsCode {
            obs_type: ObsType::CarrierPhase,
            signal: SignalCode {
                freq_band: 2,
                attribute: 'W',
            },
        };
        assert_eq!(code2.to_string(), "L2W");
    }
}
