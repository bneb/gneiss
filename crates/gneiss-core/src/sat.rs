use core::fmt;
use core::str::FromStr;

/// Represents a GNSS constellation.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
pub enum Constellation {
    Gps,
    Glonass,
    Galileo,
    Beidou,
    Sbas,
    Qzss,
    Navic,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseConstellationError;

impl fmt::Display for ParseConstellationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid constellation character")
    }
}

impl FromStr for Constellation {
    type Err = ParseConstellationError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "G" => Ok(Constellation::Gps),
            "R" => Ok(Constellation::Glonass),
            "E" => Ok(Constellation::Galileo),
            "C" => Ok(Constellation::Beidou),
            "S" => Ok(Constellation::Sbas),
            "J" => Ok(Constellation::Qzss),
            "I" => Ok(Constellation::Navic),
            _ => Err(ParseConstellationError),
        }
    }
}

impl fmt::Display for Constellation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let c = match self {
            Constellation::Gps => "G",
            Constellation::Glonass => "R",
            Constellation::Galileo => "E",
            Constellation::Beidou => "C",
            Constellation::Sbas => "S",
            Constellation::Qzss => "J",
            Constellation::Navic => "I",
        };
        write!(f, "{}", c)
    }
}

/// Identifies a specific satellite within a constellation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SatelliteId {
    pub constellation: Constellation,
    pub prn: u8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseSatelliteIdError;

impl fmt::Display for ParseSatelliteIdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid satellite id format")
    }
}

impl FromStr for SatelliteId {
    type Err = ParseSatelliteIdError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() < 2 {
            return Err(ParseSatelliteIdError);
        }
        let (const_str, prn_str) = s.split_at(1);
        let constellation =
            Constellation::from_str(const_str).map_err(|_| ParseSatelliteIdError)?;
        let prn = u8::from_str(prn_str).map_err(|_| ParseSatelliteIdError)?;
        if prn == 0 {
            return Err(ParseSatelliteIdError);
        }
        Ok(SatelliteId { constellation, prn })
    }
}

impl fmt::Display for SatelliteId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{:02}", self.constellation, self.prn)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::format;
    use alloc::string::ToString;

    #[test]
    fn test_constellation_parsing() {
        assert_eq!(Constellation::from_str("G").unwrap(), Constellation::Gps);
        assert_eq!(
            Constellation::from_str("R").unwrap(),
            Constellation::Glonass
        );
        assert_eq!(
            Constellation::from_str("E").unwrap(),
            Constellation::Galileo
        );
        assert_eq!(Constellation::from_str("C").unwrap(), Constellation::Beidou);
        assert_eq!(Constellation::from_str("S").unwrap(), Constellation::Sbas);
        assert_eq!(Constellation::from_str("J").unwrap(), Constellation::Qzss);
        assert_eq!(Constellation::from_str("I").unwrap(), Constellation::Navic);

        assert!(Constellation::from_str("X").is_err());
    }

    #[test]
    fn test_constellation_display() {
        assert_eq!(Constellation::Gps.to_string(), "G");
        assert_eq!(Constellation::Glonass.to_string(), "R");
        assert_eq!(Constellation::Galileo.to_string(), "E");
        assert_eq!(Constellation::Beidou.to_string(), "C");
        assert_eq!(Constellation::Sbas.to_string(), "S");
        assert_eq!(Constellation::Qzss.to_string(), "J");
        assert_eq!(Constellation::Navic.to_string(), "I");
    }

    #[test]
    fn test_satellite_id_parsing() {
        assert_eq!(
            SatelliteId::from_str("G01").unwrap(),
            SatelliteId {
                constellation: Constellation::Gps,
                prn: 1
            }
        );
        assert_eq!(
            SatelliteId::from_str("R24").unwrap(),
            SatelliteId {
                constellation: Constellation::Glonass,
                prn: 24
            }
        );
        assert_eq!(
            SatelliteId::from_str("E11").unwrap(),
            SatelliteId {
                constellation: Constellation::Galileo,
                prn: 11
            }
        );

        assert!(SatelliteId::from_str("G").is_err()); // Missing PRN
        assert!(SatelliteId::from_str("G00").is_err()); // Invalid PRN 0
        assert!(SatelliteId::from_str("X01").is_err()); // Unknown Constellation
        assert!(SatelliteId::from_str("G1A").is_err()); // Invalid number
    }

    #[test]
    fn test_satellite_id_display() {
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 5,
        };
        assert_eq!(sat.to_string(), "G05");

        let sat = SatelliteId {
            constellation: Constellation::Glonass,
            prn: 24,
        };
        assert_eq!(sat.to_string(), "R24");
    }

    #[test]
    fn test_parse_constellation_error_display() {
        let err = ParseConstellationError;
        assert_eq!(format!("{}", err), "invalid constellation character");
    }

    #[test]
    fn test_parse_satellite_id_error_display() {
        let err = ParseSatelliteIdError;
        assert_eq!(format!("{}", err), "invalid satellite id format");
    }

    #[test]
    fn test_satellite_id_parsing_all_constellations() {
        assert_eq!(
            SatelliteId::from_str("C15").unwrap(),
            SatelliteId {
                constellation: Constellation::Beidou,
                prn: 15
            }
        );
        assert_eq!(
            SatelliteId::from_str("S36").unwrap(),
            SatelliteId {
                constellation: Constellation::Sbas,
                prn: 36
            }
        );
        assert_eq!(
            SatelliteId::from_str("J07").unwrap(),
            SatelliteId {
                constellation: Constellation::Qzss,
                prn: 7
            }
        );
        assert_eq!(
            SatelliteId::from_str("I03").unwrap(),
            SatelliteId {
                constellation: Constellation::Navic,
                prn: 3
            }
        );
    }

    #[test]
    fn test_satellite_id_parsing_single_digit_prn() {
        // Without leading zero
        assert_eq!(
            SatelliteId::from_str("G1").unwrap(),
            SatelliteId {
                constellation: Constellation::Gps,
                prn: 1
            }
        );
        assert_eq!(
            SatelliteId::from_str("R5").unwrap(),
            SatelliteId {
                constellation: Constellation::Glonass,
                prn: 5
            }
        );
    }

    #[test]
    fn test_satellite_id_display_all_constellations() {
        let cases = [
            (Constellation::Gps, 1, "G01"),
            (Constellation::Glonass, 7, "R07"),
            (Constellation::Galileo, 12, "E12"),
            (Constellation::Beidou, 3, "C03"),
            (Constellation::Sbas, 40, "S40"),
            (Constellation::Qzss, 8, "J08"),
            (Constellation::Navic, 9, "I09"),
        ];
        for (constellation, prn, expected) in &cases {
            let sat = SatelliteId {
                constellation: *constellation,
                prn: *prn,
            };
            assert_eq!(
                sat.to_string(),
                *expected,
                "Display failed for {:?} PRN {}",
                constellation,
                prn
            );
        }
    }

    #[test]
    fn test_satellite_id_parsing_prn_overflow() {
        // u8 can hold up to 255; 256 should fail
        assert!(SatelliteId::from_str("G256").is_err());
    }

    #[test]
    fn test_satellite_id_parsing_too_short() {
        assert!(SatelliteId::from_str("").is_err());
        assert!(SatelliteId::from_str("G").is_err());
    }

    #[test]
    fn test_satellite_id_parsing_negative_prn() {
        assert!(SatelliteId::from_str("G-1").is_err());
    }

    #[test]
    fn test_satellite_id_equality() {
        let a = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let b = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let c = SatelliteId {
            constellation: Constellation::Glonass,
            prn: 1,
        };
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn test_constellation_ordering() {
        // Ensure Ord is derived and works
        assert!(Constellation::Gps < Constellation::Glonass);
        assert!(Constellation::Gps < Constellation::Galileo);
        assert!(Constellation::Galileo > Constellation::Glonass);
    }
}
