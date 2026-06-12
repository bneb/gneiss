use crate::sat::{SatelliteId, Constellation};

pub const FREQ_GPS_L1: f64 = 1575.42e6;
pub const FREQ_GPS_L2: f64 = 1227.60e6;
pub const FREQ_GPS_L5: f64 = 1176.45e6;
pub const FREQ_GAL_E5B: f64 = 1207.140e6;
pub const FREQ_BDS_B1I: f64 = 1561.098e6;
pub const FREQ_GLO_L1_NOMINAL: f64 = 1602.0e6;
pub const FREQ_GLO_L2_NOMINAL: f64 = 1246.0e6;
pub const FREQ_GLO_L1_DELTA: f64 = 0.5625e6;
pub const FREQ_GLO_L2_DELTA: f64 = 0.4375e6;

pub fn satellite_frequencies(sat: SatelliteId, freq_num: i8) -> (f64, f64) {
    match sat.constellation {
        Constellation::Gps | Constellation::Qzss => (FREQ_GPS_L1, FREQ_GPS_L2),
        Constellation::Galileo => (FREQ_GPS_L1, FREQ_GAL_E5B), // E1 shares GPS L1 freq
        Constellation::Beidou => (FREQ_BDS_B1I, FREQ_GAL_E5B), // B2I shares E5b freq
        Constellation::Glonass => {
            let f1 = FREQ_GLO_L1_NOMINAL + (freq_num as f64) * FREQ_GLO_L1_DELTA;
            let f2 = FREQ_GLO_L2_NOMINAL + (freq_num as f64) * FREQ_GLO_L2_DELTA;
            (f1, f2)
        }
        _ => (FREQ_GPS_L1, FREQ_GPS_L2)
    }
}

pub fn get_frequency(sat: SatelliteId, freq_band: u8, freq_num: i8) -> f64 {
    match freq_band {
        1 => {
            match sat.constellation {
                Constellation::Gps | Constellation::Qzss | Constellation::Galileo => FREQ_GPS_L1,
                Constellation::Beidou => FREQ_BDS_B1I,
                Constellation::Glonass => FREQ_GLO_L1_NOMINAL + (freq_num as f64) * FREQ_GLO_L1_DELTA,
                _ => FREQ_GPS_L1,
            }
        },
        2 => {
            match sat.constellation {
                Constellation::Gps | Constellation::Qzss => FREQ_GPS_L2,
                Constellation::Galileo => FREQ_GAL_E5B,
                Constellation::Beidou => FREQ_GAL_E5B, // B2I
                Constellation::Glonass => FREQ_GLO_L2_NOMINAL + (freq_num as f64) * FREQ_GLO_L2_DELTA,
                _ => FREQ_GPS_L2,
            }
        },
        5 => {
            match sat.constellation {
                Constellation::Gps | Constellation::Qzss | Constellation::Galileo => FREQ_GPS_L5,
                Constellation::Beidou => FREQ_GPS_L5,
                _ => FREQ_GPS_L5,
            }
        },
        _ => FREQ_GPS_L1,
    }
}

pub fn get_wavelength(sat: SatelliteId, freq_band: u8, freq_num: i8) -> f64 {
    let freq = get_frequency(sat, freq_band, freq_num);
    crate::constants::SPEED_OF_LIGHT_M_S / freq
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_glonass_fdma_wavelengths() {
        // Channel -4 (e.g. GLONASS PRN 6 in our dataset)
        let sat = SatelliteId { constellation: Constellation::Glonass, prn: 6 };
        let (f1, f2) = satellite_frequencies(sat, -4);
        assert_eq!(f1, FREQ_GLO_L1_NOMINAL - 4.0 * FREQ_GLO_L1_DELTA);
        assert_eq!(f2, FREQ_GLO_L2_NOMINAL - 4.0 * FREQ_GLO_L2_DELTA);
        
        let w1 = get_wavelength(sat, 1, -4);
        assert!((w1 - 0.1873).abs() < 0.01);
    }

    #[test]
    fn test_galileo_frequencies() {
        let sat = SatelliteId { constellation: Constellation::Galileo, prn: 11 };
        let (f1, f2) = satellite_frequencies(sat, 0);
        assert_eq!(f1, FREQ_GPS_L1);
        assert_eq!(f2, FREQ_GAL_E5B);
    }
}
