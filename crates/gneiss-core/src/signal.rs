use crate::sat::{SatelliteId, Constellation};

pub fn satellite_frequencies(sat: SatelliteId, freq_num: i8) -> (f64, f64) {
    match sat.constellation {
        Constellation::Gps | Constellation::Qzss => (1575.42e6, 1227.60e6), // L1, L2
        Constellation::Galileo => (1575.42e6, 1207.140e6), // E1, E5b
        Constellation::Beidou => (1561.098e6, 1207.140e6), // B1I, B2I
        Constellation::Glonass => {
            let f1 = 1602.0e6 + (freq_num as f64) * 0.5625e6;
            let f2 = 1246.0e6 + (freq_num as f64) * 0.4375e6;
            (f1, f2)
        }
        _ => (1575.42e6, 1227.60e6)
    }
}

pub fn get_wavelength(sat: SatelliteId, freq_band: u8, freq_num: i8) -> f64 {
    let (f1, f2) = satellite_frequencies(sat, freq_num);
    let freq = if freq_band == 1 { f1 } else { f2 };
    299792458.0 / freq
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_glonass_fdma_wavelengths() {
        // Channel -4 (e.g. GLONASS PRN 6 in our dataset)
        let sat = SatelliteId { constellation: Constellation::Glonass, prn: 6 };
        let (f1, f2) = satellite_frequencies(sat, -4);
        assert_eq!(f1, 1602.0e6 - 4.0 * 0.5625e6);
        assert_eq!(f2, 1246.0e6 - 4.0 * 0.4375e6);
        
        let w1 = get_wavelength(sat, 1, -4);
        assert!((w1 - 0.1873).abs() < 0.01);
    }

    #[test]
    fn test_galileo_frequencies() {
        let sat = SatelliteId { constellation: Constellation::Galileo, prn: 11 };
        let (f1, f2) = satellite_frequencies(sat, 0);
        assert_eq!(f1, 1575.42e6);
        assert_eq!(f2, 1207.140e6);
    }
}
