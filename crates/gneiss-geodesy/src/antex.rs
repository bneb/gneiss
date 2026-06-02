use alloc::string::String;
use alloc::vec::Vec;
use nalgebra::Vector3;

#[derive(Debug, Clone, PartialEq)]
pub struct PcvData {
    /// Satellite system and frequency code (e.g., 'G01' for GPS L1)
    pub freq_code: String,
    /// Phase Center Offset: North, East, Up in millimeters
    pub pco: Vector3<f64>,
    /// Phase Center Variations (independent of azimuth) in millimeters.
    /// Sampled at increments defined by zenith_grid.
    pub noazi: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AntennaModel {
    pub ant_type: String,
    pub serial_num: String,
    pub zenith_start: f64,
    pub zenith_stop: f64,
    pub zenith_inc: f64,
    pub frequencies: Vec<PcvData>,
}

impl AntennaModel {
    /// Interpolates the Phase Center Variation (PCV) in millimeters for a given frequency and zenith angle (in radians).
    pub fn interpolate_pcv(&self, freq_code: &str, zenith_rad: f64) -> Option<f64> {
        let freq = self.frequencies.iter().find(|f| f.freq_code == freq_code)?;
        let zenith_deg = zenith_rad * 180.0 / core::f64::consts::PI;
        
        if zenith_deg < self.zenith_start || zenith_deg > self.zenith_stop || self.zenith_inc <= 0.0 {
            // If outside the defined range, use the edge value or None.
            // Often antennas are only defined up to 90 deg.
            if zenith_deg > self.zenith_stop && zenith_deg <= 90.0 && freq.noazi.last().is_some() {
                 return Some(*freq.noazi.last().unwrap());
            }
            if zenith_deg < self.zenith_start && !freq.noazi.is_empty() {
                 return Some(*freq.noazi.first().unwrap());
            }
            return None;
        }
        
        let index_float = (zenith_deg - self.zenith_start) / self.zenith_inc;
        let idx0 = libm::floor(index_float) as usize;
        let idx1 = libm::ceil(index_float) as usize;
        
        if idx0 >= freq.noazi.len() {
            return None;
        }
        
        if idx0 == idx1 || idx1 >= freq.noazi.len() {
            return Some(freq.noazi[idx0]);
        }
        
        let w1 = index_float - (idx0 as f64);
        let w0 = 1.0 - w1;
        
        Some(w0 * freq.noazi[idx0] + w1 * freq.noazi[idx1])
    }
}

fn parse_antex_data_line(
    line: &str,
    m: &mut AntennaModel,
    current_freq: &mut Option<String>,
    current_pco: &mut Option<Vector3<f64>>,
) {
    let data_part = if line.len() >= 60 { &line[0..60] } else { line };
    let parts: Vec<&str> = data_part.split_whitespace().collect();

    if line.ends_with("TYPE / SERIAL NO") {
        if !parts.is_empty() {
            m.ant_type = parts[0].into();
        }
        if parts.len() >= 2 {
            m.serial_num = parts[1].into();
        }
    } else if line.ends_with("ZEN1 / ZEN2 / DZEN") {
        if parts.len() >= 3 {
            m.zenith_start = parts[0].parse().unwrap_or(0.0);
            m.zenith_stop = parts[1].parse().unwrap_or(0.0);
            m.zenith_inc = parts[2].parse().unwrap_or(0.0);
        }
    } else if line.ends_with("START OF FREQUENCY") {
        if !parts.is_empty() {
            *current_freq = Some(parts[0].into());
        }
    } else if line.ends_with("NORTH / EAST / UP") {
        if parts.len() >= 3 {
            let n: f64 = parts[0].parse().unwrap_or(0.0);
            let e: f64 = parts[1].parse().unwrap_or(0.0);
            let u: f64 = parts[2].parse().unwrap_or(0.0);
            *current_pco = Some(Vector3::new(n, e, u));
        }
    } else if line.starts_with("   NOAZI") {
        if let (Some(freq), Some(pco)) = (&current_freq, &current_pco) {
            let noazi_parts: Vec<&str> = line.split_whitespace().skip(1).collect();
            let noazi: Vec<f64> = noazi_parts.iter().map(|s| s.parse().unwrap_or(0.0)).collect();
            m.frequencies.push(PcvData {
                freq_code: freq.clone(),
                pco: *pco,
                noazi,
            });
        }
    } else if line.ends_with("END OF FREQUENCY") {
        *current_freq = None;
        *current_pco = None;
    }
}

/// Parses an ANTEX format file and returns a list of parsed Antenna Models.
pub fn parse_antex(content: &str) -> Result<Vec<AntennaModel>, String> {
    let mut models = Vec::new();
    let mut current_model: Option<AntennaModel> = None;
    let mut current_freq: Option<String> = None;
    let mut current_pco: Option<Vector3<f64>> = None;
    
    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }

        if line.ends_with("START OF ANTENNA") {
            current_model = Some(AntennaModel {
                ant_type: String::new(),
                serial_num: String::new(),
                zenith_start: 0.0,
                zenith_stop: 0.0,
                zenith_inc: 0.0,
                frequencies: Vec::new(),
            });
        } else if line.ends_with("END OF ANTENNA") {
            if let Some(m) = current_model.take() {
                models.push(m);
            }
        } else if let Some(m) = &mut current_model {
            parse_antex_data_line(line, m, &mut current_freq, &mut current_pco);
        }
    }
    
    Ok(models)
}

#[cfg(test)]
mod tests {
    use super::*;

    const MOCK_ANTEX: &str = r#"
                                                            START OF HEADER
                                                            END OF HEADER
                                                            START OF ANTENNA
IGS01/0000          NONE                                    TYPE / SERIAL NO
     0.0                                                    DAZI
     0.0  90.0   5.0                                        ZEN1 / ZEN2 / DZEN
     2                                                      # OF FREQUENCIES
G01                                                         START OF FREQUENCY
      1.00      2.00      3.00                              NORTH / EAST / UP
   NOAZI    0.00    1.00    2.00    3.00    4.00    5.00    6.00    7.00    8.00    9.00   10.00   11.00   12.00   13.00   14.00   15.00   16.00   17.00   18.00
G01                                                         END OF FREQUENCY
G02                                                         START OF FREQUENCY
      4.00      5.00      6.00                              NORTH / EAST / UP
   NOAZI    0.00   -1.00   -2.00   -3.00   -4.00   -5.00   -6.00   -7.00   -8.00   -9.00  -10.00  -11.00  -12.00  -13.00  -14.00  -15.00  -16.00  -17.00  -18.00
G02                                                         END OF FREQUENCY
                                                            END OF ANTENNA
"#;

    #[test]
    fn test_parse_antex() {
        let models = parse_antex(MOCK_ANTEX).unwrap();
        assert_eq!(models.len(), 1);
        
        let ant = &models[0];
        assert_eq!(ant.ant_type, "IGS01/0000");
        assert_eq!(ant.serial_num, "NONE");
        assert_eq!(ant.zenith_start, 0.0);
        assert_eq!(ant.zenith_stop, 90.0);
        assert_eq!(ant.zenith_inc, 5.0);
        
        assert_eq!(ant.frequencies.len(), 2);
        
        let g01 = &ant.frequencies[0];
        assert_eq!(g01.freq_code, "G01");
        assert_eq!(g01.pco, Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(g01.noazi.len(), 19); // 0 to 90 with step 5 = 19 values
        assert_eq!(g01.noazi[2], 2.0); // 10 degrees

        let g02 = &ant.frequencies[1];
        assert_eq!(g02.freq_code, "G02");
        assert_eq!(g02.pco, Vector3::new(4.0, 5.0, 6.0));
        assert_eq!(g02.noazi[2], -2.0);
    }

    #[test]
    fn test_interpolate_pcv() {
        // Assume ant model is parsed as above
        let mut ant = AntennaModel {
            ant_type: "TEST".into(),
            serial_num: "NONE".into(),
            zenith_start: 0.0,
            zenith_stop: 90.0,
            zenith_inc: 5.0,
            frequencies: Vec::new(),
        };

        // Create mock NOAZI: 0mm at 0deg, 5mm at 5deg, 10mm at 10deg.
        let mut noazi = Vec::new();
        for i in 0..=18 {
            noazi.push((i as f64) * 5.0);
        }

        ant.frequencies.push(PcvData {
            freq_code: "G01".into(),
            pco: Vector3::zeros(),
            noazi,
        });

        // 0 degrees (zenith = 0 rad) -> 0.0 mm
        assert_eq!(ant.interpolate_pcv("G01", 0.0).unwrap(), 0.0);

        // 5 degrees -> 5.0 mm
        assert_eq!(ant.interpolate_pcv("G01", 5.0 * core::f64::consts::PI / 180.0).unwrap(), 5.0);

        // 7.5 degrees (midpoint) -> 7.5 mm
        let mid_rad = 7.5 * core::f64::consts::PI / 180.0;
        let pcv = ant.interpolate_pcv("G01", mid_rad).unwrap();
        assert!((pcv - 7.5).abs() < 1e-6, "Expected 7.5, got {}", pcv);
    }
}
