use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Error as IoError};
use chrono::{DateTime, TimeZone, Utc, NaiveDate};

#[derive(Debug, Clone)]
pub struct AntennaPcv {
    pub antenna_type: String,
    pub serial_num: String,
    pub valid_from: Option<DateTime<Utc>>,
    pub valid_until: Option<DateTime<Utc>>,
    pub dzen: f64,
    pub zen1: f64,
    pub zen2: f64,
    pub dazi: f64,
    pub frequencies: HashMap<String, FrequencyPcv>,
}

#[derive(Debug, Clone)]
pub struct FrequencyPcv {
    pub frequency_code: String,
    pub pco: nalgebra::Vector3<f64>, // North, East, Up in millimeters
    pub noazi: Vec<f64>,             // Nadir-dependent (or zenith) corrections
    pub azi: Option<Vec<Vec<f64>>>,  // Azimuth-dependent corrections
}

#[derive(Debug)]
pub enum AntexError {
    Io(IoError),
    ParseError(String),
}

impl From<IoError> for AntexError {
    fn from(err: IoError) -> Self {
        AntexError::Io(err)
    }
}

pub struct AntexDatabase {
    pub antennas: Vec<AntennaPcv>,
}

impl AntexDatabase {
    pub fn parse<P: AsRef<std::path::Path>>(path: P) -> Result<Self, AntexError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        
        let mut antennas = Vec::new();
        let mut current_antenna: Option<AntennaPcv> = None;
        let mut current_frequency: Option<FrequencyPcv> = None;
        
        let mut in_antenna = false;
        
        for line_res in reader.lines() {
            let line = line_res?;
            if line.len() < 60 {
                continue;
            }
            
            let label = line[60..].trim();
            
            if label == "START OF ANTENNA" {
                in_antenna = true;
                current_antenna = Some(AntennaPcv {
                    antenna_type: String::new(),
                    serial_num: String::new(),
                    valid_from: None,
                    valid_until: None,
                    dzen: 0.0,
                    zen1: 0.0,
                    zen2: 0.0,
                    dazi: 0.0,
                    frequencies: HashMap::new(),
                });
            } else if label == "END OF ANTENNA" {
                if let Some(ant) = current_antenna.take() {
                    antennas.push(ant);
                }
                in_antenna = false;
            } else if in_antenna {
                let ant = current_antenna.as_mut().unwrap();
                match label {
                    "TYPE / SERIAL NO" => {
                        ant.antenna_type = line[0..20].trim().to_string();
                        ant.serial_num = line[20..40].trim().to_string();
                    }
                    "VALID FROM" => {
                        ant.valid_from = parse_antex_date(&line[0..60]);
                    }
                    "VALID UNTIL" => {
                        ant.valid_until = parse_antex_date(&line[0..60]);
                    }
                    "ZEN1 / ZEN2 / DZEN" => {
                        let parts: Vec<&str> = line[0..60].split_whitespace().collect();
                        if parts.len() >= 3 {
                            ant.zen1 = parts[0].parse().unwrap_or(0.0);
                            ant.zen2 = parts[1].parse().unwrap_or(0.0);
                            ant.dzen = parts[2].parse().unwrap_or(0.0);
                        }
                    }
                    "DAZI" => {
                        ant.dazi = line[0..60].trim().parse().unwrap_or(0.0);
                    }
                    "START OF FREQUENCY" => {
                        let code = line[0..60].trim().to_string();
                        current_frequency = Some(FrequencyPcv {
                            frequency_code: code.clone(),
                            pco: nalgebra::Vector3::zeros(),
                            noazi: Vec::new(),
                            azi: None,
                        });
                    }
                    "END OF FREQUENCY" => {
                        if let Some(freq) = current_frequency.take() {
                            ant.frequencies.insert(freq.frequency_code.clone(), freq);
                        }
                    }
                    "NORTH / EAST / UP" => {
                        if let Some(freq) = current_frequency.as_mut() {
                            let parts: Vec<&str> = line[0..60].split_whitespace().collect();
                            if parts.len() >= 3 {
                                let north: f64 = parts[0].parse().unwrap_or(0.0);
                                let east: f64 = parts[1].parse().unwrap_or(0.0);
                                let up: f64 = parts[2].parse().unwrap_or(0.0);
                                freq.pco = nalgebra::Vector3::new(north, east, up);
                            }
                        }
                    }
                    _ => {
                        if let Some(freq) = current_frequency.as_mut() {
                            if line.starts_with("   NOAZI") {
                                let values: Vec<f64> = line[8..60]
                                    .split_whitespace()
                                    .filter_map(|s| s.parse().ok())
                                    .collect();
                                freq.noazi = values;
                            }
                        }
                    }
                }
            }
        }
        
        Ok(AntexDatabase { antennas })
    }
    
    pub fn find_satellite(&self, prn: &str, time: DateTime<Utc>) -> Option<&AntennaPcv> {
        self.antennas.iter().find(|a| {
            if a.serial_num == prn {
                let after_from = a.valid_from.is_none_or(|from| time >= from);
                let before_until = a.valid_until.is_none_or(|until| time <= until);
                after_from && before_until
            } else {
                false
            }
        })
    }
}

fn parse_antex_date(s: &str) -> Option<DateTime<Utc>> {
    let parts: Vec<i32> = s.split_whitespace().filter_map(|x| x.parse().ok()).collect();
    if parts.len() >= 6 {
        let date = NaiveDate::from_ymd_opt(parts[0], parts[1] as u32, parts[2] as u32)?;
        let dt = date.and_hms_opt(parts[3] as u32, parts[4] as u32, parts[5] as u32)?;
        Some(Utc.from_utc_datetime(&dt))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_parse_igs14_antex() {
        let path = PathBuf::from("../../datasets/igs14.atx");
        if !path.exists() {
            return; // Skip if dataset not available
        }
        
        let db = AntexDatabase::parse(&path).unwrap();
        assert!(db.antennas.len() > 100);
        
        // Find a specific satellite (e.g., G01)
        let g01 = db.find_satellite("G01", Utc.with_ymd_and_hms(2010, 1, 1, 0, 0, 0).unwrap()).unwrap();
        assert!(g01.antenna_type.starts_with("BLOCK IIA"));
        
        let freq_g01 = g01.frequencies.get("G01").unwrap();
        assert_eq!(freq_g01.pco.x, 279.0);
        assert_eq!(freq_g01.pco.y, 0.0);
        assert_eq!(freq_g01.pco.z, 2319.5);
    }
}
