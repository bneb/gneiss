use chrono::{DateTime, NaiveDate, TimeZone, Utc};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Error as IoError};

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
    /// Index: antenna_type → position in `antennas` Vec. Built at end of parse.
    by_type: HashMap<String, usize>,
    /// Index: serial_num → list of positions for O(1) satellite lookups.
    by_sat: HashMap<String, Vec<usize>>,
}

impl AntexDatabase {
    /// Construct from a pre-built antenna list, indexing for O(1) lookup.
    pub fn new(antennas: Vec<AntennaPcv>) -> Self {
        let by_type: HashMap<String, usize> = antennas
            .iter()
            .enumerate()
            .map(|(i, a)| (a.antenna_type.clone(), i))
            .collect();
        let mut by_sat: HashMap<String, Vec<usize>> = HashMap::new();
        for (i, a) in antennas.iter().enumerate() {
            by_sat.entry(a.serial_num.clone()).or_default().push(i);
        }
        Self { antennas, by_type, by_sat }
    }

    /// O(1) lookup by antenna type. Returns None if not found.
    pub fn get_antenna(&self, antenna_type: &str) -> Option<&AntennaPcv> {
        self.by_type
            .get(antenna_type)
            .and_then(|&idx| self.antennas.get(idx))
    }
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
                let ant = current_antenna.as_mut().expect("antenna block preceded by START OF ANTENNA");
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
                            if let Some(rest) = line.strip_prefix("   NOAZI") {
                                // NOAZI records extend past the fixed 60-col
                                // label boundary (up to ~152 chars); parse
                                // the whole remainder or the grid is
                                // silently truncated to ~6 nodes.
                                let values: Vec<f64> = rest
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

        Ok(AntexDatabase::new(antennas))
    }

    pub fn find_satellite(&self, prn: &str, time: DateTime<Utc>) -> Option<&AntennaPcv> {
        // O(1) lookup via sat index, then scan only matching entries for time window
        self.by_sat.get(prn).and_then(|indices| {
            indices.iter().find_map(|&idx| {
                let a = &self.antennas[idx];
                let after_from = a.valid_from.is_none_or(|from| time >= from);
                let before_until = a.valid_until.is_none_or(|until| time <= until);
                if after_from && before_until { Some(a) } else { None }
            })
        })
    }
}

fn parse_antex_date(s: &str) -> Option<DateTime<Utc>> {
    let parts: Vec<i32> = s
        .split_whitespace()
        .filter_map(|x| x.parse().ok())
        .collect();
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
    use std::sync::atomic::{AtomicU32, Ordering};

    // ------------------------------------------------------------------
    // Helpers for building minimal ANTEX test files
    // ------------------------------------------------------------------

    /// Build a single ANTEX line: data padded to 60 chars followed by `label`.
    fn ant_line(data: &str, label: &str) -> String {
        let mut line = String::with_capacity(80);
        line.push_str(data);
        while line.len() < 60 {
            line.push(' ');
        }
        line.push_str(label);
        line.push('\n');
        line
    }

    /// Atomically-increasing counter for temp file names so parallel
    /// `cargo test` runs do not collide.
    static ANTEX_COUNTER: AtomicU32 = AtomicU32::new(0);

    /// Write `content` to a uniquely-named temp file and return the path.
    fn write_temp_antex(content: &str) -> PathBuf {
        let dir = std::env::temp_dir();
        let n = ANTEX_COUNTER.fetch_add(1, Ordering::SeqCst);
        let path = dir.join(format!("test_antex_{}.atx", n));
        std::fs::write(&path, content).unwrap();
        path
    }

    /// A minimal valid antenna block (single frequency G01).
    fn minimal_antenna_block(ant_type: &str, serial: &str, freq_code: &str) -> String {
        let mut b = String::new();
        b.push_str(&ant_line("", "START OF ANTENNA"));
        b.push_str(&ant_line(
            &format!("{:20}{:20}", ant_type, serial),
            "TYPE / SERIAL NO",
        ));
        b.push_str(&ant_line("     0.0", "DAZI"));
        b.push_str(&ant_line("     0.0  17.0   1.0", "ZEN1 / ZEN2 / DZEN"));
        b.push_str(&ant_line(
            "  2020     1    15     0     0    0.0000000",
            "VALID FROM",
        ));
        b.push_str(&ant_line(
            "  2030     1    15     0     0    0.0000000",
            "VALID UNTIL",
        ));
        b.push_str(&ant_line(&format!("   {:4}", freq_code), "START OF FREQUENCY"));
        b.push_str(&ant_line(
            "      1.00      2.00      3.00",
            "NORTH / EAST / UP",
        ));
        b.push_str(&ant_line(
            "   NOAZI    0.10    0.20    0.30",
            "",
        ));
        b.push_str(&ant_line("", "END OF FREQUENCY"));
        b.push_str(&ant_line("", "END OF ANTENNA"));
        b
    }

    // ------------------------------------------------------------------
    // Tests
    // ------------------------------------------------------------------

    #[test]
    fn test_parse_minimal_antex() {
        let content = minimal_antenna_block("TEST_ANT", "G01", "G01");
        let path = write_temp_antex(&content);
        let db = AntexDatabase::parse(&path).unwrap();
        std::fs::remove_file(&path).ok();

        assert_eq!(db.antennas.len(), 1);
        let ant = &db.antennas[0];
        assert_eq!(ant.antenna_type, "TEST_ANT");
        assert_eq!(ant.serial_num, "G01");
        assert_eq!(ant.dazi, 0.0);
        assert!((ant.zen1 - 0.0).abs() < 1e-12);
        assert!((ant.zen2 - 17.0).abs() < 1e-12);
        assert!((ant.dzen - 1.0).abs() < 1e-12);

        let freq = ant.frequencies.get("G01").unwrap();
        assert!((freq.pco.x - 1.0).abs() < 1e-12);
        assert!((freq.pco.y - 2.0).abs() < 1e-12);
        assert!((freq.pco.z - 3.0).abs() < 1e-12);
        assert_eq!(freq.noazi.len(), 3);
        assert!((freq.noazi[0] - 0.10).abs() < 1e-12);
    }

    #[test]
    fn test_parse_multiple_antennas() {
        let mut content = String::new();
        content.push_str(&minimal_antenna_block("ANT_A", "G01", "G01"));
        content.push_str(&minimal_antenna_block("ANT_B", "R02", "G01"));

        let path = write_temp_antex(&content);
        let db = AntexDatabase::parse(&path).unwrap();
        std::fs::remove_file(&path).ok();

        assert_eq!(db.antennas.len(), 2);
        assert_eq!(db.antennas[0].antenna_type, "ANT_A");
        assert_eq!(db.antennas[1].antenna_type, "ANT_B");
    }

    #[test]
    fn test_parse_invalid_path_returns_error() {
        let result = AntexDatabase::parse("/nonexistent/path/antex.atx");
        assert!(result.is_err());
        let is_io_err = matches!(result, Err(AntexError::Io(_)));
        assert!(is_io_err, "expected Io error, got unexpected result variant");
    }

    #[test]
    fn test_find_satellite_matches_serial() {
        let content = minimal_antenna_block("BLOCK_IIA", "G01", "G01");
        let path = write_temp_antex(&content);
        let db = AntexDatabase::parse(&path).unwrap();
        std::fs::remove_file(&path).ok();

        let found = db.find_satellite("G01", Utc.with_ymd_and_hms(2025, 6, 1, 0, 0, 0).unwrap());
        assert!(found.is_some());
        assert_eq!(found.unwrap().serial_num, "G01");
    }

    #[test]
    fn test_find_satellite_no_match() {
        let content = minimal_antenna_block("BLOCK_IIA", "G01", "G01");
        let path = write_temp_antex(&content);
        let db = AntexDatabase::parse(&path).unwrap();
        std::fs::remove_file(&path).ok();

        let found = db.find_satellite("R99", Utc.with_ymd_and_hms(2025, 6, 1, 0, 0, 0).unwrap());
        assert!(found.is_none());
    }

    #[test]
    fn test_parse_short_lines_skipped() {
        // Lines shorter than 60 characters are skipped.
        let short = "short line\n";
        let path = write_temp_antex(short);
        let db = AntexDatabase::parse(&path).unwrap();
        std::fs::remove_file(&path).ok();

        assert!(db.antennas.is_empty());
    }

    #[test]
    fn test_parse_frequency_multiple_frequencies() {
        let mut content = String::new();
        content.push_str(&ant_line("", "START OF ANTENNA"));
        content.push_str(&ant_line(
            &format!("{:20}{:20}", "MULTI_FREQ", "G01"),
            "TYPE / SERIAL NO",
        ));
        content.push_str(&ant_line("     0.0", "DAZI"));
        content.push_str(&ant_line("     0.0  17.0   1.0", "ZEN1 / ZEN2 / DZEN"));
        content.push_str(&ant_line(
            "  2020     1    15     0     0    0.0000000",
            "VALID FROM",
        ));
        content.push_str(&ant_line(
            "  2030     1    15     0     0    0.0000000",
            "VALID UNTIL",
        ));
        // Frequency G01
        content.push_str(&ant_line("   G01", "START OF FREQUENCY"));
        content.push_str(&ant_line(
            "     10.00     20.00     30.00",
            "NORTH / EAST / UP",
        ));
        content.push_str(&ant_line("   NOAZI    0.10    0.20", ""));
        content.push_str(&ant_line("", "END OF FREQUENCY"));
        // Frequency G02
        content.push_str(&ant_line("   G02", "START OF FREQUENCY"));
        content.push_str(&ant_line(
            "     40.00     50.00     60.00",
            "NORTH / EAST / UP",
        ));
        content.push_str(&ant_line("   NOAZI    0.30    0.40", ""));
        content.push_str(&ant_line("", "END OF FREQUENCY"));
        content.push_str(&ant_line("", "END OF ANTENNA"));

        let path = write_temp_antex(&content);
        let db = AntexDatabase::parse(&path).unwrap();
        std::fs::remove_file(&path).ok();

        assert_eq!(db.antennas[0].frequencies.len(), 2);
        let g01 = db.antennas[0].frequencies.get("G01").unwrap();
        assert!((g01.pco.z - 30.0).abs() < 1e-12);
        let g02 = db.antennas[0].frequencies.get("G02").unwrap();
        assert!((g02.pco.z - 60.0).abs() < 1e-12);
    }

    #[test]
    fn test_parse_noazi_rows_longer_than_60_cols() {
        // Real ANTEX NOAZI records run past the label column; every node
        // must survive parsing.
        let mut content = String::new();
        content.push_str(&ant_line("", "START OF ANTENNA"));
        content.push_str(&ant_line(&format!("{:<40}", "LONG_ROW"), "TYPE / SERIAL NO"));
        content.push_str(&ant_line("     0.0", "DAZI"));
        content.push_str(&ant_line("     0.0  90.0   5.0", "ZEN1 / ZEN2 / DZEN"));
        content.push_str(&ant_line("   G01", "START OF FREQUENCY"));
        content.push_str(&ant_line("      1.00      2.00      3.00", "NORTH / EAST / UP"));
        let values: Vec<String> = (0..19).map(|i| format!("{:>8.2}", i as f64 * -0.5)).collect();
        content.push_str(&format!("   NOAZI{}\n", values.join("")));
        content.push_str(&ant_line("", "END OF FREQUENCY"));
        content.push_str(&ant_line("", "END OF ANTENNA"));

        let path = write_temp_antex(&content);
        let db = AntexDatabase::parse(&path).unwrap();
        std::fs::remove_file(&path).ok();

        let noazi = &db.antennas[0].frequencies["G01"].noazi;
        assert_eq!(noazi.len(), 19);
        assert_eq!(noazi[0], 0.0);
        assert_eq!(noazi[9], -4.5);
        assert_eq!(noazi[18], -9.0);
    }

    #[test]
    fn test_parse_igs14_antex() {
        let path = PathBuf::from("../../datasets/igs14.atx");
        if !path.exists() {
            return; // Skip if dataset not available
        }

        let db = AntexDatabase::parse(&path).unwrap();
        assert!(db.antennas.len() > 100);

        // Find a specific satellite (e.g., G01)
        let g01 = db
            .find_satellite("G01", Utc.with_ymd_and_hms(2010, 1, 1, 0, 0, 0).unwrap())
            .unwrap();
        assert!(g01.antenna_type.starts_with("BLOCK IIA"));

        let freq_g01 = g01.frequencies.get("G01").unwrap();
        assert_eq!(freq_g01.pco.x, 279.0);
        assert_eq!(freq_g01.pco.y, 0.0);
        assert_eq!(freq_g01.pco.z, 2319.5);
    }
}
