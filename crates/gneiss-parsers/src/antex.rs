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

#[derive(Debug, Clone)]
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
                                let values: Vec<f64> = rest
                                    .split_whitespace()
                                    .filter_map(|s| s.parse().ok())
                                    .collect();
                                freq.noazi = values;
                            } else if ant.dazi > 0.0 {
                                let mut tokens = line.split_whitespace();
                                if let Some(_az) = tokens.next().and_then(|s| s.parse::<f64>().ok()) {
                                    let row: Vec<f64> = tokens.filter_map(|s| s.parse().ok()).collect();
                                    if !row.is_empty() {
                                        freq.azi.get_or_insert_with(Vec::new).push(row);
                                    }
                                }
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

    /// Look up satellite calibration by PRN and [`gneiss_core::time::GpsTime`].
    pub fn find_satellite_gps(&self, prn: &str, time: gneiss_core::time::GpsTime) -> Option<&AntennaPcv> {
        let unix_s = 315_964_800 + (time.week as i64) * 604_800 + (time.tow as i64);
        let dt = DateTime::from_timestamp(unix_s, 0)?;
        self.find_satellite(prn, dt)
    }

    /// Compute satellite nadir PCV in meters.
    pub fn satellite_pcv_nadir_m(&self, prn: &str, time: gneiss_core::time::GpsTime, nadir_deg: f64, freq_code: &str) -> f64 {
        self.find_satellite_gps(prn, time)
            .map_or(0.0, |ant| ant.interpolate_noazi_mm(freq_code, nadir_deg) * 1e-3)
    }
}

impl AntennaPcv {
    /// Compute nadir/zenith-dependent PCV in millimeters via linear interpolation of `noazi`.
    pub fn interpolate_noazi_mm(&self, freq_code: &str, nadir_deg: f64) -> f64 {
        let freq = match self.frequencies.get(freq_code) {
            Some(f) => f,
            None => return 0.0,
        };
        if freq.noazi.is_empty() || self.dzen <= 0.0 {
            return 0.0;
        }
        let clamped = nadir_deg.clamp(self.zen1, self.zen2);
        let idx_f = (clamped - self.zen1) / self.dzen;
        let i0 = (idx_f.floor() as usize).min(freq.noazi.len().saturating_sub(1));
        let i1 = (i0 + 1).min(freq.noazi.len().saturating_sub(1));
        let frac = idx_f - i0 as f64;
        (1.0 - frac) * freq.noazi[i0] + frac * freq.noazi[i1]
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
        let lines = [
            ("", "START OF ANTENNA"),
            (&format!("{:20}{:20}", ant_type, serial), "TYPE / SERIAL NO"),
            ("     0.0", "DAZI"),
            ("     0.0  17.0   1.0", "ZEN1 / ZEN2 / DZEN"),
            ("  2020     1    15     0     0    0.0000000", "VALID FROM"),
            ("  2030     1    15     0     0    0.0000000", "VALID UNTIL"),
            (&format!("   {:4}", freq_code), "START OF FREQUENCY"),
            ("      1.00      2.00      3.00", "NORTH / EAST / UP"),
            ("   NOAZI    0.10    0.20    0.30", ""),
            ("", "END OF FREQUENCY"),
            ("", "END OF ANTENNA"),
        ];
        lines.iter().map(|(d, l)| ant_line(d, l)).collect()
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
        let content = format!(
            "{}{}",
            minimal_antenna_block("ANT_A", "G01", "G01"),
            minimal_antenna_block("ANT_B", "R02", "G01")
        );
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
        assert!(matches!(result, Err(AntexError::Io(_))));
    }

    #[test]
    fn test_find_satellite() {
        let content = minimal_antenna_block("BLOCK_IIA", "G01", "G01");
        let path = write_temp_antex(&content);
        let db = AntexDatabase::parse(&path).unwrap();
        std::fs::remove_file(&path).ok();

        let t = Utc.with_ymd_and_hms(2025, 6, 1, 0, 0, 0).unwrap();
        let found = db.find_satellite("G01", t);
        assert_eq!(found.unwrap().serial_num, "G01");
        assert!(db.find_satellite("R99", t).is_none());
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
        let lines = [
            ("", "START OF ANTENNA"),
            (&format!("{:20}{:20}", "MULTI_FREQ", "G01"), "TYPE / SERIAL NO"),
            ("     0.0", "DAZI"),
            ("     0.0  17.0   1.0", "ZEN1 / ZEN2 / DZEN"),
            ("  2020     1    15     0     0    0.0000000", "VALID FROM"),
            ("  2030     1    15     0     0    0.0000000", "VALID UNTIL"),
            ("   G01", "START OF FREQUENCY"),
            ("     10.00     20.00     30.00", "NORTH / EAST / UP"),
            ("   NOAZI    0.10    0.20", ""),
            ("", "END OF FREQUENCY"),
            ("   G02", "START OF FREQUENCY"),
            ("     40.00     50.00     60.00", "NORTH / EAST / UP"),
            ("   NOAZI    0.30    0.40", ""),
            ("", "END OF FREQUENCY"),
            ("", "END OF ANTENNA"),
        ];
        let content: String = lines.iter().map(|(d, l)| ant_line(d, l)).collect();
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
        let values: Vec<String> = (0..19).map(|i| format!("{:>8.2}", i as f64 * -0.5)).collect();
        let noazi_line = format!("   NOAZI{}\n", values.join(""));
        let content = [
            ant_line("", "START OF ANTENNA"),
            ant_line(&format!("{:<40}", "LONG_ROW"), "TYPE / SERIAL NO"),
            ant_line("     0.0", "DAZI"),
            ant_line("     0.0  90.0   5.0", "ZEN1 / ZEN2 / DZEN"),
            ant_line("   G01", "START OF FREQUENCY"),
            ant_line("      1.00      2.00      3.00", "NORTH / EAST / UP"),
            noazi_line,
            ant_line("", "END OF FREQUENCY"),
            ant_line("", "END OF ANTENNA"),
        ].concat();

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

    #[test]
    fn test_parse_receiver_azimuth_grid() {
        let lines = [
            ("", "START OF ANTENNA"),
            (&format!("{:<40}", "RECV_WITH_AZI"), "TYPE / SERIAL NO"),
            ("    90.0", "DAZI"),
            ("     0.0  10.0  10.0", "ZEN1 / ZEN2 / DZEN"),
            ("   G01", "START OF FREQUENCY"),
            ("      1.00      2.00      3.00", "NORTH / EAST / UP"),
            ("   NOAZI    1.00    2.00", ""),
            ("     0.0    1.10    2.10", ""),
            ("    90.0    1.20    2.20", ""),
            ("   180.0    1.30    2.30", ""),
            ("   270.0    1.40    2.40", ""),
            ("   360.0    1.10    2.10", ""),
            ("", "END OF FREQUENCY"),
            ("", "END OF ANTENNA"),
        ];
        let content: String = lines.iter().map(|(d, l)| ant_line(d, l)).collect();
        let path = write_temp_antex(&content);
        let db = AntexDatabase::parse(&path).unwrap();
        std::fs::remove_file(&path).ok();

        let freq = &db.antennas[0].frequencies["G01"];
        assert_eq!(freq.noazi, vec![1.0, 2.0]);
        let azi = freq.azi.as_ref().expect("azimuth grid should be parsed when present");
        assert_eq!(azi.len(), 5); // 0, 90, 180, 270, 360
        assert_eq!(azi[0], vec![1.10, 2.10]);
        assert_eq!(azi[1], vec![1.20, 2.20]);
        assert_eq!(azi[4], vec![1.10, 2.10]);
    }

    #[test]
    fn test_satellite_nadir_pcv_interpolation() {
        let content = minimal_antenna_block("TEST_SAT", "G01", "G01");
        let path = write_temp_antex(&content);
        let db = AntexDatabase::parse(&path).unwrap();
        std::fs::remove_file(&path).ok();

        let t = gneiss_core::time::GpsTime::new(2088, 100000.0);
        let pcv_0 = db.satellite_pcv_nadir_m("G01", t, 0.0, "G01");
        assert!((pcv_0 - 0.00010).abs() < 1e-7);
        let pcv_mid = db.satellite_pcv_nadir_m("G01", t, 0.5, "G01");
        assert!((pcv_mid - 0.00015).abs() < 1e-7);
    }
}
