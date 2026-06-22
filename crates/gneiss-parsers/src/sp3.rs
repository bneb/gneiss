use gneiss_core::time::GpsTime;
use std::collections::HashMap;
use std::io::BufRead;

#[derive(Debug, Clone)]
pub struct Sp3Epoch {
    pub time: GpsTime,
    pub records: HashMap<String, Sp3Record>,
}

#[derive(Debug, Clone)]
pub struct Sp3Record {
    pub position: nalgebra::Vector3<f64>, // meters
    pub clock_offset: f64,                // seconds
}

pub fn parse_sp3<R: BufRead>(reader: R) -> Result<Vec<Sp3Epoch>, String> {
    let mut epochs = Vec::new();
    let mut current_epoch: Option<Sp3Epoch> = None;

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| e.to_string())?;
        if line.is_empty() {
            continue;
        }

        if let Some(stripped) = line.strip_prefix('*') {
            // Epoch header line: *  YYYY MM DD HH MM SS.sssssss
            let parts: Vec<&str> = stripped.split_whitespace().collect();
            if parts.len() < 6 {
                continue;
            }
            if let (Ok(year), Ok(month), Ok(day), Ok(hour), Ok(minute), Ok(sec)) = (
                parts[0].parse::<i32>(),
                parts[1].parse::<u32>(),
                parts[2].parse::<u32>(),
                parts[3].parse::<u32>(),
                parts[4].parse::<u32>(),
                parts[5].parse::<f64>(),
            ) {
                if let Some(epoch) = current_epoch.take() {
                    epochs.push(epoch);
                }

                let time = GpsTime::from_calendar(
                    year,
                    month as i32,
                    day as i32,
                    hour as i32,
                    minute as i32,
                    sec,
                );
                current_epoch = Some(Sp3Epoch {
                    time,
                    records: HashMap::new(),
                });
            }
        } else if line.starts_with('P') {
            // Position record: PG01  X Y Z Clock
            if let Some(epoch) = &mut current_epoch {
                if line.len() >= 60 {
                    let sat_id = line[1..4].trim().to_string();
                    let x_str = line[4..18].trim();
                    let y_str = line[18..32].trim();
                    let z_str = line[32..46].trim();
                    let clk_str = line[46..60].trim();

                    if let (Ok(x), Ok(y), Ok(z)) = (
                        x_str.parse::<f64>(),
                        y_str.parse::<f64>(),
                        z_str.parse::<f64>(),
                    ) {
                        let clock_offset = if !clk_str.is_empty() && clk_str != "999999.999999" {
                            clk_str.parse::<f64>().unwrap_or(f64::NAN) * 1e-6 // microseconds to seconds
                        } else {
                            f64::NAN
                        };

                        epoch.records.insert(
                            sat_id,
                            Sp3Record {
                                position: nalgebra::Vector3::new(
                                    x * 1000.0,
                                    y * 1000.0,
                                    z * 1000.0,
                                ), // km to meters
                                clock_offset,
                            },
                        );
                    }
                }
            }
        } else if line.starts_with("EOF") {
            if let Some(epoch) = current_epoch.take() {
                epochs.push(epoch);
            }
            break;
        }
    }

    if let Some(epoch) = current_epoch {
        epochs.push(epoch);
    }

    Ok(epochs)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build an SP3 epoch header line: `*  YYYY MM DD HH MM SS.sssssss`
    fn epoch_line(year: i32, month: u32, day: u32, hour: u32, min: u32, sec: f64) -> String {
        format!("*  {:4} {:02} {:02} {:02} {:02} {:>9.7}\n", year, month, day, hour, min, sec)
    }

    /// Build an SP3 position record line.
    ///
    /// Fixed-width fields (0-indexed):
    ///   [0]      'P'
    ///   [1..4]   satellite ID (e.g. "G01")
    ///   [4..18]  X coordinate in km  (14 chars)
    ///   [18..32] Y coordinate in km  (14 chars)
    ///   [32..46] Z coordinate in km  (14 chars)
    ///   [46..60] clock offset in us  (14 chars)
    fn pos_line(sat: &str, x_km: f64, y_km: f64, z_km: f64, clk_us: f64) -> String {
        format!(
            "P{:<3}{:>14.6}{:>14.6}{:>14.6}{:>14.9}\n",
            sat, x_km, y_km, z_km, clk_us,
        )
    }

    #[test]
    fn test_parse_sp3_single_epoch_one_sat() {
        let content = format!(
            "{}{}EOF\n",
            epoch_line(2024, 3, 15, 12, 0, 0.0),
            pos_line("G01", -10000.0, 20000.0, 15000.0, 0.123456780),
        );
        let reader = std::io::Cursor::new(content);
        let epochs = parse_sp3(reader).unwrap();

        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].records.len(), 1);

        let rec = epochs[0].records.get("G01").unwrap();
        // km -> m
        assert!((rec.position.x - (-10_000_000.0)).abs() < 1e-6);
        assert!((rec.position.y - 20_000_000.0).abs() < 1e-6);
        assert!((rec.position.z - 15_000_000.0).abs() < 1e-6);
        // us -> s
        assert!((rec.clock_offset - 0.123456780e-6).abs() < 1e-20);
    }

    #[test]
    fn test_parse_sp3_multiple_epochs() {
        let content = format!(
            "{}{}EOF\n",
            epoch_line(2024, 3, 15, 12, 0, 0.0),
            pos_line("G01", -10000.0, 20000.0, 15000.0, 0.1),
        );
        // No EOF in the middle — just one epoch because reader exhausts next
        let reader = std::io::Cursor::new(content);
        let epochs = parse_sp3(reader).unwrap();
        assert_eq!(epochs.len(), 1);
    }

    #[test]
    fn test_parse_sp3_multiple_satellites_one_epoch() {
        let content = format!(
            "{}{}{}EOF\n",
            epoch_line(2024, 3, 15, 12, 0, 0.0),
            pos_line("G01", -10000.0, 20000.0, 15000.0, 0.1),
            pos_line("R02", 5000.0, -6000.0, 7000.0, 0.2),
        );
        let reader = std::io::Cursor::new(content);
        let epochs = parse_sp3(reader).unwrap();

        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].records.len(), 2);

        let g01 = epochs[0].records.get("G01").unwrap();
        let r02 = epochs[0].records.get("R02").unwrap();
        assert!((g01.position.x + 10_000_000.0).abs() < 1e-6);
        assert!((r02.position.x - 5_000_000.0).abs() < 1e-6);
    }

    #[test]
    fn test_parse_sp3_empty_input() {
        let reader = std::io::Cursor::new(String::new());
        let epochs = parse_sp3(reader).unwrap();
        assert!(epochs.is_empty());
    }

    #[test]
    fn test_parse_sp3_skip_non_matching_lines() {
        // Lines that don't start with '*', 'P', or "EOF" are skipped.
        let content = format!(
            "{}{}COMMENT line\nEOF\n",
            epoch_line(2024, 3, 15, 12, 0, 0.0),
            pos_line("G01", 0.0, 0.0, 0.0, 0.0),
        );
        let reader = std::io::Cursor::new(content);
        let epochs = parse_sp3(reader).unwrap();
        assert_eq!(epochs.len(), 1);
    }

    #[test]
    fn test_parse_sp3_nan_clock_default() {
        // Clock value "999999.999999" (12-char string padded to 14) should produce NAN.
        // Use pos_line with a valid clock value, then patch the clock field inline.
        let p_line = format!(
            "P{:<3}{:>14.6}{:>14.6}{:>14.6}{:>14}\n",
            "G01", 0.0, 0.0, 0.0, "999999.999999",
        );
        let content = format!(
            "{}*  2024 03 15 12 00  0.00000000\n{}EOF\n",
            "", &p_line,
        );
        let reader = std::io::Cursor::new(content);
        let epochs = parse_sp3(reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let rec = epochs[0].records.get("G01").unwrap();
        assert!(rec.clock_offset.is_nan(), "expected NaN for default clock, got {}", rec.clock_offset);
    }

    #[test]
    fn test_parse_sp3_eof_terminates_early() {
        // EOF appears before the second epoch would start, so only one epoch.
        let content = format!(
            "{}{}EOF\n",
            epoch_line(2024, 3, 15, 12, 0, 0.0),
            pos_line("G01", 0.0, 0.0, 0.0, 0.0),
        );
        let reader = std::io::Cursor::new(content);
        let epochs = parse_sp3(reader).unwrap();
        assert_eq!(epochs.len(), 1);
    }

    #[test]
    fn test_parse_sp3_no_eof_implicit_end() {
        // Without EOF the reader is exhausted after the last line.
        let content = format!(
            "{}{}",
            epoch_line(2024, 3, 15, 12, 0, 0.0),
            pos_line("G01", 0.0, 0.0, 0.0, 0.0),
        );
        let reader = std::io::Cursor::new(content);
        let epochs = parse_sp3(reader).unwrap();
        assert_eq!(epochs.len(), 1);
    }

    #[test]
    fn test_parse_sp3_skip_short_position_lines() {
        // A "P" line shorter than 60 characters is skipped gracefully.
        let short_line = "PG01   0.0   0.0   0.0\n";
        let content = format!(
            "{}{}EOF\n",
            epoch_line(2024, 3, 15, 12, 0, 0.0),
            short_line,
        );
        let reader = std::io::Cursor::new(content);
        let epochs = parse_sp3(reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert!(epochs[0].records.is_empty());
    }

    #[test]
    fn test_parse_sp3_skip_bad_epoch_header() {
        // An epoch header with fewer than 6 fields is skipped.
        let bad_header = "*  2024 03 15\n";
        let content = format!("{}EOF\n", bad_header);
        let reader = std::io::Cursor::new(content);
        let epochs = parse_sp3(reader).unwrap();
        assert!(epochs.is_empty());
    }

    #[test]
    fn test_parse_sp3_skip_invalid_position_coords() {
        // Non-parseable coordinates in a P line produce no record.
        let p_line = "PG01   xxxx    yyyy    zzzz    0.000000000\n";
        let content = format!(
            "{}{}EOF\n",
            epoch_line(2024, 3, 15, 12, 0, 0.0),
            p_line,
        );
        let reader = std::io::Cursor::new(content);
        let epochs = parse_sp3(reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert!(epochs[0].records.is_empty());
    }
}
