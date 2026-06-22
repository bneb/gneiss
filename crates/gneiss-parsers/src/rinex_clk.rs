use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct ClockRecord {
    pub time: GpsTime,
    pub bias: f64, // in seconds
}

#[derive(Debug, Clone, Default)]
pub struct RinexClock {
    pub satellites: HashMap<SatelliteId, Vec<ClockRecord>>,
}

impl RinexClock {
    pub fn parse(content: &str) -> Self {
        let mut clk = RinexClock::default();

        for line in content.lines() {
            if line.starts_with("AS ") {
                let sat_str = &line[3..6];
                let constell = match sat_str.chars().next() {
                    Some('G') => Constellation::Gps,
                    Some('R') => Constellation::Glonass,
                    Some('E') => Constellation::Galileo,
                    Some('C') => Constellation::Beidou,
                    Some('J') => Constellation::Qzss,
                    _ => continue,
                };

                let prn = match sat_str[1..3].trim().parse::<u8>() {
                    Ok(p) => p,
                    Err(_) => continue,
                };

                let sat = SatelliteId {
                    constellation: constell,
                    prn,
                };

                // Parse time
                let year = line[8..12].trim().parse::<i32>().unwrap_or(0);
                let month = line[13..15].trim().parse::<i32>().unwrap_or(0);
                let day = line[16..18].trim().parse::<i32>().unwrap_or(0);
                let hour = line[19..21].trim().parse::<i32>().unwrap_or(0);
                let minute = line[22..24].trim().parse::<i32>().unwrap_or(0);
                let second = line[25..34].trim().parse::<f64>().unwrap_or(0.0);

                let time = GpsTime::from_calendar(year, month, day, hour, minute, second);

                // Parse number of values
                // Bias is the first value
                let bias_str = line[40..59].replace("D", "e");
                let bias = match bias_str.trim().parse::<f64>() {
                    Ok(b) => b,
                    Err(_) => continue,
                };

                clk.satellites
                    .entry(sat)
                    .or_default()
                    .push(ClockRecord { time, bias });
            }
        }

        // Sort records by time just in case
        for records in clk.satellites.values_mut() {
            records.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap());
        }

        clk
    }

    pub fn get_clock_bias(&self, sat: SatelliteId, t: GpsTime) -> Option<f64> {
        let records = self.satellites.get(&sat)?;
        if records.is_empty() {
            return None;
        }

        // Binary search for nearest or bounding interval
        let idx = match records.binary_search_by(|r| r.time.partial_cmp(&t).unwrap()) {
            Ok(i) => return Some(records[i].bias),
            Err(i) => i,
        };

        if idx == 0 {
            // Check if too far
            if (records[0].time - t).abs() > 900.0 {
                return None;
            }
            return Some(records[0].bias);
        }
        if idx >= records.len() {
            let last = records.len() - 1;
            if (t - records[last].time).abs() > 900.0 {
                return None;
            }
            return Some(records[last].bias);
        }

        let r1 = &records[idx - 1];
        let r2 = &records[idx];

        let dt = r2.time - r1.time;
        if dt == 0.0 || (t - r1.time).abs() > 900.0 || (r2.time - t).abs() > 900.0 {
            return None;
        }

        // Linear interpolation
        let bias = r1.bias + (r2.bias - r1.bias) * (t - r1.time) / dt;
        Some(bias)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::sat::Constellation;

    fn make_record(tow: f64, bias: f64) -> ClockRecord {
        ClockRecord {
            time: GpsTime::new(2000, tow),
            bias,
        }
    }

    /// Bug 24 regression test: a gap > 900 s between records must return None,
    /// not the stale r1 bias.
    #[test]
    fn test_precise_clock_gap_returns_none() {
        let mut clk = RinexClock::default();
        let sat = SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 1,
        };
        // Two records with a 1800 s gap
        let records = vec![make_record(0.0, 1.0e-7), make_record(1800.0, 2.0e-7)];
        clk.satellites.insert(sat, records);

        // Query at t = 950 s, which is 950 s past r1 (> 900 s threshold)
        let t = GpsTime::new(2000, 950.0);
        let result = clk.get_clock_bias(sat, t);
        assert!(
            result.is_none(),
            "Expected None for a clock gap > 900 s, got {:?}",
            result
        );
    }

    /// Sanity check: records within tolerance should still interpolate.
    #[test]
    fn test_precise_clock_within_tolerance_interpolates() {
        let mut clk = RinexClock::default();
        let sat = SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 2,
        };
        let records = vec![make_record(0.0, 0.0), make_record(30.0, 3.0e-9)];
        clk.satellites.insert(sat, records);

        // Midpoint should interpolate to 1.5e-9
        let t = GpsTime::new(2000, 15.0);
        let result = clk.get_clock_bias(sat, t);
        assert!(
            result.is_some(),
            "Expected Some for records within tolerance"
        );
        assert!(
            (result.unwrap() - 1.5e-9).abs() < 1e-18,
            "Expected 1.5e-9, got {:?}",
            result
        );
    }

    #[test]
    fn test_precise_clock_gap_tolerance() {
        let mut clk = RinexClock::default();
        let sat = SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 1,
        };
        // Setup records:
        // - Record 0: t = 0.0, bias = 1.0e-7
        // - Record 1: t = 1800.0, bias = 2.0e-7  (gap is 1800s > 900s)
        // - Record 2: t = 2400.0, bias = 3.0e-7  (gap is 600s <= 900s)
        let records = vec![
            make_record(0.0, 1.0e-7),
            make_record(1800.0, 2.0e-7),
            make_record(2400.0, 3.0e-7),
        ];
        clk.satellites.insert(sat, records);

        // An exact match at a timestamp (which is preceded by a gap > 900s)
        // returns the correct bias immediately and doesn't trigger the gap check.
        let result_exact = clk.get_clock_bias(sat, GpsTime::new(2000, 1800.0));
        assert_eq!(result_exact, Some(2.0e-7));

        // Non-exact matches within gaps > 900s return None.
        // Gap between record 0 (t=0.0) and record 1 (t=1800.0) is 1800s (> 900s).
        // t = 901.0 is 901.0s away from r1 (> 900s).
        let result_gap_none_1 = clk.get_clock_bias(sat, GpsTime::new(2000, 901.0));
        assert!(result_gap_none_1.is_none());

        // t = 899.0 is 901.0s away from r2 (> 900s).
        let result_gap_none_2 = clk.get_clock_bias(sat, GpsTime::new(2000, 899.0));
        assert!(result_gap_none_2.is_none());

        // Non-exact matches within gaps <= 900s correctly interpolate.
        // Gap between record 1 (t=1800.0) and record 2 (t=2400.0) is 600s (<= 900s).
        // Midpoint t = 2100.0 is 300.0s away from both r1 and r2.
        let result_interpolate = clk.get_clock_bias(sat, GpsTime::new(2000, 2100.0));
        assert!(result_interpolate.is_some());
        assert!((result_interpolate.unwrap() - 2.5e-7).abs() < 1e-18);

        // Extrapolations beyond bounds (> 900s) return None.
        // Before first record (t = -901.0)
        let result_extrap_before = clk.get_clock_bias(sat, GpsTime::new(2000, -901.0));
        assert!(result_extrap_before.is_none());

        // After last record (t = 3301.0)
        let result_extrap_after = clk.get_clock_bias(sat, GpsTime::new(2000, 3301.0));
        assert!(result_extrap_after.is_none());
    }

    // ------------------------------------------------------------------
    // Tests for RinexClock::parse (AS-format clock file parsing)
    // ------------------------------------------------------------------

    /// Helper: build a RINEX CLK "AS" line.
    ///
    /// Fixed-width column positions (0-indexed):
    ///   [0..2]   "AS"
    ///   [3..6]   satellite PRN (e.g. "G01")
    ///   [8..12]  year
    ///   [13..15] month
    ///   [16..18] day
    ///   [19..21] hour
    ///   [22..24] minute
    ///   [25..34] seconds  (9 chars, right-aligned)
    ///   [40..59] bias     (19 chars, D-exponent notation)
    fn as_line(sat: &str, year: i32, month: i32, day: i32, hour: i32, min: i32, sec: f64, bias: f64) -> String {
        let sec_fmt = format!("{:>9.6}", sec);
        // Bias with D exponent notation, padded to 19 chars.
        let bias_str = format!("{:.10e}", bias).replace('e', "D");
        format!(
            "AS {:<3}  {:4} {:02} {:02} {:02} {:02} {:<9}      {:>19}\n",
            sat, year, month, day, hour, min, sec_fmt, bias_str,
        )
    }

    #[test]
    fn test_parse_single_valid_as_line() {
        let line = as_line("G01", 2024, 3, 15, 12, 0, 0.0, 1.23456789e-7);
        let clk = RinexClock::parse(&line);

        assert_eq!(clk.satellites.len(), 1);
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let records = clk.satellites.get(&sat).unwrap();
        assert_eq!(records.len(), 1);
        assert!(
            (records[0].bias - 1.23456789e-7).abs() < 1e-18,
            "bias mismatch: got {}",
            records[0].bias
        );
    }

    #[test]
    fn test_parse_multiple_epochs_same_satellite() {
        let content = format!(
            "{}{}",
            as_line("G01", 2024, 3, 15, 12, 0, 0.0, 1.0e-7),
            as_line("G01", 2024, 3, 15, 12, 15, 0.0, 2.0e-7),
        );
        let clk = RinexClock::parse(&content);

        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let records = clk.satellites.get(&sat).unwrap();
        assert_eq!(records.len(), 2);
        assert!((records[0].bias - 1.0e-7).abs() < 1e-18);
        assert!((records[1].bias - 2.0e-7).abs() < 1e-18);
    }

    #[test]
    fn test_parse_all_constellations() {
        let content = format!(
            "{}{}{}{}{}",
            as_line("G01", 2024, 3, 15, 12, 0, 0.0, 1.0e-7),
            as_line("R02", 2024, 3, 15, 12, 0, 0.0, 2.0e-7),
            as_line("E03", 2024, 3, 15, 12, 0, 0.0, 3.0e-7),
            as_line("C04", 2024, 3, 15, 12, 0, 0.0, 4.0e-7),
            as_line("J05", 2024, 3, 15, 12, 0, 0.0, 5.0e-7),
        );
        let clk = RinexClock::parse(&content);
        assert_eq!(clk.satellites.len(), 5);

        let cases = [
            (Constellation::Gps, 1, 1.0e-7),
            (Constellation::Glonass, 2, 2.0e-7),
            (Constellation::Galileo, 3, 3.0e-7),
            (Constellation::Beidou, 4, 4.0e-7),
            (Constellation::Qzss, 5, 5.0e-7),
        ];
        for (constell, prn, expected_bias) in cases {
            let sat = SatelliteId { constellation: constell, prn };
            let recs = clk.satellites.get(&sat).unwrap_or_else(|| panic!("missing {constell:?} PRN {prn}"));
            assert!(
                (recs[0].bias - expected_bias).abs() < 1e-18,
                "bias mismatch for {constell:?} PRN {prn}: got {}",
                recs[0].bias
            );
        }
    }

    #[test]
    fn test_parse_empty_content() {
        let clk = RinexClock::parse("");
        assert!(clk.satellites.is_empty());
    }

    #[test]
    fn test_parse_non_as_lines_yield_empty() {
        let content = "COMMENT line\nanother line\n";
        let clk = RinexClock::parse(content);
        assert!(clk.satellites.is_empty());
    }

    #[test]
    fn test_parse_unknown_constellation_skipped() {
        // 'S' is not in the recognised set -> continue
        let line = as_line("S01", 2024, 3, 15, 12, 0, 0.0, 1.0e-7);
        let clk = RinexClock::parse(&line);
        assert!(clk.satellites.is_empty());
    }

    #[test]
    fn test_parse_invalid_prn_skipped() {
        // Non-numeric PRN portion ("G  " -> "  " after trimming) -> continue
        let line = "AS G   2024 03 15 12 00  0.000000      0.100000000D-06     \n";
        let clk = RinexClock::parse(line);
        assert!(clk.satellites.is_empty());
    }

    #[test]
    fn test_parse_invalid_bias_skipped() {
        // Build a line via as_line (well-padded) then overwrite the bias
        // portion with non-numeric characters so parsing fails.
        let line = as_line("G01", 2024, 3, 15, 12, 0, 0.0, 0.0);
        // The bias safely occupies [40..59] in a well-padded line.
        // Replace that span with 'x' characters.
        let mut bytes: Vec<u8> = line.into_bytes();
        for i in 40..59 {
            bytes[i] = b'x';
        }
        let content = String::from_utf8(bytes).unwrap();
        let clk = RinexClock::parse(&content);
        assert!(clk.satellites.is_empty());
    }
}
