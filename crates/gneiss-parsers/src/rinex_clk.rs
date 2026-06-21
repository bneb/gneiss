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
}
