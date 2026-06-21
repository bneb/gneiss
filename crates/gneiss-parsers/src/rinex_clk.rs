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
            Ok(i) => i,
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
        // Bug 24 fix: gap too large → stale clock, return None rather than
        // propagating a drifting bias into the EKF.
        if dt == 0.0 || (t - r1.time).abs() > 900.0 {
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
        let records = vec![
            make_record(0.0, 1.0e-7),
            make_record(1800.0, 2.0e-7),
        ];
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
        let records = vec![
            make_record(0.0, 0.0),
            make_record(30.0, 3.0e-9),
        ];
        clk.satellites.insert(sat, records);

        // Midpoint should interpolate to 1.5e-9
        let t = GpsTime::new(2000, 15.0);
        let result = clk.get_clock_bias(sat, t);
        assert!(result.is_some(), "Expected Some for records within tolerance");
        assert!(
            (result.unwrap() - 1.5e-9).abs() < 1e-18,
            "Expected 1.5e-9, got {:?}",
            result
        );
    }
}
