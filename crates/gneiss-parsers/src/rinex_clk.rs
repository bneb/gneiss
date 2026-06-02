use std::collections::HashMap;
use gneiss_core::time::GpsTime;
use gneiss_core::sat::{SatelliteId, Constellation};

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
                
                let sat = SatelliteId { constellation: constell, prn };
                
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
                
                clk.satellites.entry(sat).or_default().push(ClockRecord {
                    time,
                    bias,
                });
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
        if records.is_empty() { return None; }
        
        // Binary search for nearest or bounding interval
        let idx = match records.binary_search_by(|r| r.time.partial_cmp(&t).unwrap()) {
            Ok(i) => i,
            Err(i) => i,
        };
        
        if idx == 0 {
            // Check if too far
            if (records[0].time - t).abs() > 60.0 { return None; }
            return Some(records[0].bias);
        }
        if idx >= records.len() {
            let last = records.len() - 1;
            if (t - records[last].time).abs() > 60.0 { return None; }
            return Some(records[last].bias);
        }
        
        let r1 = &records[idx - 1];
        let r2 = &records[idx];
        
        let dt = r2.time - r1.time;
        if dt == 0.0 || (t - r1.time).abs() > 60.0 {
            return Some(r1.bias);
        }
        
        // Linear interpolation
        let bias = r1.bias + (r2.bias - r1.bias) * (t - r1.time) / dt;
        Some(bias)
    }
}
