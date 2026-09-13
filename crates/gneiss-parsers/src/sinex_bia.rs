use gneiss_core::obs::ObsCode;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use std::collections::HashMap;
use std::io::BufRead;
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BiasType {
    Osb,
    Dcb,
}

impl FromStr for BiasType {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim() {
            "OSB" => Ok(BiasType::Osb),
            "DCB" => Ok(BiasType::Dcb),
            _ => Err(()),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BiasRecord {
    pub bias_type: BiasType,
    pub sat: SatelliteId,
    pub station: Option<String>,
    pub obs1: ObsCode,
    pub obs2: Option<ObsCode>,
    pub start_time: GpsTime,
    pub end_time: GpsTime,
    pub unit: String,
    pub value: f64,
    pub std_dev: f64,
}

#[derive(Debug, Clone, Default)]
pub struct SinexBias {
    pub records: Vec<BiasRecord>,
    by_sat_obs: HashMap<(SatelliteId, ObsCode), Vec<usize>>,
}

fn doy_to_month_day(y: i32, doy: i32) -> Option<(i32, i32)> {
    let mut dim = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    if y % 4 == 0 && (y % 100 != 0 || y % 400 == 0) { dim[1] = 29; }
    let mut d = doy;
    for (i, &days) in dim.iter().enumerate() {
        if d <= days { return Some((i as i32 + 1, d)); }
        d -= days;
    }
    None
}

fn parse_yds(s: &str) -> Option<GpsTime> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 3 { return None; }
    let year = parts[0].parse::<i32>().ok()?;
    let doy = parts[1].parse::<i32>().ok()?;
    let sec = parts[2].parse::<f64>().ok()?;
    let y = if year < 100 { if year >= 80 { year + 1900 } else { year + 2000 } } else { year };
    let (m, d) = doy_to_month_day(y, doy)?;
    Some(GpsTime::from_calendar(y, m, d, 0, 0, sec))
}

fn parse_sat_and_obs(line: &str) -> Option<(SatelliteId, ObsCode, Option<ObsCode>)> {
    let prn_str = &line[11..14];
    if prn_str.trim().is_empty() { return None; }
    let constell = match prn_str.chars().next()? {
        'G' => Constellation::Gps,
        'R' => Constellation::Glonass,
        'E' => Constellation::Galileo,
        'C' => Constellation::Beidou,
        'J' => Constellation::Qzss,
        _ => return None,
    };
    let prn = prn_str[1..3].trim().parse::<u8>().ok()?;
    let obs1 = ObsCode::from_str(line[25..29].trim()).ok()?;
    let obs2_str = line[30..34].trim();
    let obs2 = if obs2_str.is_empty() { None } else { ObsCode::from_str(obs2_str).ok() };
    Some((SatelliteId { constellation: constell, prn }, obs1, obs2))
}

fn parse_bias_record(line: &str) -> Option<BiasRecord> {
    let bias_type = BiasType::from_str(&line[0..4]).ok()?;
    let (sat, obs1, obs2) = parse_sat_and_obs(line)?;
    let station = match line[15..24].trim() {
        "" => None,
        s => Some(s.to_string()),
    };
    let start_time = parse_yds(line[35..49].trim())?;
    let end_time = parse_yds(line[50..64].trim())?;
    let unit = line[65..69].trim().to_string();
    let value = line[70..91].trim().parse::<f64>().ok()?;
    let std_dev = if line.len() >= 102 {
        line[92..103].trim().parse::<f64>().unwrap_or(0.0)
    } else {
        0.0
    };
    Some(BiasRecord { bias_type, sat, station, obs1, obs2, start_time, end_time, unit, value, std_dev })
}

impl SinexBias {
    pub fn parse<R: BufRead>(mut reader: R) -> Result<Self, String> {
        let mut records = Vec::new();
        let mut in_solution = false;
        let mut buf = Vec::new();
        while let Ok(bytes_read) = reader.read_until(b'\n', &mut buf) {
            if bytes_read == 0 { break; }
            let line = String::from_utf8_lossy(&buf).trim_end().to_string();
            buf.clear();
            if line.starts_with("+BIAS/SOLUTION") {
                in_solution = true;
            } else if line.starts_with("-BIAS/SOLUTION") {
                in_solution = false;
            } else if in_solution && !line.starts_with('*') && line.len() >= 90 {
                if let Some(rec) = parse_bias_record(&line) {
                    records.push(rec);
                }
            }
        }
        Ok(SinexBias::new(records))
    }

    pub fn new(records: Vec<BiasRecord>) -> Self {
        let mut by_sat_obs: HashMap<(SatelliteId, ObsCode), Vec<usize>> = HashMap::new();
        for (i, r) in records.iter().enumerate() {
            by_sat_obs.entry((r.sat, r.obs1)).or_default().push(i);
        }
        Self { records, by_sat_obs }
    }

    pub fn get_bias(&self, sat: SatelliteId, obs: ObsCode, t: GpsTime) -> Option<f64> {
        if let Some(val) = self.get_exact_bias(sat, obs, t) {
            return Some(val);
        }
        for &code_str in fallback_codes(sat.constellation, &obs.to_string()) {
            if let Ok(code) = ObsCode::from_str(code_str) {
                if let Some(val) = self.get_exact_bias(sat, code, t) {
                    return Some(val);
                }
            }
        }
        None
    }

    pub fn get_exact_bias(&self, sat: SatelliteId, obs: ObsCode, t: GpsTime) -> Option<f64> {
        if let Some(indices) = self.by_sat_obs.get(&(sat, obs)) {
            for &idx in indices {
                let rec = &self.records[idx];
                if rec.bias_type == BiasType::Osb && t >= rec.start_time && t <= rec.end_time {
                    return Some(rec.value);
                }
            }
        }
        None
    }

    pub fn lookup_bias(&self, sat: &SatelliteId, obs_code: &str, time: GpsTime) -> Option<f64> {
        let code = ObsCode::from_str(obs_code).ok()?;
        self.get_bias(*sat, code, time)
    }

    pub fn lookup_bias_m(&self, sat: SatelliteId, obs: ObsCode, t: GpsTime) -> Option<f64> {
        let val = self.get_bias(sat, obs, t)?;
        let unit = self.by_sat_obs.get(&(sat, obs))
            .and_then(|idxs| idxs.first())
            .map(|&idx| self.records[idx].unit.as_str())
            .unwrap_or("ns");
        match unit {
            "m" => Some(val),
            _ => Some(val * gneiss_core::constants::SPEED_OF_LIGHT_M_S * 1e-9),
        }
    }

    pub fn wide_lane_satellite_bias(&self, sat: &SatelliteId, time: GpsTime) -> Option<f64> {
        let (f1, f2, p1_s, p2_s, l1_s, l2_s) = sat_bias_bands(sat.constellation);
        let p1 = ObsCode::from_str(p1_s).ok()?;
        let p2 = ObsCode::from_str(p2_s).ok()?;
        let l1 = ObsCode::from_str(l1_s).ok()?;
        let l2 = ObsCode::from_str(l2_s).ok()?;
        self.compute_wl_bias(*sat, f1, f2, p1, p2, l1, l2, time)
    }

    pub fn narrow_lane_satellite_bias(&self, sat: &SatelliteId, time: GpsTime) -> Option<f64> {
        let (f1, f2, _, _, l1_s, l2_s) = sat_bias_bands(sat.constellation);
        let l1 = ObsCode::from_str(l1_s).ok()?;
        let l2 = ObsCode::from_str(l2_s).ok()?;
        self.compute_nl_bias(*sat, f1, f2, l1, l2, time)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn compute_wl_bias(
        &self, sat: SatelliteId, f1: f64, f2: f64,
        p1: ObsCode, p2: ObsCode, l1: ObsCode, l2: ObsCode, t: GpsTime,
    ) -> Option<f64> {
        let d_p1 = self.lookup_bias_m(sat, p1, t)?;
        let d_p2 = self.lookup_bias_m(sat, p2, t)?;
        let d_l1 = self.lookup_bias_m(sat, l1, t)?;
        let d_l2 = self.lookup_bias_m(sat, l2, t)?;
        let wl_lambda = gneiss_core::constants::SPEED_OF_LIGHT_M_S / (f1 - f2);
        let term_phase = (f1 * d_l1 - f2 * d_l2) / (f1 - f2);
        let term_code = (f1 * d_p1 + f2 * d_p2) / (f1 + f2);
        Some((term_phase - term_code) / wl_lambda)
    }

    pub fn compute_nl_bias(
        &self, sat: SatelliteId, f1: f64, f2: f64,
        l1: ObsCode, l2: ObsCode, t: GpsTime,
    ) -> Option<f64> {
        let d_l1 = self.lookup_bias_m(sat, l1, t)?;
        let d_l2 = self.lookup_bias_m(sat, l2, t)?;
        let nl_lambda = gneiss_core::constants::SPEED_OF_LIGHT_M_S / (f1 + f2);
        let if_phase = (f1 * f1 * d_l1 - f2 * f2 * d_l2) / (f1 * f1 - f2 * f2);
        Some(if_phase / nl_lambda)
    }
}

fn sat_bias_bands(c: Constellation) -> (f64, f64, &'static str, &'static str, &'static str, &'static str) {
    match c {
        Constellation::Galileo => (1575.42e6, 1176.45e6, "C1C", "C5Q", "L1C", "L5Q"),
        Constellation::Beidou => (1561.098e6, 1207.14e6, "C2I", "C7I", "L2I", "L7I"),
        Constellation::Qzss => (1575.42e6, 1227.60e6, "C1C", "C2L", "L1C", "L2L"),
        _ => (1575.42e6, 1227.60e6, "C1W", "C2W", "L1W", "L2W"),
    }
}

fn fallback_codes(constellation: Constellation, obs_str: &str) -> &'static [&'static str] {
    match obs_str {
        "C1C" | "C1X" => if constellation == Constellation::Gps { &["C1W"] } else { &["C1C"] },
        "L1C" | "L1X" => if constellation == Constellation::Gps { &["L1W"] } else { &["L1C"] },
        "C2X" | "C2L" | "C2S" => if constellation == Constellation::Gps { &["C2C", "C2W"] } else { &["C2C"] },
        "L2X" | "L2L" | "L2S" => if constellation == Constellation::Gps { &["L2C", "L2W"] } else { &["L2C"] },
        "C5X" | "C5I" => &["C5Q"],
        "L5X" | "L5I" => &["L5Q"],
        "C7X" | "C7I" => &["C7Q"],
        "L7X" | "L7I" => &["L7Q"],
        "C8X" | "C8I" => &["C8Q"],
        "L8X" | "L8I" => &["L8Q"],
        _ => &[],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_sinex_bias_parsing() {
        let content = r#"%=BIA 1.00
+BIAS/SOLUTION
*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___
 OSB  G002 G02           C1C       2021:123:00000 2021:123:86400 ns           -2.5186000000    0.000000
 DCB  G002 G02           C1W  C2W  2021:123:00000 2021:123:86400 ns            1.2345678901    0.000000
-BIAS/SOLUTION
"#;

        let cursor = Cursor::new(content);
        let bias = SinexBias::parse(cursor).unwrap();

        assert_eq!(bias.records.len(), 2);

        let r1 = &bias.records[0];
        assert_eq!(r1.bias_type, BiasType::Osb);
        assert_eq!(r1.sat.constellation, Constellation::Gps);
        assert_eq!(r1.sat.prn, 2);
        assert_eq!(r1.obs1, ObsCode::from_str("C1C").unwrap());
        assert_eq!(r1.obs2, None);
        assert_eq!(r1.unit, "ns");
        assert_eq!(r1.value, -2.5186);

        let r2 = &bias.records[1];
        assert_eq!(r2.bias_type, BiasType::Dcb);
        assert_eq!(r2.obs2, Some(ObsCode::from_str("C2W").unwrap()));
        assert_eq!(r2.value, 1.2345678901);

        // test get_bias
        let t = parse_yds("2021:123:43200").unwrap();
        let val = bias.get_bias(
            SatelliteId {
                constellation: Constellation::Gps,
                prn: 2,
            },
            ObsCode::from_str("C1C").unwrap(),
            t,
        );
        assert_eq!(val, Some(-2.5186));
    }
}
#[cfg(test)]
mod fallback_tests {
    use super::*;
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;
    use std::io::Cursor;
    use std::str::FromStr;

    #[test]
    fn test_sinex_bias_fallback() {
        let content = r#"%=BIA 1.00
+BIAS/SOLUTION
*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___
 OSB  G002 G02           C1W       2021:123:00000 2021:123:86400 ns            1.0000000000    0.000000
 OSB  G002 G02           L1W       2021:123:00000 2021:123:86400 ns            2.0000000000    0.000000
 OSB  G002 G02           C2W       2021:123:00000 2021:123:86400 ns            3.0000000000    0.000000
 OSB  G002 G02           L2W       2021:123:00000 2021:123:86400 ns            4.0000000000    0.000000
 OSB  E002 E02           C2C       2021:123:00000 2021:123:86400 ns            5.0000000000    0.000000
 OSB  E002 E02           L2C       2021:123:00000 2021:123:86400 ns            6.0000000000    0.000000
-BIAS/SOLUTION
"#;
        let cursor = Cursor::new(content);
        let bias = SinexBias::parse(cursor).unwrap();
        let t = GpsTime::new(2156, 129600.0); // 2021:123 at 12:00:00

        let g02 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };
        let e02 = SatelliteId {
            constellation: Constellation::Galileo,
            prn: 2,
        };

        // Test GPS Fallbacks
        assert_eq!(
            bias.get_bias(g02, ObsCode::from_str("C1C").unwrap(), t),
            Some(1.0)
        );
        assert_eq!(
            bias.get_bias(g02, ObsCode::from_str("L1C").unwrap(), t),
            Some(2.0)
        );
        assert_eq!(
            bias.get_bias(g02, ObsCode::from_str("C2L").unwrap(), t),
            Some(3.0)
        );
        assert_eq!(
            bias.get_bias(g02, ObsCode::from_str("L2X").unwrap(), t),
            Some(4.0)
        );

        // Test non-GPS Fallbacks
        assert_eq!(
            bias.get_bias(e02, ObsCode::from_str("C2X").unwrap(), t),
            Some(5.0)
        );
        assert_eq!(
            bias.get_bias(e02, ObsCode::from_str("L2S").unwrap(), t),
            Some(6.0)
        );

        // Ensure no fallback if direct is present
        assert_eq!(
            bias.get_exact_bias(g02, ObsCode::from_str("C1C").unwrap(), t),
            None
        );

        // Test lookup_bias
        assert_eq!(bias.lookup_bias(&g02, "C1W", t), Some(1.0));
        assert_eq!(bias.lookup_bias(&g02, "INVALID", t), None);

        // Test wide_lane_satellite_bias and narrow_lane_satellite_bias
        let wl = bias.wide_lane_satellite_bias(&g02, t);
        assert!(wl.is_some());
        let nl = bias.narrow_lane_satellite_bias(&g02, t);
        assert!(nl.is_some());
    }
}
