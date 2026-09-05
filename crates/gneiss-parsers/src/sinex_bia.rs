use gneiss_core::obs::ObsCode;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
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
}

fn parse_yds(s: &str) -> Option<GpsTime> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 3 {
        return None;
    }
    let year = parts[0].parse::<i32>().ok()?;
    let doy = parts[1].parse::<i32>().ok()?;
    let sec = parts[2].parse::<f64>().ok()?;

    let y = if year < 100 {
        if year >= 80 {
            year + 1900
        } else {
            year + 2000
        }
    } else {
        year
    };

    let mut days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    if y % 4 == 0 && (y % 100 != 0 || y % 400 == 0) {
        days_in_month[1] = 29;
    }

    let mut d = doy;
    let mut m = 0;
    for (i, &dim) in days_in_month.iter().enumerate() {
        if d <= dim {
            m = i as i32 + 1;
            break;
        }
        d -= dim;
    }

    if m == 0 {
        return None;
    }

    Some(GpsTime::from_calendar(y, m, d, 0, 0, sec))
}

impl SinexBias {
    pub fn parse<R: BufRead>(mut reader: R) -> Result<Self, String> {
        let mut records = Vec::new();
        let mut in_solution = false;

        let mut buf = Vec::new();
        while let Ok(bytes_read) = reader.read_until(b'\n', &mut buf) {
            if bytes_read == 0 {
                break;
            }
            let line = String::from_utf8_lossy(&buf).trim_end().to_string();
            buf.clear();

            if line.starts_with("+BIAS/SOLUTION") {
                in_solution = true;
                continue;
            } else if line.starts_with("-BIAS/SOLUTION") {
                in_solution = false;
                continue;
            }

            if in_solution && !line.starts_with('*') && line.len() >= 90 {
                let bias_type_str = &line[0..4];
                let bias_type = match BiasType::from_str(bias_type_str) {
                    Ok(b) => b,
                    Err(_) => continue,
                };

                let prn_str = &line[11..14];
                if prn_str.trim().is_empty() {
                    // Could be a station-only bias
                    continue;
                }

                let constell = match prn_str.chars().next() {
                    Some('G') => Constellation::Gps,
                    Some('R') => Constellation::Glonass,
                    Some('E') => Constellation::Galileo,
                    Some('C') => Constellation::Beidou,
                    Some('J') => Constellation::Qzss,
                    _ => continue,
                };

                let prn = match prn_str[1..3].trim().parse::<u8>() {
                    Ok(p) => p,
                    Err(_) => continue,
                };
                let sat = SatelliteId {
                    constellation: constell,
                    prn,
                };

                let station_str = line[15..24].trim();
                let station = if station_str.is_empty() {
                    None
                } else {
                    Some(station_str.to_string())
                };

                let obs1_str = line[25..29].trim();
                let obs1 = match ObsCode::from_str(obs1_str) {
                    Ok(o) => o,
                    Err(_) => continue,
                };

                let obs2_str = line[30..34].trim();
                let obs2 = if obs2_str.is_empty() {
                    None
                } else {
                    ObsCode::from_str(obs2_str).ok()
                };

                let start_time_str = line[35..49].trim();
                let start_time = match parse_yds(start_time_str) {
                    Some(t) => t,
                    None => continue,
                };

                let end_time_str = line[50..64].trim();
                let end_time = match parse_yds(end_time_str) {
                    Some(t) => t,
                    None => continue,
                };

                let unit = line[65..69].trim().to_string();
                let value = match line[70..91].trim().parse::<f64>() {
                    Ok(v) => v,
                    Err(_) => continue,
                };

                let std_dev = if line.len() >= 102 {
                    line[92..103].trim().parse::<f64>().unwrap_or(0.0)
                } else {
                    0.0
                };

                records.push(BiasRecord {
                    bias_type,
                    sat,
                    station,
                    obs1,
                    obs2,
                    start_time,
                    end_time,
                    unit,
                    value,
                    std_dev,
                });
            }
        }

        Ok(SinexBias { records })
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

impl SinexBias {

    pub fn get_exact_bias(&self, sat: SatelliteId, obs: ObsCode, t: GpsTime) -> Option<f64> {
        for rec in &self.records {
            if rec.bias_type == BiasType::Osb && rec.sat == sat && rec.obs1 == obs
                && t >= rec.start_time && t <= rec.end_time {
                    return Some(rec.value);
                }
        }
        None
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
    }
}
