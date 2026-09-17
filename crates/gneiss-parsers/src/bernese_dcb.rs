use std::io::BufRead;
use std::str::FromStr;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use gneiss_core::obs::ObsCode;
use crate::sinex_bia::{BiasRecord, BiasType};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DcbDiffType {
    P1C1,
    P2C2,
    P1P2,
    P1C2,
}

pub fn parse_dcb_diff_type(line: &str) -> Option<DcbDiffType> {
    if line.contains("(P1-C1)") {
        Some(DcbDiffType::P1C1)
    } else if line.contains("(P2-C2)") {
        Some(DcbDiffType::P2C2)
    } else if line.contains("(P1-P2)") {
        Some(DcbDiffType::P1P2)
    } else if line.contains("(P1-C2)") {
        Some(DcbDiffType::P1C2)
    } else {
        None
    }
}

pub fn parse_sat_id(s: &str) -> Option<SatelliteId> {
    let s = s.trim();
    if s.len() < 2 {
        return None;
    }
    let constellation = match s.chars().next()? {
        'G' => Constellation::Gps,
        'R' => Constellation::Glonass,
        'E' => Constellation::Galileo,
        'C' => Constellation::Beidou,
        'J' => Constellation::Qzss,
        _ => return None,
    };
    let prn = s[1..].parse::<u8>().ok()?;
    Some(SatelliteId { constellation, prn })
}

fn parse_ymd_hms(s: &str) -> Option<GpsTime> {
    let parts: Vec<&str> = s.split_whitespace().collect();
    if parts.len() < 6 {
        return None;
    }
    let y = parts[0].parse::<i32>().ok()?;
    let m = parts[1].parse::<i32>().ok()?;
    let d = parts[2].parse::<i32>().ok()?;
    let hh = parts[3].parse::<i32>().ok()?;
    let mm = parts[4].parse::<i32>().ok()?;
    let ss = parts[5].parse::<f64>().ok()?;
    Some(GpsTime::from_calendar(y, m, d, hh, mm, ss))
}

fn map_diff_type_to_bias(
    diff_type: DcbDiffType,
    val_ns: f64,
) -> Option<(BiasType, ObsCode, Option<ObsCode>, f64)> {
    match diff_type {
        DcbDiffType::P1C1 => {
            let obs1 = ObsCode::from_str("C1C").ok()?;
            Some((BiasType::Osb, obs1, None, -val_ns))
        }
        DcbDiffType::P2C2 => {
            let obs1 = ObsCode::from_str("C2C").ok()?;
            Some((BiasType::Osb, obs1, None, -val_ns))
        }
        DcbDiffType::P1P2 => {
            let obs1 = ObsCode::from_str("C1W").ok()?;
            let obs2 = ObsCode::from_str("C2W").ok()?;
            Some((BiasType::Dcb, obs1, Some(obs2), val_ns))
        }
        DcbDiffType::P1C2 => {
            let obs1 = ObsCode::from_str("C2W").ok()?;
            Some((BiasType::Osb, obs1, None, -val_ns))
        }
    }
}

pub fn parse_dcb_sat_line(
    line: &str,
    diff_type: DcbDiffType,
    default_span: (GpsTime, GpsTime),
) -> Option<BiasRecord> {
    if line.len() < 35 {
        return None;
    }
    let station_field = line.get(3..20)?;
    if !station_field.trim().is_empty() {
        return None;
    }
    let sat_field = line.get(0..3)?;
    let sat = parse_sat_id(sat_field)?;
    let val_ns = line.get(20..35)?.trim().parse::<f64>().ok()?;
    let rms_end = 48.min(line.len());
    let std_dev = line.get(35..rms_end).and_then(|s| s.trim().parse::<f64>().ok()).unwrap_or(0.0);
    let (bias_type, obs1, obs2, value) = map_diff_type_to_bias(diff_type, val_ns)?;

    let (start_time, end_time) = if line.len() >= 88 {
        let s1 = line.get(48..68).and_then(parse_ymd_hms).unwrap_or(default_span.0);
        let s2 = line.get(68..88).and_then(parse_ymd_hms).unwrap_or(default_span.1);
        (s1, s2)
    } else {
        default_span
    };

    Some(BiasRecord {
        bias_type,
        sat,
        station: None,
        obs1,
        obs2,
        start_time,
        end_time,
        unit: "ns".to_string(),
        value,
        std_dev,
    })
}

pub fn parse_bernese_dcb<R: BufRead>(mut reader: R) -> Result<Vec<BiasRecord>, String> {
    let mut records = Vec::new();
    let mut current_diff_type: Option<DcbDiffType> = None;
    let default_span = (
        GpsTime::from_calendar(2000, 1, 1, 0, 0, 0.0),
        GpsTime::from_calendar(2050, 1, 1, 0, 0, 0.0),
    );
    let mut line_buf = String::new();
    while reader.read_line(&mut line_buf).map_err(|e| e.to_string())? > 0 {
        let trimmed = line_buf.trim_end();
        if trimmed.starts_with("DIFFERENTIAL") {
            current_diff_type = parse_dcb_diff_type(trimmed);
        } else if let Some(diff_type) = current_diff_type {
            if let Some(rec) = parse_dcb_sat_line(trimmed, diff_type, default_span) {
                records.push(rec);
            }
        }
        line_buf.clear();
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_diff_type() {
        assert_eq!(parse_dcb_diff_type("DIFFERENTIAL (P1-C1) CODE BIASES"), Some(DcbDiffType::P1C1));
        assert_eq!(parse_dcb_diff_type("DIFFERENTIAL (P2-C2) CODE BIASES"), Some(DcbDiffType::P2C2));
        assert_eq!(parse_dcb_diff_type("DIFFERENTIAL (P1-P2) CODE BIASES"), Some(DcbDiffType::P1P2));
        assert_eq!(parse_dcb_diff_type("HEADER LINE"), None);
    }

    #[test]
    fn test_parse_sat_id() {
        let sat_g = parse_sat_id("G01").unwrap();
        assert_eq!(sat_g.constellation, Constellation::Gps);
        assert_eq!(sat_g.prn, 1);

        let sat_r = parse_sat_id("R16").unwrap();
        assert_eq!(sat_r.constellation, Constellation::Glonass);
        assert_eq!(sat_r.prn, 16);

        assert!(parse_sat_id("XYZ").is_none());
        assert!(parse_sat_id("G").is_none());
    }

    #[test]
    fn test_parse_dcb_sat_and_station_filtering() {
        let default_span = (
            GpsTime::from_calendar(2020, 1, 1, 0, 0, 0.0),
            GpsTime::from_calendar(2021, 1, 1, 0, 0, 0.0),
        );
        let sat_line = "R16                           2.738       0.032";
        let rec = parse_dcb_sat_line(sat_line, DcbDiffType::P2C2, default_span).unwrap();
        assert_eq!(rec.sat, SatelliteId { constellation: Constellation::Glonass, prn: 16 });
        assert_eq!(rec.obs1.to_string(), "C2C");
        assert!((rec.value - (-2.738)).abs() < 1e-6);
        assert!((rec.std_dev - 0.032).abs() < 1e-6);

        let station_line = "R     ABMF                    0.361       0.005";
        assert!(parse_dcb_sat_line(station_line, DcbDiffType::P2C2, default_span).is_none());
    }

    #[test]
    fn test_parse_bernese_dcb_buffer() {
        let text = "CODE'S MONTHLY GNSS P2-C2 DCB SOLUTION\n\
                    DIFFERENTIAL (P2-C2) CODE BIASES FOR SATELLITES AND RECEIVERS:\n\
                    PRN / STATION NAME        VALUE (NS)  RMS (NS)\n\
                    G01                           1.352       0.004\n\
                    R07                          -0.472       0.006\n\
                    R16                           2.738       0.032\n\
                    G     ABMF                   -1.413       0.005\n";
        let recs = parse_bernese_dcb(Cursor::new(text)).unwrap();
        assert_eq!(recs.len(), 3);
        assert_eq!(recs[0].sat.prn, 1);
        assert!((recs[0].value - (-1.352)).abs() < 1e-6);
        assert_eq!(recs[1].sat.prn, 7);
        assert!((recs[1].value - 0.472).abs() < 1e-6);
        assert_eq!(recs[2].sat.prn, 16);
        assert!((recs[2].value - (-2.738)).abs() < 1e-6);
    }

    #[test]
    fn test_sinex_bias_append_bernese_dcb() {
        let text = "CODE'S MONTHLY GNSS P2-C2 DCB SOLUTION\n\
                    DIFFERENTIAL (P2-C2) CODE BIASES FOR SATELLITES AND RECEIVERS:\n\
                    PRN / STATION NAME        VALUE (NS)  RMS (NS)\n\
                    R16                           2.738       0.032\n";
        let mut bia = crate::sinex_bia::SinexBias::new(vec![]);
        let count = bia.load_bernese_dcb(Cursor::new(text)).unwrap();
        assert_eq!(count, 1);
        let sat_r16 = SatelliteId { constellation: Constellation::Glonass, prn: 16 };
        let obs_c2c = ObsCode::from_str("C2C").unwrap();
        let t = GpsTime::from_calendar(2020, 12, 24, 12, 0, 0.0);
        let bias = bia.get_bias(sat_r16, obs_c2c, t);
        assert_eq!(bias, Some(-2.738));
        let bias_m = bia.lookup_bias_m(sat_r16, obs_c2c, t).unwrap();
        let expected_m = -2.738e-9 * gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        assert!((bias_m - expected_m).abs() < 1e-6);
    }
}
