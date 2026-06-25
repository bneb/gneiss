//! IONEX (IONosphere Map Exchange) format parser.
//!
//! Parses IGS IONEX v1.0 files containing global TEC (Total Electron Content) maps.
//! These maps are used as ionospheric priors for UDUC PPP — replacing the Klobuchar
//! broadcast model (~3 m std) with CODE/IGS gridded TEC (~5 cm accuracy).
//!
//! Format reference: IONEX v1.0, Werner Gurtner, 1998.
//! Sample: ftp://ftp.aiub.unibe.ch/CODE/2019/CODG3350.19I.Z

use std::io::BufRead;

/// A single TEC map at a given epoch.
#[derive(Debug, Clone)]
pub struct TecMap {
    /// GPS time of this map
    pub time: gneiss_core::time::GpsTime,
    /// VTEC values in TECU (10^16 el/m²), shape: [nlat][nlon]
    /// Row-major: lat index varies slowest, lon index fastest.
    pub tec: Vec<Vec<f64>>,
}

/// Parsed IONEX grid covering a 24h window.
#[derive(Debug, Clone)]
pub struct IonexGrid {
    /// Latitude grid: start (degrees north, positive), end, step
    pub lat1: f64,
    pub lat2: f64,
    pub dlat: f64,
    /// Longitude grid: start, end, step (degrees east)
    pub lon1: f64,
    pub lon2: f64,
    pub dlon: f64,
    /// Ionospheric shell height in km (typically 450 km)
    pub height_km: f64,
    /// RMS map (optional, same shape as tec maps)
    pub rms_maps: Vec<Vec<Vec<f64>>>,
    /// TEC maps ordered by time
    pub tec_maps: Vec<TecMap>,
    /// Exponent applied to TEC values (e.g., -1 means divide by 10)
    pub exponent: i32,
}

/// Parse an IONEX file from a buffered reader.
///
/// Returns `Ok(IonexGrid)` on success, or `Err(String)` with a human-readable
/// description of the parse failure.
pub fn parse_ionex<R: BufRead>(reader: R) -> Result<IonexGrid, String> {
    let mut lines = Vec::new();
    for line in reader.lines() {
        lines.push(line.map_err(|e| format!("I/O error: {}", e))?);
    }

    let mut grid = IonexGrid {
        lat1: 87.5,
        lat2: -87.5,
        dlat: -2.5,
        lon1: -180.0,
        lon2: 180.0,
        dlon: 5.0,
        height_km: 450.0,
        rms_maps: Vec::new(),
        tec_maps: Vec::new(),
        exponent: -1,
    };

    let mut in_tec = false;
    let mut in_rms = false;
    let mut current_time: Option<gneiss_core::time::GpsTime> = None;
    let mut current_tec: Vec<Vec<f64>> = Vec::new();
    let mut current_rms: Vec<Vec<f64>> = Vec::new();
    let mut _lat_count: usize = 0;
    let mut lon_count: usize = 0;

    let mut i = 0;
    while i < lines.len() {
        let line = &lines[i];

        if line.len() < 60 {
            i += 1;
            continue;
        }

        let label = line[60..].trim();

        match label {
            "EXPONENT" => {
                grid.exponent = line[0..6].trim().parse::<i32>().unwrap_or(-1);
            }
            "HGT1 / HGT2 / DHGT" => {
                let parts: Vec<f64> = line[0..20]
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if !parts.is_empty() {
                    grid.height_km = parts[0];
                }
            }
            "LAT1 / LAT2 / DLAT" => {
                let parts: Vec<f64> = line[0..20]
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if parts.len() >= 3 {
                    grid.lat1 = parts[0];
                    grid.lat2 = parts[1];
                    grid.dlat = parts[2];
                }
            }
            "LON1 / LON2 / DLON" => {
                let parts: Vec<f64> = line[0..20]
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if parts.len() >= 3 {
                    grid.lon1 = parts[0];
                    grid.lon2 = parts[1];
                    grid.dlon = parts[2];
                }
            }
            "START OF TEC MAP" => {
                in_tec = true;
                in_rms = false;
                current_tec = Vec::new();
                current_rms = Vec::new();
                _lat_count = 0;
                // Compute expected longitude count. When lon2 <= lon1, the grid
                // wraps around 360° (e.g., LON1=0, LON2=0, DLON=60 → 6 values).
                lon_count = if grid.lon2 > grid.lon1 {
                    ((grid.lon2 - grid.lon1) / grid.dlon.abs()).round() as usize + 1
                } else {
                    (360.0 / grid.dlon.abs()).round() as usize
                };
                current_time = None;
            }
            "START OF RMS MAP" => {
                in_rms = true;
                in_tec = false;
                current_rms = Vec::new();
                _lat_count = 0;
            }
            "END OF TEC MAP" | "END OF RMS MAP" => {
                if in_tec && current_time.is_some() {
                    grid.tec_maps.push(TecMap {
                        time: current_time.take().unwrap(),
                        tec: std::mem::take(&mut current_tec),
                    });
                }
                if in_rms {
                    grid.rms_maps.push(std::mem::take(&mut current_rms));
                }
                in_tec = false;
                in_rms = false;
            }
            "EPOCH OF CURRENT MAP" => {
                if in_tec && current_time.is_none() {
                    current_time = parse_ionex_epoch(&line[0..36]);
                }
            }
            "LAT/LON1/LON2/DLON/H" => {
                // Data descriptor line: e.g. "   87.5-180.0 180.0   5.0 450.0"
                // The TEC values follow on the NEXT line(s).  The descriptor
                // line itself does NOT contain TEC data — only metadata.
                if in_tec || in_rms {
                    let mut tec_row = Vec::with_capacity(lon_count);
                    let mut col = 0;

                    // Read TEC values from subsequent lines (skip the descriptor)
                    loop {
                        i += 1;
                        if i >= lines.len() {
                            break;
                        }
                        let data_line = &lines[i];
                        // If this line has an IONEX label, it's not TEC data
                        if data_line.len() >= 60 {
                            let lbl = data_line[60..].trim();
                            if lbl.len() > 2 && !lbl.starts_with(' ') {
                                // This is a new IONEX record — back up and stop
                                i -= 1;
                                break;
                            }
                        }
                        // Parse up to 12 values (5 chars each) from columns 0-59
                        let end = data_line.len().min(60);
                        for j in (0..end).step_by(5) {
                            if j + 5 <= end && col < lon_count {
                                let val_str = data_line[j..j + 5].trim();
                                if let Ok(val) = val_str.parse::<f64>() {
                                    tec_row.push(val);
                                    col += 1;
                                } else if !val_str.is_empty() {
                                    tec_row.push(0.0);
                                    col += 1;
                                }
                            }
                        }
                        if col >= lon_count {
                            break;
                        }
                    }

                    if in_tec {
                        current_tec.push(tec_row);
                    } else {
                        current_rms.push(tec_row);
                    }
                    _lat_count += 1;
                }
            }
            "END OF FILE" => {
                break;
            }
            _ => {}
        }
        i += 1;
    }

    if grid.tec_maps.is_empty() {
        return Err("No TEC maps found in IONEX file".into());
    }

    Ok(grid)
}

/// Parse an IONEX epoch line: "  YYYY  MM  DD  HH  MM  SS"
fn parse_ionex_epoch(s: &str) -> Option<gneiss_core::time::GpsTime> {
    let parts: Vec<&str> = s.split_whitespace().collect();
    if parts.len() >= 6 {
        let year = parts[0].parse::<i32>().ok()?;
        let month = parts[1].parse::<i32>().ok()?;
        let day = parts[2].parse::<i32>().ok()?;
        let hour = parts[3].parse::<i32>().ok()?;
        let min = parts[4].parse::<i32>().ok()?;
        let sec = parts[5].parse::<f64>().ok()?;
        Some(gneiss_core::time::GpsTime::from_calendar(
            year, month, day, hour, min, sec,
        ))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_ionex_minimal() {
        // Simple test with compact 2×3 grid (2 lats, 3 lons)
        let data = r#"     1.0            IONOSPHERE MAPS     GNSS                IONEX VERSION / TYPE
SAMPLE              TEST AGENCY          05-DEC-19 20:20     PGM / RUN BY / DATE
  2019    12     1     0     0     0                        EPOCH OF FIRST MAP
  2019    12     1     2     0     0                        EPOCH OF LAST MAP
   450.0 450.0   0.0                                        HGT1 / HGT2 / DHGT
    60.0 -60.0 -60.0                                        LAT1 / LAT2 / DLAT
     0.0   0.0  60.0                                        LON1 / LON2 / DLON
    -1                                                      EXPONENT
     1                                                      START OF TEC MAP
  2019    12     1     0     0     0                        EPOCH OF CURRENT MAP
    60.0   0.0   0.0  60.0 450.0                            LAT/LON1/LON2/DLON/H
    5    6    7
   -60.0   0.0   0.0  60.0 450.0                            LAT/LON1/LON2/DLON/H
  105  106  107
     1                                                      END OF TEC MAP
     1                                                      END OF FILE
"#;
        let grid = parse_ionex(data.as_bytes()).unwrap();
        assert_eq!(grid.tec_maps.len(), 1);
        assert_eq!(grid.tec_maps[0].tec.len(), 2); // 2 lat bands
        assert_eq!(grid.tec_maps[0].tec[0].len(), 3); // 3 longitudes
        assert!((grid.tec_maps[0].tec[0][0] - 5.0).abs() < 0.01);
        assert!((grid.tec_maps[0].tec[1][0] - 105.0).abs() < 0.01);
        assert_eq!(grid.height_km, 450.0);
        assert_eq!(grid.exponent, -1);
    }

    #[test]
    fn test_parse_ionex_real_file() {
        // Parse the actual CODE IONEX file if available
        let path = std::path::Path::new("datasets/igs/codg3350.19i");
        if path.exists() {
            let file = std::fs::File::open(path).unwrap();
            let reader = std::io::BufReader::new(file);
            let grid = parse_ionex(reader).unwrap();
            assert!(grid.tec_maps.len() >= 1);
            // Verify grid dimensions
            let nlat = ((grid.lat1 - grid.lat2).abs() / grid.dlat.abs()).round() as usize + 1;
            let nlon = ((grid.lon2 - grid.lon1) / grid.dlon.abs()).round() as usize;
            for map in &grid.tec_maps {
                assert_eq!(map.tec.len(), nlat);
                for row in &map.tec {
                    assert_eq!(row.len(), nlon);
                }
            }
        }
    }
}
