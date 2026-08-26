//! BLQ Ocean Tide Loading (OTL) File Parser.
//!
//! Parses standard IERS format BLQ files (from Onsala Space Observatory /
//! H.-G. Scherneck OTL provider) containing 11 tidal constituent parameters
//! (M2, S2, N2, K2, K1, O1, P1, Q1, Mf, Mm, Ssa) per station:
//!
//! 1. Radial/Up amplitude (m)
//! 2. West amplitude (m)
//! 3. South amplitude (m)
//! 4. Radial/Up phase (degrees)
//! 5. West phase (degrees)
//! 6. South phase (degrees)

use core::f64::consts::PI;
use gneiss_core::tides::{OceanTideParams, OTL_CONSTITUENT_COUNT};
use std::collections::HashMap;
use std::io::BufRead;

const DEG_TO_RAD: f64 = PI / 180.0;

/// Database of station Ocean Tide Loading parameters parsed from BLQ.
#[derive(Debug, Clone, Default)]
pub struct BlqDatabase {
    pub stations: HashMap<String, OceanTideParams>,
}

impl BlqDatabase {
    /// Retrieve OTL parameters for a station (case-insensitive 4-character ID).
    #[must_use]
    pub fn get(&self, station: &str) -> Option<&OceanTideParams> {
        self.stations.get(&station.to_ascii_uppercase())
    }
}

/// Parse a BLQ file reader into a [`BlqDatabase`].
pub fn parse_blq<R: BufRead>(reader: R) -> Result<BlqDatabase, &'static str> {
    let mut db = BlqDatabase::default();
    let mut lines = reader.lines();
    
    while let Some(line_res) = lines.next() {
        let line = line_res.map_err(|_| "Failed to read BLQ line")?;
        let trimmed = line.trim();
        if trimmed.starts_with("$$") || trimmed.is_empty() {
            continue;
        }

        let station_name = parse_station_header(trimmed);
        if let Some(name) = station_name {
            let params = parse_station_block(&mut lines)?;
            db.stations.insert(name, params);
        }
    }

    Ok(db)
}

fn parse_station_header(line: &str) -> Option<String> {
    let first_token = line.split_whitespace().next()?;
    if first_token.len() >= 3 && !first_token.starts_with('$') {
        Some(first_token.to_ascii_uppercase())
    } else {
        None
    }
}

fn parse_station_block<I>(lines: &mut I) -> Result<OceanTideParams, &'static str>
where
    I: Iterator<Item = std::io::Result<String>>,
{
    let mut rows: Vec<[f64; OTL_CONSTITUENT_COUNT]> = Vec::with_capacity(6);

    while rows.len() < 6 {
        let next_line = match lines.next() {
            Some(Ok(l)) => l,
            _ => return Err("Unexpected EOF while parsing 6 BLQ data rows"),
        };
        let t = next_line.trim();
        if t.starts_with("$$") || t.is_empty() {
            continue;
        }

        let vals = parse_11_floats(t)?;
        rows.push(vals);
    }

    let mut amp_radial_m = [0.0; OTL_CONSTITUENT_COUNT];
    let mut amp_west_m = [0.0; OTL_CONSTITUENT_COUNT];
    let mut amp_south_m = [0.0; OTL_CONSTITUENT_COUNT];
    let mut ph_radial_rad = [0.0; OTL_CONSTITUENT_COUNT];
    let mut ph_west_rad = [0.0; OTL_CONSTITUENT_COUNT];
    let mut ph_south_rad = [0.0; OTL_CONSTITUENT_COUNT];

    for i in 0..OTL_CONSTITUENT_COUNT {
        amp_radial_m[i] = rows[0][i];
        amp_west_m[i] = rows[1][i];
        amp_south_m[i] = rows[2][i];
        ph_radial_rad[i] = rows[3][i] * DEG_TO_RAD;
        ph_west_rad[i] = rows[4][i] * DEG_TO_RAD;
        ph_south_rad[i] = rows[5][i] * DEG_TO_RAD;
    }

    Ok(OceanTideParams {
        amp_radial_m,
        amp_west_m,
        amp_south_m,
        ph_radial_rad,
        ph_west_rad,
        ph_south_rad,
    })
}

fn parse_11_floats(line: &str) -> Result<[f64; OTL_CONSTITUENT_COUNT], &'static str> {
    let mut arr = [0.0; OTL_CONSTITUENT_COUNT];
    let mut count = 0;

    for token in line.split_whitespace() {
        if count >= OTL_CONSTITUENT_COUNT {
            break;
        }
        let val = token.parse::<f64>().map_err(|_| "Failed to parse float in BLQ row")?;
        arr[count] = val;
        count += 1;
    }

    if count == OTL_CONSTITUENT_COUNT {
        Ok(arr)
    } else {
        Err("Row does not contain 11 constituent values")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_BLQ: &str = r#"
$$ Ocean loading displacement parameters
$$ Column order: M2 S2 N2 K2 K1 O1 P1 Q1 Mf Mm Ssa
$$
  P224
$$
  .01524  .00512  .00341  .00140  .02103  .01452  .00690  .00301  .00110  .00060  .00030
  .00612  .00201  .00120  .00054  .00841  .00582  .00271  .00120  .00041  .00021  .00010
  .00411  .00150  .00092  .00041  .00520  .00391  .00170  .00080  .00030  .00015  .00008
  295.40  300.20  278.10  302.50  210.30  205.10  210.80  200.20   15.00   12.00    5.00
  112.30  118.50   95.20  120.10   45.20   40.10   46.00   35.20  340.00  335.00  320.00
   85.40   90.20   70.10   92.00   25.10   20.50   26.00   15.00  310.00  305.00  290.00
$$
  OHLN
$$
  .01820  .00620  .00410  .00170  .02340  .01620  .00780  .00340  .00120  .00070  .00035
  .00710  .00240  .00140  .00065  .00950  .00660  .00310  .00135  .00045  .00025  .00012
  .00480  .00180  .00110  .00050  .00610  .00460  .00200  .00095  .00035  .00018  .00009
  296.10  301.00  279.00  303.20  211.00  206.00  211.50  201.00   15.50   12.50    5.50
  113.00  119.00   96.00  121.00   46.00   41.00   47.00   36.00  341.00  336.00  321.00
   86.00   91.00   71.00   93.00   26.00   21.00   27.00   16.00  311.00  306.00  291.00
"#;

    #[test]
    fn test_parse_blq_multiple_stations() {
        let cursor = std::io::Cursor::new(SAMPLE_BLQ);
        let db = parse_blq(cursor).expect("must parse sample BLQ");
        assert_eq!(db.stations.len(), 2);

        let p224 = db.get("P224").expect("P224 must exist");
        assert!((p224.amp_radial_m[0] - 0.01524).abs() < 1e-6);
        assert!((p224.amp_west_m[0] - 0.00612).abs() < 1e-6);
        assert!((p224.amp_south_m[0] - 0.00411).abs() < 1e-6);
        // Phase check: 295.4 deg in rad
        let exp_ph_rad = 295.40 * DEG_TO_RAD;
        assert!((p224.ph_radial_rad[0] - exp_ph_rad).abs() < 1e-6);

        let ohln = db.get("ohln").expect("case-insensitive lookup must find OHLN");
        assert!((ohln.amp_radial_m[0] - 0.01820).abs() < 1e-6);
    }
}
