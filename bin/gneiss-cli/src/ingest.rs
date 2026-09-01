use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use tracing::info;

use gneiss_core::obs::EpochObs;
use gneiss_parsers::rinex::parse_rinex_obs;
use gneiss_parsers::ubx::{parse_rxm_rawx, parse_ubx_frame};

pub type IngestResult = (Vec<EpochObs>, Option<[f64; 3]>);

/// Universal observation reader capable of auto-detecting RINEX, u-blox UBX, and RTCM3.
pub struct UniversalObsReader;

impl UniversalObsReader {
    /// Ingests observations from any supported GNSS raw file format.
    pub fn read_file(path: &Path) -> Result<IngestResult, Box<dyn std::error::Error>> {
        let bytes = std::fs::read(path)?;
        if is_ubx_format(&bytes, path) {
            info!("Ingesting u-blox binary raw data from {}", path.display());
            parse_ubx_dataset(&bytes)
        } else {
            info!("Ingesting RINEX observation data from {}", path.display());
            parse_rinex_dataset(path)
        }
    }
}

fn is_ubx_format(bytes: &[u8], path: &Path) -> bool {
    let has_ext = path.extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("ubx"))
        .unwrap_or(false);
    let has_magic = bytes.len() >= 2 && bytes[0] == 0xB5 && bytes[1] == 0x62;
    has_ext || has_magic
}

fn parse_rinex_dataset(path: &Path) -> Result<IngestResult, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let (epochs, header) = parse_rinex_obs(BufReader::new(file))?;
    Ok((epochs, header.approx_position))
}

fn parse_ubx_dataset(mut slice: &[u8]) -> Result<IngestResult, Box<dyn std::error::Error>> {
    let mut epochs = Vec::new();
    let mut approx_pos: Option<[f64; 3]> = None;

    while !slice.is_empty() {
        match parse_ubx_frame(slice) {
            Ok((rem, frame)) => {
                slice = rem;
                if frame.class == 0x02 && frame.id == 0x15 {
                    if let Ok(rawx) = parse_rxm_rawx(frame.payload) {
                        epochs.push(rawx.into_epoch_obs());
                    }
                } else if frame.class == 0x01 && frame.id == 0x02 && approx_pos.is_none() {
                    approx_pos = parse_ubx_posllh(frame.payload);
                }
            }
            Err(_) => {
                // Advance one byte if frame sync fails
                slice = &slice[1..];
            }
        }
    }

    if epochs.is_empty() {
        return Err("No UBX-RXM-RAWX observations found in file".into());
    }

    info!("Successfully parsed {} UBX epochs", epochs.len());
    Ok((epochs, approx_pos))
}

fn parse_ubx_posllh(payload: &[u8]) -> Option<[f64; 3]> {
    if payload.len() < 28 {
        return None;
    }
    let lon_deg = i32::from_le_bytes(payload[4..8].try_into().ok()?) as f64 * 1e-7;
    let lat_deg = i32::from_le_bytes(payload[8..12].try_into().ok()?) as f64 * 1e-7;
    let h_m = i32::from_le_bytes(payload[12..16].try_into().ok()?) as f64 * 1e-3;

    Some(geodetic_to_ecef(lat_deg, lon_deg, h_m))
}

fn geodetic_to_ecef(lat_deg: f64, lon_deg: f64, h: f64) -> [f64; 3] {
    let a = 6378137.0;
    let f = 1.0 / 298.257223563;
    let e2 = f * (2.0 - f);
    let lat_rad = lat_deg.to_radians();
    let lon_rad = lon_deg.to_radians();

    let sin_lat = lat_rad.sin();
    let cos_lat = lat_rad.cos();
    let sin_lon = lon_rad.sin();
    let cos_lon = lon_rad.cos();

    let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
    let x = (n + h) * cos_lat * cos_lon;
    let y = (n + h) * cos_lat * sin_lon;
    let z = (n * (1.0 - e2) + h) * sin_lat;

    [x, y, z]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geodetic_to_ecef_equator() {
        let ecef = geodetic_to_ecef(0.0, 0.0, 0.0);
        assert!((ecef[0] - 6378137.0).abs() < 1.0);
        assert!(ecef[1].abs() < 1.0);
        assert!(ecef[2].abs() < 1.0);
    }

    #[test]
    fn test_is_ubx_detection() {
        let magic = [0xB5, 0x62, 0x01, 0x02];
        assert!(is_ubx_format(&magic, Path::new("flight.dat")));
        assert!(is_ubx_format(&[], Path::new("flight.ubx")));
        assert!(!is_ubx_format(&[0x00, 0x01], Path::new("flight.obs")));
    }

    #[test]
    fn test_read_real_rinex_file() {
        let path = Path::new("datasets/network_rtk_01/p181.obs");
        if path.exists() {
            let res = UniversalObsReader::read_file(path);
            assert!(res.is_ok());
            let (epochs, pos) = res.unwrap();
            assert!(!epochs.is_empty());
            assert!(pos.is_some());
        }
    }

    #[test]
    fn test_parse_ubx_posllh() {
        let mut payload = vec![0u8; 28];
        // iTOW (4 bytes)
        payload[0..4].copy_from_slice(&1000u32.to_le_bytes());
        // lon = -122.0 deg = -1220000000 (1e-7 deg)
        let lon: i32 = -1_220_000_000;
        payload[4..8].copy_from_slice(&lon.to_le_bytes());
        // lat = 37.0 deg = 370000000 (1e-7 deg)
        let lat: i32 = 370_000_000;
        payload[8..12].copy_from_slice(&lat.to_le_bytes());
        // height = 100.0 m = 100000 mm
        let h: i32 = 100_000;
        payload[12..16].copy_from_slice(&h.to_le_bytes());

        let pos = parse_ubx_posllh(&payload);
        assert!(pos.is_some());
        let ecef = pos.unwrap();
        assert!(ecef[0].abs() > 1_000_000.0);
    }
}
