//! NRCan (Natural Resources Canada) BYN binary geoid grid parser.

use alloc::vec::Vec;
use super::GeoidGrid;

/// Parse an NRCan `.byn` binary geoid grid (CGG2013, HTv2).
///
/// Header: 80 bytes little-endian.
/// Bounds and steps are stored in milliarcseconds (0.001" or mas).
/// 1 degree = 3,600,000 milliarcseconds.
pub fn parse_byn(bytes: &[u8]) -> Result<GeoidGrid, &'static str> {
    if bytes.len() < 80 {
        return Err("BYN file too small (less than 80-byte header)");
    }

    let min_lat_mas = i32::from_le_bytes(bytes[8..12].try_into().map_err(|_| "invalid min_lat")?);
    let max_lat_mas = i32::from_le_bytes(bytes[12..16].try_into().map_err(|_| "invalid max_lat")?);
    let min_lon_mas = i32::from_le_bytes(bytes[16..20].try_into().map_err(|_| "invalid min_lon")?);
    let max_lon_mas = i32::from_le_bytes(bytes[20..24].try_into().map_err(|_| "invalid max_lon")?);
    let dlat_mas = i32::from_le_bytes(bytes[24..28].try_into().map_err(|_| "invalid dlat")?);
    let dlon_mas = i32::from_le_bytes(bytes[28..32].try_into().map_err(|_| "invalid dlon")?);
    let data_type = i16::from_le_bytes(bytes[42..44].try_into().map_err(|_| "invalid data_type")?);
    let scale_factor = f32::from_le_bytes(bytes[44..48].try_into().map_err(|_| "invalid scale_factor")?);

    let scale = if scale_factor > 0.0 { scale_factor } else { 0.001 };
    let mas_to_deg = 1.0 / 3_600_000.0;

    let min_lat = min_lat_mas as f64 * mas_to_deg;
    let max_lat = max_lat_mas as f64 * mas_to_deg;
    let min_lon = min_lon_mas as f64 * mas_to_deg;
    let max_lon = max_lon_mas as f64 * mas_to_deg;
    let dlat = dlat_mas as f64 * mas_to_deg;
    let dlon = dlon_mas as f64 * mas_to_deg;

    if dlat <= 0.0 || dlon <= 0.0 || max_lat <= min_lat || max_lon <= min_lon {
        return Err("invalid BYN grid bounds or step size");
    }

    let n_lat = ((max_lat_mas - min_lat_mas) / dlat_mas + 1) as usize;
    let n_lon = ((max_lon_mas - min_lon_mas) / dlon_mas + 1) as usize;

    let mut raw_rows: Vec<Vec<f32>> = Vec::with_capacity(n_lat);
    let mut offset = 80;

    for _r in 0..n_lat {
        let mut row = Vec::with_capacity(n_lon);
        for _c in 0..n_lon {
            let val = match data_type {
                1 => {
                    if bytes.len() < offset + 2 { return Err("truncated BYN i16 data"); }
                    let raw = i16::from_le_bytes(bytes[offset..offset + 2].try_into().map_err(|_| "read i16")?);
                    offset += 2;
                    raw as f32 * scale
                }
                2 => {
                    if bytes.len() < offset + 4 { return Err("truncated BYN i32 data"); }
                    let raw = i32::from_le_bytes(bytes[offset..offset + 4].try_into().map_err(|_| "read i32")?);
                    offset += 4;
                    raw as f32 * scale
                }
                _ => {
                    if bytes.len() < offset + 4 { return Err("truncated BYN f32 data"); }
                    let raw = f32::from_le_bytes(bytes[offset..offset + 4].try_into().map_err(|_| "read f32")?);
                    offset += 4;
                    raw
                }
            };
            row.push(val);
        }
        raw_rows.push(row);
    }

    // NRCan BYN rows are ordered North-to-South; reverse to South-to-North for GeoidGrid
    let mut grid = Vec::with_capacity(n_lat * n_lon);
    for row in raw_rows.into_iter().rev() {
        grid.extend(row);
    }

    GeoidGrid::new(min_lat, max_lat, min_lon, max_lon, dlat, dlon, n_lat, n_lon, grid)
        .ok_or("failed to construct GeoidGrid from BYN")
}
