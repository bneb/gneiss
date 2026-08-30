//! NOAA VDatum / PROJ GTX binary geoid grid parser.

use alloc::vec::Vec;
use super::GeoidGrid;

/// Parse a NOAA VDatum / PROJ `.gtx` binary grid.
///
/// Format: 40-byte big-endian header (min_lat, min_lon, dlat, dlon, n_lat, n_lon)
/// followed by n_lat * n_lon big-endian f32 elevation values.
pub fn parse_gtx(bytes: &[u8]) -> Result<GeoidGrid, &'static str> {
    if bytes.len() < 40 {
        return Err("GTX file too small (less than 40-byte header)");
    }

    let min_lat = f64::from_be_bytes(bytes[0..8].try_into().map_err(|_| "invalid min_lat")?);
    let min_lon = f64::from_be_bytes(bytes[8..16].try_into().map_err(|_| "invalid min_lon")?);
    let dlat = f64::from_be_bytes(bytes[16..24].try_into().map_err(|_| "invalid dlat")?);
    let dlon = f64::from_be_bytes(bytes[24..32].try_into().map_err(|_| "invalid dlon")?);
    let n_lat = u32::from_be_bytes(bytes[32..36].try_into().map_err(|_| "invalid n_lat")?) as usize;
    let n_lon = u32::from_be_bytes(bytes[36..40].try_into().map_err(|_| "invalid n_lon")?) as usize;

    if n_lat < 2 || n_lon < 2 || dlat <= 0.0 || dlon <= 0.0 {
        return Err("invalid GTX grid dimensions or step size");
    }

    let expected_data_len = n_lat * n_lon * 4;
    if bytes.len() < 40 + expected_data_len {
        return Err("truncated GTX grid data");
    }

    let mut grid = Vec::with_capacity(n_lat * n_lon);
    let mut offset = 40;
    for _ in 0..(n_lat * n_lon) {
        let val = f32::from_be_bytes(bytes[offset..offset + 4].try_into().map_err(|_| "read float error")?);
        grid.push(val);
        offset += 4;
    }

    let max_lat = min_lat + (n_lat - 1) as f64 * dlat;
    let max_lon = min_lon + (n_lon - 1) as f64 * dlon;

    GeoidGrid::new(min_lat, max_lat, min_lon, max_lon, dlat, dlon, n_lat, n_lon, grid)
        .ok_or("failed to construct GeoidGrid from GTX")
}
