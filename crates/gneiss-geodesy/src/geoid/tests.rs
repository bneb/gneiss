#![allow(clippy::unwrap_used)]

use super::*;
use alloc::vec;
use alloc::vec::Vec;

fn create_sample_grid() -> GeoidGrid {
    let grid = alloc::vec![
        30.0, 31.0, 32.0,
        32.0, 33.0, 34.0,
        34.0, 35.0, 36.0,
    ];
    GeoidGrid::new(30.0, 32.0, 130.0, 132.0, 1.0, 1.0, 3, 3, grid).expect("valid grid")
}

#[test]
fn test_geoid_grid_exact_points_and_boundaries() {
    let grid = create_sample_grid();
    let n_sw = grid.undulation(30.0, 130.0).expect("undulation");
    assert!((n_sw - 30.0).abs() < 1e-6);

    let n_se = grid.undulation(30.0, 132.0).expect("undulation");
    assert!((n_se - 32.0).abs() < 1e-6);

    let n_nw = grid.undulation(32.0, 130.0).expect("undulation");
    assert!((n_nw - 34.0).abs() < 1e-6);

    let n_ne = grid.undulation(32.0, 132.0).expect("undulation");
    assert!((n_ne - 36.0).abs() < 1e-6);

    let n_mid = grid.undulation(31.0, 131.0).expect("undulation");
    assert!((n_mid - 33.0).abs() < 1e-6);
}

#[test]
fn test_gtx_roundtrip_parsing() {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&40.0_f64.to_be_bytes()); // min_lat
    bytes.extend_from_slice(&(-100.0_f64).to_be_bytes()); // min_lon
    bytes.extend_from_slice(&0.5_f64.to_be_bytes()); // dlat
    bytes.extend_from_slice(&0.5_f64.to_be_bytes()); // dlon
    bytes.extend_from_slice(&2_u32.to_be_bytes()); // n_lat
    bytes.extend_from_slice(&2_u32.to_be_bytes()); // n_lon

    // 4 float values: [10.0, 20.0, 30.0, 40.0]
    bytes.extend_from_slice(&10.0_f32.to_be_bytes());
    bytes.extend_from_slice(&20.0_f32.to_be_bytes());
    bytes.extend_from_slice(&30.0_f32.to_be_bytes());
    bytes.extend_from_slice(&40.0_f32.to_be_bytes());

    let grid = GeoidGrid::from_gtx_bytes(&bytes).expect("parse gtx");
    assert_eq!(grid.n_lat, 2);
    assert_eq!(grid.n_lon, 2);
    let n = grid.undulation(40.0, -100.0).expect("undulation");
    assert!((n - 10.0).abs() < 1e-5);
}

#[test]
fn test_byn_roundtrip_parsing() {
    let mut bytes = vec![0u8; 80]; // 80-byte header
    let min_lat_mas: i32 = 40 * 3_600_000;
    let max_lat_mas: i32 = 41 * 3_600_000;
    let min_lon_mas: i32 = -100 * 3_600_000;
    let max_lon_mas: i32 = -99 * 3_600_000;
    let dlat_mas: i32 = 1 * 3_600_000;
    let dlon_mas: i32 = 1 * 3_600_000;

    bytes[8..12].copy_from_slice(&min_lat_mas.to_le_bytes());
    bytes[12..16].copy_from_slice(&max_lat_mas.to_le_bytes());
    bytes[16..20].copy_from_slice(&min_lon_mas.to_le_bytes());
    bytes[20..24].copy_from_slice(&max_lon_mas.to_le_bytes());
    bytes[24..28].copy_from_slice(&dlat_mas.to_le_bytes());
    bytes[28..32].copy_from_slice(&dlon_mas.to_le_bytes());
    bytes[42..44].copy_from_slice(&1_i16.to_le_bytes()); // i16
    bytes[44..48].copy_from_slice(&0.001_f32.to_le_bytes()); // scale factor (mm -> m)

    // Data: 2 rows x 2 cols. Row 0 (North, lat=41): [30_000, 25_000], Row 1 (South, lat=40): [10_000, 20_000]
    bytes.extend_from_slice(&30_000_i16.to_le_bytes());
    bytes.extend_from_slice(&25_000_i16.to_le_bytes());
    bytes.extend_from_slice(&10_000_i16.to_le_bytes());
    bytes.extend_from_slice(&20_000_i16.to_le_bytes());

    let grid = GeoidGrid::from_byn_bytes(&bytes).expect("parse byn");
    assert_eq!(grid.n_lat, 2);
    assert_eq!(grid.n_lon, 2);
    let n_sw = grid.undulation(40.0, -100.0).expect("undulation");
    assert!((n_sw - 10.0).abs() < 1e-3);
    let n_nw = grid.undulation(41.0, -100.0).expect("undulation");
    assert!((n_nw - 30.0).abs() < 1e-3);
}
