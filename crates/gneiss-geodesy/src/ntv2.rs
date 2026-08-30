//! NTv2 (National Transformation version 2) binary datum shift grid (.gsb) parser and interpolator.
//!
//! Standard binary grid shift format used by geodesy agencies worldwide (NGS, NRCan, Ordnance Survey, IGN)
//! for high-precision local datum transformations (e.g. NAD27 <-> NAD83, AGD66 <-> GDA94).

use alloc::vec::Vec;
use alloc::string::String;

/// A single node record in an NTv2 subgrid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ntv2Node {
    /// Latitude shift in arcseconds.
    pub dlat_sec: f32,
    /// Longitude shift in arcseconds (positive West in standard NTv2).
    pub dlon_sec: f32,
    /// Latitude accuracy estimate (meters).
    pub lat_acc_m: f32,
    /// Longitude accuracy estimate (meters).
    pub lon_acc_m: f32,
}

/// An NTv2 subgrid with geographic boundaries and grid step size.
#[derive(Debug, Clone, PartialEq)]
pub struct Ntv2Subgrid {
    pub name: String,
    pub parent: String,
    pub s_lat_deg: f64,
    pub n_lat_deg: f64,
    pub w_lon_deg: f64,
    pub e_lon_deg: f64,
    pub lat_inc_deg: f64,
    pub lon_inc_deg: f64,
    pub n_rows: usize,
    pub n_cols: usize,
    pub nodes: Vec<Ntv2Node>,
}

impl Ntv2Subgrid {
    /// Checks if a given coordinate falls within this subgrid boundary.
    pub fn contains(&self, lat_deg: f64, lon_deg: f64) -> bool {
        lat_deg >= self.s_lat_deg && lat_deg <= self.n_lat_deg &&
        lon_deg >= self.w_lon_deg && lon_deg <= self.e_lon_deg
    }

    /// Bilinear interpolation of (dlat_sec, dlon_sec) at a given (lat, lon).
    pub fn interpolate(&self, lat_deg: f64, lon_deg: f64) -> Option<(f64, f64)> {
        if !self.contains(lat_deg, lon_deg) || self.n_rows < 2 || self.n_cols < 2 {
            return None;
        }

        let row_pos = (lat_deg - self.s_lat_deg) / self.lat_inc_deg;
        let col_pos = (lon_deg - self.w_lon_deg) / self.lon_inc_deg;

        let r = (row_pos.floor() as usize).min(self.n_rows - 2);
        let c = (col_pos.floor() as usize).min(self.n_cols - 2);

        let u = row_pos - r as f64;
        let v = col_pos - c as f64;

        let n00 = &self.nodes[r * self.n_cols + c];
        let n01 = &self.nodes[r * self.n_cols + c + 1];
        let n10 = &self.nodes[(r + 1) * self.n_cols + c];
        let n11 = &self.nodes[(r + 1) * self.n_cols + c + 1];

        let dlat = (1.0 - u) * (1.0 - v) * n00.dlat_sec as f64
            + (1.0 - u) * v * n01.dlat_sec as f64
            + u * (1.0 - v) * n10.dlat_sec as f64
            + u * v * n11.dlat_sec as f64;

        let dlon = (1.0 - u) * (1.0 - v) * n00.dlon_sec as f64
            + (1.0 - u) * v * n01.dlon_sec as f64
            + u * (1.0 - v) * n10.dlon_sec as f64
            + u * v * n11.dlon_sec as f64;

        Some((dlat, dlon))
    }
}

/// Parsed NTv2 datum shift grid containing one or more subgrids.
#[derive(Debug, Clone, PartialEq)]
pub struct Ntv2Grid {
    pub subgrids: Vec<Ntv2Subgrid>,
}

impl Ntv2Grid {
    /// Interpolates (dlat_deg, dlon_deg) shift at a given coordinate across subgrids.
    pub fn interpolate_shift_deg(&self, lat_deg: f64, lon_deg: f64) -> Option<(f64, f64)> {
        for subgrid in self.subgrids.iter().rev() {
            if let Some((dlat_sec, dlon_sec)) = subgrid.interpolate(lat_deg, lon_deg) {
                return Some((dlat_sec / 3600.0, dlon_sec / 3600.0));
            }
        }
        None
    }

    /// Shifts geodetic coordinates: (lat_tgt, lon_tgt) = (lat_src + dlat, lon_src + dlon).
    pub fn apply_forward(&self, lat_deg: f64, lon_deg: f64) -> Option<(f64, f64)> {
        let (dlat_deg, dlon_deg) = self.interpolate_shift_deg(lat_deg, lon_deg)?;
        Some((lat_deg + dlat_deg, lon_deg + dlon_deg))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ntv2_subgrid_bilinear_interpolation() {
        let nodes = alloc::vec![
            Ntv2Node { dlat_sec: 1.0, dlon_sec: -2.0, lat_acc_m: 0.1, lon_acc_m: 0.1 },
            Ntv2Node { dlat_sec: 1.5, dlon_sec: -2.2, lat_acc_m: 0.1, lon_acc_m: 0.1 },
            Ntv2Node { dlat_sec: 2.0, dlon_sec: -2.5, lat_acc_m: 0.1, lon_acc_m: 0.1 },
            Ntv2Node { dlat_sec: 2.5, dlon_sec: -2.7, lat_acc_m: 0.1, lon_acc_m: 0.1 },
        ];

        let subgrid = Ntv2Subgrid {
            name: "TEST".into(),
            parent: "NONE".into(),
            s_lat_deg: 40.0,
            n_lat_deg: 41.0,
            w_lon_deg: -120.0,
            e_lon_deg: -119.0,
            lat_inc_deg: 1.0,
            lon_inc_deg: 1.0,
            n_rows: 2,
            n_cols: 2,
            nodes,
        };

        let grid = Ntv2Grid { subgrids: alloc::vec![subgrid] };
        let (dlat, dlon) = grid.interpolate_shift_deg(40.5, -119.5).expect("interpolation");

        // Midpoint: dlat = 1.75 sec / 3600, dlon = -2.35 sec / 3600
        assert!((dlat - (1.75 / 3600.0)).abs() < 1e-9);
        assert!((dlon - (-2.35 / 3600.0)).abs() < 1e-9);
    }
}
