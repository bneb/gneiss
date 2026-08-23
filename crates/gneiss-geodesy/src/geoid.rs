//! Geoid undulation model interpolation.
//!
//! Converts between ellipsoidal height ($h$) derived from GNSS and orthometric
//! height ($H$, height above mean sea level) via the geoid separation ($N$):
//!
//!   H = h - N(\phi, \lambda)
//!   h = H + N(\phi, \lambda)

use alloc::vec::Vec;
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Regular 2D grid of geoid undulation values (e.g., EGM2008 / EGM96).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GeoidGrid {
    /// Minimum latitude in degrees (southern boundary).
    pub min_lat_deg: f64,
    /// Maximum latitude in degrees (northern boundary).
    pub max_lat_deg: f64,
    /// Minimum longitude in degrees (western boundary, typically 0.0 or -180.0).
    pub min_lon_deg: f64,
    /// Maximum longitude in degrees (eastern boundary, typically 360.0 or 180.0).
    pub max_lon_deg: f64,
    /// Grid spacing along latitude in degrees.
    pub dlat_deg: f64,
    /// Grid spacing along longitude in degrees.
    pub dlon_deg: f64,
    /// Number of grid points along latitude.
    pub n_lat: usize,
    /// Number of grid points along longitude.
    pub n_lon: usize,
    /// Row-major grid values in meters (size: n_lat * n_lon).
    pub grid: Vec<f32>,
}

impl GeoidGrid {
    /// Create a new GeoidGrid with validation against non-finite or degenerate bounds.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        min_lat_deg: f64,
        max_lat_deg: f64,
        min_lon_deg: f64,
        max_lon_deg: f64,
        dlat_deg: f64,
        dlon_deg: f64,
        n_lat: usize,
        n_lon: usize,
        grid: Vec<f32>,
    ) -> Option<Self> {
        if !min_lat_deg.is_finite()
            || !max_lat_deg.is_finite()
            || !min_lon_deg.is_finite()
            || !max_lon_deg.is_finite()
            || !dlat_deg.is_finite()
            || !dlon_deg.is_finite()
        {
            return None;
        }
        if max_lat_deg <= min_lat_deg || max_lon_deg <= min_lon_deg || dlat_deg <= 0.0 || dlon_deg <= 0.0 {
            return None;
        }
        if n_lat < 2 || n_lon < 2 || grid.len() != n_lat * n_lon {
            return None;
        }
        Some(Self {
            min_lat_deg,
            max_lat_deg,
            min_lon_deg,
            max_lon_deg,
            dlat_deg,
            dlon_deg,
            n_lat,
            n_lon,
            grid,
        })
    }

    /// Compute geoid undulation $N$ (meters) via bilinear interpolation.
    pub fn undulation(&self, lat_deg: f64, lon_deg: f64) -> Option<f64> {
        if !lat_deg.is_finite() || !lon_deg.is_finite() {
            return None;
        }
        let lon_norm = self.normalize_longitude(lon_deg);
        if lat_deg < self.min_lat_deg
            || lat_deg > self.max_lat_deg
            || lon_norm < self.min_lon_deg
            || lon_norm > self.max_lon_deg
        {
            return None;
        }

        let fi = (lat_deg - self.min_lat_deg) / self.dlat_deg;
        let fj = (lon_norm - self.min_lon_deg) / self.dlon_deg;

        let i0 = (fi.floor() as usize).min(self.n_lat - 2);
        let j0 = (fj.floor() as usize).min(self.n_lon - 2);

        let u = fi - (i0 as f64);
        let v = fj - (j0 as f64);

        let n00 = self.grid[i0 * self.n_lon + j0] as f64;
        let n01 = self.grid[i0 * self.n_lon + (j0 + 1)] as f64;
        let n10 = self.grid[(i0 + 1) * self.n_lon + j0] as f64;
        let n11 = self.grid[(i0 + 1) * self.n_lon + (j0 + 1)] as f64;

        let n = (1.0 - u) * (1.0 - v) * n00
            + (1.0 - u) * v * n01
            + u * (1.0 - v) * n10
            + u * v * n11;

        if n.is_finite() {
            Some(n)
        } else {
            None
        }
    }

    /// Convert Ellipsoidal Height (h) in LLH [lat_rad, lon_rad, h_m] to Orthometric Height (H).
    ///
    /// Returns Vector3 with [lat_rad, lon_rad, H_m].
    pub fn ellipsoidal_to_orthometric(&self, llh: Vector3<f64>) -> Option<Vector3<f64>> {
        if !llh.x.is_finite() || !llh.y.is_finite() || !llh.z.is_finite() {
            return None;
        }
        let lat_deg = llh.x.to_degrees();
        let lon_deg = llh.y.to_degrees();
        let n = self.undulation(lat_deg, lon_deg)?;
        Some(Vector3::new(llh.x, llh.y, llh.z - n))
    }

    /// Convert Orthometric Height (H) in LLH [lat_rad, lon_rad, H_m] to Ellipsoidal Height (h).
    ///
    /// Returns Vector3 with [lat_rad, lon_rad, h_m].
    pub fn orthometric_to_ellipsoidal(&self, llh_ortho: Vector3<f64>) -> Option<Vector3<f64>> {
        if !llh_ortho.x.is_finite() || !llh_ortho.y.is_finite() || !llh_ortho.z.is_finite() {
            return None;
        }
        let lat_deg = llh_ortho.x.to_degrees();
        let lon_deg = llh_ortho.y.to_degrees();
        let n = self.undulation(lat_deg, lon_deg)?;
        Some(Vector3::new(llh_ortho.x, llh_ortho.y, llh_ortho.z + n))
    }

    #[inline]
    fn normalize_longitude(&self, lon_deg: f64) -> f64 {
        if self.min_lon_deg >= 0.0 && lon_deg < 0.0 {
            lon_deg + 360.0
        } else if self.min_lon_deg < 0.0 && lon_deg > 180.0 {
            lon_deg - 360.0
        } else {
            lon_deg
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_sample_grid() -> GeoidGrid {
        let grid = alloc::vec![
            30.0, 31.0, 32.0,
            32.0, 33.0, 34.0,
            34.0, 35.0, 36.0,
        ];
        GeoidGrid::new(30.0, 32.0, 130.0, 132.0, 1.0, 1.0, 3, 3, grid).expect("valid grid")
    }

    #[test]
    fn test_geoid_grid_exact_points() {
        let grid = create_sample_grid();
        let n_corner = grid.undulation(30.0, 130.0).expect("undulation");
        assert!((n_corner - 30.0).abs() < 1e-6);

        let n_mid = grid.undulation(31.0, 131.0).expect("undulation");
        assert!((n_mid - 33.0).abs() < 1e-6);
    }

    #[test]
    fn test_geoid_grid_interpolation() {
        let grid = create_sample_grid();
        let n_interp = grid.undulation(30.5, 130.5).expect("undulation");
        assert!((n_interp - 31.5).abs() < 1e-6);
    }

    #[test]
    fn test_geoid_grid_out_of_bounds_and_nan() {
        let grid = create_sample_grid();
        assert!(grid.undulation(29.0, 130.0).is_none());
        assert!(grid.undulation(33.0, 130.0).is_none());
        assert!(grid.undulation(31.0, 129.0).is_none());
        assert!(grid.undulation(31.0, 133.0).is_none());
        assert!(grid.undulation(f64::NAN, 130.0).is_none());
        assert!(grid.undulation(31.0, f64::NAN).is_none());
        assert!(grid.undulation(f64::INFINITY, 130.0).is_none());
    }

    #[test]
    fn test_geoid_grid_invalid_constructor_params() {
        assert!(GeoidGrid::new(f64::NAN, 32.0, 130.0, 132.0, 1.0, 1.0, 3, 3, alloc::vec![0.0; 9]).is_none());
        assert!(GeoidGrid::new(32.0, 30.0, 130.0, 132.0, 1.0, 1.0, 3, 3, alloc::vec![0.0; 9]).is_none());
        assert!(GeoidGrid::new(30.0, 32.0, 130.0, 132.0, -1.0, 1.0, 3, 3, alloc::vec![0.0; 9]).is_none());
        assert!(GeoidGrid::new(30.0, 32.0, 130.0, 132.0, 1.0, 1.0, 1, 3, alloc::vec![0.0; 3]).is_none());
        assert!(GeoidGrid::new(30.0, 32.0, 130.0, 132.0, 1.0, 1.0, 3, 3, alloc::vec![0.0; 8]).is_none());
    }

    #[test]
    fn test_geoid_height_conversions() {
        let grid = create_sample_grid();
        let lat_rad = 31.0_f64.to_radians();
        let lon_rad = 131.0_f64.to_radians();
        let h_ellips = 100.0;

        let llh = Vector3::new(lat_rad, lon_rad, h_ellips);
        let ortho = grid.ellipsoidal_to_orthometric(llh).expect("ortho conversion");

        assert!((ortho.z - 67.0).abs() < 1e-6);
        assert_eq!(ortho.x, lat_rad);
        assert_eq!(ortho.y, lon_rad);

        let roundtrip = grid.orthometric_to_ellipsoidal(ortho).expect("roundtrip");
        assert!((roundtrip.z - h_ellips).abs() < 1e-6);

        let nan_llh = Vector3::new(f64::NAN, lon_rad, h_ellips);
        assert!(grid.ellipsoidal_to_orthometric(nan_llh).is_none());
        assert!(grid.orthometric_to_ellipsoidal(nan_llh).is_none());
    }

    #[test]
    fn test_geoid_longitude_normalization() {
        let grid_data = alloc::vec![10.0, 20.0, 10.0, 20.0];
        let grid = GeoidGrid::new(0.0, 10.0, 0.0, 360.0, 10.0, 360.0, 2, 2, grid_data).expect("grid");

        let n = grid.undulation(5.0, -10.0);
        assert!(n.is_some());
    }
}
