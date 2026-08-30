//! Geoid undulation model interpolation and binary grid ingestion.
//!
//! Converts between ellipsoidal height ($h$) derived from GNSS and orthometric
//! height ($H$, height above mean sea level) via the geoid separation ($N$):
//!
//!   H = h - N(\phi, \lambda)
//!   h = H + N(\phi, \lambda)

pub mod byn;
pub mod gtx;
#[cfg(test)]
mod tests;

use alloc::vec::Vec;
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Regular 2D grid of geoid undulation values (e.g., EGM2008, GEOID18, CGG2013).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GeoidGrid {
    pub min_lat_deg: f64,
    pub max_lat_deg: f64,
    pub min_lon_deg: f64,
    pub max_lon_deg: f64,
    pub dlat_deg: f64,
    pub dlon_deg: f64,
    pub n_lat: usize,
    pub n_lon: usize,
    pub grid: Vec<f32>,
}

impl GeoidGrid {
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

    /// Ingest a NOAA VDatum / PROJ `.gtx` binary grid.
    pub fn from_gtx_bytes(bytes: &[u8]) -> Result<Self, &'static str> {
        gtx::parse_gtx(bytes)
    }

    /// Ingest an NRCan `.byn` binary grid.
    pub fn from_byn_bytes(bytes: &[u8]) -> Result<Self, &'static str> {
        byn::parse_byn(bytes)
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

    pub fn ellipsoidal_to_orthometric(&self, llh: Vector3<f64>) -> Option<Vector3<f64>> {
        if !llh.x.is_finite() || !llh.y.is_finite() || !llh.z.is_finite() {
            return None;
        }
        let lat_deg = llh.x.to_degrees();
        let lon_deg = llh.y.to_degrees();
        let n = self.undulation(lat_deg, lon_deg)?;
        Some(Vector3::new(llh.x, llh.y, llh.z - n))
    }

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
        if self.min_lon_deg >= 0.0 {
            if lon_deg < 0.0 {
                lon_deg + 360.0
            } else if lon_deg > 360.0 {
                lon_deg - 360.0
            } else {
                lon_deg
            }
        } else if self.min_lon_deg < 0.0 && lon_deg > 180.0 {
            lon_deg - 360.0
        } else {
            lon_deg
        }
    }
}
