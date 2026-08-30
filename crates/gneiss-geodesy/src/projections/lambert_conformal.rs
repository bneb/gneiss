//! Lambert Conformal Conic (LCC) 2-Parallel Projection.
//!
//! Standard mapping projection for aeronautical charts, national topographic mapping,
//! and US State Plane Coordinate Systems (SPCS).

use core::f64::consts::PI;
use super::transverse_mercator::Ellipsoid;

/// Lambert Conformal Conic Projection parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LambertConformalConic {
    pub ellipsoid: Ellipsoid,
    pub lat1_rad: f64,
    pub lat2_rad: f64,
    pub lat0_rad: f64,
    pub lon0_rad: f64,
    pub false_easting: f64,
    pub false_northing: f64,
    n: f64,
    f_const: f64,
    rho0: f64,
}

impl LambertConformalConic {
    /// Creates a 2-parallel Lambert Conformal Conic projection.
    pub fn new(
        ellipsoid: Ellipsoid,
        lat1_rad: f64,
        lat2_rad: f64,
        lat0_rad: f64,
        lon0_rad: f64,
        false_easting: f64,
        false_northing: f64,
    ) -> Option<Self> {
        let a = ellipsoid.a;
        let f = ellipsoid.f;
        let e2 = f * (2.0 - f);
        let e = libm::sqrt(e2);

        let m1 = libm::cos(lat1_rad) / libm::sqrt(1.0 - e2 * libm::sin(lat1_rad) * libm::sin(lat1_rad));
        let m2 = libm::cos(lat2_rad) / libm::sqrt(1.0 - e2 * libm::sin(lat2_rad) * libm::sin(lat2_rad));

        let t0 = Self::calc_t(lat0_rad, e);
        let t1 = Self::calc_t(lat1_rad, e);
        let t2 = Self::calc_t(lat2_rad, e);

        let n = if (lat1_rad - lat2_rad).abs() < 1e-10 {
            libm::sin(lat1_rad)
        } else {
            (libm::log(m1) - libm::log(m2)) / (libm::log(t1) - libm::log(t2))
        };

        if n.abs() < 1e-12 {
            return None;
        }

        let f_const = m1 / (n * libm::pow(t1, n));
        let rho0 = a * f_const * libm::pow(t0, n);

        Some(Self {
            ellipsoid,
            lat1_rad,
            lat2_rad,
            lat0_rad,
            lon0_rad,
            false_easting,
            false_northing,
            n,
            f_const,
            rho0,
        })
    }

    fn calc_t(lat: f64, e: f64) -> f64 {
        let s = libm::sin(lat);
        let es = e * s;
        libm::tan(PI / 4.0 - lat / 2.0) / libm::pow((1.0 - es) / (1.0 + es), e / 2.0)
    }

    /// Forward projection: geodetic (lat, lon in radians) -> Easting, Northing (meters).
    pub fn forward(&self, lat_rad: f64, lon_rad: f64) -> (f64, f64) {
        let a = self.ellipsoid.a;
        let e = libm::sqrt(self.ellipsoid.f * (2.0 - self.ellipsoid.f));
        let t = Self::calc_t(lat_rad, e);
        let rho = a * self.f_const * libm::pow(t, self.n);
        let theta = self.n * (lon_rad - self.lon0_rad);

        let easting = self.false_easting + rho * libm::sin(theta);
        let northing = self.false_northing + self.rho0 - rho * libm::cos(theta);

        (easting, northing)
    }

    /// Inverse projection: Easting, Northing (meters) -> geodetic (lat, lon in radians).
    pub fn inverse(&self, easting: f64, northing: f64) -> (f64, f64) {
        let a = self.ellipsoid.a;
        let e = libm::sqrt(self.ellipsoid.f * (2.0 - self.ellipsoid.f));
        let x_p = easting - self.false_easting;
        let y_p = self.rho0 - (northing - self.false_northing);

        let rho = libm::sqrt(x_p * x_p + y_p * y_p) * libm::copysign(1.0, self.n);
        let theta = libm::atan2(x_p, y_p);

        let t = libm::pow(rho / (a * self.f_const), 1.0 / self.n);
        let lon = self.lon0_rad + theta / self.n;

        let mut lat = PI / 2.0 - 2.0 * libm::atan(t);
        for _ in 0..6 {
            let s = libm::sin(lat);
            let es = e * s;
            let factor = libm::pow((1.0 - es) / (1.0 + es), e / 2.0);
            lat = PI / 2.0 - 2.0 * libm::atan(t * factor);
        }

        (lat, lon)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lambert_conformal_conic_roundtrip() {
        // France Lambert-93 parameters: Lat1=44, Lat2=49, Lat0=46.5, Lon0=3.0, X0=700000, Y0=6600000
        let lcc = LambertConformalConic::new(
            Ellipsoid::GRS80,
            44.0 * PI / 180.0,
            49.0 * PI / 180.0,
            46.5 * PI / 180.0,
            3.0 * PI / 180.0,
            700_000.0,
            6_600_000.0,
        ).expect("LCC projection");

        let lat_rad = 48.8566 * PI / 180.0; // Paris
        let lon_rad = 2.3522 * PI / 180.0;

        let (e, n) = lcc.forward(lat_rad, lon_rad);
        assert!((e - 652_473.0).abs() < 100.0, "Easting ~ 652.4 km");
        assert!((n - 6_862_076.0).abs() < 100.0, "Northing ~ 6862.0 km");

        let (lat_back, lon_back) = lcc.inverse(e, n);
        assert!((lat_back - lat_rad).abs() < 1e-9, "Lat roundtrip sub-mm");
        assert!((lon_back - lon_rad).abs() < 1e-9, "Lon roundtrip sub-mm");
    }

    #[test]
    fn test_lcc_wide_grid_roundtrip() {
        let lcc = LambertConformalConic::new(
            Ellipsoid::GRS80,
            33.0 * PI / 180.0,
            45.0 * PI / 180.0,
            39.0 * PI / 180.0,
            -96.0 * PI / 180.0,
            0.0,
            0.0,
        ).expect("LCC US continental");

        let lats = [30.0, 35.0, 40.0, 48.0];
        let lons = [-110.0, -100.0, -96.0, -90.0, -80.0];

        for &lat_deg in &lats {
            let lat = lat_deg * PI / 180.0;
            for &lon_deg in &lons {
                let lon = lon_deg * PI / 180.0;
                let (e, n) = lcc.forward(lat, lon);
                let (lat_inv, lon_inv) = lcc.inverse(e, n);

                assert!(
                    (lat_inv - lat).abs() < 1e-10,
                    "Lat error for ({}, {}) -> ({}, {})",
                    lat_deg, lon_deg, lat_inv, lat
                );
                assert!(
                    (lon_inv - lon).abs() < 1e-10,
                    "Lon error for ({}, {}) -> ({}, {})",
                    lat_deg, lon_deg, lon_inv, lon
                );
            }
        }
    }
}
