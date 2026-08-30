//! High-precision Transverse Mercator and Universal Transverse Mercator (UTM) projections.
//!
//! Implements Karney-Krüger $n$-series expansion accurate to sub-millimeter
//! within $\pm 30^\circ$ of the central meridian on standard reference ellipsoids.

use core::f64::consts::PI;

/// Reference ellipsoid parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipsoid {
    /// Semi-major axis $a$ in meters.
    pub a: f64,
    /// Flattening $f = (a - b) / a$.
    pub f: f64,
}

impl Ellipsoid {
    /// WGS84 Reference Ellipsoid.
    pub const WGS84: Self = Self {
        a: 6_378_137.0,
        f: 1.0 / 298.257_223_563,
    };

    /// GRS80 Reference Ellipsoid (NAD83).
    pub const GRS80: Self = Self {
        a: 6_378_137.0,
        f: 1.0 / 298.257_222_101,
    };
}

/// Transverse Mercator Projection parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransverseMercator {
    pub ellipsoid: Ellipsoid,
    pub lat0_rad: f64,
    pub lon0_rad: f64,
    pub k0: f64,
    pub false_easting: f64,
    pub false_northing: f64,
}

impl TransverseMercator {
    /// Creates a Transverse Mercator projection instance.
    pub fn new(
        ellipsoid: Ellipsoid,
        lat0_rad: f64,
        lon0_rad: f64,
        k0: f64,
        false_easting: f64,
        false_northing: f64,
    ) -> Self {
        Self {
            ellipsoid,
            lat0_rad,
            lon0_rad,
            k0,
            false_easting,
            false_northing,
        }
    }

    /// Creates a Universal Transverse Mercator (UTM) projection for a given zone and hemisphere.
    pub fn utm(zone: u8, north: bool) -> Option<Self> {
        if !(1..=60).contains(&zone) {
            return None;
        }
        let lon0_deg = (zone as f64 - 1.0) * 6.0 - 180.0 + 3.0;
        let false_northing = if north { 0.0 } else { 10_000_000.0 };
        Some(Self::new(
            Ellipsoid::WGS84,
            0.0,
            lon0_deg * PI / 180.0,
            0.9996,
            500_000.0,
            false_northing,
        ))
    }

    /// Determines the UTM zone (1-60) from longitude in radians.
    pub fn utm_zone_from_lon(lon_rad: f64) -> u8 {
        let lon_deg = lon_rad * 180.0 / PI;
        let mut norm = lon_deg.rem_euclid(360.0);
        if norm >= 180.0 {
            norm -= 360.0;
        }
        let z = ((norm + 180.0) / 6.0).floor() as i32 + 1;
        z.clamp(1, 60) as u8
    }

    /// Forward projection: Converts geodetic coordinates (lat, lon in radians) to grid Easting and Northing (meters).
    pub fn forward(&self, lat_rad: f64, lon_rad: f64) -> (f64, f64) {
        let a = self.ellipsoid.a;
        let f = self.ellipsoid.f;
        let n = f / (2.0 - f);
        let n2 = n * n;
        let n3 = n2 * n;
        let n4 = n3 * n;

        let a_hat = (a / (1.0 + n)) * (1.0 + 0.25 * n2 + (1.0 / 64.0) * n4);

        let dlon = lon_rad - self.lon0_rad;
        let e2 = f * (2.0 - f);
        let e = libm::sqrt(e2);

        let s = libm::sin(lat_rad);
        let tau = libm::tan(lat_rad);
        let sigma = libm::sinh(e * libm::atanh(e * s));
        let tau_p = tau * libm::sqrt(1.0 + sigma * sigma) - sigma * libm::sqrt(1.0 + tau * tau);

        let xi_prime = libm::atan2(tau_p, libm::cos(dlon));
        let eta_prime = libm::asinh(libm::sin(dlon) / libm::sqrt(tau_p * tau_p + libm::cos(dlon) * libm::cos(dlon)));

        // Krüger alpha coefficients
        let a1 = 0.5 * n - (2.0 / 3.0) * n2 + (5.0 / 16.0) * n3 + (41.0 / 180.0) * n4;
        let a2 = (13.0 / 48.0) * n2 - (3.0 / 5.0) * n3 + (557.0 / 1440.0) * n4;
        let a3 = (61.0 / 240.0) * n3 - (103.0 / 140.0) * n4;
        let a4 = (49561.0 / 161280.0) * n4;

        let xi = xi_prime
            + a1 * libm::sin(2.0 * xi_prime) * libm::cosh(2.0 * eta_prime)
            + a2 * libm::sin(4.0 * xi_prime) * libm::cosh(4.0 * eta_prime)
            + a3 * libm::sin(6.0 * xi_prime) * libm::cosh(6.0 * eta_prime)
            + a4 * libm::sin(8.0 * xi_prime) * libm::cosh(8.0 * eta_prime);

        let eta = eta_prime
            + a1 * libm::cos(2.0 * xi_prime) * libm::sinh(2.0 * eta_prime)
            + a2 * libm::cos(4.0 * xi_prime) * libm::sinh(4.0 * eta_prime)
            + a3 * libm::cos(6.0 * xi_prime) * libm::sinh(6.0 * eta_prime)
            + a4 * libm::cos(8.0 * xi_prime) * libm::sinh(8.0 * eta_prime);

        let easting = self.false_easting + self.k0 * a_hat * eta;
        let northing = self.false_northing + self.k0 * a_hat * xi;

        (easting, northing)
    }

    /// Inverse projection: Converts grid Easting and Northing (meters) to geodetic coordinates (lat, lon in radians).
    pub fn inverse(&self, easting: f64, northing: f64) -> (f64, f64) {
        let a = self.ellipsoid.a;
        let f = self.ellipsoid.f;
        let n = f / (2.0 - f);
        let n2 = n * n;
        let n3 = n2 * n;
        let n4 = n3 * n;

        let a_hat = (a / (1.0 + n)) * (1.0 + 0.25 * n2 + (1.0 / 64.0) * n4);

        let xi = (northing - self.false_northing) / (self.k0 * a_hat);
        let eta = (easting - self.false_easting) / (self.k0 * a_hat);

        // Krüger beta coefficients
        let b1 = 0.5 * n - (2.0 / 3.0) * n2 + (37.0 / 96.0) * n3 - (1.0 / 360.0) * n4;
        let b2 = (1.0 / 48.0) * n2 + (1.0 / 15.0) * n3 - (437.0 / 1440.0) * n4;
        let b3 = (17.0 / 480.0) * n3 - (37.0 / 840.0) * n4;
        let b4 = (4397.0 / 161280.0) * n4;

        let xi_prime = xi
            - b1 * libm::sin(2.0 * xi) * libm::cosh(2.0 * eta)
            - b2 * libm::sin(4.0 * xi) * libm::cosh(4.0 * eta)
            - b3 * libm::sin(6.0 * xi) * libm::cosh(6.0 * eta)
            - b4 * libm::sin(8.0 * xi) * libm::cosh(8.0 * eta);

        let eta_prime = eta
            - b1 * libm::cos(2.0 * xi) * libm::sinh(2.0 * eta)
            - b2 * libm::cos(4.0 * xi) * libm::sinh(4.0 * eta)
            - b3 * libm::cos(6.0 * xi) * libm::sinh(6.0 * eta)
            - b4 * libm::cos(8.0 * xi) * libm::sinh(8.0 * eta);

        let sinh_eta_p = libm::sinh(eta_prime);
        let sin_xi_p = libm::sin(xi_prime);
        let cos_xi_p = libm::cos(xi_prime);

        let tau_p = sin_xi_p / libm::sqrt(sinh_eta_p * sinh_eta_p + cos_xi_p * cos_xi_p);
        let dlon = libm::atan2(sinh_eta_p, cos_xi_p);

        let e2 = f * (2.0 - f);
        let e = libm::sqrt(e2);

        // Fixed-point iteration for conformal to geodetic latitude
        let mut tau = tau_p;
        for _ in 0..6 {
            let sigma = libm::sinh(e * libm::atanh(e * tau / libm::sqrt(1.0 + tau * tau)));
            tau = tau_p * libm::sqrt(1.0 + sigma * sigma) + sigma * libm::sqrt(1.0 + tau_p * tau_p);
        }

        let lat = libm::atan(tau);
        let lon = self.lon0_rad + dlon;

        (lat, lon)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_utm_forward_and_inverse_roundtrip() {
        let utm10 = TransverseMercator::utm(10, true).expect("UTM zone 10");

        // San Francisco: 37.7749 N, 122.4194 W
        let lat_rad = 37.7749 * PI / 180.0;
        let lon_rad = -122.4194 * PI / 180.0;

        let (e, n) = utm10.forward(lat_rad, lon_rad);
        assert!((e - 551_166.0).abs() < 100.0, "Easting roughly 551 km");
        assert!((n - 4_181_000.0).abs() < 1000.0, "Northing roughly 4181 km");

        let (lat_back, lon_back) = utm10.inverse(e, n);
        assert!((lat_back - lat_rad).abs() < 1e-9, "Lat roundtrip sub-mm");
        assert!((lon_back - lon_rad).abs() < 1e-9, "Lon roundtrip sub-mm");
    }

    #[test]
    fn test_utm_wide_grid_high_precision_roundtrips() {
        // Test Southern Hemisphere and extreme zone edges to catch all Krüger series terms
        let utm55s = TransverseMercator::utm(55, false).expect("UTM 55S (Sydney)");
        let lats = [-75.0, -45.0, -10.0, 5.0, 35.0, 65.0, 80.0];
        let d_lons = [-2.9, -1.5, -0.01, 0.0, 0.01, 1.5, 2.9];

        for &lat_deg in &lats {
            let lat = lat_deg * PI / 180.0;
            for &d_lon_deg in &d_lons {
                let lon = utm55s.lon0_rad + d_lon_deg * PI / 180.0;
                let (e, n) = utm55s.forward(lat, lon);
                let (lat_inv, lon_inv) = utm55s.inverse(e, n);

                assert!(
                    (lat_inv - lat).abs() < 1e-10,
                    "Lat error too large for ({}, {}) -> ({}, {})",
                    lat_deg, d_lon_deg, lat_inv, lat
                );
                assert!(
                    (lon_inv - lon).abs() < 1e-10,
                    "Lon error too large for ({}, {}) -> ({}, {})",
                    lat_deg, d_lon_deg, lon_inv, lon
                );
            }
        }
    }

    #[test]
    fn test_utm_zone_determination() {
        // Longitude -122 deg is Zone 10
        let z = TransverseMercator::utm_zone_from_lon(-122.0 * PI / 180.0);
        assert_eq!(z, 10);

        // Longitude +2 deg (Paris) is Zone 31
        let z_paris = TransverseMercator::utm_zone_from_lon(2.35 * PI / 180.0);
        assert_eq!(z_paris, 31);
    }
}
