#[cfg(test)]
mod tests {
    use super::super::coords::ecef_to_llh;
    use nalgebra::Vector3;
    

    #[test]
    fn test_ecef_to_llh_precision() {
        // WGS84 ECEF: [crate::constants::WGS84_SEMI_MAJOR_AXIS_M, 0, 0] (Equator, Lon 0, Height 0)
        let pos = Vector3::new(crate::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let llh = ecef_to_llh(pos);
        assert!(llh.x.abs() < 1e-12); // Lat 0
        assert!(llh.y.abs() < 1e-12); // Lon 0
        assert!(llh.z.abs() < 1e-4);  // Height 0
    }

    #[test]
    fn test_ecef_llh_round_trip_poles() {
        use super::super::coords::{ecef_to_llh, llh_to_ecef};
        use nalgebra::Vector3;
        use core::f64::consts::FRAC_PI_2;

        // North pole at sea level
        let north_pole_llh = Vector3::new(FRAC_PI_2, 0.0, 0.0);
        let ecef = llh_to_ecef(north_pole_llh);
        let llh_back = ecef_to_llh(ecef);
        assert!((llh_back.x - FRAC_PI_2).abs() < 1e-10, "Latitude error at North Pole: {}", (llh_back.x - FRAC_PI_2).abs());
        assert!(llh_back.z.abs() < 1e-3, "Height error at North Pole: {}", llh_back.z.abs());

        // South pole at sea level
        let south_pole_llh = Vector3::new(-FRAC_PI_2, 0.0, 0.0);
        let ecef = llh_to_ecef(south_pole_llh);
        let llh_back = ecef_to_llh(ecef);
        assert!((llh_back.x - (-FRAC_PI_2)).abs() < 1e-10, "Latitude error at South Pole: {}", (llh_back.x + FRAC_PI_2).abs());
        assert!(llh_back.z.abs() < 1e-3, "Height error at South Pole: {}", llh_back.z.abs());

        // North pole at 1000m altitude
        let north_pole_high = Vector3::new(FRAC_PI_2, 0.0, 1000.0);
        let ecef = llh_to_ecef(north_pole_high);
        let llh_back = ecef_to_llh(ecef);
        assert!((llh_back.x - FRAC_PI_2).abs() < 1e-10);
        assert!((llh_back.z - 1000.0).abs() < 1e-3);
    }

    #[test]
    fn test_ecef_llh_round_trip_high_altitude() {
        use super::super::coords::{ecef_to_llh, llh_to_ecef};
        use nalgebra::Vector3;

        // Geostationary orbit altitude (~35786 km)
        let geo_llh = Vector3::new(0.0, 0.0, 35786000.0);
        let ecef = llh_to_ecef(geo_llh);
        let llh_back = ecef_to_llh(ecef);
        assert!(llh_back.x.abs() < 1e-10, "Lat error at GEO: {}", llh_back.x);
        assert!(llh_back.y.abs() < 1e-10, "Lon error at GEO: {}", llh_back.y);
        assert!((llh_back.z - 35786000.0).abs() < 1e-3, "Height error at GEO: {}", (llh_back.z - 35786000.0).abs());

        // GPS orbit altitude (~20200 km) at 45° latitude
        let gps_llh = Vector3::new(core::f64::consts::FRAC_PI_4, core::f64::consts::FRAC_PI_4, 20200000.0);
        let ecef = llh_to_ecef(gps_llh);
        let llh_back = ecef_to_llh(ecef);
        assert!((llh_back.x - core::f64::consts::FRAC_PI_4).abs() < 1e-10);
        assert!((llh_back.y - core::f64::consts::FRAC_PI_4).abs() < 1e-10);
        assert!((llh_back.z - 20200000.0).abs() < 1e-3);
    }

    #[test]
    fn test_ecef_llh_round_trip_parametric() {
        use super::super::coords::{ecef_to_llh, llh_to_ecef};
        use nalgebra::Vector3;

        // Test 20 points spanning the full range of lat/lon/height
        let test_points = [
            (0.0, 0.0, 0.0),
            (0.5, 1.0, 100.0),
            (-0.5, -1.0, 500.0),
            (1.0, 3.0, 0.0),
            (-1.0, -3.0, 0.0),
            (0.1, 0.1, 8848.0),    // Everest altitude
            (-0.6, 2.5, -50.0),    // Dead Sea (below sea level)
            (0.8, -2.0, 10000.0),
            (1.5, 0.5, 0.0),       // Near pole
            (-1.5, -0.5, 0.0),     // Near south pole
            (0.3, -1.5, 35000.0),  // Aircraft altitude
            (-0.3, 1.5, 400000.0), // ISS altitude
            (0.7, 0.0, 100.0),
            (0.0, core::f64::consts::PI, 0.0),      // Near date line
            (0.0, -core::f64::consts::PI, 0.0),     // Near date line (other side)
            (0.01, 0.01, 0.0),     // Near equator/prime meridian
            (-0.01, -0.01, 0.0),
            (1.4, 2.0, 5000.0),    // High latitude
            (-1.4, -2.0, 5000.0),
            (0.785, 1.571, 1000.0),// 45deg lat, 90deg lon
        ];

        for (lat, lon, height) in test_points {
            let llh = Vector3::new(lat, lon, height);
            let ecef = llh_to_ecef(llh);
            let llh_back = ecef_to_llh(ecef);

            assert!((llh_back.x - lat).abs() < 1e-9,
                "Lat round-trip failed for ({}, {}, {}): error = {}", lat, lon, height, (llh_back.x - lat).abs());

            // Handle longitude wrap-around for ±π
            let lon_err = ((llh_back.y - lon) + core::f64::consts::PI).rem_euclid(2.0 * core::f64::consts::PI) - core::f64::consts::PI;
            assert!(lon_err.abs() < 1e-9,
                "Lon round-trip failed for ({}, {}, {}): error = {}", lat, lon, height, lon_err.abs());

            assert!((llh_back.z - height).abs() < 1e-3,
                "Height round-trip failed for ({}, {}, {}): error = {}m", lat, lon, height, (llh_back.z - height).abs());
        }
    }
}
