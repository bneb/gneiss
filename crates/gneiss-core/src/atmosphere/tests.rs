#![allow(clippy::unwrap_used)]

use super::*;
use super::mapping::*;
use alloc::vec;

    #[test]
    fn test_tropo_delay() {
        let params = TropoParams::default();
        let delay = AtmosphereModel::tropo_saastamoinen(&params, 0.5, 100.0);
        assert!(
            (delay - 4.94).abs() < 0.05,
            "Delay should be approximately 4.94m, got {}",
            delay
        );
    }

    /// Bug 6 regression test: the GMF mapping function uses fully normalized
    /// associated Legendre functions.  Verify a few known values.
    ///
    /// For P̄_{0,0}(x) the normalization gives sqrt(1) * 1 = 1.
    /// For P̄_{1,0}(x) the norm is sqrt(3) and P_{1,0}(x) = x → P̄_{1,0}(x) = sqrt(3)*x.
    /// For P̄_{1,1}(x) the norm is sqrt(3) and P_{1,1}(x) = sqrt(1-x²) → P̄_{1,1}(x) = sqrt(3*(1-x²)).
    /// For P̄_{2,0}(x) the norm is sqrt(5) and P_{2,0}(x) = (3x²-1)/2 → P̄_{2,0} = sqrt(5)*(3x²-1)/2.
    #[test]
    fn test_legendre_normalization() {
        let x = 0.5_f64; // sin(lat) for lat = 30°

        // P̄_{0,0}(x) = 1
        let p00 = _legendre_norm(0, 0, x);
        assert!((p00 - 1.0).abs() < 1e-12, "P̄_00 = 1, got {p00}");

        // P̄_{1,0}(x) = sqrt(3) * x
        let p10 = _legendre_norm(1, 0, x);
        let expected_p10 = (3.0_f64).sqrt() * x;
        assert!(
            (p10 - expected_p10).abs() < 1e-12,
            "P̄_10 = sqrt(3)*x = {expected_p10}, got {p10}"
        );

        // P̄_{1,1}(x) = sqrt(3) * sqrt(1-x²)
        let p11 = _legendre_norm(1, 1, x);
        let expected_p11 = (3.0_f64).sqrt() * (1.0 - x * x).sqrt();
        assert!(
            (p11 - expected_p11).abs() < 1e-12,
            "P̄_11 = sqrt(3*(1-x²)) = {expected_p11}, got {p11}"
        );

        // P̄_{2,0}(x) = sqrt(5) * (3x²-1)/2
        let p20 = _legendre_norm(2, 0, x);
        let expected_p20 = (5.0_f64).sqrt() * (3.0 * x * x - 1.0) / 2.0;
        assert!(
            (p20 - expected_p20).abs() < 1e-12,
            "P̄_20 = sqrt(5)*(3x²-1)/2 = {expected_p20}, got {p20}"
        );
    }

    /// Bug 5 regression test: GMF must produce different mapping factors for
    /// different longitudes (confirming the cos(m*lon) term is active).
    #[test]
    fn test_gmf_longitude_variation() {
        let t = GpsTime::new(2000, 100000.0);
        let lat = 0.6_f64; // ~34°N
        let el = 0.3_f64; // ~17° elevation
        let h = 100.0_f64;

        // Same position but different longitudes
        let pos_lon0 = Vector3::new(lat, 0.0, h);
        let pos_lon90 = Vector3::new(lat, core::f64::consts::FRAC_PI_2, h);
        let pos_lon180 = Vector3::new(lat, core::f64::consts::PI, h);

        let (mh0, mw0) = gmf_impl(pos_lon0, el, t);
        let (mh90, mw90) = gmf_impl(pos_lon90, el, t);
        let (mh180, mw180) = gmf_impl(pos_lon180, el, t);

        // The mapping factors must vary with longitude (spherical harmonic terms include cos(m*lon))
        // m=0 terms are longitude-independent but m≥1 terms are not.
        let h_range = (mh0 - mh90).abs().max((mh0 - mh180).abs());
        let w_range = (mw0 - mw90).abs().max((mw0 - mw180).abs());
        assert!(
            h_range > 1e-6,
            "GMF dry mapping factor must vary with longitude (h_range={h_range})"
        );
        assert!(
            w_range > 1e-8,
            "GMF wet mapping factor must vary with longitude (w_range={w_range})"
        );

        // All values must be > 1 (mapping factors are always ≥ 1 in the valid range)
        assert!(mh0 > 1.0, "GMF m_h must be > 1, got {mh0}");
        assert!(mw0 > 1.0, "GMF m_w must be > 1, got {mw0}");
    }

    /// Bug 23 regression test: Klobuchar model must be evaluated at the IPP,
    /// not at the receiver position.  The fix activates the azimuth parameter,
    /// so delays at different azimuths (but same elevation) must differ.
    #[test]
    fn test_klobuchar_ipp_uses_azimuth() {
        use crate::atmosphere::{AtmosphereModel, KlobucharParams};
        let params = KlobucharParams {
            alpha: [3.82e-8, 1.49e-8, -1.79e-7, 0.0],
            beta: [1.43e5, 0.0, -3.28e5, 1.13e5],
        };
        let pos_llh = Vector3::new(0.6, 0.3, 100.0); // ~34°N, ~17°E
        let el = 0.4_f64; // ~23° elevation
        let t = GpsTime::new(2000, 50000.0);

        // Azimuth north vs. south — IPP moves in opposite latitude directions,
        // so the Klobuchar geomagnetic latitude and hence the delay differ.
        let delay_north = AtmosphereModel::iono_klobuchar(&params, pos_llh, 0.0, el, t);
        let delay_south =
            AtmosphereModel::iono_klobuchar(&params, pos_llh, core::f64::consts::PI, el, t);

        // The delays must differ because the IPP geomagnetic latitude differs.
        assert!(
            (delay_north - delay_south).abs() > 1e-4,
            "Klobuchar delay must differ for opposite azimuths (north={delay_north:.6}, south={delay_south:.6})"
        );

        // Both delays must be non-negative (Klobuchar is always ≥ 0)
        assert!(
            delay_north >= 0.0,
            "Klobuchar delay must be ≥ 0, got {delay_north}"
        );
        assert!(
            delay_south >= 0.0,
            "Klobuchar delay must be ≥ 0, got {delay_south}"
        );
    }

    fn make_test_tec_grid() -> Vec<Vec<f64>> {
        let mut grid = Vec::with_capacity(3);
        for i in 0..3 {
            let lat = 30.0 + (i as f64) * 5.0;
            let mut row = Vec::with_capacity(4);
            for j in 0..4 {
                let lon = (j as f64 - 1.0) * 5.0;
                row.push(lat + lon / 10.0);
            }
            grid.push(row);
        }
        grid
    }

    #[test]
    fn test_ionex_bilinear_grid_center() {
        let grid = make_test_tec_grid();
        let tec_maps: Vec<(GpsTime, &Vec<Vec<f64>>)> = vec![(GpsTime::new(2000, 0.0), &grid)];
        let delay = AtmosphereModel::iono_ionex(
            &tec_maps, 30.0, 40.0, 5.0, -5.0, 10.0, 5.0, 350.0,
            Vector3::new(35.0_f64.to_radians(), 7.5_f64.to_radians(), 100.0),
            0.0, core::f64::consts::FRAC_PI_2,
            GpsTime::new(2000, 0.0),
        );
        assert!((delay - 0.1624 * 35.75).abs() < 0.001);
    }

    #[test]
    fn test_ionex_elevation_mapping() {
        let grid = make_test_tec_grid();
        let tec_maps: Vec<(GpsTime, &Vec<Vec<f64>>)> = vec![(GpsTime::new(2000, 0.0), &grid)];
        let d90 = AtmosphereModel::iono_ionex(
            &tec_maps, 30.0, 40.0, 5.0, -5.0, 10.0, 5.0, 350.0,
            Vector3::new(30.0_f64.to_radians(), 0.0_f64.to_radians(), 100.0),
            0.0, core::f64::consts::FRAC_PI_2,
            GpsTime::new(2000, 0.0),
        );
        let d45 = AtmosphereModel::iono_ionex(
            &tec_maps, 30.0, 40.0, 5.0, -5.0, 10.0, 5.0, 350.0,
            Vector3::new(30.0_f64.to_radians(), 0.0_f64.to_radians(), 100.0),
            0.0, 45.0_f64.to_radians(),
            GpsTime::new(2000, 0.0),
        );
        let ratio = d45 / d90.max(1e-9);
        assert!(ratio > 1.3 && ratio < 1.5, "MF ratio should be ~1.4, got {:.3}", ratio);
    }

    #[test]
    fn test_ionex_temporal_interpolation() {
        let grid1 = make_test_tec_grid();
        let mut grid2 = make_test_tec_grid();
        for row in &mut grid2 { for v in row { *v += 10.0; } }
        let tec_maps: Vec<(GpsTime, &Vec<Vec<f64>>)> = vec![
            (GpsTime::new(2000, 0.0), &grid1),
            (GpsTime::new(2000, 7200.0), &grid2),
        ];
        let d1 = AtmosphereModel::iono_ionex(
            &tec_maps, 30.0, 40.0, 5.0, -5.0, 10.0, 5.0, 350.0,
            Vector3::new(30.0_f64.to_radians(), 0.0_f64.to_radians(), 100.0),
            0.0, core::f64::consts::FRAC_PI_2,
            GpsTime::new(2000, 0.0),
        );
        let d2 = AtmosphereModel::iono_ionex(
            &tec_maps, 30.0, 40.0, 5.0, -5.0, 10.0, 5.0, 350.0,
            Vector3::new(30.0_f64.to_radians(), 0.0_f64.to_radians(), 100.0),
            0.0, core::f64::consts::FRAC_PI_2,
            GpsTime::new(2000, 7200.0),
        );
        let dmid = AtmosphereModel::iono_ionex(
            &tec_maps, 30.0, 40.0, 5.0, -5.0, 10.0, 5.0, 350.0,
            Vector3::new(30.0_f64.to_radians(), 0.0_f64.to_radians(), 100.0),
            0.0, core::f64::consts::FRAC_PI_2,
            GpsTime::new(2000, 3600.0),
        );
        assert!((dmid - (d1 + d2) / 2.0).abs() < 0.001);
    }

    #[test]
    fn test_ionex_zero_elevation_returns_zero() {
        let grid = make_test_tec_grid();
        let tec_maps = vec![(GpsTime::new(2000, 0.0), &grid)];
        let delay = AtmosphereModel::iono_ionex(
            &tec_maps, 30.0, 40.0, 5.0, -5.0, 10.0, 5.0, 350.0,
            Vector3::new(35.0_f64.to_radians(), 0.0_f64.to_radians(), 100.0),
            0.0, 0.0,
            GpsTime::new(2000, 0.0),
        );
        assert_eq!(delay, 0.0);
    }

    #[test]
    fn test_ionex_empty_maps_returns_zero() {
        let tec_maps: Vec<(GpsTime, &Vec<Vec<f64>>)> = vec![];
        let delay = AtmosphereModel::iono_ionex(
            &tec_maps, 30.0, 40.0, 5.0, -5.0, 10.0, 5.0, 350.0,
            Vector3::new(35.0_f64.to_radians(), 0.0_f64.to_radians(), 100.0),
            0.0, core::f64::consts::FRAC_PI_2,
            GpsTime::new(2000, 0.0),
        );
        assert_eq!(delay, 0.0);
    }

