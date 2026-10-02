//! Double-difference corrections, the geometry path and the RO-AZI grid.

#![allow(clippy::unwrap_used)]

use super::*;

    #[test]
    fn test_azimuth_elevation_2d_interpolation() {
        let azi_data = vec![
            vec![1.0, 2.0, 3.0], // Azimuth 0 deg: zen 0, 10, 20
            vec![2.0, 4.0, 6.0], // Azimuth 90 deg
            vec![3.0, 6.0, 9.0], // Azimuth 180 deg
            vec![2.0, 4.0, 6.0], // Azimuth 270 deg
        ];
        let pcv = ReceiverPcv {
            ant_type: "AZI_TEST".into(),
            radome: "NONE".into(),
            pco_neu_mm: nalgebra::Vector3::zeros(),
            pcv_grid_mm: vec![1.0, 2.0, 3.0],
            azi_grid_mm: Some(azi_data),
            dazi_deg: 90.0,
            zen_start_deg: 0.0,
            zen_step_deg: 10.0,
        };

        // Exact node: az 0, zen 0 -> 1.0
        assert_near(pcv.interpolate_az_zen(0.0, 0.0), 1.0, TOL, "node 0,0");
        // Exact node: az 90, zen 10 -> 4.0
        assert_near(pcv.interpolate_az_zen(90.0, 10.0), 4.0, TOL, "node 90,10");
        // Bilinear interpolation: az 45, zen 5 -> 0.5 * (1.5) + 0.5 * (3.0) = 2.25
        assert_near(pcv.interpolate_az_zen(45.0, 5.0), 2.25, TOL, "bilinear midpoint");
        // Azimuth wrap: az 360 -> az 0
        assert_near(pcv.interpolate_az_zen(360.0, 10.0), 2.0, TOL, "wrap 360");
    }

    /// Build a synthetic calibration with an explicit azimuth grid.
    fn synth_azi_pcv(azi: Vec<Vec<f64>>, dazi: f64, dzen: f64, noazi: Vec<f64>) -> ReceiverPcv {
        ReceiverPcv {
            ant_type: "AZI".into(),
            radome: "NONE".into(),
            pco_neu_mm: nalgebra::Vector3::zeros(),
            pcv_grid_mm: noazi,
            azi_grid_mm: Some(azi),
            dazi_deg: dazi,
            zen_start_deg: 0.0,
            zen_step_deg: dzen,
        }
    }

    // --- Requirement: the geometry path uses the station ECEF position. ---

    /// `dd_correction_from_geometry` used to hand the *geodetic* station
    /// vector to `az_el`'s `pos_ecef` slot, so every line of sight was
    /// taken from the wrong origin. Reference case, hand-derived:
    ///
    /// Station at lat 45 deg, lon 0 deg, height 0 -> ECEF (4.5e6, 0, 4.5e6).
    /// At lon 0 / lat 45 the local up is (cos45, 0, sin45) and local
    /// north is (-sin45, 0, cos45), so
    ///   * an overhead satellite sits at zen 0 deg,
    ///   * a satellite 20 000 km due north sits exactly on the horizon,
    ///     zen 90 deg.
    ///
    /// Rover PCV rises 1 mm per 5 deg of zenith (0, 1, 2, ... 18 mm over
    /// 0..90 deg, so the table reaches the horizon without clamping) and
    /// the base PCV is flat, so the double difference is
    ///   (0 - 18) - (0 - 0) = -18 mm = -0.018 m.
    fn linear_rover() -> ReceiverPcv {
        synth_pcv((0..19).map(|i| i as f64).collect(), 0.0, 5.0)
    }

    #[test]
    fn geometry_path_uses_the_station_ecef_position() {
        use gneiss_core::coords::llh_to_ecef;
        let llh = nalgebra::Vector3::new(45f64.to_radians(), 0.0, 0.0);
        let ecef = llh_to_ecef(llh);
        let (up, north) = (
            nalgebra::Vector3::new(45f64.to_radians().cos(), 0.0, 45f64.to_radians().sin()),
            nalgebra::Vector3::new(-45f64.to_radians().sin(), 0.0, 45f64.to_radians().cos()),
        );
        let overhead = ecef + up * 2.0e7;
        let horizon = ecef + north * 2.0e7;

        let rov = linear_rover();
        let bas = synth_pcv(vec![0.0; 19], 0.0, 5.0);
        // The elevation-only model over the same zeniths is the reference.
        let want = ReceiverPcv::dd_correction_m(&rov, &bas, 0.0, 90.0, 0.0, 90.0);
        assert_near(want, -0.018, TOL, "hand-derived DD (m)");
        let got = ReceiverPcv::dd_correction_from_geometry(
            &rov, &bas, llh, llh, overhead, horizon,
        );
        // Feeding the wrong origin instead put the horizon look-up at
        // atan(6378137 / 2.0e7) = 17.70 deg of elevation, i.e. zen 72.30,
        // which on this table is 14.46 mm rather than 18 mm.
        assert_near(got, -0.018, 1e-6, "geometry DD must match the same zeniths");
    }

    #[test]
    fn geometry_path_agrees_with_explicit_zeniths_on_real_antennas() {
        use gneiss_core::coords::llh_to_ecef;
        let Some(db) = load_db() else { return };
        let rov = ReceiverPcv::from_antex(&db, "ASH701945B_M", "SCIT").expect("ASH");
        let bas = ReceiverPcv::from_antex(&db, "LEIAR20", "LEIM").expect("LEIAR");
        let llh = nalgebra::Vector3::new(47.0f64.to_radians(), 8.0f64.to_radians(), 500.0);
        let ecef = llh_to_ecef(llh);
        let lat = 47.0f64.to_radians();
        let lon = 8.0f64.to_radians();
        let up = nalgebra::Vector3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin());
        let east = nalgebra::Vector3::new(-lon.sin(), lon.cos(), 0.0);
        let sat = ecef + up * 2.1e7 + east * 1.0e6;
        let r#ref_sat = ecef + up * 2.1e7 - east * 1.0e6;
        let geo = ReceiverPcv::dd_correction_from_geometry(&rov, &bas, llh, llh, sat, r#ref_sat);
        // Same pair, zeniths computed the long way round with the correct
        // ENU projection; both satellites share an elevation, so the
        // explicit form must give the same double difference.
        let d = sat - ecef;
        let horiz = d.x.hypot(d.y);
        let zen = (90.0 - d.z.atan2(horiz).to_degrees()).max(0.0);
        let explicit = ReceiverPcv::dd_correction_m(&rov, &bas, zen, zen, zen, zen);
        assert!(explicit.abs() < 1e-12, "symmetric pair cancels");
        assert_near(geo, explicit, 1e-9, "geometry vs explicit zeniths");
    }

    // --- Requirement: the 2D grid cannot be indexed out of bounds. ---

    #[test]
    fn ragged_azimuth_grid_falls_back_to_noazi_instead_of_panicking() {
        // Row 1 is one column short of row 0. Indexing `azi[az_next][zen_i]`
        // with a zenith chosen from row 0's length reads past row 1.
        let ragged = synth_azi_pcv(
            vec![vec![10.0, 99.0, 99.0], vec![20.0], vec![30.0], vec![20.0]],
            90.0,
            10.0,
            vec![1.0, 2.0, 3.0],
        );
        // Fallback is the NOAZI table: 2 + 0.5 * (3 - 2) = 2.5 mm at 15 deg.
        assert_near(ragged.interpolate_az_zen(45.0, 15.0), 2.5, TOL, "ragged grid falls back");
        // A well-formed grid of the same shape still uses 2D: at az 45
        // (rows 0,1) and zen 15 (columns 1,2) that is the mean of the
        // four corner nodes, 0.25 * (11 + 12 + 21 + 22) = 16.5 mm.
        let square = synth_azi_pcv(
            vec![
                vec![10.0, 11.0, 12.0],
                vec![20.0, 21.0, 22.0],
                vec![30.0, 31.0, 32.0],
                vec![20.0, 21.0, 22.0],
            ],
            90.0,
            10.0,
            vec![1.0, 2.0, 3.0],
        );
        assert_near(square.interpolate_az_zen(45.0, 15.0), 16.5, TOL, "square grid stays 2D");
    }

    #[test]
    fn empty_or_degenerate_grids_never_panic() {
        let empty_rows = synth_azi_pcv(Vec::new(), 90.0, 10.0, vec![4.0, 5.0]);
        assert_near(empty_rows.interpolate_az_zen(10.0, 5.0), 4.5, TOL, "no rows -> NOAZI");
        let empty_cols = synth_azi_pcv(vec![Vec::new(); 4], 90.0, 10.0, vec![4.0, 5.0]);
        assert_near(empty_cols.interpolate_az_zen(10.0, 5.0), 4.5, TOL, "no columns -> NOAZI");
        let no_dazi = synth_azi_pcv(vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0]], 0.0, 10.0, vec![7.0, 8.0]);
        assert_near(no_dazi.interpolate_az_zen(10.0, 5.0), 7.5, TOL, "dazi 0 -> NOAZI");
        let degenerate = synth_pcv(Vec::new(), 0.0, 5.0);
        assert_eq!(degenerate.interpolate(30.0), 0.0, "empty NOAZI -> 0");
    }

    #[test]
    fn az_zen_grid_reproduces_every_tabulated_node() {
        // d(azi)/d(zen) = 1 mm/deg^2 on both axes, so node (k, j) of the
        // grid is exactly 10*k + j. Sampling must return it bit-exactly,
        // which is what the bilinear weights reduce to at frac == 0.
        let (n_azi, n_zen, dazi, dzen) = (36usize, 19usize, 10.0, 5.0);
        let grid: Vec<Vec<f64>> = (0..n_azi)
            .map(|k| (0..n_zen).map(|j| (10.0 * k as f64) + j as f64).collect())
            .collect();
        let pcv = synth_azi_pcv(grid, dazi, dzen, vec![0.0; n_zen]);
        for k in 0..n_azi {
            for j in 0..n_zen {
                let want = (10.0 * k as f64) + j as f64;
                assert_near(
                    pcv.interpolate_az_zen(k as f64 * dazi, j as f64 * dzen),
                    want,
                    1e-9,
                    &format!("node az{k} zen{j}"),
                );
            }
        }
        // Azimuth wrap: 370 deg is 10 deg, and -10 deg is 350 deg.
        assert_near(pcv.interpolate_az_zen(370.0, 0.0), 10.0, 1e-9, "wrap +370");
        assert_near(pcv.interpolate_az_zen(-10.0, 0.0), 350.0, 1e-9, "wrap -10");
    }

    // --- Invariants that must hold for any valid NOAZI table. ---

    #[test]
    fn interpolation_is_continuous_clamped_and_range_bounded() {
        let mut grids: Vec<ReceiverPcv> = vec![
            synth_pcv(vec![1.0, 3.0, -5.0], 10.0, 5.0),
            synth_pcv(vec![-9.34, -9.18, 0.0, 3.03], 0.0, 5.0),
            synth_pcv(vec![0.0; 1], 0.0, 5.0),
        ];
        if let Some(db) = load_db() {
            grids.push(ReceiverPcv::from_antex(&db, "ASH701945B_M", "SCIT").expect("ASH"));
            grids.push(ReceiverPcv::from_antex(&db, "TRM59800.80", "SCIT").expect("TRM"));
        }
        for pcv in grids {
            let (lo, hi) = grid_range(&pcv);
            let last = (pcv.pcv_grid_mm.len() - 1) as f64;
            let top = pcv.zen_start_deg + last * pcv.zen_step_deg;
            for i in 0..4001 {
                let zen = pcv.zen_start_deg - 20.0 + (top + 40.0 - pcv.zen_start_deg) * (i as f64 / 4000.0);
                let v = pcv.interpolate(zen);
                assert!(v >= lo - 1e-12 && v <= hi + 1e-12, "{} @ {zen}: {v} outside [{lo},{hi}]", pcv.ant_type);
            }
            // Node continuity: the one-sided limits meet the node value.
            for (i, &node) in pcv.pcv_grid_mm.iter().enumerate() {
                let z = pcv.zen_start_deg + i as f64 * pcv.zen_step_deg;
                let eps = pcv.zen_step_deg * 1e-9;
                let left = pcv.interpolate(z - eps);
                let right = pcv.interpolate(z + eps);
                assert!((left - node).abs() < 1e-6 && (right - node).abs() < 1e-6,
                    "{} node {i} @{z}: left {left} right {right} node {node}", pcv.ant_type);
            }
            assert_eq!(pcv.interpolate(pcv.zen_start_deg - 1.0), pcv.pcv_grid_mm[0]);
            assert_eq!(pcv.interpolate(top + 1.0), *pcv.pcv_grid_mm.last().unwrap());
        }
    }

    fn grid_range(pcv: &ReceiverPcv) -> (f64, f64) {
        let lo = pcv.pcv_grid_mm.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = pcv.pcv_grid_mm.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (lo, hi)
    }

    #[test]
    fn dd_correction_is_exactly_zero_for_one_shared_calibration() {
        let pcv = synth_pcv(vec![1.0, 3.0, -5.0], 10.0, 5.0);
        // Same object, same satellite pair seen from both places: the two
        // differences are the same subtraction, so the result is exactly
        // 0.0 - 0.0 rather than merely small.
        assert_eq!(ReceiverPcv::dd_correction_m(&pcv, &pcv, 75.0, 22.5, 75.0, 22.5), 0.0);
        assert_eq!(ReceiverPcv::dd_correction_m(&pcv, &pcv, 13.7, 44.4, 13.7, 44.4), 0.0);
        // Identical table, different instances: same cancellation.
        let clone = pcv.clone();
        assert_eq!(ReceiverPcv::dd_correction_m(&pcv, &clone, 3.0, 19.0, 3.0, 19.0), 0.0);
    }

    #[test]
    fn pco_projects_antex_north_east_up_onto_the_line_of_sight() {
        // ANTEX "NORTH / EAST / UP" order, so x=north, y=east, z=up.
        let pcv = ReceiverPcv {
            ant_type: "PCO".into(),
            radome: "NONE".into(),
            pco_neu_mm: nalgebra::Vector3::new(20.0, 7.0, 100.0),
            pcv_grid_mm: vec![0.0],
            azi_grid_mm: None,
            dazi_deg: 0.0,
            zen_start_deg: 0.0,
            zen_step_deg: 5.0,
        };
        // Zenith: the LOS unit vector is (0, 0, 1) -> up component only.
        assert_near(pcv.pco_correction_m(0.0, 90.0), 0.100, TOL, "up at zenith");
        // Azimuth 0 is north, elevation 0 is the horizon: LOS = (1, 0, 0).
        assert_near(pcv.pco_correction_m(0.0, 0.0), 0.020, TOL, "north at horizon");
        // Azimuth 90 is east: LOS = (0, 1, 0).
        assert_near(pcv.pco_correction_m(90.0, 0.0), 0.007, TOL, "east at horizon");
        // Swapping north and east would swap the first two expectations.
        assert!((pcv.pco_correction_m(0.0, 0.0) - pcv.pco_correction_m(90.0, 0.0)).abs() > 1e-4);
    }

    #[test]
    fn real_ro_azi_grid_is_reproduced_node_for_node() {
        // LEIAR20 LEIM ships a RO-AZI table in igs14.atx: 73 azimuth rows
        // (0..360 step 5 deg) by 19 zenith columns (0..90 step 5 deg).
        let Some(db) = load_db() else { return };
        let pcv = ReceiverPcv::from_antex(&db, "LEIAR20", "LEIM").expect("LEIAR20");
        let grid = pcv.azi_grid_mm.as_ref().expect("LEIAR20 has RO-AZI");
        assert_eq!(grid.len(), 73, "azimuth rows");
        assert_eq!(pcv.dazi_deg, 5.0);
        for row in grid.iter() {
            assert_eq!(row.len(), 19, "zenith columns must be uniform");
        }
        for (k, row) in grid.iter().enumerate() {
            for (j, &node) in row.iter().enumerate() {
                let az = k as f64 * pcv.dazi_deg;
                let zen = j as f64 * pcv.zen_step_deg;
                assert_near(pcv.interpolate_az_zen(az, zen), node, 1e-9, &format!("az{k} zen{j}"));
            }
        }
        // Hand-derived bilinear midpoint of four real nodes: az 2.5 deg is
        // halfway between rows 0 (az 0) and 1 (az 5), zen 7.5 deg halfway
        // between columns 1 (zen 5) and 2 (zen 10).
        let want = 0.25 * (grid[0][1] + grid[0][2] + grid[1][1] + grid[1][2]);
        assert_near(pcv.interpolate_az_zen(2.5, 7.5), want, 1e-9, "real bilinear midpoint");
        // The zenith clamp holds both end columns of the real table.
        assert_eq!(pcv.interpolate_az_zen(0.0, 95.0), grid[0][18]);
        assert_eq!(pcv.interpolate_az_zen(0.0, -30.0), grid[0][0]);
    }

    #[test]
    fn noazi_only_antenna_ignores_the_azimuth_argument() {
        // ASH701945B_M SCIT declares DAZI 0 and ships no RO-AZI table.
        let Some(db) = load_db() else { return };
        let pcv = ReceiverPcv::from_antex(&db, "ASH701945B_M", "SCIT").expect("ASH");
        assert!(pcv.azi_grid_mm.is_none());
        for zen in [0.0, 17.5, 45.0, 80.0] {
            let north = pcv.interpolate_az_zen(0.0, zen);
            assert_near(pcv.interpolate_az_zen(217.3, zen), north, TOL, "azimuth invariance");
            assert_near(north, pcv.interpolate(zen), TOL, "NOAZI fallback");
        }
    }

    #[test]
    fn test_standalone_receiver_pco_pcv() {
        let pcv = ReceiverPcv {
            ant_type: "STANDALONE_TEST".into(),
            radome: "NONE".into(),
            pco_neu_mm: nalgebra::Vector3::new(10.0, 20.0, 30.0),
            pcv_grid_mm: vec![1.0, 2.0, 3.0],
            azi_grid_mm: None,
            dazi_deg: 0.0,
            zen_start_deg: 0.0,
            zen_step_deg: 10.0,
        };

        let pco_zen = pcv.pco_correction_m(0.0, 90.0);
        assert_near(pco_zen, 0.030, TOL, "pco at zenith");
        let pcv_zen = pcv.pcv_correction_m(0.0, 90.0);
        assert_near(pcv_zen, 0.001, TOL, "pcv at zenith");
        assert_near(pcv.total_correction_m(0.0, 90.0), 0.031, TOL, "total at zenith");
    }
