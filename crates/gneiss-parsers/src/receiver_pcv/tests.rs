#![allow(clippy::unwrap_used)]

    use super::*;
    use crate::antex::{AntennaPcv, FrequencyPcv};
    use std::collections::HashMap;

    /// Locate the cached igs14.atx across checkout layouts: the main tree
    /// keeps `datasets/` at the repo root, while worktrees expose the
    /// shared store through a nested `datasets/datasets` symlink. Returns
    /// None when the file is not cached locally (tests then skip).
    fn igs14_path() -> Option<std::path::PathBuf> {
        [
            "../../datasets/igs14.atx",
            "../../datasets/datasets/igs14.atx",
        ]
        .iter()
        .map(|rel| std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(rel))
        .find(|p| p.exists())
    }

    /// Tight tolerance for values that are exact grid nodes (fp noise only).
    const TOL: f64 = 1e-9;

    // --- Ground truth: datasets/igs14.atx G01 NOAZI rows (coordinator-verified). ---

    /// ASH701945B_M SCIT: zen 0..=80 deg, step 5 deg.
    const ASH_G01_GRID: [f64; 17] = [
        0.00, -0.53, -1.51, -2.84, -4.24, -5.83, -7.18, -8.31, -9.10, -9.34,
        -9.18, -8.42, -7.16, -5.48, -3.31, -0.44, 3.03,
    ];
    /// LEIAR20 LEIM: zen 0..=90 deg, step 5 deg.
    const LEIAR_G01_GRID: [f64; 19] = [
        0.00, -0.06, -0.23, -0.55, -1.02, -1.61, -2.26, -2.90, -3.42, -3.77,
        -3.92, -3.85, -3.53, -2.88, -1.77, -0.05, 2.33, 5.28, 8.43,
    ];
    /// TRM59800.00 SCIT and TRM59800.80 SCIT (identical tables).
    const TRM_G01_GRID: [f64; 19] = [
        0.00, -0.40, -1.51, -3.05, -4.69, -6.17, -7.38, -8.27, -8.85, -9.07,
        -8.82, -7.97, -6.45, -4.33, -1.74, 1.26, 4.82, 9.34, 15.24,
    ];

    /// `None` when igs14.atx isn't cached locally -- callers skip (matches
    /// `igs14_path`'s own documented contract, which this used to violate
    /// by panicking instead: every real-data test here hard-failed on a
    /// fresh checkout with no `datasets/` present, e.g. a CI runner).
    fn load_db() -> Option<AntexDatabase> {
        let path = igs14_path()?;
        Some(AntexDatabase::parse(path).expect("igs14.atx must parse"))
    }

    fn assert_near(actual: f64, expected: f64, tol: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= tol,
            "{what}: got {actual:.9}, expected {expected:.9} (tol {tol})"
        );
    }

    fn assert_grid_matches(pcv: &ReceiverPcv, expected: &[f64]) {
        let label = format!("{} {}", pcv.ant_type, pcv.radome);
        assert_eq!(pcv.pcv_grid_mm.len(), expected.len(), "{label}: grid length");
        for (i, (got, want)) in pcv.pcv_grid_mm.iter().zip(expected).enumerate() {
            assert_eq!(got, want, "{label} node {i}");
        }
    }

    // --- Synthetic-database fixtures for edge cases igs14.atx cannot express. ---

    fn freq(code: &str, neu: [f64; 3], noazi: Vec<f64>) -> (String, FrequencyPcv) {
        (
            code.to_string(),
            FrequencyPcv {
                frequency_code: code.to_string(),
                pco: nalgebra::Vector3::new(neu[0], neu[1], neu[2]),
                noazi,
                azi: None,
            },
        )
    }

    fn synth_db(ant_type: &str, dzen: f64, freqs: Vec<(String, FrequencyPcv)>) -> AntexDatabase {
        let mut map = HashMap::new();
        for (code, freq) in freqs {
            map.insert(code, freq);
        }
        AntexDatabase::new(vec![AntennaPcv {
            antenna_type: ant_type.to_string(),
            serial_num: String::new(),
            valid_from: None,
            valid_until: None,
            dzen,
            zen1: 10.0,
            zen2: 20.0,
            dazi: 0.0,
            frequencies: map,
        }])
    }

    /// Three-node fixture with non-zero start, used for interpolation math.
    fn synth_pcv(grid: Vec<f64>, start: f64, step: f64) -> ReceiverPcv {
        ReceiverPcv {
            ant_type: "SYNTH".to_string(),
            radome: "NONE".to_string(),
            pco_neu_mm: nalgebra::Vector3::zeros(),
            pcv_grid_mm: grid,
            azi_grid_mm: None,
            dazi_deg: 0.0,
            zen_start_deg: start,
            zen_step_deg: step,
        }
    }

    // --- Requirement 1: known antennas load with expected PCO/PCV values. ---

    #[test]
    fn loads_ash701945b_scit_with_ground_truth() {
        let Some(db) = load_db() else { return };
        let pcv = ReceiverPcv::from_antex(&db, "ASH701945B_M", "SCIT")
            .expect("ASH701945B_M SCIT must exist in igs14.atx");
        assert_eq!(pcv.ant_type, "ASH701945B_M");
        assert_eq!(pcv.radome, "SCIT");
        assert_near(pcv.pco_neu_mm[0], 0.76, TOL, "PCO north");
        assert_near(pcv.pco_neu_mm[1], -0.43, TOL, "PCO east");
        assert_near(pcv.pco_neu_mm[2], 89.24, TOL, "PCO up");
        assert_eq!(pcv.zen_start_deg, 0.0);
        assert_eq!(pcv.zen_step_deg, 5.0);
        assert_grid_matches(&pcv, &ASH_G01_GRID);
        assert_near(pcv.interpolate(10.0), -1.51, TOL, "PCV @zen10");
        assert_near(pcv.interpolate(30.0), -7.18, TOL, "PCV @zen30");
        assert_near(pcv.interpolate(60.0), -7.16, TOL, "PCV @zen60");
    }

    #[test]
    fn loads_leiar20_leim_with_ground_truth() {
        let Some(db) = load_db() else { return };
        let pcv = ReceiverPcv::from_antex(&db, "LEIAR20", "LEIM")
            .expect("LEIAR20 LEIM must exist in igs14.atx");
        assert_eq!(pcv.ant_type, "LEIAR20");
        assert_eq!(pcv.radome, "LEIM");
        assert_near(pcv.pco_neu_mm[0], 0.66, TOL, "PCO north");
        assert_near(pcv.pco_neu_mm[1], 0.08, TOL, "PCO east");
        assert_near(pcv.pco_neu_mm[2], 124.58, TOL, "PCO up");
        assert_eq!(pcv.zen_start_deg, 0.0);
        assert_eq!(pcv.zen_step_deg, 5.0);
        assert_grid_matches(&pcv, &LEIAR_G01_GRID);
        assert_near(pcv.interpolate(10.0), -0.23, TOL, "PCV @zen10");
        assert_near(pcv.interpolate(30.0), -2.26, TOL, "PCV @zen30");
        assert_near(pcv.interpolate(60.0), -3.53, TOL, "PCV @zen60");
    }

    #[test]
    fn loads_trm59800_scit_variants_with_ground_truth() {
        let Some(db) = load_db() else { return };
        for family in ["TRM59800.00", "TRM59800.80"] {
            let pcv = ReceiverPcv::from_antex(&db, family, "SCIT")
                .unwrap_or_else(|| panic!("{family} SCIT must exist in igs14.atx"));
            assert_eq!(pcv.ant_type, family);
            assert_eq!(pcv.radome, "SCIT");
            assert_near(pcv.pco_neu_mm[0], 1.32, TOL, "PCO north");
            assert_near(pcv.pco_neu_mm[1], 0.88, TOL, "PCO east");
            assert_near(pcv.pco_neu_mm[2], 84.85, TOL, "PCO up");
            assert_grid_matches(&pcv, &TRM_G01_GRID);
            assert_near(pcv.interpolate(10.0), -1.51, TOL, "PCV @zen10");
            assert_near(pcv.interpolate(30.0), -7.38, TOL, "PCV @zen30");
            assert_near(pcv.interpolate(60.0), -6.45, TOL, "PCV @zen60");
        }
    }

    // --- Requirement 2: interpolation exact at grid nodes. ---

    #[test]
    fn interpolation_is_exact_at_every_grid_node() {
        let Some(db) = load_db() else { return };
        let cases = [
            ("ASH701945B_M", "SCIT"), ("LEIAR20", "LEIM"),
            ("TRM59800.00", "SCIT"), ("TRM59800.80", "SCIT"),
        ];
        for (family, radome) in cases {
            let pcv = ReceiverPcv::from_antex(&db, family, radome)
                .unwrap_or_else(|| panic!("{family} {radome} must load"));
            for (i, &expected) in pcv.pcv_grid_mm.iter().enumerate() {
                let zen = pcv.zen_start_deg + i as f64 * pcv.zen_step_deg;
                assert_eq!(pcv.interpolate(zen), expected, "{family} {radome} node {i} @ {zen}");
            }
        }
    }

    // --- Requirement 3: linear formula verified by hand calculation. ---

    #[test]
    fn midpoint_interpolation_matches_hand_calculation() {
        // Nodes: z=10 -> 1.0 mm, z=15 -> 3.0 mm, z=20 -> -5.0 mm.
        let pcv = synth_pcv(vec![1.0, 3.0, -5.0], 10.0, 5.0);
        // Halfway 12.5 deg: 1.0*0.5 + 3.0*0.5 = 2.0 mm.
        assert_eq!(pcv.interpolate(12.5), 2.0, "midpoint of first pair");
        // Quarter point 11.25 deg: 1.0*0.75 + 3.0*0.25 = 1.5 mm.
        assert_eq!(pcv.interpolate(11.25), 1.5, "quarter point");
        // Halfway 17.5 deg: 3.0*0.5 + (-5.0)*0.5 = -1.0 mm.
        assert_eq!(pcv.interpolate(17.5), -1.0, "midpoint of second pair");
    }

    #[test]
    fn fractional_interpolation_on_real_nodes_matches_hand_calc() {
        let Some(db) = load_db() else { return };
        let leiar = ReceiverPcv::from_antex(&db, "LEIAR20", "LEIM").expect("LEIAR20");
        // Between -2.26 mm (z=30) and -2.90 mm (z=35): midpoint -2.58 mm.
        assert_near(leiar.interpolate(32.5), -2.58, TOL, "LEIAR20 @32.5 deg");
        let trm = ReceiverPcv::from_antex(&db, "TRM59800.00", "SCIT").expect("TRM");
        // Between -6.17 mm (z=25) and -7.38 mm (z=30): midpoint -6.775 mm.
        assert_near(trm.interpolate(27.5), -6.775, TOL, "TRM @27.5 deg");
    }

    // --- Requirement 4: clamping outside the grid range. ---

    #[test]
    fn clamps_outside_grid_range() {
        let pcv = synth_pcv(vec![1.0, 3.0, -5.0], 10.0, 5.0);
        assert_eq!(pcv.interpolate(-40.0), 1.0, "far below clamps to node 0");
        assert_eq!(pcv.interpolate(9.999), 1.0, "just below start clamps");
        assert_eq!(pcv.interpolate(1234.0), -5.0, "far above clamps to last");
        let Some(db) = load_db() else { return };
        let ash = ReceiverPcv::from_antex(&db, "ASH701945B_M", "SCIT").expect("ASH");
        assert_eq!(ash.interpolate(-1.0), 0.00, "below ASH grid");
        assert_eq!(ash.interpolate(85.0), 3.03, "ASH grid ends at 80 deg");
    }

    #[test]
    fn single_node_grid_clamps_to_that_value_everywhere() {
        let db = synth_db("ONE_NODE", 5.0, vec![freq("G01", [0.0; 3], vec![-7.5])]);
        let pcv = ReceiverPcv::from_antex(&db, "ONE_NODE", "NONE").expect("one node");
        assert_eq!(pcv.interpolate(0.0), -7.5);
        assert_eq!(pcv.interpolate(42.0), -7.5);
        assert_eq!(pcv.interpolate(88.0), -7.5);
    }

    // --- Requirement 5: same-type antennas give ~0 DD correction. ---

    #[test]
    fn dd_correction_same_type_is_zero() {
        let Some(db) = load_db() else { return };
        let rov = ReceiverPcv::from_antex(&db, "TRM59800.80", "SCIT").expect("rover");
        let bas = ReceiverPcv::from_antex(&db, "TRM59800.80", "SCIT").expect("base");
        // Both stations observe the SAME satellite pair, so the sat/ref
        // zeniths match across the arguments; identical tables then cancel.
        let corr = ReceiverPcv::dd_correction_m(&rov, &bas, 75.0, 22.5, 75.0, 22.5);
        assert!(
            corr.abs() < 1e-12,
            "identical antennas must cancel, got {corr}"
        );
        // Long-baseline site parallax (millidegree zenith shifts between
        // stations observing one satellite pair) leaves only a slope-times-
        // delta residual: ~0.71 mm/deg * 0.002 deg ≈ 0.8 um here, orders
        // below the mm-scale cross-type signature the correction targets.
        let parallax =
            ReceiverPcv::dd_correction_m(&rov, &bas, 75.002, 22.498, 75.0, 22.5);
        assert!(
            parallax.abs() < 1e-6,
            "parallax residual must stay sub-micrometre, got {parallax}"
        );
        // Aliasing one instance through both arguments cancels exactly.
        let alias = ReceiverPcv::dd_correction_m(&rov, &rov, 13.7, 44.4, 13.7, 44.4);
        assert_eq!(alias, 0.0);
    }

    // --- Requirement 6: cross-type correction non-zero and millimetre-scale. ---

    #[test]
    fn dd_correction_cross_type_hand_calculated() {
        let Some(db) = load_db() else { return };
        let leiar = ReceiverPcv::from_antex(&db, "LEIAR20", "LEIM").expect("LEIAR20");
        let ash = ReceiverPcv::from_antex(&db, "ASH701945B_M", "SCIT").expect("ASH");
        // All four angles sit on 5-deg nodes, so this is pure arithmetic:
        //   LEIAR20: (-2.26) - (-0.23) = -2.03 mm
        //   ASH:     (-7.18) - (-1.51) = -5.67 mm
        //   DD = (-2.03) - (-5.67) = +3.64 mm = +0.00364 m
        let corr = ReceiverPcv::dd_correction_m(&leiar, &ash, 30.0, 10.0, 30.0, 10.0);
        assert_near(corr, 0.00364, TOL, "cross-type DD correction (m)");
        assert!(corr.abs() > 1e-6, "cross-type correction must be non-zero");
        assert!(corr.abs() < 0.05, "correction must stay millimetre-scale");

        // Asymmetric rover/base geometry catches argument-order mix-ups:
        //   TRM59800.80: (-6.45) - (-3.05) = -3.40 mm
        //   ASH:         (-7.16) - (-2.84) = -4.32 mm
        //   DD = -3.40 - (-4.32) = +0.92 mm = +0.00092 m
        let trm = ReceiverPcv::from_antex(&db, "TRM59800.80", "SCIT").expect("TRM");
        let asym =
            ReceiverPcv::dd_correction_m(&trm, &ash, 60.0, 15.0, 60.0, 15.0);
        assert_near(asym, 0.00092, TOL, "asymmetric cross-type DD (m)");

        // Swapping satellite/reference zeniths negates the correction.
        let swapped =
            ReceiverPcv::dd_correction_m(&leiar, &ash, 10.0, 30.0, 10.0, 30.0);
        assert_near(swapped, -0.00364, TOL, "swap negates sign");
    }

    // --- Requirement 7: graceful None for missing antennas. ---

    #[test]
    fn missing_antenna_returns_none() {
        let Some(db) = load_db() else { return };
        assert!(
            ReceiverPcv::from_antex(&db, "NO_SUCH_ANTENNA_XY", "NONE").is_none(),
            "unknown family must return None"
        );
        assert!(
            ReceiverPcv::from_antex(&db, "TRM59800.80", "LEIM").is_none(),
            "known family with wrong radome must return None"
        );
        assert!(
            ReceiverPcv::from_antex(&db, "LEIAR20", "SCIT").is_none(),
            "known family paired with foreign radome must return None"
        );
    }

    #[test]
    fn none_radome_and_single_token_entries_resolve() {
        let Some(db) = load_db() else { return };
        let none_entry =
            ReceiverPcv::from_antex(&db, "ASH701945B_M", "NONE").expect("NONE entry");
        assert_eq!(none_entry.radome, "NONE");
        assert_near(none_entry.pco_neu_mm[0], 0.07, TOL, "NONE PCO north");
        assert_near(none_entry.pco_neu_mm[1], -0.69, TOL, "NONE PCO east");
        assert_near(none_entry.pco_neu_mm[2], 89.92, TOL, "NONE PCO up");

        // Single-token TYPE/SERIAL entries behave as unradomed ("NONE").
        let db2 =
            synth_db("AOAD/M_T", 5.0, vec![freq("G01", [1.0, 2.0, 3.0], vec![0.5])]);
        let got =
            ReceiverPcv::from_antex(&db2, "AOAD/M_T", "NONE").expect("token-less entry");
        assert_eq!(got.radome, "NONE");
        assert!(ReceiverPcv::from_antex(&db2, "AOAD/M_T", "SCIT").is_none());
        assert!(
            ReceiverPcv::from_antex(&db2, "AOAD/M_T", "").is_some(),
            "empty radome string resolves unradomed entries"
        );
    }

    #[test]
    fn missing_l1_frequency_returns_none() {
        let db = synth_db("ONLY_L2", 5.0, vec![freq("G02", [0.0; 3], vec![1.0])]);
        assert!(ReceiverPcv::from_antex(&db, "ONLY_L2", "NONE").is_none());
    }

    #[test]
    fn legacy_g1_frequency_code_fallback_loads() {
        let db = synth_db(
            "LEGACY",
            5.0,
            vec![freq("G1", [4.0, 5.0, 6.0], vec![0.25, 0.75])],
        );
        let pcv =
            ReceiverPcv::from_antex(&db, "LEGACY", "NONE").expect("G1 fallback");
        assert_near(pcv.pco_neu_mm[2], 6.0, TOL, "PCO up via G1 fallback");
        assert_eq!(pcv.pcv_grid_mm.len(), 2);
        // Zenith metadata comes from the block, not assumed defaults.
        assert_eq!(pcv.zen_start_deg, 10.0);
        assert_eq!(pcv.zen_step_deg, 5.0);
    }

    #[test]
    fn degenerate_or_empty_grids_return_none() {
        let zero_step =
            synth_db("ZERO_STEP", 0.0, vec![freq("G01", [0.0; 3], vec![1.0, 2.0])]);
        assert!(
            ReceiverPcv::from_antex(&zero_step, "ZERO_STEP", "NONE").is_none(),
            "dzen == 0 would divide by zero during interpolation"
        );
        let negative =
            synth_db("NEG_STEP", -5.0, vec![freq("G01", [0.0; 3], vec![1.0, 2.0])]);
        assert!(
            ReceiverPcv::from_antex(&negative, "NEG_STEP", "NONE").is_none(),
            "negative dzen is meaningless"
        );
        let empty =
            synth_db("EMPTY_GRID", 5.0, vec![freq("G01", [0.0; 3], vec![])]);
        assert!(
            ReceiverPcv::from_antex(&empty, "EMPTY_GRID", "NONE").is_none(),
            "empty NOAZI table carries no calibration"
        );
    }

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
