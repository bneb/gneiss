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

mod corrections;
mod loading;

