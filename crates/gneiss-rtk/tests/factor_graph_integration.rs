use gneiss_rtk::factor_graph::gnss_factors::{CarrierPhaseFactor, PseudorangeFactor};
use gneiss_rtk::factor_graph::{FactorGraphOptimizer, PriorFactor};
use nalgebra::{DMatrix, DVector, Vector3};

// ---------- Tests with PriorFactor constraints (observable) ----------

/// PriorFactor + PseudorangeFactor: the prior constrains position so the
/// optimizer converges to a known solution.
#[test]
fn test_prior_with_pseudorange_convergence() {
    let sat_pos: Vector3<f64> = Vector3::new(2.0, 3.0, 6.0);
    let dist: f64 = 7.0;
    let truth_dt: f64 = 0.1;
    let pr: f64 = dist + truth_dt;

    let mut optimizer = FactorGraphOptimizer::new();

    // Strong prior on position at origin — weight = 1000 (variance = 0.001)
    let prior_info: DMatrix<f64> = DMatrix::from_diagonal(&DVector::from_vec(vec![
        1000.0, 1000.0, 1000.0, 0.001,
    ]));
    optimizer.add_factor(Box::new(PriorFactor {
        information: prior_info,
    }));

    // Pseudorange measurement with weaker weight
    optimizer.add_factor(Box::new(PseudorangeFactor {
        sat_pos,
        measured_pr: pr,
        variance: 1.0,
        sat_clock_bias: 0.0,
        tropo_dry_delay: 0.0,
        map_wet: 0.0,
        index_x: 0,
        index_y: 1,
        index_z: 2,
        index_dt: 3,
        index_zwd: None,
        robust_threshold: 3.0,
    }));

    let initial = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
    let (optimized, _cov) = optimizer.optimize(&initial, 30, 1e-6);

    // Prior dominates position, PR helps determine clock
    assert!((optimized[3] - truth_dt).abs() < 1e-3, "dt should be ~0.1, got {}", optimized[3]);
}

/// Two GNSS factors: ensure the optimizer build_normal_equations does not
/// panic with mixed factor types.
#[test]
fn test_mixed_factor_types_no_panic() {
    let mut optimizer = FactorGraphOptimizer::new();

    optimizer.add_factor(Box::new(PseudorangeFactor {
        sat_pos: Vector3::new(2.0, 3.0, 6.0),
        measured_pr: 7.0,
        variance: 1.0,
        sat_clock_bias: 0.0,
        tropo_dry_delay: 0.0,
        map_wet: 0.0,
        index_x: 0,
        index_y: 1,
        index_z: 2,
        index_dt: 3,
        index_zwd: None,
        robust_threshold: 3.0,
    }));

    optimizer.add_factor(Box::new(CarrierPhaseFactor {
        sat_pos: Vector3::new(2.0, 3.0, 6.0),
        measured_cp: 7.48,
        variance: 0.01,
        sat_clock_bias: 0.0,
        tropo_dry_delay: 0.0,
        map_wet: 0.0,
        wavelength: 0.19,
        index_x: 0,
        index_y: 1,
        index_z: 2,
        index_dt: 3,
        index_zwd: None,
        index_amb: 4,
        robust_threshold: 3.0,
    }));

    // Must not panic during optimization
    let initial = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 0.0]);
    let (_optimized, _cov) = optimizer.optimize(&initial, 10, 1e-6);
    // Test passes if no panic
}

// ---------- Multi-satellite tests (overdetermined) ----------

/// Four pseudorange factors with good geometry should converge to the true position.
#[test]
fn test_four_satellite_position_solution() {
    let mut optimizer = FactorGraphOptimizer::new();

    // 4 satellites in tetrahedral-ish configuration, far enough for good geometry
    let sats: [Vector3<f64>; 4] = [
        Vector3::new(20000000.0, 0.0, 0.0),
        Vector3::new(0.0, 20000000.0, 0.0),
        Vector3::new(0.0, 0.0, 20000000.0),
        Vector3::new(14142135.6, 14142135.6, 0.0), // sqrt(2)*1e7 on x,y
    ];
    let truth_pos: Vector3<f64> = Vector3::new(100.0, 200.0, 300.0);
    let truth_dt: f64 = 1e-4;

    for sat_pos in &sats {
        let dx: f64 = sat_pos.x - truth_pos.x;
        let dy: f64 = sat_pos.y - truth_pos.y;
        let dz: f64 = sat_pos.z - truth_pos.z;
        let dist: f64 = (dx * dx + dy * dy + dz * dz).sqrt();
        optimizer.add_factor(Box::new(PseudorangeFactor {
            sat_pos: *sat_pos,
            measured_pr: dist + truth_dt,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            robust_threshold: 3.0,
        }));
    }

    let initial = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
    let (optimized, _cov) = optimizer.optimize(&initial, 30, 1e-6);

    assert!((optimized[0] - truth_pos.x).abs() < 1.0, "x should be ~100, got {}", optimized[0]);
    assert!((optimized[1] - truth_pos.y).abs() < 1.0, "y should be ~200, got {}", optimized[1]);
    assert!((optimized[2] - truth_pos.z).abs() < 1.0, "z should be ~300, got {}", optimized[2]);
    assert!((optimized[3] - truth_dt).abs() < 1e-3, "dt should be ~1e-4, got {}", optimized[3]);
}

/// Five satellites + carrier phase should solve for position, clock, and ambiguity.
#[test]
fn test_multi_satellite_with_carrier_phase() {
    let mut optimizer = FactorGraphOptimizer::new();

    let sats: [Vector3<f64>; 5] = [
        Vector3::new(20000000.0, 0.0, 0.0),
        Vector3::new(0.0, 20000000.0, 0.0),
        Vector3::new(0.0, 0.0, 20000000.0),
        Vector3::new(14142135.6, 14142135.6, 0.0),
        Vector3::new(-10000000.0, 10000000.0, 14142135.6),
    ];

    let truth_pos: Vector3<f64> = Vector3::new(100.0, 200.0, 300.0);
    let truth_dt: f64 = 1e-4;
    let truth_amb: f64 = 2.0;
    let wavelength: f64 = 0.19;

    for sat_pos in &sats {
        let dx: f64 = sat_pos.x - truth_pos.x;
        let dy: f64 = sat_pos.y - truth_pos.y;
        let dz: f64 = sat_pos.z - truth_pos.z;
        let dist: f64 = (dx * dx + dy * dy + dz * dz).sqrt();

        optimizer.add_factor(Box::new(PseudorangeFactor {
            sat_pos: *sat_pos,
            measured_pr: dist + truth_dt,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            robust_threshold: 3.0,
        }));

        optimizer.add_factor(Box::new(CarrierPhaseFactor {
            sat_pos: *sat_pos,
            measured_cp: dist + truth_dt + truth_amb * wavelength,
            variance: 0.01,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            wavelength,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            index_amb: 4,
            robust_threshold: 3.0,
        }));
    }

    let initial = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 0.0]);
    let (optimized, _cov) = optimizer.optimize(&initial, 50, 1e-6);

    // With 10 measurements (5 PR + 5 CP) for 5 states, should converge well
    assert!((optimized[0] - truth_pos.x).abs() < 1.0, "x should be ~100, got {}", optimized[0]);
    assert!((optimized[1] - truth_pos.y).abs() < 1.0, "y should be ~200, got {}", optimized[1]);
    assert!((optimized[2] - truth_pos.z).abs() < 1.0, "z should be ~300, got {}", optimized[2]);
    assert!((optimized[3] - truth_dt).abs() < 1e-3, "dt should be ~1e-4, got {}", optimized[3]);
    assert!((optimized[4] - truth_amb).abs() < 1.0, "amb should be ~2.0, got {}", optimized[4]);
}

/// Robust outlier rejection works with multiple identical satellites
/// (one good, one outlier). The good measurements dominate.
#[test]
fn test_pseudorange_outlier_rejection() {
    let mut optimizer = FactorGraphOptimizer::new();

    // 4 good satellite measurements
    let sats: [Vector3<f64>; 4] = [
        Vector3::new(20000000.0, 0.0, 0.0),
        Vector3::new(0.0, 20000000.0, 0.0),
        Vector3::new(0.0, 0.0, 20000000.0),
        Vector3::new(14142135.6, 14142135.6, 0.0),
    ];
    let truth_pos: Vector3<f64> = Vector3::new(100.0, 200.0, 300.0);
    let truth_dt: f64 = 1e-4;

    for sat_pos in &sats {
        let dx: f64 = sat_pos.x - truth_pos.x;
        let dy: f64 = sat_pos.y - truth_pos.y;
        let dz: f64 = sat_pos.z - truth_pos.z;
        let dist: f64 = (dx * dx + dy * dy + dz * dz).sqrt();
        optimizer.add_factor(Box::new(PseudorangeFactor {
            sat_pos: *sat_pos,
            measured_pr: dist + truth_dt,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            index_x: 0,
            index_y: 1,
            index_z: 2,
            index_dt: 3,
            index_zwd: None,
            robust_threshold: 3.0,
        }));
    }

    // Add an outlier satellite measurement with 1000m error
    let outlier_pos: Vector3<f64> = Vector3::new(-10000000.0, 10000000.0, 14142135.6);
    let outlier_dist: f64 = ((outlier_pos.x - truth_pos.x).powi(2)
        + (outlier_pos.y - truth_pos.y).powi(2)
        + (outlier_pos.z - truth_pos.z).powi(2))
    .sqrt();
    optimizer.add_factor(Box::new(PseudorangeFactor {
        sat_pos: outlier_pos,
        measured_pr: outlier_dist + truth_dt + 1000.0,
        variance: 1.0,
        sat_clock_bias: 0.0,
        tropo_dry_delay: 0.0,
        map_wet: 0.0,
        index_x: 0,
        index_y: 1,
        index_z: 2,
        index_dt: 3,
        index_zwd: None,
        robust_threshold: 3.0,
    }));

    let initial = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
    let (optimized, _cov) = optimizer.optimize(&initial, 30, 1e-6);

    // Even with one outlier, the 4 good measurements should keep solution near truth
    assert!((optimized[0] - truth_pos.x).abs() < 10.0, "x should stay near 100 despite outlier, got {}", optimized[0]);
    assert!((optimized[1] - truth_pos.y).abs() < 10.0, "y should stay near 200 despite outlier, got {}", optimized[1]);
    assert!((optimized[2] - truth_pos.z).abs() < 10.0, "z should stay near 300 despite outlier, got {}", optimized[2]);
}

// ---------- Covariance / structural tests ----------

/// Verify final covariance is symmetric positive-definite with GNSS factors.
#[test]
fn test_final_covariance_with_gnss_factors() {
    let mut optimizer = FactorGraphOptimizer::new();
    optimizer.add_factor(Box::new(PseudorangeFactor {
        sat_pos: Vector3::new(2.0, 3.0, 6.0),
        measured_pr: 7.0,
        variance: 1.0,
        sat_clock_bias: 0.0,
        tropo_dry_delay: 0.0,
        map_wet: 0.0,
        index_x: 0,
        index_y: 1,
        index_z: 2,
        index_dt: 3,
        index_zwd: None,
        robust_threshold: 3.0,
    }));

    let initial = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
    let (_optimized, cov) = optimizer.optimize(&initial, 10, 1e-6);

    for i in 0..cov.nrows() {
        assert!(
            cov[(i, i)] > 0.0,
            "Covariance diagonal [{}] should be positive, got {}",
            i,
            cov[(i, i)]
        );
    }
    for i in 0..cov.nrows() {
        for j in 0..cov.ncols() {
            assert!(
                (cov[(i, j)] - cov[(j, i)]).abs() < 1e-10,
                "Covariance not symmetric at [{},{}]: {} vs {}",
                i,
                j,
                cov[(i, j)],
                cov[(j, i)]
            );
        }
    }
}

/// Test that PriorFactor alone with the optimizer converges to zero.
#[test]
fn test_prior_factor_only_optimization() {
    let mut optimizer = FactorGraphOptimizer::new();
    let prior_info: DMatrix<f64> = DMatrix::from_diagonal(&DVector::from_vec(vec![2.0, 2.0]));
    optimizer.add_factor(Box::new(PriorFactor {
        information: prior_info,
    }));

    let initial_state = DVector::from_vec(vec![5.0, -3.0]);
    let (optimized, cov) = optimizer.optimize(&initial_state, 10, 1e-6);

    assert!(
        optimized[0].abs() < 1e-3,
        "Prior-only should converge to 0. Got {}",
        optimized[0]
    );
    assert!(
        optimized[1].abs() < 1e-3,
        "Prior-only should converge to 0. Got {}",
        optimized[1]
    );
    assert!(cov[(0, 0)] > 0.0, "Covariance diagonal should be positive");
}
