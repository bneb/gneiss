use crate::engine::measurement::build_dense_covariance_matrix;
use gneiss_core::sat::{Constellation, SatelliteId};

#[test]
fn test_dd_covariance_off_diagonals() {
    let sat1 = SatelliteId {
        constellation: Constellation::Gps,
        prn: 1,
    };
    let sat2 = SatelliteId {
        constellation: Constellation::Gps,
        prn: 2,
    };

    let r_diagonals = vec![1.0, 2.0];
    let meas_types = vec![
        (sat1, 0, 0.5), // G01, L1 PR, ref_var = 0.5
        (sat2, 0, 0.5), // G02, L1 PR, ref_var = 0.5
    ];

    let r_mat = build_dense_covariance_matrix(&r_diagonals, &meas_types);

    assert_eq!(r_mat[(0, 0)], 1.0);
    assert_eq!(r_mat[(1, 1)], 2.0);
    assert_eq!(r_mat[(0, 1)], 0.5);
    assert_eq!(r_mat[(1, 0)], 0.5);
}

#[test]
fn test_dd_covariance_frequency_independence() {
    let sat1 = SatelliteId {
        constellation: Constellation::Gps,
        prn: 1,
    };
    let sat2 = SatelliteId {
        constellation: Constellation::Gps,
        prn: 2,
    };

    let r_diagonals = vec![1.0, 2.0];
    let meas_types = vec![
        (sat1, 0, 0.5), // G01, L1 PR
        (sat2, 4, 0.5), // G02, L2 PR
    ];

    let r_mat = build_dense_covariance_matrix(&r_diagonals, &meas_types);

    assert_eq!(r_mat[(0, 0)], 1.0);
    assert_eq!(r_mat[(1, 1)], 2.0);
    // Should be zero because types (0 vs 4) are different
    assert_eq!(r_mat[(0, 1)], 0.0);
    assert_eq!(r_mat[(1, 0)], 0.0);
}

#[test]
fn test_dd_covariance_constellation_independence() {
    let sat1 = SatelliteId {
        constellation: Constellation::Gps,
        prn: 1,
    };
    let sat2 = SatelliteId {
        constellation: Constellation::Galileo,
        prn: 2,
    };

    let r_diagonals = vec![1.0, 2.0];
    let meas_types = vec![
        (sat1, 0, 0.5), // G01, L1 PR
        (sat2, 0, 0.5), // E02, L1 PR
    ];

    let r_mat = build_dense_covariance_matrix(&r_diagonals, &meas_types);

    assert_eq!(r_mat[(0, 0)], 1.0);
    assert_eq!(r_mat[(1, 1)], 2.0);
    // Should be zero because constellations are different
    assert_eq!(r_mat[(0, 1)], 0.0);
    assert_eq!(r_mat[(1, 0)], 0.0);
}
