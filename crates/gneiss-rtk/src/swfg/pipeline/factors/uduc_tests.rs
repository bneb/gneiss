use super::*;
use nalgebra::DVector;

#[test]
fn test_uduc_pseudorange_zero_residual_at_truth() {
    let obs = CorrectedObservation {
        satellite: 1,
        constellation_id: 0,
        pr_l1: 20_000_000.0,
        pr_l2: None,
        cp_l1: None,
        cp_l1_lli: None,
        cp_l2: None,
        doppler: 0.0,
        snr_dbhz: 45.0,
        sat_pos_ecef: Vector3::new(26_560_000.0, 0.0, 0.0),
        sat_clock_m: 100.0,
        f1: 1575.42e6,
        f2: 1227.60e6,
        freq_num: 0,
        elevation_rad: 1.0,
        tropo_dry_m: 2.30,
        tropo_map_wet: 1.50,
        iono_l1_m: 4.50,
        variance_m2: 1.0,
        cp_variance_m2: 0.0001,
    };

    let norm = GeodeticNormalizations {
        sat_relativity_m: -8.859,
        shapiro_delay_m: 0.012,
        sat_pco_ecef_m: Vector3::new(1.5, 0.0, 0.0),
        solid_earth_tide_m: Vector3::new(0.15, 0.0, 0.0),
        ocean_tide_loading_m: Vector3::new(0.02, 0.0, 0.0),
    };

    let p_id = VariableId::new(1);
    let c_id = VariableId::new(2);
    let z_id = VariableId::new(3);
    let i_id = VariableId::new(4);

    let mut graph_vars = std::collections::BTreeMap::new();
    graph_vars.insert(p_id, crate::swfg::variables::VariableNode {
        id: p_id,
        kind: crate::swfg::variables::VariableKind::Pose { epoch: 0 },
        value: DVector::from_vec(vec![6_378_137.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    });
    graph_vars.insert(c_id, crate::swfg::variables::VariableNode {
        id: c_id,
        kind: crate::swfg::variables::VariableKind::ClockBias { epoch: 0, constellation_id: 0 },
        value: DVector::from_vec(vec![50.0]),
    });
    graph_vars.insert(z_id, crate::swfg::variables::VariableNode {
        id: z_id,
        kind: crate::swfg::variables::VariableKind::TropoZwd { epoch: 0 },
        value: DVector::from_vec(vec![0.20]),
    });
    graph_vars.insert(i_id, crate::swfg::variables::VariableNode {
        id: i_id,
        kind: crate::swfg::variables::VariableKind::IonosphereSlant { constellation_id: 0, epoch: 0, satellite: 1 },
        value: DVector::from_vec(vec![4.50]),
    });

    let values = VariableValues::build(&graph_vars);

    let rx_pos = Vector3::new(6_378_137.0 + 0.15 + 0.02, 0.0, 0.0);
    let sat_apc = Vector3::new(26_560_000.0 + 1.5, 0.0, 0.0);
    let geom_range = (sat_apc - rx_pos).norm();
    let expected_pr = geom_range + 50.0 - 100.0 - 8.859 + 0.012 + (2.30 + 1.50 * 0.20) + 4.50;

    let factor = UducPseudorangeFactor::new(obs, 0, expected_pr, p_id, c_id, Some(z_id), Some(i_id), 1.0, norm);
    let res = factor.residual(&values);
    assert!(res[0].abs() < 1e-10, "Residual was {:.10}m", res[0]);

    let jac = factor.jacobian(&values);
    assert_eq!(jac[(0, 0)], 1.0);  // Satellite along +X: d(res)/dx = +1.0
    assert_eq!(jac[(0, 6)], -1.0); // d(res)/d(clock) = -1.0
    assert_eq!(jac[(0, 7)], -1.50); // d(res)/d(ZWD) = -1.50
    assert_eq!(jac[(0, 8)], -1.0);  // d(res)/d(Iono) = -1.0
}

#[test]
fn test_uduc_carrier_phase_zero_residual_at_truth() {
    let obs = CorrectedObservation {
        satellite: 1,
        constellation_id: 0,
        pr_l1: 20_000_000.0,
        pr_l2: None,
        cp_l1: Some(100_000_000.0),
        cp_l1_lli: None,
        cp_l2: None,
        doppler: 0.0,
        snr_dbhz: 45.0,
        sat_pos_ecef: Vector3::new(26_560_000.0, 0.0, 0.0),
        sat_clock_m: 100.0,
        f1: 1575.42e6,
        f2: 1227.60e6,
        freq_num: 0,
        elevation_rad: 1.0,
        tropo_dry_m: 2.30,
        tropo_map_wet: 1.50,
        iono_l1_m: 4.50,
        variance_m2: 1.0,
        cp_variance_m2: 0.0001,
    };

    let norm = GeodeticNormalizations {
        sat_relativity_m: -8.859,
        shapiro_delay_m: 0.012,
        sat_pco_ecef_m: Vector3::new(1.5, 0.0, 0.0),
        solid_earth_tide_m: Vector3::new(0.15, 0.0, 0.0),
        ocean_tide_loading_m: Vector3::new(0.02, 0.0, 0.0),
    };

    let p_id = VariableId::new(1);
    let c_id = VariableId::new(2);
    let z_id = VariableId::new(3);
    let i_id = VariableId::new(4);
    let a_id = VariableId::new(5);

    let mut graph_vars = std::collections::BTreeMap::new();
    graph_vars.insert(p_id, crate::swfg::variables::VariableNode {
        id: p_id,
        kind: crate::swfg::variables::VariableKind::Pose { epoch: 0 },
        value: DVector::from_vec(vec![6_378_137.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    });
    graph_vars.insert(c_id, crate::swfg::variables::VariableNode {
        id: c_id,
        kind: crate::swfg::variables::VariableKind::ClockBias { epoch: 0, constellation_id: 0 },
        value: DVector::from_vec(vec![50.0]),
    });
    graph_vars.insert(z_id, crate::swfg::variables::VariableNode {
        id: z_id,
        kind: crate::swfg::variables::VariableKind::TropoZwd { epoch: 0 },
        value: DVector::from_vec(vec![0.20]),
    });
    graph_vars.insert(i_id, crate::swfg::variables::VariableNode {
        id: i_id,
        kind: crate::swfg::variables::VariableKind::IonosphereSlant { constellation_id: 0, epoch: 0, satellite: 1 },
        value: DVector::from_vec(vec![4.50]),
    });
    graph_vars.insert(a_id, crate::swfg::variables::VariableNode {
        id: a_id,
        kind: crate::swfg::variables::VariableKind::Ambiguity { constellation_id: 0, satellite: 1, frequency: 1, arc: 0 },
        value: DVector::from_vec(vec![12345.0]),
    });

    let values = VariableValues::build(&graph_vars);

    let wavelength_l1 = 299_792_458.0 / 1575.42e6;
    let rx_pos = Vector3::new(6_378_137.0 + 0.15 + 0.02, 0.0, 0.0);
    let sat_apc = Vector3::new(26_560_000.0 + 1.5, 0.0, 0.0);
    let geom_range = (sat_apc - rx_pos).norm();
    let windup_m = 0.05;
    let expected_cp_m = geom_range + 50.0 - 100.0 - 8.859 + 0.012 + (2.30 + 1.50 * 0.20) - 4.50 + wavelength_l1 * 12345.0;
    let raw_cp_cycles = (expected_cp_m + windup_m) / wavelength_l1;

    let factor = UducCarrierPhaseFactor::new(
        obs, 0, raw_cp_cycles, p_id, c_id, Some(z_id), Some(i_id), a_id, wavelength_l1, 1.0, windup_m, norm,
    );

    let res = factor.residual(&values);
    assert!(res[0].abs() < 1e-10, "Residual was {:.10}m", res[0]);

    let jac = factor.jacobian(&values);
    assert_eq!(jac[(0, 0)], 1.0);
    assert_eq!(jac[(0, 6)], -1.0);
    assert_eq!(jac[(0, 7)], -1.50);
    assert_eq!(jac[(0, 8)], 1.0); // Phase advance (+1.0 in Jacobian)
    assert!((jac[(0, 9)] - (-wavelength_l1)).abs() < 1e-10);
}
