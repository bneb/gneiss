use super::*;
use crate::swfg::variables::{VariableKind, VariableNode};
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

// ------------------------------------------------------------------ helpers

pub(crate) fn base_obs(variance: f64, cp_variance: f64) -> CorrectedObservation {
CorrectedObservation {
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
variance_m2: variance,
cp_variance_m2: cp_variance,
}
}

pub(crate) fn base_norm() -> GeodeticNormalizations {
GeodeticNormalizations {
sat_relativity_m: -8.859,
shapiro_delay_m: 0.012,
sat_pco_ecef_m: Vector3::new(1.5, 0.0, 0.0),
solid_earth_tide_m: Vector3::new(0.15, 0.0, 0.0),
ocean_tide_loading_m: Vector3::new(0.02, 0.0, 0.0),
}
}

pub(crate) fn l1_wavelength() -> f64 {
gneiss_core::constants::SPEED_OF_LIGHT_M_S / 1575.42e6
}

/// Four scalars in the packing order pose(6) | clock(1) | zwd(1) | iono(1).
pub(crate) fn uduc_values(pose: [f64; 6], clock: f64, zwd: f64, iono: f64) -> VariableValues {
let mut vars = std::collections::BTreeMap::new();
vars.insert(
VariableId::new(1),
VariableNode { id: VariableId::new(1), kind: VariableKind::Pose { epoch: 0 }, value: DVector::from_vec(pose.to_vec()) },
);
vars.insert(
VariableId::new(2),
VariableNode { id: VariableId::new(2), kind: VariableKind::ClockBias { epoch: 0, constellation_id: 0 }, value: DVector::from_element(1, clock) },
);
vars.insert(
VariableId::new(3),
VariableNode { id: VariableId::new(3), kind: VariableKind::TropoZwd { epoch: 0 }, value: DVector::from_element(1, zwd) },
);
vars.insert(
VariableId::new(4),
VariableNode { id: VariableId::new(4), kind: VariableKind::IonosphereSlant { constellation_id: 0, epoch: 0, satellite: 1 }, value: DVector::from_element(1, iono) },
);
VariableValues::build(&vars)
}

// ------------------------------------------------------------- pseudorange

#[test]
fn uduc_pr_variable_order_follows_the_constructor() {
let f = UducPseudorangeFactor::new(
base_obs(1.0, 1e-4), 0, 0.0, VariableId::new(1), VariableId::new(2),
Some(VariableId::new(3)), Some(VariableId::new(4)), 1.0, base_norm(),
);
assert_eq!(f.variables(), &[VariableId::new(1), VariableId::new(2), VariableId::new(3), VariableId::new(4)]);
// Optional states omitted: the list shrinks but the order stays pose, clock.
let no_opt = UducPseudorangeFactor::new(
base_obs(1.0, 1e-4), 0, 0.0, VariableId::new(1), VariableId::new(2), None, None, 1.0, base_norm(),
);
assert_eq!(no_opt.variables(), &[VariableId::new(1), VariableId::new(2)]);
// Ionosphere only.
let iono_only = UducPseudorangeFactor::new(
base_obs(1.0, 1e-4), 0, 0.0, VariableId::new(1), VariableId::new(2), None, Some(VariableId::new(4)), 1.0, base_norm(),
);
assert_eq!(iono_only.variables(), &[VariableId::new(1), VariableId::new(2), VariableId::new(4)]);
}

#[test]
fn uduc_pr_information_is_inverse_variance_with_a_1e_4_floor() {
// variance_m2.max(1e-4): 1 m^2 -> W = 1, and a degenerate 0 -> W = 1e4.
let mk = |v: f64| UducPseudorangeFactor::new(
base_obs(v, 1e-4), 0, 0.0, VariableId::new(1), VariableId::new(2), None, None, 1.0, base_norm(),
);
assert_eq!(mk(1.0).information()[(0, 0)], 1.0);
assert_eq!(mk(0.25).information()[(0, 0)], 4.0);
assert_eq!(mk(0.0).information()[(0, 0)], 1e4);
assert_eq!(mk(-3.0).information()[(0, 0)], 1e4);
}

#[test]
fn uduc_pr_iono_scaling_follows_the_frequency_ratio() {
// On L2 the ionospheric code delay is (f1/f2)^2 times the L1 value:
//   1575.42 / 1227.60 = 1.2833333... (exactly 77/60 in MHz terms)
//   (77/60)^2 = 5929/3600 = 1.64694444...
let mu = (1575.42e6 / 1227.60e6f64).powi(2);
assert!((mu - 5929.0 / 3600.0).abs() < 1e-12, "mu = {mu}");
let values = uduc_values([6_378_137.0, 0., 0., 0., 0., 0.], 0.0, 0.0, 4.0);
let f = UducPseudorangeFactor::new(
base_obs(1.0, 1e-4), 0, 0.0, VariableId::new(1), VariableId::new(2), None, Some(VariableId::new(4)), mu, base_norm(),
);
// With no clock and no zwd, r = raw_pr - (range + mu*iono), so the residual
// moves by exactly 4*mu when the raw observation is moved by 4*mu.
let r0 = f.residual(&values)[0];
let shifted = UducPseudorangeFactor::new(
base_obs(1.0, 1e-4), 0, 4.0 * mu, VariableId::new(1), VariableId::new(2), None, Some(VariableId::new(4)), mu, base_norm(),
);
assert!((shifted.residual(&values)[0] - r0 - 4.0 * mu).abs() < 1e-9);
// And the Jacobian entry is -mu for a code measurement: r = raw - modeled
// and the modeled code carries +mu * I, so dr/dI = -mu.
let j = f.jacobian(&values);
let (s_iono, _) = values.index_of(VariableId::new(4)).expect("iono offset");
assert_eq!(s_iono, 8, "pose(6) | clock(1) | zwd(1) | iono(1)");
assert!((j[(0, s_iono)] + mu).abs() < 1e-15, "J iono = {}", j[(0, s_iono)]);
}

#[test]
fn uduc_pr_clock_and_zwd_signs_are_negative() {
// modeled_pr = range + c*dt - clk + ... + dry + map*zwd + mu*iono and
// r = raw - modeled, so every state enters the residual with a minus sign.
let values = uduc_values([6_378_137.0, 0., 0., 0., 0., 0.], 0.0, 0.0, 0.0);
let f = UducPseudorangeFactor::new(
base_obs(1.0, 1e-4), 0, 0.0, VariableId::new(1), VariableId::new(2), Some(VariableId::new(3)), None, 1.0, base_norm(),
);
let j = f.jacobian(&values);
let (s_clk, _) = values.index_of(VariableId::new(2)).expect("clock offset");
let (s_zwd, _) = values.index_of(VariableId::new(3)).expect("zwd offset");
assert_eq!(j[(0, s_clk)], -1.0, "clock");
assert_eq!(j[(0, s_zwd)], -1.50, "zwd through the 1.5 wet mapping function");
// A +1 m clock bias must lower the residual by exactly 1 m.
let biased = uduc_values([6_378_137.0, 0., 0., 0., 0., 0.], 1.0, 0.0, 0.0);
assert!((f.residual(&biased)[0] - f.residual(&values)[0] + 1.0).abs() < 1e-9);
}

#[test]
fn uduc_pr_jacobian_range_is_the_unit_line_of_sight() {
// The fixture satellite sits on the +X axis and the receiver at the origin,
// so the unit line of sight is exactly (1, 0, 0): dr/dx = +1, the rest zero.
let values = uduc_values([0.0, 0.0, 0.0, 0., 0., 0.], 0.0, 0.0, 0.0);
let f = UducPseudorangeFactor::new(
base_obs(1.0, 1e-4), 0, 0.0, VariableId::new(1), VariableId::new(2), None, None, 1.0, GeodeticNormalizations::default(),
);
let j = f.jacobian(&values);
assert_eq!(j[(0, 0)], 1.0);
assert_eq!(j[(0, 1)], 0.0);
assert_eq!(j[(0, 2)], 0.0);
for c in 3..6 {
assert_eq!(j[(0, c)], 0.0, "attitude must not enter a range-only pseudorange");
}
}

#[test]
fn uduc_pr_degenerates_gracefully_when_the_pose_is_missing() {
// With no pose the receiver collapses to the ECEF origin and every optional
// state falls back to zero, so the residual is still well defined:
//   r = raw - (|sat_apc| - sat_clock + rel + shapiro + tropo_dry)
// sat_apc = (26 560 000 + 1.5, 0, 0) (no PCO/earth tide when norm is default).
let mut vars = std::collections::BTreeMap::new();
vars.insert(
VariableId::new(2),
VariableNode { id: VariableId::new(2), kind: VariableKind::ClockBias { epoch: 0, constellation_id: 0 }, value: DVector::from_element(1, 0.0) },
);
let values = VariableValues::build(&vars);
let f = UducPseudorangeFactor::new(
base_obs(1.0, 1e-4), 0, 0.0, VariableId::new(1), VariableId::new(2), None, None, 1.0, GeodeticNormalizations::default(),
);
// `raw_pr_m` is the constructor argument (0.0 here), not `obs.pr_l1`.
let modeled: f64 = 26_560_000.0 - 100.0 + 2.30;
let want = -modeled;
let got = f.residual(&values)[0];
assert!((got + modeled).abs() < 1e-6, "residual {got}, want {want}");
// The Jacobian keeps the clock column and the position columns are simply
// left at zero because the pose offset cannot be resolved.
let j = f.jacobian(&values);
assert_eq!((j.nrows(), j.ncols()), (1, 1));
assert_eq!(j[(0, 0)], -1.0);
}

// ----------------------------------------------------------- carrier phase
