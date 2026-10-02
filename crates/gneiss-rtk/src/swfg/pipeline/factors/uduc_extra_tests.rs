use super::uduc_tests::{base_norm, base_obs, l1_wavelength, uduc_values};
use super::*;
use crate::swfg::variables::{VariableKind, VariableNode};

#[test]
fn uduc_cp_windup_is_subtracted_from_the_observation() {
// Repo convention (AGENTS.md): phase windup is always *subtracted*.
// observed_cp_m = cycles * lambda - windup, so flipping the windup sign must
// move the residual by exactly 2 * windup.
let lambda = l1_wavelength();
let values = uduc_values([6_378_137.0, 0., 0., 0., 0., 0.], 0.0, 0.0, 0.0);
let mk = |windup: f64| UducCarrierPhaseFactor::new(
base_obs(1.0, 1e-4), 0, 100_000_000.0, VariableId::new(1), VariableId::new(2),
Some(VariableId::new(3)), Some(VariableId::new(4)), VariableId::new(5), lambda, 1.0, windup, base_norm(),
);
// A larger windup removes more range from the observation, so the residual
// must fall: r(-0.05) - r(+0.05) = +0.10 exactly.
let plus = mk(0.05).residual(&values)[0];
let minus = mk(-0.05).residual(&values)[0];
// The residuals here are O(1e6) m (1e8 cycles of L1 carrier minus a modelled
// 2e4 km range), so f64 carries ~1e-9 m of absolute rounding at that magnitude.
// 1e-6 m still catches a flipped sign, which would show up as a 0.2 m error.
assert!((minus - plus - 0.10).abs() < 1e-6, "plus={plus} minus={minus}");
// Zero windup sits exactly halfway.
let zero = mk(0.0).residual(&values)[0];
assert!((zero - (plus + minus) / 2.0).abs() < 1e-6);
}

#[test]
fn uduc_cp_iono_enters_with_the_opposite_sign_to_pseudorange() {
// Code is delayed by +mu*I, carrier phase advances by -mu*I, so the two
// factors must disagree on the sign of the ionosphere column.
let lambda = l1_wavelength();
let values = uduc_values([6_378_137.0, 0., 0., 0., 0., 0.], 0.0, 0.0, 4.0);
let cp = UducCarrierPhaseFactor::new(
base_obs(1.0, 1e-4), 0, 100_000_000.0, VariableId::new(1), VariableId::new(2), None, Some(VariableId::new(4)),
VariableId::new(5), lambda, 2.0, 0.0, base_norm(),
);
let pr = UducPseudorangeFactor::new(
base_obs(1.0, 1e-4), 0, 0.0, VariableId::new(1), VariableId::new(2), None, Some(VariableId::new(4)), 2.0, base_norm(),
);
let (s_iono, _) = values.index_of(VariableId::new(4)).expect("iono offset");
assert_eq!(cp.jacobian(&values)[(0, s_iono)], 2.0);
assert_eq!(pr.jacobian(&values)[(0, s_iono)], -2.0);
// A 4 m L1 ionosphere delay raises the modelled code by mu*I = 8 m and lowers
// the modelled carrier by the same 8 m; r = observed - modelled, so the code
// residual FALLS by 8 m while the carrier residual RISES by 8 m.
let with_iono = uduc_values([6_378_137.0, 0., 0., 0., 0., 0.], 0.0, 0.0, 4.0);
let no_iono = uduc_values([6_378_137.0, 0., 0., 0., 0., 0.], 0.0, 0.0, 0.0);
assert!((cp.residual(&with_iono)[0] - cp.residual(&no_iono)[0] - 8.0).abs() < 1e-9);
assert!((pr.residual(&with_iono)[0] - pr.residual(&no_iono)[0] + 8.0).abs() < 1e-9);
}

#[test]
fn uduc_cp_ambiguity_scales_by_the_wavelength() {
let lambda = l1_wavelength();
let f = UducCarrierPhaseFactor::new(
base_obs(1.0, 1e-4), 0, 100_000_000.0, VariableId::new(1), VariableId::new(2), None, None,
VariableId::new(5), lambda, 1.0, 0.0, base_norm(),
);
// The ambiguity is not in `values`, so it reads as 0; the Jacobian still
// reports -lambda for its column once the variable exists.
let mut vars = std::collections::BTreeMap::new();
vars.insert(
VariableId::new(1),
VariableNode { id: VariableId::new(1), kind: VariableKind::Pose { epoch: 0 }, value: DVector::from_vec(vec![6_378_137.0, 0., 0., 0., 0., 0.]) },
);
vars.insert(
VariableId::new(5),
VariableNode { id: VariableId::new(5), kind: VariableKind::Ambiguity { constellation_id: 0, satellite: 1, frequency: 1, arc: 0 }, value: DVector::from_element(1, 10.0) },
);
let with_amb = VariableValues::build(&vars);
let n0 = UducCarrierPhaseFactor::new(
base_obs(1.0, 1e-4), 0, 100_000_000.0, VariableId::new(1), VariableId::new(2), None, None,
VariableId::new(9), lambda, 1.0, 0.0, base_norm(),
).residual(&with_amb)[0];
let n10 = f.residual(&with_amb)[0];
// modeled_cp carries +lambda * N and r = observed - modeled, so ten extra
// cycles on the ambiguity lower the residual by exactly 10 * lambda.
assert!((n0 - n10 - 10.0 * lambda).abs() < 1e-6, "n0={n0} n10={n10}");
let (s_amb2, dim_amb) = with_amb.index_of(VariableId::new(5)).expect("amb offset");
assert_eq!(dim_amb, 1);
assert_eq!(s_amb2 + dim_amb, with_amb.total_dim(), "the ambiguity is the last, 1-column block");
let (s_amb, _) = with_amb.index_of(VariableId::new(5)).expect("amb offset");
assert!((f.jacobian(&with_amb)[(0, s_amb)] + lambda).abs() < 1e-15);
}

#[test]
fn uduc_cp_information_uses_the_carrier_variance_with_a_1e_6_floor() {
// Note the carrier floor is 1e-6 m^2, ten times tighter than the pseudorange
// factor's 1e-4 m^2 floor.
let mk = |v: f64| UducCarrierPhaseFactor::new(
base_obs(1.0, v), 0, 100_000_000.0, VariableId::new(1), VariableId::new(2), None, None,
VariableId::new(5), l1_wavelength(), 1.0, 0.0, base_norm(),
);
assert_eq!(mk(1e-4).information()[(0, 0)], 1e4);
assert_eq!(mk(4e-4).information()[(0, 0)], 2500.0);
assert_eq!(mk(0.0).information()[(0, 0)], 1e6);
}

#[test]
fn uduc_cp_degenerates_gracefully_when_the_ambiguity_is_missing() {
// No ambiguity variable: it reads as zero, the residual is still well defined,
// and the Jacobian drops the ambiguity column entirely.
let mut vars = std::collections::BTreeMap::new();
vars.insert(
VariableId::new(1),
VariableNode { id: VariableId::new(1), kind: VariableKind::Pose { epoch: 0 }, value: DVector::from_vec(vec![6_378_137.0, 0., 0., 0., 0., 0.]) },
);
vars.insert(
VariableId::new(2),
VariableNode { id: VariableId::new(2), kind: VariableKind::ClockBias { epoch: 0, constellation_id: 0 }, value: DVector::from_element(1, 0.0) },
);
let values = VariableValues::build(&vars);
let f = UducCarrierPhaseFactor::new(
base_obs(1.0, 1e-4), 0, 100_000_000.0, VariableId::new(1), VariableId::new(2), None, None,
VariableId::new(5), l1_wavelength(), 1.0, 0.0, base_norm(),
);
assert!(f.residual(&values)[0].is_finite());
let j = f.jacobian(&values);
let (s_clk, _) = values.index_of(VariableId::new(2)).expect("clock offset");
assert_eq!((j.nrows(), j.ncols()), (1, 7));
assert_eq!(j[(0, s_clk)], -1.0, "the clock column survives");
}

// --------------------------------------------------- slant iono random walk

#[test]
fn slant_iono_variance_is_the_process_noise_times_dt() {
// new() floors dt at 0.1 s so a zero-length epoch cannot divide by zero, and
// the variance itself at 1e-5 m^2 so a zero q_I cannot blow the information up.
let f = SlantIonoRandomWalkFactor::new(VariableId::new(1), VariableId::new(2), 10.0, 2e-4);
assert!((f.variance - 2e-3).abs() < 1e-15, "variance = {}", f.variance);
assert!((f.information()[(0, 0)] - 500.0).abs() < 1e-9);
// dt below the 0.1 s floor is clamped.
let floored = SlantIonoRandomWalkFactor::new(VariableId::new(1), VariableId::new(2), 0.01, 2e-4);
assert!((floored.variance - 2e-5).abs() < 1e-15, "variance = {}", floored.variance);
assert!((floored.information()[(0, 0)] - 5e4).abs() < 1e-6);
// A zero process noise floors at 1e-5 m^2 rather than producing an infinity.
let zero = SlantIonoRandomWalkFactor::new(VariableId::new(1), VariableId::new(2), 10.0, 0.0);
assert_eq!(zero.variance, 1e-5);
assert!((zero.information()[(0, 0)] - 1e5).abs() < 1e-6);
assert_eq!(zero.variables(), &[VariableId::new(1), VariableId::new(2)]);
}

#[test]
fn slant_iono_residual_is_the_epoch_to_epoch_change() {
let mut vars = std::collections::BTreeMap::new();
vars.insert(
VariableId::new(1),
VariableNode { id: VariableId::new(1), kind: VariableKind::IonosphereSlant { constellation_id: 0, epoch: 0, satellite: 1 }, value: DVector::from_element(1, 4.0) },
);
vars.insert(
VariableId::new(2),
VariableNode { id: VariableId::new(2), kind: VariableKind::IonosphereSlant { constellation_id: 0, epoch: 1, satellite: 1 }, value: DVector::from_element(1, 6.5) },
);
let values = VariableValues::build(&vars);
let f = SlantIonoRandomWalkFactor::new(VariableId::new(1), VariableId::new(2), 10.0, 2e-4);
assert!((f.residual(&values)[0] - 2.5).abs() < 1e-15);
let j = f.jacobian(&values);
assert_eq!((j.nrows(), j.ncols()), (1, 2));
assert_eq!(j[(0, 0)], -1.0);
assert_eq!(j[(0, 1)], 1.0);
}

#[test]
fn slant_iono_degrades_to_zero_when_a_state_is_missing() {
// Only the current epoch is in the graph; the previous reads as zero.
let mut vars = std::collections::BTreeMap::new();
vars.insert(
VariableId::new(2),
VariableNode { id: VariableId::new(2), kind: VariableKind::IonosphereSlant { constellation_id: 0, epoch: 1, satellite: 1 }, value: DVector::from_element(1, 6.5) },
);
let values = VariableValues::build(&vars);
let f = SlantIonoRandomWalkFactor::new(VariableId::new(1), VariableId::new(2), 10.0, 2e-4);
assert!((f.residual(&values)[0] - 6.5).abs() < 1e-15);
let j = f.jacobian(&values);
assert_eq!((j.nrows(), j.ncols()), (1, 1));
assert_eq!(j[(0, 0)], 1.0, "the surviving state keeps its +1 column");
}
