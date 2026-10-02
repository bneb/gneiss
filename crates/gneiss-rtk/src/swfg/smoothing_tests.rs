use super::*;
use crate::swfg::factor::{compute_robust_error, Factor};
use nalgebra::DMatrix;
use crate::swfg::graph::EstimationGraph;
use gneiss_core::ephemeris::{Ephemeris, GpsEphemeris};
use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};

/// Batch cost `E(x) = sum_i rho(||r_i||)`, the quantity the LM solver minimises.
/// Written from the definition in `factor.rs` using the public loss helpers, so
/// it can check the *result* of a solve without re-running the solver.
pub(crate) fn batch_cost(graph: &EstimationGraph) -> f64 {
let values = VariableValues::build(&graph.variables);
graph.factors.iter().map(|f| f.as_ref()) .map(|f: &dyn Factor| {
let r = f.residual(&values);
let w = f.information();
let raw = (&r.transpose() * &w * &r)[(0, 0)];
match f.robust_threshold() {
None => raw,
Some(k) => {
let n = r.norm();
if n < 1e-12 { raw } else { raw * compute_robust_error(n, k, f.use_cauchy()) / (n * n) }
}
}
}).sum()
}

fn make_gps_eph(sat: SatelliteId, toc: GpsTime, m0: f64) -> Ephemeris {
Ephemeris::Gps(GpsEphemeris {
sat,
toe: toc,
toc,
af0: 0.0,
af1: 0.0,
af2: 0.0,
crs: 0.0,
crc: 0.0,
cuc: 0.0,
cus: 0.0,
cic: 0.0,
cis: 0.0,
m0,
e: 0.001,
sqrt_a: 26_560_000.0_f64.sqrt(),
delta_n: 0.0,
omega0: 0.0,
omega_dot: 0.0,
i0: 0.96,
idot: 0.0,
omega: 0.0,
tgd: 0.0,
iode: 0,
iodc: 0,
})
}

/// Eight GPS ephemerides at the epoch time, with mean anomalies spread evenly
/// over one revolution so the constellation has usable geometry (the same
/// construction `swfg/engine/tests.rs` uses for its end-to-end epoch test).
pub(crate) fn constellation_n(t: GpsTime, n: u8) -> Vec<Ephemeris> {
(1..=n)
.map(|prn| {
let sat = SatelliteId { constellation: Constellation::Gps, prn };
let m0 = (prn as f64 - 1.0) * std::f64::consts::TAU / n as f64;
make_gps_eph(sat, t, m0)
})
.collect()
}

pub(crate) fn constellation(t: GpsTime) -> Vec<Ephemeris> {
constellation_n(t, 8)
}

pub(crate) fn l1_wavelength() -> f64 {
gneiss_core::constants::SPEED_OF_LIGHT_M_S / 1575.42e6
}

pub(crate) fn make_test_obs(time: GpsTime) -> EpochObs {
make_obs_n(time, 8)
}

/// As `make_obs_n`, with a constant added to every pseudorange (the carrier
/// phase follows so the two stay consistent to within one cycle of `delta`).
pub(crate) fn make_obs_offset(time: GpsTime, n: u8, offset_m: f64) -> EpochObs {
let lambda = l1_wavelength();
let mut obs = EpochObs { time, satellites: Vec::new() };
for prn in 1..=n {
let sat = SatelliteId { constellation: Constellation::Gps, prn };
let range = 20_000_000.0 + (prn as f64) * 100.0 + offset_m;
obs.satellites.push(SatObs {
sat,
observations: vec![
Observation {
code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 1, attribute: 'C' } },
value: range,
lock_time: None,
lli: None,
},
Observation {
code: ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band: 1, attribute: 'C' } },
value: range / lambda,
lock_time: None,
lli: None,
},
],
});
}
obs
}

pub(crate) fn make_obs_n(time: GpsTime, n: u8) -> EpochObs {
let lambda = l1_wavelength();
let mut obs = EpochObs { time, satellites: Vec::new() };
for prn in 1..=n {
let sat = SatelliteId { constellation: Constellation::Gps, prn };
let range = 20_000_000.0 + (prn as f64) * 100.0;
obs.satellites.push(SatObs {
sat,
observations: vec![
Observation {
code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 1, attribute: 'C' } },
value: range,
lock_time: None,
lli: None,
},
Observation {
code: ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band: 1, attribute: 'C' } },
value: range / lambda,
lock_time: None,
lli: None,
},
],
});
}
obs
}

#[test]
fn test_batch_smoothing_scaffolding() {
let config = BatchSmoothingConfig::default();
let mut batch = BatchFactorGraph::new(config, Vec::new());
let t0 = GpsTime::new(2200, 100.0);
let obs0 = make_test_obs(t0);
let p0 = batch.add_epoch(0, &obs0, None, None);
assert_eq!(batch.pose_ids.len(), 1);
assert_eq!(batch.pose_ids[0], p0);
}

#[test]
fn default_config_matches_its_documented_values() {
let c = BatchSmoothingConfig::default();
assert_eq!(c.max_iterations, 20);
assert_eq!(c.convergence_tol, 1e-4);
assert!(c.enable_ar);
assert_eq!(c.ar_ratio_threshold, 2.0);
}

#[test]
fn new_graph_starts_empty() {
let batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), Vec::new());
assert!(batch.epochs.is_empty());
assert!(batch.pose_ids.is_empty());
assert!(batch.graph.variables.is_empty());
assert!(batch.graph.factors.is_empty());
assert!(batch.graph.validate_graph_structure().is_ok());
}

/// Epoch 0 is anchored so the batch has an absolute datum; later epochs are
/// carried by their measurement factors and the smoothness chain alone.
#[test]
fn only_the_first_epoch_is_anchored_by_a_prior() {
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), constellation(GpsTime::new(2200, 100.0)));
for k in 0..3u32 {
batch.add_epoch(k, &make_test_obs(GpsTime::new(2200, 100.0 + k as f64)), None, None);
}
let priors = batch.graph.factors.iter().filter(|f| f.variables().len() == 1
&& f.residual(&VariableValues::build(&batch.graph.variables)).len() == 6).count();
assert!(priors >= 1, "epoch 0 must be anchored");
assert_eq!(priors, 1, "only epoch 0 may carry a 6-DOF anchor");
}

#[test]
fn add_epoch_seeds_pose_clock_and_zwd_without_a_base() {
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), Vec::new());
let t = GpsTime::new(2200, 100.0);
let pose = batch.add_epoch(0, &make_test_obs(t), None, None);

let node = batch.graph.variables.get(&pose).expect("pose inserted");
assert!(matches!(node.kind, VariableKind::Pose { epoch: 0 }));
assert_eq!(node.value.len(), 6);

let clocks: Vec<_> = batch.graph.variables.values()
.filter(|n| matches!(n.kind, VariableKind::ClockBias { .. })).collect();
assert_eq!(clocks.len(), 1, "one receiver clock per epoch");
assert!(batch.graph.variables.values().any(|n| matches!(n.kind, VariableKind::TropoZwd { epoch: 0 })),
"an SPP-style epoch without a base must estimate its own zenith wet delay");
}

#[test]
fn add_epoch_with_a_base_does_not_create_a_troposphere_state() {
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), constellation(GpsTime::new(2200, 100.0)));
let t = GpsTime::new(2200, 100.0);
let base_pos = Vector3::new(-3963427.0, 3350882.0, 3694866.0);
let base_obs = make_test_obs(t);
batch.add_epoch(0, &base_obs, Some((&base_obs, base_pos)), None);
assert!(!batch.graph.variables.values().any(|n| matches!(n.kind, VariableKind::TropoZwd { .. })),
"a double-differenced epoch takes its troposphere from the base station");
}

#[test]
fn add_epoch_honours_the_supplied_initial_position() {
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), Vec::new());
let truth = Vector3::new(6_378_137.0, -100.0, 4_500_000.0);
let pose = batch.add_epoch(0, &make_test_obs(GpsTime::new(2200, 100.0)), None, Some(truth));
let v = batch.graph.variables.get(&pose).unwrap().value.clone();
assert_eq!(v[0], truth.x);
assert_eq!(v[1], truth.y);
assert_eq!(v[2], truth.z);
assert_eq!(v.rows(3, 3).into_owned(), DVector::zeros(3), "attitude always starts level");
}

#[test]
fn add_epoch_falls_back_to_the_default_datum() {
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), Vec::new());
let pose = batch.add_epoch(0, &make_test_obs(GpsTime::new(2200, 100.0)), None, None);
let v = batch.graph.variables.get(&pose).unwrap().value.clone();
assert_eq!(v[0], -3963427.0);
assert_eq!(v[1], 3350882.0);
assert_eq!(v[2], 3694866.0);
}

#[test]
fn add_epoch_appends_epochs_and_poses_in_order() {
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), Vec::new());
let mut times = Vec::new();
let mut ids = Vec::new();
for k in 0..4u32 {
let t = GpsTime::new(2200, 100.0 + k as f64);
times.push(t);
ids.push(batch.add_epoch(k, &make_test_obs(t), None, None));
}
assert_eq!(batch.epochs, times);
assert_eq!(batch.pose_ids, ids);
assert!(ids.windows(2).all(|w| w[0] < w[1]), "pose ids must be distinct and increasing");
assert_eq!(batch.graph.variables.values().filter(|n| matches!(n.kind, VariableKind::Pose { .. })).count(), 4);
}

#[test]
fn ensure_ambiguity_reuses_a_variable_for_the_same_key() {
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), Vec::new());
let a = batch.ensure_ambiguity(0, 7, 1, 0);
let b = batch.ensure_ambiguity(0, 7, 1, 0);
assert_eq!(a, b, "the same arc/satellite/frequency must map to one variable");
let c = batch.ensure_ambiguity(0, 8, 1, 0);
let d = batch.ensure_ambiguity(0, 7, 2, 0);
let e = batch.ensure_ambiguity(0, 7, 1, 1);
let f = batch.ensure_ambiguity(1, 7, 1, 0);
for other in [c, d, e, f] {
assert_ne!(other, a, "a different key must mint a new variable");
}
assert_eq!(batch.graph.variables.values().filter(|n| matches!(n.kind, VariableKind::Ambiguity { .. })).count(), 5);
}

// ------------------------------------------------------- smoothness factors

#[test]
fn smoothness_variance_is_velocity_times_dt_squared() {
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), Vec::new());
for k in 0..4u32 {
batch.add_epoch(k, &make_test_obs(GpsTime::new(2200, 100.0 + k as f64)), None, None);
}
let before = batch.graph.factors.len();
batch.add_smoothness_factors(10.0, 1.0);
assert_eq!(batch.graph.factors.len(), before + 3, "n epochs get n-1 links");
let links: Vec<_> = batch.graph.factors.iter().skip(before).collect();
assert_eq!(links.len(), 3);
for (i, link) in links.iter().enumerate() {
let vars = link.variables();
assert_eq!(vars.len(), 2);
assert_eq!(vars[0], batch.pose_ids[i]);
assert_eq!(vars[1], batch.pose_ids[i + 1]);
// W = I(6) / (v*dt)^2 = I(6) / 100.
assert_eq!(link.information(), DMatrix::<f64>::identity(6, 6) * 0.01);
}
}

#[test]
fn smoothness_variance_scales_with_the_square_of_the_speed_limit() {
let mut a = BatchFactorGraph::new(BatchSmoothingConfig::default(), Vec::new());
let mut b = BatchFactorGraph::new(BatchSmoothingConfig::default(), Vec::new());
for batch in [&mut a, &mut b] {
for k in 0..2u32 {
batch.add_epoch(k, &make_test_obs(GpsTime::new(2200, 100.0 + k as f64)), None, None);
}
}
a.add_smoothness_factors(5.0, 1.0);
b.add_smoothness_factors(5.0, 2.0);
// dt = 2 s doubles the allowed step, so the variance quadruples and the
// information drops to a quarter: 1/(5*2)^2 vs 1/(5*1)^2.
let ia = a.graph.factors.last().unwrap().information();
let ib = b.graph.factors.last().unwrap().information();
assert_eq!(ia, DMatrix::<f64>::identity(6, 6) / 25.0);
assert_eq!(ib, DMatrix::<f64>::identity(6, 6) / 100.0);
assert_eq!(ia[(0, 0)], 4.0 * ib[(0, 0)]);
}

#[test]
fn smoothness_on_one_epoch_links_nothing() {
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), Vec::new());
batch.add_epoch(0, &make_test_obs(GpsTime::new(2200, 100.0)), None, None);
let before = batch.graph.factors.len();
batch.add_smoothness_factors(10.0, 1.0);
assert_eq!(batch.graph.factors.len(), before);
}

#[test]
fn smoothness_on_an_empty_graph_is_a_no_op() {
// `0..(len - 1)` must not underflow when no epoch has been added yet.
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), Vec::new());
assert!(batch.pose_ids.is_empty());
batch.add_smoothness_factors(10.0, 1.0);
assert!(batch.graph.factors.is_empty());
}

#[test]
fn smoothness_links_are_relative_pose_factors() {
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), Vec::new());
for k in 0..3u32 {
batch.add_epoch(k, &make_test_obs(GpsTime::new(2200, 100.0 + k as f64)), None, None);
}
batch.add_smoothness_factors(10.0, 1.0);
// Every new factor must be a two-pose RelativePoseFactor.
let new_factors: Vec<&Box<dyn Factor>> = batch.graph.factors.iter()
.filter(|f| f.variables().len() == 2).collect();
assert_eq!(new_factors.len(), 2);
for f in new_factors {
assert_eq!(f.information().nrows(), 6);
let vars = VariableValues::build(&batch.graph.variables);
assert_eq!(f.residual(&vars).len(), 6);
}
}
