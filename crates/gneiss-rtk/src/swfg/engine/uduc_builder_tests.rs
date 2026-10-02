//! Golden-vector tests for the UDUC factor builder.
//!
//! Geometry is chosen so every hand computation is exact: the receiver sits
//! on the WGS84 equator / prime meridian and the satellite on the +X axis at
//! GPS orbit radius, so the ECEF line of sight is
//! `RANGE = 26 560 000 - 6 378 137 = 20 181 863 m` exactly, with the receiver
//! clock, satellite clock, relativity, Shapiro and troposphere terms all
//! zeroed.  Every expected value below is a literal worked out by hand in the
//! test that uses it.
#![allow(clippy::unwrap_used)]

use super::*;
use crate::swfg::config::{EngineConfig, RtkConfig};
use crate::swfg::variables::VariableValues;
use gneiss_core::time::GpsTime;

pub(crate) const RX_X: f64 = 6_378_137.0;
pub(crate) const SAT_X: f64 = 26_560_000.0;
/// `|sat - rx|` on the +X axis: 26 560 000 - 6 378 137 = 20 181 863 m.
pub(crate) const RANGE: f64 = SAT_X - RX_X;
pub(crate) const F1: f64 = 1575.42e6;
pub(crate) const F2: f64 = 1227.60e6;
/// (f1/f2)^2 = (1575.42/1227.60)^2.  1575.42/1227.60 = 77/60 exactly, since
/// 77 x 1227.60 = 94 525.2 and 60 x 1575.42 = 94 525.2, so gamma = 5929/3600.
pub(crate) const GAMMA: f64 = 5929.0 / 3600.0;
/// Slant ionosphere delay used throughout: 10.000 m.
pub(crate) const IONO: f64 = 10.0;
/// Wavelengths from the defining relation lambda = c/f.
pub(crate) const LAMBDA1: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / F1;
pub(crate) const LAMBDA2: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / F2;

pub(crate) fn rx_pos() -> Vector3<f64> {
    Vector3::new(RX_X, 0.0, 0.0)
}

fn rover_time() -> GpsTime {
    GpsTime::new(2000, 100.0)
}

/// Carrier phase in cycles for an integer ambiguity of exactly `n`, i.e.
/// `cp * lambda = pr + n * lambda`.
pub(crate) fn cp_cycles(pr_m: f64, lambda: f64, n: f64) -> f64 {
    (pr_m + n * lambda) / lambda
}

/// Satellite on the +X axis, code range set to `RANGE + IONO` so the L1
/// pseudorange model is exactly satisfied.
///
/// Note: with the satellite on the +X axis the receiver-dipole projection
/// `d_rx = north - k (k . north) + k x east` is identically zero (k = -x_hat,
/// north = z_hat, east = y_hat gives z_hat + (-x_hat) x y_hat = 0), so the Wu
/// model returns a windup of exactly 0 for this geometry.  Tests that need a
/// non-zero windup use [`obs_at`] with the satellite on the +Y axis.
pub(crate) fn obs(satellite: u16, elevation_rad: f64) -> CorrectedObservation {
    obs_at(satellite, elevation_rad, Vector3::new(SAT_X, 0.0, 0.0))
}

/// As [`obs`], with the satellite placed at `sat_pos`.  `pr_l1` is always
/// `RANGE + IONO`; tests that move the satellite only assert on the ambiguity
/// seeding, which does not depend on the geometry.
pub(crate) fn obs_at(satellite: u16, elevation_rad: f64, sat_pos: Vector3<f64>) -> CorrectedObservation {
    CorrectedObservation {
        satellite,
        constellation_id: 0,
        pr_l1: RANGE + IONO,
        pr_l2: None,
        cp_l1: None,
        cp_l1_lli: None,
        cp_l2: None,
        doppler: 0.0,
        snr_dbhz: 45.0,
        sat_pos_ecef: sat_pos,
        sat_clock_m: 0.0,
        f1: F1,
        f2: F2,
        freq_num: 0,
        elevation_rad,
        tropo_dry_m: 0.0,
        tropo_map_wet: 0.0,
        iono_l1_m: IONO,
        variance_m2: 1.0,
        cp_variance_m2: 1.0,
    }
}

pub(crate) struct Graph {
    pub(crate) solver: SlidingWindowSolver,
    pose: VariableId,
    zwd: VariableId,
}

impl Graph {
    /// A graph with one clock-bias state per epoch in `clock_epochs`, one
    /// pose, and one zenith wet delay, all seeded to zero except the pose.
    pub(crate)     fn new(clock_epochs: &[u32]) -> Self {
        let mut solver = SlidingWindowSolver::new(&EngineConfig::Rtk(RtkConfig::default()));
        let pose = solver.graph.add_variable(VariableKind::Pose { epoch: 0 });
        for &e in clock_epochs {
            let id = solver.graph.add_variable(VariableKind::ClockBias { epoch: e, constellation_id: 0 });
            solver.graph.set_value(id, &[0.0]);
        }
        let zwd = solver.graph.add_variable(VariableKind::TropoZwd { epoch: 0 });
        solver.graph.set_value(zwd, &[0.0]);
        solver.graph.set_value(pose, &[RX_X, 0.0, 0.0, 0.0, 0.0, 0.0]);
        Self { solver, pose, zwd }
    }

    pub(crate)     fn values(&self) -> VariableValues {
        VariableValues::build(&self.solver.graph.variables)
    }

    /// All ambiguity states as `(id, frequency * 1000 + arc)`, ascending by id.
    pub(crate)     fn ambiguities(&self) -> Vec<(VariableId, u32)> {
        self.solver
            .graph
            .variables
            .iter()
            .filter_map(|(id, n)| match n.kind {
                VariableKind::Ambiguity { frequency, arc, .. } => Some((*id, frequency as u32 * 1000 + arc)),
                _ => None,
            })
            .collect()
    }

    /// The slant ionosphere state created for `epoch`.  Must be matched on
    /// the epoch: a previous-epoch state also lives in the graph and has a
    /// lower `VariableId`, so a plain "first IonosphereSlant" lookup returns
    /// the wrong variable.
    pub(crate)     fn iono_var(&self, epoch: u32) -> Option<VariableId> {
        self.solver
            .graph
            .variables
            .iter()
            .find(|(_, n)| matches!(n.kind, VariableKind::IonosphereSlant { epoch: e, .. } if e == epoch))
            .map(|(id, _)| *id)
    }

    pub(crate)     fn n_slant(&self) -> usize {
        self.solver
            .graph
            .variables
            .values()
            .filter(|n| matches!(n.kind, VariableKind::IonosphereSlant { .. }))
            .count()
    }
}

/// Split the factor list by role, identified from the variable count: a
/// pseudorange factor spans [pose, clock, zwd, iono] (4 variables) and a
/// carrier factor adds the ambiguity (5).
pub(crate) fn factor_indices(g: &Graph) -> (Vec<usize>, Vec<usize>) {
    let mut pr = Vec::new();
    let mut cp = Vec::new();
    for (i, f) in g.solver.graph.factors.iter().enumerate() {
        match f.variables().len() {
            4 => pr.push(i),
            5 => cp.push(i),
            _ => {}
        }
    }
    (pr, cp)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run(
    g: &mut Graph,
    corrected: &[CorrectedObservation],
    epoch: u32,
    prev_epoch: Option<u32>,
    dt_sec: f64,
    slips: &mut HashMap<(u8, u16), u32>,
    trackers: &mut HashMap<(u8, u16), gneiss_geodesy::windup::PhaseWindupTracker>,
) {
    let mut ctx = UducFactorContext {
        epoch,
        pose_id: g.pose,
        zwd_id: Some(g.zwd),
        prev_epoch,
        dt_sec,
        norm: GeodeticNormalizations::default(),
        slip_counts: slips,
        windup_trackers: trackers,
        rover_time: rover_time(),
        rx_pos: rx_pos(),
    };
    build_uduc_factors(&mut g.solver, corrected, &mut ctx);
}

/// Phase windup for the standard geometry, from an INDEPENDENTLY constructed
/// tracker.
///
/// The Wu windup model is an *input* to the builder, not the thing under test.
/// A fresh `PhaseWindupTracker` fed the same satellite and Sun positions
/// supplies the reference value.  The basis vectors are written out literally
/// (a receiver on the equator / prime meridian gives up = +X, north = +Z,
/// east = +Y) rather than taken from `local_enu_basis`, so a basis error
/// cannot mask itself here.
pub(crate) fn reference_windup_rad(sat_pos: &Vector3<f64>) -> f64 {
    let (sun, _) = gneiss_geodesy::tides::solar_lunar_positions(rover_time().tow, rover_time().week);
    let mut tracker = gneiss_geodesy::windup::PhaseWindupTracker::new();
    tracker.update(
        sat_pos,
        &sun,
        &rx_pos(),
        &Vector3::new(1.0, 0.0, 0.0),
        &Vector3::new(0.0, 0.0, 1.0),
        &Vector3::new(0.0, 1.0, 0.0),
    )
}

// ---- local ENU basis --------------------------------------------------------

#[test]
fn enu_basis_on_equator_and_prime_meridian_is_exact() {
    // lat = lon = 0  =>  up    = (cos0cos0, cos0sin0, sin0) = (1, 0, 0)
    //                    north = (-sin0cos0, -sin0sin0, cos0) = (0, 0, 1)
    //                    east  = (-sin0, cos0, 0)             = (0, 1, 0)
    let (up, north, east) = local_enu_basis(rx_pos());
    assert!((up - Vector3::new(1.0, 0.0, 0.0)).norm() < 1e-12, "up = {up:?}");
    assert!((north - Vector3::new(0.0, 0.0, 1.0)).norm() < 1e-12, "north = {north:?}");
    assert!((east - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12, "east = {east:?}");
}

#[test]
fn enu_basis_is_orthonormal_and_left_handed_everywhere() {
    // With rows (up, north, east), R^T R must be the identity.  Hand expansion
    // of the cross product (s = sin lat, c = cos lat, l = lon):
    //   up x north = (c sin l, -c cos l, 0) = -(-sin l, cos l, 0) = -east
    // so the triple is LEFT handed: det(R) = -1, not +1.
    let positions = [
        gneiss_core::coords::llh_to_ecef(Vector3::new(0.0, 0.0, 0.0)),
        gneiss_core::coords::llh_to_ecef(Vector3::new(45.0_f64.to_radians(), 0.0, 500.0)),
        gneiss_core::coords::llh_to_ecef(Vector3::new(-30.0_f64.to_radians(), 2.1, -300.0)),
        gneiss_core::coords::llh_to_ecef(Vector3::new(89.0_f64.to_radians(), -3.0, 1.0e4)),
    ];
    for pos in positions {
        let (up, north, east) = local_enu_basis(pos);
        assert!((up.norm() - 1.0).abs() < 1e-12, "up not unit at {pos:?}");
        assert!((north.norm() - 1.0).abs() < 1e-12, "north not unit at {pos:?}");
        assert!((east.norm() - 1.0).abs() < 1e-12, "east not unit at {pos:?}");
        let rows = [up, north, east];
        for r in 0..3 {
            for c in 0..3 {
                let dot: f64 = rows[r].iter().zip(rows[c].iter()).map(|(a, b)| a * b).sum();
                let want = if r == c { 1.0 } else { 0.0 };
                assert!((dot - want).abs() < 1e-12, "R^T R [{r},{c}] = {dot} at {pos:?}");
            }
        }
        let det = up.dot(&north.cross(&east));
        assert!((det + 1.0).abs() < 1e-12, "det(up, north, east) = {det} at {pos:?}, must be -1");
    }
}

#[test]
fn enu_up_is_the_geodetic_normal_not_the_geocentric_direction() {
    // At geodetic latitude 45 deg the geocentric latitude is
    //   atan((1 - e^2) tan 45 deg) = atan(0.9933056) = 44.8082 deg,
    // i.e. 0.1918 deg lower, so `up` must differ from normalize(rx) by
    // ~0.192 deg.  Getting this wrong biases the windup dipole projection.
    let lat = 45.0_f64.to_radians();
    let pos = gneiss_core::coords::llh_to_ecef(Vector3::new(lat, 0.0, 0.0));
    let (up, _, _) = local_enu_basis(pos);
    let want = Vector3::new(std::f64::consts::FRAC_1_SQRT_2, 0.0, std::f64::consts::FRAC_1_SQRT_2);
    assert!((up - want).norm() < 1e-12, "up = {up:?}, want {want:?}");
    let separation_deg = pos.normalize().angle(&up).to_degrees();
    assert!(
        (0.19..0.20).contains(&separation_deg),
        "geodetic vs geocentric separation was {separation_deg} deg, expected ~0.192 deg"
    );
}

// ---- pseudorange factors ----------------------------------------------------

#[test]
fn pseudorange_only_epoch_builds_one_factor_and_seeds_the_slant_iono() {
    let mut g = Graph::new(&[0]);
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    run(&mut g, &[obs(1, 1.0)], 0, None, 1.0, &mut slips, &mut trackers);

    let (pr, cp) = factor_indices(&g);
    assert_eq!(pr.len(), 1, "expected one pseudorange factor, got {pr:?}");
    assert!(cp.is_empty(), "no carrier phase was supplied, got {cp:?}");
    // pr_l1 = RANGE + IONO and the model is
    //   range + c*dt - sat_clock + rel + shapiro + tropo + 1.0 * iono
    //     = 20181863 + 0 - 0 + 0 + 0 + 0 + 10 = 20 181 873
    let v = g.values();
    let r = g.solver.graph.factors[pr[0]].residual(&v);
    assert!(r[0].abs() < 1e-9, "pseudorange residual = {}", r[0]);

    let iono = g.iono_var(0).expect("slant ionosphere state must exist");
    assert!(
        (g.solver.graph.variables[&iono].value[0] - IONO).abs() < 1e-15,
        "slant iono must be seeded from the observation"
    );
    assert!(g.ambiguities().is_empty(), "no carrier phase means no ambiguity");
}

#[test]
fn l2_pseudorange_is_scaled_by_the_squared_frequency_ratio() {
    // The L2 model adds gamma * iono with gamma = 5929/3600, so
    //   5929/3600 * 10 = 16.469 444 4... m
    // and pr_l2 = 20 181 879.469 444 4 gives a zero residual.  A builder that
    // used gamma = 1 would leave 6.469 m; one using f1/f2 unsquared, 3.909 m.
    let mut g = Graph::new(&[0]);
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    let mut o = obs(2, 1.0);
    o.pr_l2 = Some(RANGE + GAMMA * IONO);
    run(&mut g, &[o], 0, None, 1.0, &mut slips, &mut trackers);

    let (pr, _) = factor_indices(&g);
    assert_eq!(pr.len(), 2, "L1 and L2 pseudorange factors expected");
    let v = g.values();
    for (i, &idx) in pr.iter().enumerate() {
        let r = g.solver.graph.factors[idx].residual(&v)[0];
        assert!(r.abs() < 1e-6, "pseudorange residual {i} = {r}");
    }
}

#[test]
fn observation_without_a_clock_bias_for_this_epoch_is_skipped() {
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    // The clock state exists for epoch 5 but the builder is asked for epoch 7.
    let mut g = Graph::new(&[0]);
    run(&mut g, &[obs(3, 1.0)], 7, None, 1.0, &mut slips, &mut trackers);
    assert_eq!(g.solver.graph.n_factors(), 0, "no factors without a clock state");
    assert!(g.iono_var(7).is_none(), "no state may be created for a skipped observation");

    // Same graph, matching epoch: now it must build.
    let mut g = Graph::new(&[7]);
    run(&mut g, &[obs(3, 1.0)], 7, None, 1.0, &mut slips, &mut trackers);
    assert_eq!(g.solver.graph.n_factors(), 1);
}
