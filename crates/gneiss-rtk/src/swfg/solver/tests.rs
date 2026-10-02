#![allow(clippy::unwrap_used)]

use super::*;
use crate::swfg::config::RtkConfig;
use crate::swfg::factor::{Factor, PriorFactor};
use crate::swfg::graph::MarginalPriorFactor;

#[test]
fn solver_creates_epoch_variables() {
    let config = EngineConfig::Rtk(RtkConfig::default());
    let mut solver = SlidingWindowSolver::new(&config);
    let vars = solver.create_epoch_variables(0, 10, &[0], false, false);
    assert_eq!(vars.len(), 3);
    assert_eq!(solver.graph.n_variables(), 3);
    assert_eq!(solver.graph.total_dim(), 6 + 1 + 1);
}

#[test]
fn solver_ensure_ambiguity_is_idempotent() {
    let config = EngineConfig::Rtk(RtkConfig::default());
    let mut solver = SlidingWindowSolver::new(&config);
    let id1 = solver.ensure_ambiguity(0, 1, 1, 0);
    let id2 = solver.ensure_ambiguity(0, 1, 1, 0);
    assert_eq!(id1, id2, "second call should return existing ambiguity");
    assert_eq!(solver.graph.n_variables(), 1);
}

#[test]
fn solver_extract_state_vector_matches_variable_values() {
    let config = EngineConfig::Rtk(RtkConfig::default());
    let mut solver = SlidingWindowSolver::new(&config);
    solver.create_epoch_variables(0, 1, &[0], false, false);
    let state = solver.extract_state_vector();
    let vars = VariableValues::build(&solver.graph.variables);
    assert_eq!(state.len(), vars.total_dim());
}

#[test]
fn solver_lm_converges_with_prior_factor() {
    use crate::swfg::factor::PriorFactor;
    let config = EngineConfig::Rtk(crate::swfg::config::RtkConfig::default());
    let mut solver = SlidingWindowSolver::new(&config);
    let var_ids = solver.create_epoch_variables(0, 1, &[0], false, false);
    let pose_id = var_ids[0];

    for &id in &var_ids {
        let dim = solver.graph.variables[&id].kind.dim().size();
        let mu = if id == pose_id {
            DVector::from_vec(vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0])
        } else {
            DVector::zeros(dim)
        };
        solver.graph.add_factor(Box::new(PriorFactor {
            variable: id,
            mu,
            information: DMatrix::identity(dim, dim),
        }));
    }

    solver.graph.set_value(pose_id, &[1.0, 2.0, 3.0, 0.0, 0.0, 0.0]);

    let result = solver.solve();
    assert!(result.is_ok(), "solve failed: {:?}", result.err());
    let _state = result.unwrap();

    let vals = VariableValues::build(&solver.graph.variables);
    let pose = vals.get(pose_id).unwrap();
    assert!((pose[0] - 1.0).abs() < 0.01);
    assert!((pose[1] - 2.0).abs() < 0.01);
    assert!((pose[2] - 3.0).abs() < 0.01);
}

// ---- configuration and epoch variables --------------------------------------

fn rtk(window: usize) -> EngineConfig {
    EngineConfig::Rtk(RtkConfig { window_size: window, ..RtkConfig::default() })
}

#[test]
fn window_size_is_taken_from_the_active_configuration() {
    // SPP is a single-epoch filter; every other mode uses its configured
    // window, floored at one epoch.
    assert_eq!(SlidingWindowSolver::new(&rtk(7)).window_size(), 7);
    assert_eq!(
        SlidingWindowSolver::new(&EngineConfig::Ppp(crate::swfg::config::PppConfig { window_size: 4, ..Default::default() })).window_size(),
        4
    );
    assert_eq!(
        SlidingWindowSolver::new(&EngineConfig::Spp(crate::swfg::config::SppConfig::default())).window_size(),
        1
    );
    assert_eq!(SlidingWindowSolver::new(&rtk(0)).window_size(), 1, "window size is floored at one epoch");
}

#[test]
fn an_imu_epoch_creates_velocity_and_exactly_one_session_bias() {
    let mut solver = SlidingWindowSolver::new(&rtk(10));
    let first = solver.create_epoch_variables(0, 8, &[0, 6], true, false);
    // pose, velocity, two clock biases, zwd, imu bias
    assert_eq!(first.len(), 6, "expected pose, velocity, 2 clocks, zwd and bias");
    let kinds: Vec<VariableKind> = first.iter().map(|&id| solver.graph.variables[&id].kind).collect();
    assert!(kinds.iter().any(|k| matches!(k, VariableKind::ImuBias)));
    assert_eq!(
        solver.graph.variables.values().filter(|n| matches!(n.kind, VariableKind::ImuBias)).count(),
        1,
        "the bias is per session, not per epoch"
    );
    // The orientation prior is skipped when an IMU supplies the attitude.
    assert!(
        !solver
            .graph
            .factors
            .iter()
            .any(|f| f.variables() == [first[0]]),
        "GNSS-only orientation prior must not be added alongside the IMU"
    );

    // A second epoch reuses the same bias state and does not re-add its prior.
    let n_before = solver.graph.n_factors();
    let second = solver.create_epoch_variables(1, 8, &[0, 6], true, false);
    assert_eq!(second.len(), 5, "no new bias variable in epoch 1");
    assert_eq!(solver.graph.n_factors(), n_before);
}

#[test]
fn double_difference_mode_creates_no_clock_or_zwd_state() {
    let mut solver = SlidingWindowSolver::new(&rtk(10));
    let vars = solver.create_epoch_variables(0, 8, &[0, 6], false, true);
    assert_eq!(vars.len(), 1, "double differencing removes receiver clocks and ZWD");
    assert!(matches!(solver.graph.variables[&vars[0]].kind, VariableKind::Pose { .. }));

    // Without IMU the unobserved orientation is pinned by a 1e8 prior.
    let prior = solver.graph.factors.iter().find(|f| f.variables() == [vars[0]]).expect("orientation prior");
    let info = prior.information();
    for i in 0..3 {
        assert_eq!(info[(i, i)], 0.0, "position must not be pinned");
    }
    for i in 3..6 {
        assert_eq!(info[(i, i)], 1e8, "orientation is unobservable without IMU");
    }
}

#[test]
fn session_wide_states_are_created_exactly_once() {
    let mut solver = SlidingWindowSolver::new(&rtk(10));
    let pose_a = solver.ensure_static_pose();
    let pose_b = solver.ensure_static_pose();
    assert_eq!(pose_a, pose_b);
    assert_eq!(solver.graph.variables.values().filter(|n| matches!(n.kind, VariableKind::StaticPose)).count(), 1);

    let ifb_a = solver.ensure_ifb_glonass();
    let ifb_b = solver.ensure_ifb_glonass();
    assert_eq!(ifb_a, ifb_b);
    assert_eq!(solver.graph.variables.values().filter(|n| matches!(n.kind, VariableKind::IfbGlonass)).count(), 1);
}

#[test]
fn dd_ambiguity_is_created_once_and_seeded_with_its_float_value() {
    let mut solver = SlidingWindowSolver::new(&rtk(10));
    let (first, is_new) = solver.ensure_dd_ambiguity(0, 5, 1, 1, 0, 12.5);
    assert!(is_new, "the first call creates the state");
    assert!((solver.graph.variables[&first].value[0] - 12.5).abs() < 1e-15, "seed value must be the float estimate");

    let (again, is_new) = solver.ensure_dd_ambiguity(0, 5, 1, 1, 0, 99.0);
    assert!(!is_new && again == first, "a second call must reuse the state and keep its value");
    assert!((solver.graph.variables[&first].value[0] - 12.5).abs() < 1e-15, "value must not be re-seeded");

    // A different reference satellite is a different ambiguity.
    let (other, is_new) = solver.ensure_dd_ambiguity(0, 5, 2, 1, 0, 12.5);
    assert!(is_new && other != first, "(sat, ref) pairs are distinct ambiguities");
}

// ---- state update -----------------------------------------------------------

#[test]
fn apply_delta_translates_position_and_rotates_attitude_on_so3() {
    let mut solver = SlidingWindowSolver::new(&rtk(10));
    let pose = solver.graph.add_variable(VariableKind::Pose { epoch: 0 });
    let zwd = solver.graph.add_variable(VariableKind::TropoZwd { epoch: 0 });
    solver.graph.set_value(pose, &[1.0, 2.0, 3.0, 0.0, 0.0, 0.0]);
    solver.graph.set_value(zwd, &[0.25]);

    let theta = std::f64::consts::FRAC_PI_2;
    solver.apply_delta(&DVector::from_vec(vec![1.0, -1.0, 0.5, 0.0, 0.0, theta, 0.75]));

    let p = solver.graph.variables[&pose].value.clone();
    assert!((p[0] - 2.0).abs() < 1e-12 && (p[1] - 1.0).abs() < 1e-12 && (p[2] - 3.5).abs() < 1e-12, "position {p:?}");
    // Starting from the identity, composing exp(theta * z_hat) leaves the
    // rotation vector exactly (0, 0, theta).
    assert!(p[3].abs() < 1e-12 && p[4].abs() < 1e-12, "yaw-only rotation {p:?}");
    assert!((p[5] - theta).abs() < 1e-12, "rot vec = {}", p[5]);
    assert!((solver.graph.variables[&zwd].value[0] - 1.0).abs() < 1e-12, "scalars are updated linearly");
}

#[test]
fn extract_state_vector_packs_variables_in_id_order() {
    let mut solver = SlidingWindowSolver::new(&rtk(10));
    let a = solver.graph.add_variable(VariableKind::TropoZwd { epoch: 0 });
    let b = solver.graph.add_variable(VariableKind::Velocity { epoch: 0 });
    solver.graph.set_value(a, &[7.0]);
    solver.graph.set_value(b, &[1.0, 2.0, 3.0]);
    let state = solver.extract_state_vector();
    assert_eq!(state.as_slice(), &[7.0, 1.0, 2.0, 3.0], "packing must follow VariableId order");
}

#[test]
fn epoch_counter_advances_monotonically() {
    let mut solver = SlidingWindowSolver::new(&rtk(10));
    assert_eq!(solver.current_epoch(), 0);
    solver.advance_epoch();
    solver.advance_epoch();
    assert_eq!(solver.current_epoch(), 2);
}

// ---- least squares ----------------------------------------------------------

/// A prior pulling `var` towards `mu` with information `w`, sized to the
/// variable as it exists in `solver`.
fn weighted_prior(solver: &SlidingWindowSolver, var: VariableId, mu: f64, w: f64) -> Box<dyn Factor> {
    let dim = solver.graph.variables[&var].value.len();
    Box::new(PriorFactor {
        variable: var,
        mu: DVector::from_element(dim, mu),
        information: DMatrix::identity(dim, dim) * w,
    })
}

#[test]
fn lm_recovers_the_weighted_least_squares_minimum() {
    // minimise (x - 3)^2 + 0.25 (x - 7)^2:
    //   d/dx = 2(x - 3) + 0.5(x - 7) = 0  ->  x = 9.5 / 2.5 = 3.8
    let mut solver = SlidingWindowSolver::new(&rtk(10));
    let x = solver.graph.add_variable(VariableKind::TropoZwd { epoch: 0 });
    solver.graph.add_factor(weighted_prior(&solver, x, 3.0, 1.0));
    solver.graph.add_factor(weighted_prior(&solver, x, 7.0, 0.25));
    solver.graph.set_value(x, &[0.0]);

    solver.solve().expect("solve");
    let got = solver.graph.variables[&x].value[0];
    assert!((got - 3.8).abs() < 1e-6, "weighted least-squares minimum is 3.8, got {got}");
}

#[test]
fn solve_refuses_a_graph_with_an_orphan_variable() {
    let mut solver = SlidingWindowSolver::new(&rtk(10));
    let x = solver.graph.add_variable(VariableKind::TropoZwd { epoch: 0 });
    solver.graph.add_factor(weighted_prior(&solver, x, 1.0, 1.0));
    let _lonely = solver.graph.add_variable(VariableKind::TropoZwd { epoch: 1 });
    match solver.solve() {
        Err(SolveError::OrphanVariable(msg)) => assert!(msg.contains("Orphan"), "message was {msg}"),
        other => panic!("expected OrphanVariable, got {other:?}"),
    }
}

#[test]
fn total_error_is_the_weighted_sum_of_factor_costs() {
    // At x = 5: (5 - 3)^2 + 0.25 (5 - 7)^2 = 4 + 1 = 5.
    let mut solver = SlidingWindowSolver::new(&rtk(10));
    let x = solver.graph.add_variable(VariableKind::TropoZwd { epoch: 0 });
    solver.graph.add_factor(weighted_prior(&solver, x, 3.0, 1.0));
    solver.graph.add_factor(weighted_prior(&solver, x, 7.0, 0.25));
    solver.graph.set_value(x, &[5.0]);
    let err = solver.evaluate_total_error(&VariableValues::build(&solver.graph.variables));
    assert!((err - 5.0).abs() < 1e-12, "total error = {err}");
}

#[test]
fn normal_equations_add_the_marginal_prior_hessian_and_gradient() {
    // One 6-DOF pose with a marginal prior H = 2 I, g = 0, x0 = (1, 0, ...).
    // With the state at the origin the offset is d = x - x0 = (-1, 0, ...), so
    //   Jtr = g + H d = (2 * -1, 0, ...) = (-2, 0, ...)
    //   cost = 1/2 d' H d + g' d = 1/2 * 2 * 1 = 1
    let mut solver = SlidingWindowSolver::new(&rtk(10));
    let pose = solver.graph.add_variable(VariableKind::Pose { epoch: 1 });
    solver.graph.add_factor(weighted_prior(&solver, pose, 0.0, 0.0)); // keeps the graph connected, W = 0
    solver.graph.marginal_prior = Some(MarginalPriorFactor {
        variables: vec![(pose, 6)],
        hessian: DMatrix::identity(6, 6) * 2.0,
        gradient: DVector::zeros(6),
        x0: DVector::from_vec(vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    });

    let values = VariableValues::build(&solver.graph.variables);
    let (jtj, jtr, cost) = solver.build_normal_equations(&values);
    for i in 0..6 {
        assert!((jtj[(i, i)] - 2.0).abs() < 1e-12, "H[{i},{i}] = {}", jtj[(i, i)]);
    }
    assert!((jtr[0] + 2.0).abs() < 1e-12, "Jtr[0] = {}", jtr[0]);
    assert!(jtr.rows(1, 5).iter().all(|v| v.abs() < 1e-12), "only the offset axis moves");
    assert!((cost - 1.0).abs() < 1e-12, "marginal prior cost = {cost}");

    // evaluate_total_error must price the same prior the same way.
    assert!((solver.evaluate_total_error(&values) - 1.0).abs() < 1e-12);
}

#[test]
fn marginalizing_an_empty_window_is_a_noop() {
    let mut solver = SlidingWindowSolver::new(&rtk(10));
    assert!(solver.marginalize_oldest_epoch(0).is_ok());
    assert_eq!(solver.graph.n_factors(), 0);
    assert!(solver.graph.marginal_prior.is_none());
}
