#![allow(clippy::unwrap_used)]

use super::*;
use crate::swfg::config::RtkConfig;

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
    let id1 = solver.ensure_ambiguity(1, 1, 0);
    let id2 = solver.ensure_ambiguity(1, 1, 0);
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
