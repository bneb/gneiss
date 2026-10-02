//! Sprint-hunt tests driving `execute_ar_step` and the post-fit slip detector
//! end-to-end through a real sliding-window factor graph.
#![allow(clippy::unwrap_used)]

use std::collections::HashMap;

use nalgebra::{DVector, Vector3};

use crate::swfg::ar_integration::{collect_ambiguity_variables, extract_ambiguity_state};
use crate::swfg::config::EngineConfig;
use crate::swfg::engine::ar_handler::{
    check_postfit_cycle_slips, execute_ar_step, reset_ppp_tracker,
};
use crate::swfg::engine::builder::CpMeasurementRecord;
use crate::swfg::factor::PriorFactor;
use crate::swfg::pipeline::dd_factors::DdCarrierPhaseFactor;
use crate::swfg::solver::SlidingWindowSolver;
use crate::swfg::variables::{VariableId, VariableKind, VariableValues};

use super::sprint_common::{BASE, LAMBDA_L1, TRUE_POS};

// ===========================================================================
// execute_ar_step through a real factor graph
// ===========================================================================

/// Build a graph whose double-difference carrier phases are exactly consistent
/// with integer ambiguities at `TRUE_POS`.
///
/// Hand arithmetic: with pose = TRUE_POS and ambiguity = N,
///   dd_cp_obs = (|p_i - TRUE| - |p_ref - TRUE|) - base_dd + lambda * N
/// makes every `DdCarrierPhaseFactor` residual exactly zero at (TRUE, N).
///
/// Each ambiguity additionally carries a weak prior centred `prior_bias_cycles`
/// away from N, which is what a code-derived float looks like before AR. That
/// keeps the system determined (8 phase rows alone cannot constrain 14
/// unknowns) AND guarantees the float ambiguities are genuinely fractional, so
/// the test below cannot pass vacuously.
pub(super) fn rtk_graph(
    n_sats: usize,
    prior_bias_cycles: f64,
    pose_offset_m: f64,
) -> (SlidingWindowSolver, VariableId, Vec<VariableId>) {
    let mut solver = SlidingWindowSolver::new(&EngineConfig::default());
    let pose = solver.graph.add_variable(VariableKind::Pose { epoch: 0 });
    solver.graph.set_value(
        pose,
        &[TRUE_POS.x + pose_offset_m, TRUE_POS.y, TRUE_POS.z, 0.0, 0.0, 0.0],
    );

    let mut sats = Vec::new();
    for i in 0..n_sats {
        let az = 2.0 * std::f64::consts::PI * (i as f64) / (n_sats as f64);
        sats.push(
            TRUE_POS + Vector3::new(az.cos(), az.sin(), 0.0).normalize() * 2.0e7
                + Vector3::new(0.0, 0.0, 1.0e7),
        );
    }
    let ref_pos = sats[0];

    let mut amb_ids = Vec::new();
    for (i, &p_i) in sats.iter().enumerate().skip(1) {
        let rover_dd = (p_i - TRUE_POS).norm() - (ref_pos - TRUE_POS).norm();
        let base_dd = (p_i - BASE).norm() - (ref_pos - BASE).norm();
        let n = 100_000.0 + 1000.0 * (i as f64);
        let amb = solver.ensure_dd_ambiguity(0, i as u16, 0, 1, 0, n);
        solver.graph.add_factor(Box::new(DdCarrierPhaseFactor {
            var_pose: pose,
            var_amb: amb.0,
            dd_cp_obs_m: rover_dd - base_dd + LAMBDA_L1 * n,
            sat_pos: p_i,
            ref_pos,
            base_pos: BASE,
            base_dd_range: base_dd,
            lambda: LAMBDA_L1,
            variance_m2: 1e-6,
            elevation_rad: 1.2,
            ref_elevation_rad: 1.2,
            is_new_amb: amb.1,
            variables: vec![pose, amb.0],
        }));
        // Weak prior: 1 sigma = 0.1 cycle from the true integer.
        // NOTE: PriorFactor::new's third argument is a VARIANCE, not an
        // information (factor.rs:62-65).
        solver.graph.add_factor(Box::new(PriorFactor::new(
            amb.0,
            DVector::from_element(1, n + prior_bias_cycles),
            0.01,
        )));
        amb_ids.push(amb.0);
    }
    // Pose prior centred on the offset (variance 9 m^2 => information 1/9).
    solver.graph.add_factor(Box::new(PriorFactor::new(
        pose,
        DVector::from_vec(vec![
            TRUE_POS.x + pose_offset_m,
            TRUE_POS.y,
            TRUE_POS.z,
            0.0,
            0.0,
            0.0,
        ]),
        9.0,
    )));
    (solver, pose, amb_ids)
}

#[test]
fn execute_ar_step_snaps_the_pose_to_the_true_integer_consistent_position() {
    reset_ppp_tracker();
    let (mut solver, pose, amb_ids) = rtk_graph(9, 0.4, 0.3);
    solver.solve().expect("float solve must succeed");
    let before = VariableValues::build(&solver.graph.variables).get(pose).unwrap()[0];
    assert!(
        (before - (TRUE_POS.x + 0.3)).abs() < 1.0,
        "float must still sit near the offset, got {before}"
    );
    // The float ambiguities must genuinely be fractional, otherwise the AR
    // assertions below would be satisfied without AR running at all.
    for id in &amb_ids {
        let v = VariableValues::build(&solver.graph.variables).get(*id).unwrap()[0];
        assert!(
            (v - v.round()).abs() > 1e-3,
            "float ambiguity {v} is already an integer; the fixture is vacuous"
        );
    }
    // Sanity: the ambiguities must be inside the selection window that
    // `execute_ar_step` applies, otherwise this test would prove nothing.
    let values = VariableValues::build(&solver.graph.variables);
    let (jtj, _, _) = solver.build_normal_equations(&values);
    let ids = collect_ambiguity_variables(&solver.graph);
    let (_, cov) = extract_ambiguity_state(&solver.graph, &jtj, &ids);
        assert!(ids.len() >= amb_ids.len());
    for i in 0..ids.len() {
        assert!(
            cov[(i, i)] > 1e-6 && cov[(i, i)] < 2.0,
            "ambiguity {i} variance {} is outside the (1e-6, 2.0) AR window; \
             the Cauchy robust weight is suppressing the phase factors",
            cov[(i, i)]
        );
    }

    execute_ar_step(&mut solver, pose, BASE, true, 1);

    // A successful DD-AR step is defined by the ambiguities being EXACT
    // integers afterwards. This is the strongest available statement and it is
    // independent of where the pose lands.
    let after = VariableValues::build(&solver.graph.variables);
    for id in &amb_ids {
        let v = after.get(*id).unwrap()[0];
        assert!(
            (v - v.round()).abs() < 1e-6,
            "ambiguity {id:?} = {v} is not an integer after a successful AR fix"
        );
    }
    let pos = after.get(pose).unwrap()[0];
    assert!(
        (pos - TRUE_POS.x).abs() < 0.2,
        "AR did not snap the pose to the true position: {pos} vs {}",
        TRUE_POS.x
    );
}

/// The rejection path pops factors off the end of the factor vector and
/// re-solves. Any off-by-one there corrupts the window for every later epoch,
/// so the graph must still be structurally valid and solvable afterwards.
#[test]
fn graph_still_solves_and_validates_after_repeated_ar_steps() {
    reset_ppp_tracker();
    let (mut solver, pose, _) = rtk_graph(9, 0.4, 0.3);
    solver.solve().expect("initial solve");
    for epoch in 1..=5u32 {
        execute_ar_step(&mut solver, pose, BASE, true, epoch);
        solver
            .solve()
            .unwrap_or_else(|e| panic!("graph failed to solve after AR at epoch {epoch}: {e:?}"));
        solver
            .graph
            .validate_graph_structure()
            .unwrap_or_else(|e| panic!("graph structure invalid after AR at epoch {epoch}: {e}"));
    }
}

/// A second AR step on an already-fixed graph must be a no-op in the pose: the
/// ambiguities are integers, LAMBDA returns the same set, and the pose cannot
/// move again. Any drift here means the fixed solution is not a fixed point.
#[test]
fn a_second_ar_step_on_an_already_fixed_graph_is_a_pose_fixed_point() {
    reset_ppp_tracker();
    let (mut solver, pose, _) = rtk_graph(9, 0.4, 0.3);
    solver.solve().expect("float solve");
    execute_ar_step(&mut solver, pose, BASE, true, 1);
    let first = VariableValues::build(&solver.graph.variables).get(pose).unwrap()[0];
    execute_ar_step(&mut solver, pose, BASE, true, 2);
    let second = VariableValues::build(&solver.graph.variables).get(pose).unwrap()[0];
    assert!(
        (second - first).abs() < 1e-9,
        "AR is not idempotent: pose moved {first} -> {second}"
    );
}

/// The integration-level zero-false-fix test.
///
/// Eight ambiguities are given priors centred exactly halfway between two
/// integers (N + 0.5) with a 0.1-cycle standard deviation, and no carrier
/// factor constrains them further. Hand arithmetic: for each ambiguity the two
/// nearest integers are at +-0.5 cycles, i.e. 5 sigma away and equally good, so
/// LAMBDA's best and second-best candidates are exactly equidistant and the
/// ratio is 1.0. `execute_ar_step` uses `min_ratio = max(3.0, ffrt_threshold)`,
/// so the step MUST be refused and the graph MUST come out byte-identical.
#[test]
fn ar_step_refuses_a_fix_the_float_cannot_support_and_leaves_the_graph_intact() {
    reset_ppp_tracker();
    let mut solver = SlidingWindowSolver::new(&EngineConfig::default());
    let pose = solver.graph.add_variable(VariableKind::Pose { epoch: 0 });
    solver
        .graph
        .set_value(pose, &[TRUE_POS.x, TRUE_POS.y, TRUE_POS.z, 0.0, 0.0, 0.0]);
    let mut amb_ids = Vec::new();
    for i in 1..=8u16 {
        let n = 100_000.0 + 1000.0 * (i as f64);
        let id = solver
            .ensure_dd_ambiguity(0, i, 0, 1, 0, n + 0.5)
            .0;
        solver.graph.add_factor(Box::new(PriorFactor::new(
            id,
            DVector::from_element(1, n + 0.5),
            0.01,
        )));
        amb_ids.push(id);
    }
    solver
        .graph
        .add_factor(Box::new(PriorFactor::new(pose, DVector::from_element(6, 0.0), 9.0)));
    solver.solve().expect("solve");

    let n_before = solver.graph.n_factors();
    let vals_before = VariableValues::build(&solver.graph.variables);
    let pos_before = vals_before.get(pose).unwrap()[0];

    execute_ar_step(&mut solver, pose, BASE, true, 1);

    let vals_after = VariableValues::build(&solver.graph.variables);
    assert_eq!(
        solver.graph.n_factors(),
        n_before,
        "a refused fix must not inject or remove any factor"
    );
    assert!((vals_after.get(pose).unwrap()[0] - pos_before).abs() < 1e-12, "pose moved");
    for id in &amb_ids {
        let v = vals_after.get(*id).unwrap()[0];
        assert!(
            (v - (v + 0.5).round()).abs() > 1e-3,
            "ambiguity {id:?} moved to {v}; the float was never fixed"
        );
    }
    solver.solve().expect("graph must still solve after a refused AR step");
}

// ===========================================================================
// check_postfit_cycle_slips
// ===========================================================================

/// The detector's model must agree with `DdCarrierPhaseFactor::residual`:
/// `predicted = (|p_i - pos| - |p_ref - pos|) - base_dd + lambda * amb`.
fn slip_record(dd_cp_m: f64) -> (SlidingWindowSolver, VariableId, CpMeasurementRecord) {
    let p_i = TRUE_POS + Vector3::new(1.0e7, 5.0e6, 8.0e6);
    let p_ref = TRUE_POS + Vector3::new(-6.0e6, 1.1e7, 1.4e7);
    let base_dd = (p_i - BASE).norm() - (p_ref - BASE).norm();
    let geom = (p_i - TRUE_POS).norm() - (p_ref - TRUE_POS).norm() - base_dd;
    let true_amb = 100.0_f64;

    let mut solver = SlidingWindowSolver::new(&EngineConfig::default());
    let pose = solver.graph.add_variable(VariableKind::Pose { epoch: 0 });
    solver
        .graph
        .set_value(pose, &[TRUE_POS.x, TRUE_POS.y, TRUE_POS.z, 0.0, 0.0, 0.0]);
    let amb = solver.ensure_dd_ambiguity(0, 5, 1, 2, 0, true_amb);
    solver.graph.add_factor(Box::new(PriorFactor::new(
        amb.0,
        DVector::from_element(1, true_amb),
        1e8,
    )));
    solver
        .graph
        .add_factor(Box::new(PriorFactor::new(pose, DVector::from_element(6, 0.0), 1.0)));

    let rec = CpMeasurementRecord {
        sat: 5,
        constellation_id: 0,
        var_amb: amb.0,
        dd_cp_m: geom + LAMBDA_L1 * true_amb + dd_cp_m,
        base_dd,
        sat_pos: p_i,
        ref_pos: p_ref,
        lambda: LAMBDA_L1,
    };
    (solver, amb.0, rec)
}

#[test]
fn postfit_slip_check_is_silent_on_a_zero_residual() {
    let (mut solver, amb, rec) = slip_record(0.0);
    let mut counts = HashMap::new();
    check_postfit_cycle_slips(&mut solver, TRUE_POS, &[rec], &mut counts);
    assert!(counts.is_empty(), "a zero residual must not flag a slip: {counts:?}");
    assert!(solver.graph.variables.contains_key(&amb), "ambiguity must survive");
}

#[test]
fn postfit_slip_check_respects_its_half_metre_threshold() {
    let (mut solver, amb, rec) = slip_record(0.49);
    let mut counts = HashMap::new();
    check_postfit_cycle_slips(&mut solver, TRUE_POS, &[rec], &mut counts);
    assert!(counts.is_empty(), "0.49 m is below the 0.5 m gate");
    assert!(solver.graph.variables.contains_key(&amb));

    let (mut solver, amb, rec) = slip_record(0.51);
    let mut counts = HashMap::new();
    check_postfit_cycle_slips(&mut solver, TRUE_POS, &[rec], &mut counts);
    assert_eq!(counts.get(&(0, 5)).copied(), Some(1), "0.51 m must flag exactly one slip");
    assert!(!solver.graph.variables.contains_key(&amb), "slipped ambiguity must be dropped");
}

#[test]
fn postfit_slip_check_catches_a_ten_metre_residual() {
    let (mut solver, amb, rec) = slip_record(10.0);
    let mut counts = HashMap::new();
    check_postfit_cycle_slips(&mut solver, TRUE_POS, &[rec], &mut counts);
    assert_eq!(counts.get(&(0, 5)).copied(), Some(1), "slip count {counts:?}");
    assert!(!solver.graph.variables.contains_key(&amb));
}
