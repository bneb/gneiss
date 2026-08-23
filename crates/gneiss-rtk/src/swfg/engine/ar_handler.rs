//! Ambiguity Resolution (AR) and post-fit cycle slip handling for SWFG engine.

use std::collections::HashMap;
use nalgebra::Vector3;

use crate::swfg::ar_integration::{
    attempt_ar_fix, collect_ambiguity_variables, extract_ambiguity_state, inject_fixed_priors,
    validate_fix_geometry, ArResult,
};
use crate::swfg::engine::builder::CpMeasurementRecord;
use crate::swfg::solver::SlidingWindowSolver;
use crate::swfg::variables::{VariableId, VariableValues};

/// Execute LAMBDA + PAR ambiguity resolution step on active factor graph.
pub fn execute_ar_step(
    solver: &mut SlidingWindowSolver,
    pose_id: VariableId,
    init_pos: Vector3<f64>,
    is_rtk: bool,
    epoch: u32,
) {
    if (!is_rtk && !epoch.is_multiple_of(5)) || epoch == 0 {
        return;
    }
    let all_amb_ids = collect_ambiguity_variables(&solver.graph);
    if all_amb_ids.len() < 4 {
        return;
    }

    let values = VariableValues::build(&solver.graph.variables);
    let (jtj, _, _) = solver.build_normal_equations(&values);
    let (float_amb_all, amb_cov_all) =
        extract_ambiguity_state(&solver.graph, &jtj, &all_amb_ids);

    let mut indexed_vars: Vec<(usize, f64)> = (0..all_amb_ids.len())
        .map(|i| (i, amb_cov_all[(i, i)]))
        .filter(|(_, var)| *var > 1e-6 && *var < 2.0)
        .collect();
    indexed_vars.sort_by(|a, b| a.1.total_cmp(&b.1));

    let selected_indices: Vec<usize> = indexed_vars.iter().take(8).map(|(i, _)| *i).collect();
    if selected_indices.len() < 4 {
        return;
    }

    let amb_ids: Vec<_> = selected_indices.iter().map(|&i| all_amb_ids[i]).collect();
    let float_amb = nalgebra::DVector::from_iterator(
        amb_ids.len(),
        selected_indices.iter().map(|&i| float_amb_all[i]),
    );
    let mut amb_cov = nalgebra::DMatrix::zeros(amb_ids.len(), amb_ids.len());
    for (r_idx, &r) in selected_indices.iter().enumerate() {
        for (c_idx, &c) in selected_indices.iter().enumerate() {
            amb_cov[(r_idx, c_idx)] = amb_cov_all[(r, c)];
        }
    }

    let ffrt_threshold = crate::ambiguity::ffrt::calculate_threshold(amb_ids.len(), 0.001);
    let min_ratio = 3.0_f64.max(ffrt_threshold);
    let ar_result = attempt_ar_fix(&float_amb, &amb_cov, min_ratio);

    if let ArResult::Fixed { integers, .. } = ar_result {
        apply_validated_fix(solver, pose_id, init_pos, &amb_ids, &integers);
    }
}

/// Applies fixed priors and validates the resulting position jump.
fn apply_validated_fix(
    solver: &mut SlidingWindowSolver,
    pose_id: VariableId,
    init_pos: Vector3<f64>,
    amb_ids: &[VariableId],
    integers: &nalgebra::DVector<f64>,
) {
    let vals_before = VariableValues::build(&solver.graph.variables);
    let float_pos = vals_before
        .get(pose_id)
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .unwrap_or(init_pos);

    let float_amb_values: Vec<f64> = amb_ids
        .iter()
        .map(|id| solver.graph.variables.get(id).map_or(0.0, |n| n.value[0]))
        .collect();

    for (i, &amb_id) in amb_ids.iter().enumerate() {
        solver.graph.set_value(amb_id, &[integers[i]]);
    }
    inject_fixed_priors(&mut solver.graph, amb_ids, integers);

    let _ = solver.solve();

    let vals_fixed = VariableValues::build(&solver.graph.variables);
    let pose_fixed = vals_fixed
        .get(pose_id)
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .unwrap_or(init_pos);

    let is_valid = validate_fix_geometry(
        &solver.graph,
        &[float_pos.x, float_pos.y, float_pos.z],
        &[pose_fixed.x, pose_fixed.y, pose_fixed.z],
        5.0,
    );

    if !is_valid {
        for (i, &amb_id) in amb_ids.iter().enumerate() {
            solver.graph.set_value(amb_id, &[float_amb_values[i]]);
        }
        let n_fixed = amb_ids.len();
        for _ in 0..n_fixed {
            solver.graph.factors.pop();
        }
        let _ = solver.solve();
    }
}

/// Checks post-fit residuals of carrier-phase measurements and flags cycle slips.
pub fn check_postfit_cycle_slips(
    solver: &mut SlidingWindowSolver,
    pos_ecef: Vector3<f64>,
    cp_records: &[CpMeasurementRecord],
    slip_counts: &mut HashMap<u16, u32>,
) {
    let vals = VariableValues::build(&solver.graph.variables);
    for rec in cp_records {
        let r_sat = (rec.sat_pos - pos_ecef).norm();
        let r_ref = (rec.ref_pos - pos_ecef).norm();
        let expected_dd_cp = (r_sat - r_ref) - rec.base_dd;
        if let Some(amb_val) = vals.get(rec.var_amb) {
            let predicted = expected_dd_cp + rec.lambda * amb_val[0];
            let residual = rec.dd_cp_m - predicted;
            if residual.abs() > 0.5 {
                *slip_counts.entry(rec.sat).or_insert(0) += 1;
                solver.graph.remove_variable(rec.var_amb);
            }
        }
    }
}
