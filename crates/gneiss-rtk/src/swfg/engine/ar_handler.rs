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

thread_local! {
    static PPP_MW_TRACKER: std::cell::RefCell<Option<crate::ambiguity::ppp_ar::PppMwTracker>> = const { std::cell::RefCell::new(None) };
}

pub fn record_mw_sample(constellation_id: u8, satellite: u16, mw_cycles: f64, epoch: u32, slip: bool) {
    PPP_MW_TRACKER.with(|c| {
        let mut guard = c.borrow_mut();
        let tracker = guard.get_or_insert_with(crate::ambiguity::ppp_ar::PppMwTracker::new);
        tracker.update((constellation_id, satellite), mw_cycles, epoch, slip);
    });
}

pub fn reset_ppp_tracker() {
    PPP_MW_TRACKER.with(|c| *c.borrow_mut() = Some(crate::ambiguity::ppp_ar::PppMwTracker::new()));
    super::epoch::reset_epoch_tracker();
}

/// Execute LAMBDA + PAR ambiguity resolution step on active factor graph.
pub fn execute_ar_step(
    solver: &mut SlidingWindowSolver,
    pose_id: VariableId,
    init_pos: Vector3<f64>,
    is_rtk: bool,
    epoch: u32,
) {
    if epoch == 0 {
        reset_ppp_tracker();
    }
    let all_amb_ids = collect_ambiguity_variables(&solver.graph);
    if all_amb_ids.len() < 4 {
        return;
    }

    if !is_rtk {
        let min_epoch = std::env::var("PPP_AR_MIN_EPOCH").ok().and_then(|v| v.parse().ok()).unwrap_or(1800);
        if epoch >= min_epoch {
            execute_ppp_ar_step(solver, pose_id, init_pos, &all_amb_ids);
        }
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
        apply_validated_fix(solver, pose_id, init_pos, &amb_ids, &integers, 5.0);
    }
}

#[derive(Debug, Clone)]
pub struct SdAmbiguityConstraintFactor {
    pub vars: [VariableId; 2],
    pub fixed_diff: f64,
    pub information: nalgebra::DMatrix<f64>,
}

impl SdAmbiguityConstraintFactor {
    pub fn new(cand_amb: VariableId, ref_amb: VariableId, fixed_diff: f64, weight: f64) -> Self {
        Self {
            vars: [cand_amb, ref_amb],
            fixed_diff,
            information: nalgebra::DMatrix::from_element(1, 1, weight),
        }
    }
}

impl crate::swfg::factor::Factor for SdAmbiguityConstraintFactor {
    fn variables(&self) -> &[VariableId] {
        &self.vars
    }

    fn residual(&self, values: &VariableValues) -> nalgebra::DVector<f64> {
        let v_cand = values.get(self.vars[0]).map_or(0.0, |v| v[0]);
        let v_ref = values.get(self.vars[1]).map_or(0.0, |v| v[0]);
        nalgebra::DVector::from_element(1, (v_cand - v_ref) - self.fixed_diff)
    }

    fn jacobian(&self, values: &VariableValues) -> nalgebra::DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = nalgebra::DMatrix::zeros(1, total_dim);
        if let Some((start, _)) = values.index_of(self.vars[0]) {
            j[(0, start)] = 1.0;
        }
        if let Some((start, _)) = values.index_of(self.vars[1]) {
            j[(0, start)] = -1.0;
        }
        j
    }

    fn information(&self) -> nalgebra::DMatrix<f64> {
        self.information.clone()
    }
}

fn apply_sd_constraints(
    graph: &mut crate::swfg::graph::EstimationGraph,
    amb_ids: &[VariableId],
    fixed_undiff: &nalgebra::DVector<f64>,
) -> usize {
    let n_cand = amb_ids.len() - 1;
    let ref_id = amb_ids[0];
    let ref_val = graph.variables.get(&ref_id).map_or(0.0, |n| n.value[0]);
    let mut added = 0;
    for k in 0..n_cand {
        let cand_id = amb_ids[k + 1];
        let already_has = graph.factors.iter().any(|f| {
            let v = f.variables();
            v.len() == 2 && v[0] == cand_id && v[1] == ref_id
        });
        if !already_has {
            let diff = fixed_undiff[k + 1] - fixed_undiff[0];
            graph.set_value(cand_id, &[ref_val + diff]);
            let factor = SdAmbiguityConstraintFactor::new(cand_id, ref_id, diff, 1e8);
            graph.add_factor(Box::new(factor));
            added += 1;
        }
    }
    added
}

fn revert_sd_constraints(
    solver: &mut SlidingWindowSolver,
    amb_ids: &[VariableId],
    float_values: &[f64],
    n_added: usize,
) {
    for (i, &amb_id) in amb_ids.iter().enumerate() {
        solver.graph.set_value(amb_id, &[float_values[i]]);
    }
    for _ in 0..n_added {
        solver.graph.factors.pop();
    }
    let _ = solver.solve();
}

fn apply_validated_ppp_fix(
    solver: &mut SlidingWindowSolver,
    pose_id: VariableId,
    init_pos: Vector3<f64>,
    amb_ids: &[VariableId],
    fixed_undiff: &nalgebra::DVector<f64>,
    max_jump_m: f64,
) -> bool {
    let vals_before = VariableValues::build(&solver.graph.variables);
    let float_pos = vals_before.get(pose_id).map(|p| Vector3::new(p[0], p[1], p[2])).unwrap_or(init_pos);
    let float_amb_values: Vec<f64> = amb_ids.iter()
        .map(|id| solver.graph.variables.get(id).map_or(0.0, |n| n.value[0]))
        .collect();

    let n_added = apply_sd_constraints(&mut solver.graph, amb_ids, fixed_undiff);
    if n_added == 0 { return true; }
    let _ = solver.solve();

    let vals_fixed = VariableValues::build(&solver.graph.variables);
    let pose_fixed = vals_fixed.get(pose_id).map(|p| Vector3::new(p[0], p[1], p[2])).unwrap_or(init_pos);

    if std::env::var("PPP_AR_DEBUG").is_ok() {
        println!("DEBUG PPP-AR FIX: float_pos={:?} pose_fixed={:?} jump={:.4}m max={:.1}m",
            float_pos, pose_fixed, (pose_fixed - float_pos).norm(), max_jump_m);
    }

    let is_valid = validate_fix_geometry(
        &solver.graph,
        &[float_pos.x, float_pos.y, float_pos.z],
        &[pose_fixed.x, pose_fixed.y, pose_fixed.z],
        max_jump_m,
    );

    if !is_valid {
        revert_sd_constraints(solver, amb_ids, &float_amb_values, n_added);
        false
    } else {
        if std::env::var("PPP_AR_DEBUG").is_ok() {
            println!("DEBUG PPP-AR FIX ACCEPTED: jump={:.4}m", (pose_fixed - float_pos).norm());
        }
        true
    }
}

fn build_par_subset(
    sub_ids: &[VariableId],
    float_vec: &[f64],
    cov_sub: &nalgebra::DMatrix<f64>,
    wl_integers: &[i32],
    drop_idx: usize,
) -> (Vec<VariableId>, Vec<f64>, nalgebra::DMatrix<f64>, Vec<i32>) {
    let mut par_ids = vec![sub_ids[0]];
    let mut par_float = vec![float_vec[0]];
    let mut par_wl = Vec::new();
    let mut kept_indices = vec![0];
    for (k, &wl) in wl_integers.iter().enumerate() {
        if k != drop_idx && k + 1 < sub_ids.len() && k + 1 < float_vec.len() {
            par_ids.push(sub_ids[k + 1]);
            par_float.push(float_vec[k + 1]);
            par_wl.push(wl);
            kept_indices.push(k + 1);
        }
    }
    let m = kept_indices.len();
    let mut par_cov = nalgebra::DMatrix::zeros(m, m);
    for (r_i, &r) in kept_indices.iter().enumerate() {
        for (c_i, &c) in kept_indices.iter().enumerate() {
            par_cov[(r_i, c_i)] = cov_sub[(r, c)];
        }
    }
    (par_ids, par_float, par_cov, par_wl)
}

fn solve_ppp_ar_with_par(
    float_vec: &[f64],
    cov_sub: &nalgebra::DMatrix<f64>,
    f1: f64,
    f2: f64,
    wl_integers: &[i32],
    sub_ids: &[VariableId],
) -> Option<(Vec<VariableId>, nalgebra::DVector<f64>)> {
    if let Some((fixed_undiff, ratio)) = crate::ambiguity::ppp_ar::PppArSolver::resolve_sd_with_fixed_wl(
        float_vec, cov_sub, f1, f2, wl_integers, 2.0,
    ) {
        if std::env::var("PPP_AR_DEBUG").is_ok() {
            println!("DEBUG PPP-AR: full set fixed with ratio={:.3}", ratio);
        }
        return Some((sub_ids.to_vec(), fixed_undiff));
    }
    let n_cand = wl_integers.len();
    if n_cand < 4 { return None; }
    for drop_idx in 0..n_cand {
        let (par_ids, par_float, par_cov, par_wl) = build_par_subset(
            sub_ids, float_vec, cov_sub, wl_integers, drop_idx,
        );
        if let Some((fixed_undiff, ratio)) = crate::ambiguity::ppp_ar::PppArSolver::resolve_sd_with_fixed_wl(
            &par_float, &par_cov, f1, f2, &par_wl, 2.0,
        ) {
            if std::env::var("PPP_AR_DEBUG").is_ok() {
                println!("DEBUG PPP-AR: PAR subset drop {} fixed with ratio={:.3}", drop_idx, ratio);
            }
            return Some((par_ids, fixed_undiff));
        }
    }
    None
}

fn collect_sub_wl(
    ref_item: (usize, VariableId, u16),
    items: &[(usize, VariableId, u16)],
    fixed_wl: &[((u8, u16), i32)],
) -> (Vec<VariableId>, Vec<usize>, Vec<i32>) {
    let mut sub_amb_ids = vec![ref_item.1];
    let mut sub_indices = vec![ref_item.0];
    let mut wl_integers = Vec::new();
    for &(cand_key, n_wl) in fixed_wl {
        if sub_amb_ids.len() >= 8 { break; }
        if let Some(item) = items.iter().find(|s| s.2 == cand_key.1) {
            sub_amb_ids.push(item.1);
            sub_indices.push(item.0);
            wl_integers.push(n_wl);
        }
    }
    (sub_amb_ids, sub_indices, wl_integers)
}

struct AmbiguityState<'a> {
    ids: &'a [VariableId],
    float_vec: &'a nalgebra::DVector<f64>,
    cov: &'a nalgebra::DMatrix<f64>,
}

fn try_fix_constellation(
    solver: &mut SlidingWindowSolver,
    pose_id: VariableId,
    init_pos: Vector3<f64>,
    tracker: &crate::ambiguity::ppp_ar::PppMwTracker,
    constell: u8,
    items: &[(usize, VariableId, u16)],
    state: &AmbiguityState<'_>,
) {
    if items.len() < 4 { return; }
    let ref_opt = items.iter().find(|item| {
        tracker.get_smoothed((constell, item.2)).is_some_and(|(_, c)| c >= 5)
    });
    let ref_item = match ref_opt {
        Some(&it) => it,
        None => return,
    };
    let ref_sat = (constell, ref_item.2);
    let cand_keys: Vec<(u8, u16)> = items.iter().filter(|s| s.2 != ref_item.2).map(|s| (constell, s.2)).collect();
    let fixed_wl = tracker.fix_sd_wide_lane_subset(ref_sat, &cand_keys, 5, 0.40);
    if fixed_wl.len() < 3 { return; }

    let (sub_amb_ids, sub_indices, wl_integers) = collect_sub_wl(ref_item, items, &fixed_wl);
    if sub_amb_ids.len() < 4 { return; }

    let (sub_ids, float_vec, cov_sub) = extract_sub_ambiguity_system(
        state.ids, state.float_vec, state.cov, &sub_indices,
    );
    let (f1, f2) = sat_frequencies(constell);
    if let Some((fix_ids, fixed_undiff)) = solve_ppp_ar_with_par(
        &float_vec, &cov_sub, f1, f2, &wl_integers, &sub_ids,
    ) {
        apply_validated_ppp_fix(solver, pose_id, init_pos, &fix_ids, &fixed_undiff, 5.0);
    }
}

fn execute_ppp_ar_step(
    solver: &mut SlidingWindowSolver,
    pose_id: VariableId,
    init_pos: Vector3<f64>,
    all_amb_ids: &[VariableId],
) {
    let values = VariableValues::build(&solver.graph.variables);
    let (jtj, _, _) = solver.build_normal_equations(&values);
    let (float_amb_all, amb_cov_all) = extract_ambiguity_state(&solver.graph, &jtj, all_amb_ids);

    let state = AmbiguityState { ids: all_amb_ids, float_vec: &float_amb_all, cov: &amb_cov_all };
    let groups = group_ambiguities_for_ppp(solver, all_amb_ids, &amb_cov_all);
    PPP_MW_TRACKER.with(|c| {
        let guard = c.borrow();
        if let Some(tracker) = guard.as_ref() {
            for (constell, items) in &groups {
                if *constell != 0 {
                    continue;
                }
                try_fix_constellation(solver, pose_id, init_pos, tracker, *constell, items, &state);
            }
        }
    });
}

fn group_ambiguities_for_ppp(
    solver: &SlidingWindowSolver,
    all_amb_ids: &[VariableId],
    amb_cov_all: &nalgebra::DMatrix<f64>,
) -> HashMap<u8, Vec<(usize, VariableId, u16)>> {
    let mut by_group: HashMap<u8, Vec<(usize, VariableId, u16)>> = HashMap::new();
    for (i, &amb_id) in all_amb_ids.iter().enumerate() {
        let Some(node) = solver.graph.variables.get(&amb_id) else { continue };
        let crate::swfg::variables::VariableKind::Ambiguity { constellation_id, satellite, .. } = node.kind else { continue };
        let var = amb_cov_all[(i, i)];
        if var > 1e-6 && var < 10.0 {
            by_group.entry(constellation_id).or_default().push((i, amb_id, satellite));
        }
    }
    for list in by_group.values_mut() {
        list.sort_by(|a, b| amb_cov_all[(a.0, a.0)].total_cmp(&amb_cov_all[(b.0, b.0)]));
    }
    by_group
}

fn extract_sub_ambiguity_system(
    all_amb_ids: &[VariableId],
    float_all: &nalgebra::DVector<f64>,
    cov_all: &nalgebra::DMatrix<f64>,
    indices: &[usize],
) -> (Vec<VariableId>, Vec<f64>, nalgebra::DMatrix<f64>) {
    let m = indices.len();
    let sub_ids: Vec<VariableId> = indices.iter().map(|&i| all_amb_ids[i]).collect();
    let float_vec: Vec<f64> = indices.iter().map(|&i| float_all[i]).collect();
    let mut cov_sub = nalgebra::DMatrix::zeros(m, m);
    for (r_i, &r) in indices.iter().enumerate() {
        for (c_i, &c) in indices.iter().enumerate() {
            cov_sub[(r_i, c_i)] = cov_all[(r, c)];
        }
    }
    (sub_ids, float_vec, cov_sub)
}

fn sat_frequencies(constell: u8) -> (f64, f64) {
    match constell {
        2 => (1575.42e6, 1207.14e6),
        3 => (1561.098e6, 1207.14e6),
        _ => (1575.42e6, 1227.60e6),
    }
}

/// Applies fixed priors and validates the resulting position jump.
fn apply_validated_fix(
    solver: &mut SlidingWindowSolver,
    pose_id: VariableId,
    init_pos: Vector3<f64>,
    amb_ids: &[VariableId],
    integers: &nalgebra::DVector<f64>,
    max_jump_m: f64,
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

    if std::env::var("PPP_AR_DEBUG").is_ok() {
        println!("DEBUG PPP-AR FIX: float_pos={:?} pose_fixed={:?} jump={:.4}m",
            float_pos, pose_fixed, (pose_fixed - float_pos).norm());
    }

    let is_valid = validate_fix_geometry(
        &solver.graph,
        &[float_pos.x, float_pos.y, float_pos.z],
        &[pose_fixed.x, pose_fixed.y, pose_fixed.z],
        max_jump_m,
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
    slip_counts: &mut HashMap<(u8, u16), u32>,
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
                *slip_counts.entry((rec.constellation_id, rec.sat)).or_insert(0) += 1;
                solver.graph.remove_variable(rec.var_amb);
            }
        }
    }
}
