//! Schur complement marginalization for the sliding-window factor graph.
//!
//! When the window is full and a new epoch arrives, the oldest epoch's
//! per-epoch variables (Pose, Velocity, ClockBias, TropoZwd) are marginalized
//! out.  Variables that persist across epochs (ambiguities, IMU bias, IFB)
//! are kept in the active set.  The marginalized information is condensed
//! into a dense `MarginalPriorFactor` on the remaining variables.
//!
//! Mathematical detail:
//!   Partition the full Hessian H and gradient g into:
//!     H = [A  B; B^T  C]    g = [g_a; g_b]
//!   where A = marginalized-marginalized block, C = remaining-remaining.
//!   Schur complement of A:  S = C - B^T A^{-1} B
//!   Reduced gradient:        g_rem = g_b - B^T A^{-1} g_a
//!   The prior factor on remaining variables encodes (S, g_rem).

use nalgebra::{DMatrix, DVector};

use crate::swfg::graph::{EstimationGraph, MarginalPriorFactor};
use crate::swfg::variables::{VariableId, VariableValues};

/// Identify which variables should be marginalized (per-epoch variables
/// from the specified epoch) vs kept (cross-epoch variables).
pub fn partition_variables(
    graph: &EstimationGraph,
    epoch_to_marginalize: u32,
) -> (Vec<VariableId>, Vec<VariableId>) {
    let mut explicit_marg = std::collections::HashSet::new();

    for (id, node) in graph.variables.iter() {
        let should_marginalize = node.kind.is_per_epoch()
            && node.kind.epoch() == Some(epoch_to_marginalize);
        if should_marginalize {
            explicit_marg.insert(*id);
        }
    }

    let mut unconsumed_factors = Vec::new();
    for (i, factor) in graph.factors.iter().enumerate() {
        if !factor.variables().iter().any(|v| explicit_marg.contains(v)) {
            unconsumed_factors.push(i);
        }
    }

    let mut marginalized = Vec::new();
    let mut kept = Vec::new();

    for id in graph.variables.keys() {
        if explicit_marg.contains(id) {
            marginalized.push(*id);
            continue;
        }

        let is_in_unconsumed_factor = unconsumed_factors.iter().any(|&f_idx| {
            let vars = graph.factors[f_idx].variables();
            vars.contains(id) && vars.len() > 1
        });

        if !is_in_unconsumed_factor {
            // This variable is no longer observed by any active factor.
            // Since the old marginal prior is fully absorbed into the new one during marginalization,
            // we can safely marginalize out this dead variable to keep the state vector small.
            marginalized.push(*id);
        } else {
            kept.push(*id);
        }
    }

    (marginalized, kept)
}

/// `kept_ids` are the variables to keep.
///
/// Returns a `MarginalPriorFactor` that encodes all information from the
/// marginalized variables about the kept variables, or `None` if the
/// marginalized block is singular.
pub fn marginalize(
    graph: &EstimationGraph,
    hessian: &DMatrix<f64>,
    gradient: &DVector<f64>,
    marginalized_ids: &[VariableId],
    kept_ids: &[VariableId],
) -> Option<MarginalPriorFactor> {
    if marginalized_ids.is_empty() || kept_ids.is_empty() {
        return None;
    }

    let values = VariableValues::build(&graph.variables);

    // Build index sets: which rows/cols of the full Hessian belong to
    // marginalized vs kept variables.
    let marg_indices: Vec<usize> = marginalized_ids
        .iter()
        .filter_map(|id| values.index_of(*id))
        .flat_map(|(start, dim)| (start..start + dim).collect::<Vec<_>>())
        .collect();

    let kept_indices: Vec<usize> = kept_ids
        .iter()
        .filter_map(|id| values.index_of(*id))
        .flat_map(|(start, dim)| (start..start + dim).collect::<Vec<_>>())
        .collect();

    if marg_indices.is_empty() || kept_indices.is_empty() {
        return None;
    }

    // Extract blocks A (marginalized-marginalized), B (marginalized-kept), C (kept-kept)
    let n_marg = marg_indices.len();
    let n_kept = kept_indices.len();
    let mut a = DMatrix::zeros(n_marg, n_marg);
    let mut b = DMatrix::zeros(n_marg, n_kept);
    let mut c = DMatrix::zeros(n_kept, n_kept);

    for (i, &mi) in marg_indices.iter().enumerate() {
        for (j, &mj) in marg_indices.iter().enumerate() {
            a[(i, j)] = hessian[(mi, mj)];
        }
        for (j, &kj) in kept_indices.iter().enumerate() {
            b[(i, j)] = hessian[(mi, kj)];
        }
    }
    for (i, &ki) in kept_indices.iter().enumerate() {
        for (j, &kj) in kept_indices.iter().enumerate() {
            c[(i, j)] = hessian[(ki, kj)];
        }
    }

    // Extract gradient blocks and x0 (linearization point)
    let mut g_a = DVector::zeros(n_marg);
    let mut g_b = DVector::zeros(n_kept);
    let mut x0 = DVector::zeros(n_kept);
    for (i, &mi) in marg_indices.iter().enumerate() {
        g_a[i] = gradient[mi];
    }
    
    // Kept indices can span multiple variables. We just pull the current values 
    // from the graph and pack them into x0.
    // Wait, the variables in `kept_ids` map to `kept_indices`.
    // Let's do it cleanly by iterating over kept_ids.
    let mut offset = 0;
    for &id in kept_ids {
        if let Some(node) = graph.variables.get(&id) {
            let dim = node.value.len();
            for i in 0..dim {
                x0[offset + i] = node.value[i];
            }
            offset += dim;
        }
    }
    
    for (i, &ki) in kept_indices.iter().enumerate() {
        g_b[i] = gradient[ki];
    }

    // Schur complement: S = C - B^T A^{-1} B
    // Use SVD pseudo-inverse for numerical stability
    let svd = a.clone().svd(true, true);
    let mut a_inv = DMatrix::zeros(a.nrows(), a.ncols());
    for (i, &sigma) in svd.singular_values.iter().enumerate() {
        if sigma > 1e-10 {
            let u_col = svd.u.as_ref()?.column(i);
            let v_col = svd.v_t.as_ref()?.row(i).transpose();
            a_inv += (v_col * u_col.transpose()) / sigma;
        }
    }

    let s = &c - &b.transpose() * &a_inv * &b;
    let g_rem = &g_b - &b.transpose() * &a_inv * &g_a;

    let mut prior_variables = Vec::new();
    for &id in kept_ids {
        if let Some(var) = graph.variables.get(&id) {
            prior_variables.push((id, var.value.len()));
        }
    }

    Some(MarginalPriorFactor {
        variables: prior_variables,
        hessian: s,
        gradient: g_rem,
        x0,
    })
}

/// Add the marginal prior to the graph and remove marginalized variables.
pub fn apply_marginalization(
    graph: &mut EstimationGraph,
    prior: MarginalPriorFactor,
    marginalized_ids: &[VariableId],
) {
    graph.marginal_prior = Some(prior);
    for &id in marginalized_ids {
        graph.remove_variable(id);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;
    use crate::swfg::graph::EstimationGraph;
    use crate::swfg::variables::{VariableId, VariableKind, VariableValues};

    #[test]
    fn test_partition_variables() {
        let mut graph = EstimationGraph::new();
        
        let pose0 = graph.add_variable(VariableKind::Pose { epoch: 0 });
        let vel0 = graph.add_variable(VariableKind::Velocity { epoch: 0 });
        let pose1 = graph.add_variable(VariableKind::Pose { epoch: 1 });
        let amb = graph.add_variable(VariableKind::DdAmbiguity { constellation_id: 1, satellite: 1, ref_satellite: 0, frequency: 1, arc: 0 });
        
        // Add a factor connecting pose1 and amb to keep them in unconsumed factors
        let factor = crate::swfg::pipeline::DdCarrierPhaseFactor {
            var_pose: pose1,
            var_amb: amb,
            dd_cp_obs_m: 10.0,
            sat_pos: Vector3::new(100.0, 0.0, 0.0),
            ref_pos: Vector3::new(0.0, 100.0, 0.0),
            base_pos: Vector3::zeros(),
            base_dd_range: 100.0,
            lambda: 0.19,
            variance_m2: 1e-4,
            elevation_rad: 1.0,
            ref_elevation_rad: 1.0,
            is_new_amb: false,
            variables: vec![pose1, amb],
        };
        graph.add_factor(Box::new(factor));

        let (marg_ids, kept_ids) = partition_variables(&graph, 0);
        
        assert_eq!(marg_ids.len(), 2);
        assert!(marg_ids.contains(&pose0));
        assert!(marg_ids.contains(&vel0));
        
        assert_eq!(kept_ids.len(), 2);
        assert!(kept_ids.contains(&pose1));
        assert!(kept_ids.contains(&amb));
    }

    #[test]
    fn test_marginalize_completely_disconnected() {
        let mut graph = EstimationGraph::new();
        let _p0 = graph.add_variable(VariableKind::Pose { epoch: 0 });
        let p1 = graph.add_variable(VariableKind::Pose { epoch: 1 });
        let amb = graph.add_variable(VariableKind::DdAmbiguity { constellation_id: 1, satellite: 1, ref_satellite: 0, frequency: 1, arc: 0 });

        let factor = crate::swfg::pipeline::DdCarrierPhaseFactor {
            var_pose: p1,
            var_amb: amb,
            dd_cp_obs_m: 10.0,
            sat_pos: Vector3::new(100.0, 0.0, 0.0),
            ref_pos: Vector3::new(0.0, 100.0, 0.0),
            base_pos: Vector3::zeros(),
            base_dd_range: 100.0,
            lambda: 0.19,
            variance_m2: 1e-4,
            elevation_rad: 1.0,
            ref_elevation_rad: 1.0,
            is_new_amb: false,
            variables: vec![p1, amb],
        };
        graph.add_factor(Box::new(factor));

        let (marg_ids, kept_ids) = partition_variables(&graph, 0);
        let values = VariableValues::build(&graph.variables);
        let total_dim = values.total_dim();
        let hessian = DMatrix::identity(total_dim, total_dim);
        let gradient = DVector::zeros(total_dim);

        let prior = marginalize(&graph, &hessian, &gradient, &marg_ids, &kept_ids);
        assert!(prior.is_some());
        let p = prior.unwrap();
        assert_eq!(p.variables.len(), kept_ids.len());
    }

    #[test]
    fn marginalize_empty_marg_returns_none() {
        let graph = EstimationGraph::new();
        let hessian = DMatrix::identity(1, 1);
        let gradient = DVector::zeros(1);
        let result = marginalize(&graph, &hessian, &gradient, &[], &[VariableId::new(0)]);
        assert!(result.is_none());
    }

    #[test]
    fn apply_marginalization_removes_variables() {
        let mut graph = EstimationGraph::new();
        let _p0 = graph.add_variable(VariableKind::Pose { epoch: 0 });
        let p1 = graph.add_variable(VariableKind::Pose { epoch: 1 });
        let (marg_ids, _kept_ids) = partition_variables(&graph, 0);

        let prior = MarginalPriorFactor {
            variables: vec![(p1, 6)],
            x0: DVector::zeros(6),
            hessian: DMatrix::identity(6, 6),
            gradient: DVector::zeros(6),
        };
        let n_vars_before = graph.n_variables();
        apply_marginalization(&mut graph, prior, &marg_ids);
        assert!(graph.n_variables() < n_vars_before);
        assert!(graph.marginal_prior.is_some());
    }
}
