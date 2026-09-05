//! Integer ambiguity resolution integrated into the sliding-window factor graph.
//!
//! Two-step process:
//!   1. Float solve — run LM to convergence with continuous ambiguities
//!   2. Extract the marginal covariance of the ambiguity variables
//!   3. Call LAMBDA to find the best integer candidate
//!   4. If the ratio test passes, inject FixedAmbiguityPriorFactors
//!      (near-infinite information → hard constraint)
//!   5. Partial LM re-optimize to snap Pose/Velocity to the fixed integers

use nalgebra::{DMatrix, DVector};

use crate::swfg::factor::PriorFactor;
use crate::swfg::graph::EstimationGraph;
use crate::swfg::variables::{VariableId, VariableKind, VariableValues};

/// Result of attempting to fix integer ambiguities.
#[derive(Debug, Clone)]
pub enum ArResult {
    /// Fix succeeded.  Contains the fixed integer values and the ratio.
    Fixed {
        /// Fixed integer ambiguities (cycles), indexed parallel to `amb_vars`.
        integers: DVector<f64>,
        /// LAMBDA ratio test value.
        ratio: f64,
        /// Number of ambiguities that were fixed.
        n_fixed: usize,
    },
    /// Fix failed — keep float solution.
    Float {
        reason: String,
        best_ratio: f64,
    },
}

/// Extract the float ambiguity values and their marginal covariance from
/// the graph after a full LM solve.
///
/// The marginal covariance is the block of the full Hessian inverse
/// corresponding to the ambiguity variables.
pub fn extract_ambiguity_state(
    graph: &EstimationGraph,
    hessian: &DMatrix<f64>, // J^T W J from the converged float solution
    amb_var_ids: &[VariableId],
) -> (DVector<f64>, DMatrix<f64>) {
    // Build mapping: VariableId → (start, dim) in the full state
    let values = VariableValues::build(&graph.variables);

    // Collect indices of ambiguity variables
    let amb_indices: Vec<(usize, usize)> = amb_var_ids
        .iter()
        .filter_map(|id| values.index_of(*id))
        .collect();

    let n_amb = amb_indices.len();
    let mut float_amb = DVector::zeros(n_amb);
    let mut amb_cov = DMatrix::zeros(n_amb, n_amb);

    if n_amb == 0 {
        return (float_amb, amb_cov);
    }

    // Extract ambiguity values
    for (i, (start, _dim)) in amb_indices.iter().enumerate() {
        float_amb[i] = values.state()[*start];
    }

    // Extract marginal covariance using QR solve on standard basis vectors for each ambiguity variable.
    // Damping (1e-6 * I) guards against rank-deficiency in unconstrained gauge parameters.
    let damped = hessian + DMatrix::identity(hessian.nrows(), hessian.ncols()) * 1e-6;
    let qr = damped.qr();
    for (i, (si, _)) in amb_indices.iter().enumerate() {
        let mut e = DVector::zeros(hessian.nrows());
        e[*si] = 1.0;
        if let Some(col) = qr.solve(&e) {
            for (j, (sj, _)) in amb_indices.iter().enumerate() {
                amb_cov[(i, j)] = col[*sj];
            }
        }
    }

    (float_amb, amb_cov)
}

/// Attempt to fix integer ambiguities using LAMBDA.
///
/// Returns `ArResult::Fixed` with the integer solution if the ratio test
/// passes, or `ArResult::Float` with the reason for failure.
pub fn attempt_ar_fix(
    float_amb: &DVector<f64>,
    amb_cov: &DMatrix<f64>,
    min_ratio: f64,
) -> ArResult {
    let n = float_amb.len();
    if n < 3 {
        return ArResult::Float {
            reason: format!("too few ambiguities: {} < 3", n),
            best_ratio: 0.0,
        };
    }

    // Call LAMBDA — uses the existing implementation from the ambiguity crate.
    let result = match crate::ambiguity::lambda::resolve_lambda(float_amb, amb_cov) {
        Ok(res) => res,
        Err(e) => {
            return ArResult::Float {
                reason: format!("LAMBDA failed: {}", e),
                best_ratio: 0.0,
            };
        }
    };

    if result.ratio >= min_ratio {
        ArResult::Fixed {
            integers: result.best_integers,
            ratio: result.ratio,
            n_fixed: n,
        }
    } else {
        ArResult::Float {
            reason: format!("ratio {:.2} < {:.2}", result.ratio, min_ratio),
            best_ratio: result.ratio,
        }
    }
}

/// Partial Ambiguity Resolution (PAR): If full AR ratio test fails,
/// attempt to fix a subset of the highest-quality float ambiguities.
pub fn attempt_partial_ar_fix(
    float_amb: &DVector<f64>,
    amb_cov: &DMatrix<f64>,
    min_ratio: f64,
) -> (ArResult, Vec<usize>) {
    let full_res = attempt_ar_fix(float_amb, amb_cov, min_ratio);
    let n = float_amb.len();
    if matches!(full_res, ArResult::Fixed { .. }) || n <= 4 {
        let indices: Vec<usize> = (0..n).collect();
        return (full_res, indices);
    }

    // Sort ambiguity indices by variance (worst/highest variance first)
    let mut indices_by_var: Vec<(usize, f64)> = (0..n)
        .map(|i| (i, amb_cov[(i, i)]))
        .collect();
    indices_by_var.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut current_subset: Vec<usize> = (0..n).collect();
    for (drop_idx, _) in indices_by_var {
        if current_subset.len() <= 4 {
            break;
        }
        current_subset.retain(|&idx| idx != drop_idx);

        let k = current_subset.len();
        let mut sub_float = DVector::zeros(k);
        let mut sub_cov = DMatrix::zeros(k, k);
        for (r_out, &r_in) in current_subset.iter().enumerate() {
            sub_float[r_out] = float_amb[r_in];
            for (c_out, &c_in) in current_subset.iter().enumerate() {
                sub_cov[(r_out, c_out)] = amb_cov[(r_in, c_in)];
            }
        }

        let sub_res = attempt_ar_fix(&sub_float, &sub_cov, min_ratio + 0.1);
        if let ArResult::Fixed { integers, ratio, n_fixed: _ } = sub_res {
            return (
                ArResult::Fixed {
                    integers,
                    ratio,
                    n_fixed: k,
                },
                current_subset,
            );
        }
    }

    (full_res, (0..n).collect())
}

/// Inject fixed-ambiguity prior factors into the graph.
///
/// Each factor is a strong quadratic prior centered at the fixed integer
/// value.  The information is set high enough (1e8) that the optimizer
/// treats it as a hard constraint.
pub fn inject_fixed_priors(
    graph: &mut EstimationGraph,
    amb_var_ids: &[VariableId],
    fixed_integers: &DVector<f64>,
) {
    for (i, &amb_id) in amb_var_ids.iter().enumerate() {
        let factor = PriorFactor {
            variable: amb_id,
            mu: DVector::from_element(1, fixed_integers[i]),
            information: DMatrix::from_element(1, 1, 1e8), // near-infinite → hard constraint
        };
        graph.add_factor(Box::new(factor));
    }
}

/// Validate that the fixed solution is consistent with the pseudorange
/// evidence.  If the fixed position differs from the float position by
/// more than `max_jump_m`, reject the fix.
pub fn validate_fix_geometry(
    _graph: &EstimationGraph,
    _float_position: &[f64],
    _fixed_position: &[f64],
    max_jump_m: f64,
) -> bool {
    let dx = _float_position[0] - _fixed_position[0];
    let dy = _float_position[1] - _fixed_position[1];
    let dz = _float_position[2] - _fixed_position[2];
    let jump = (dx * dx + dy * dy + dz * dz).sqrt();
    jump <= max_jump_m
}

/// Filter active ambiguity variables from the graph's variable set.
pub fn collect_ambiguity_variables(graph: &EstimationGraph) -> Vec<VariableId> {
    let has_factors = !graph.factors.is_empty();
    let active_vars: std::collections::HashSet<VariableId> = if has_factors {
        graph
            .factors
            .iter()
            .flat_map(|f| f.variables())
            .copied()
            .collect()
    } else {
        std::collections::HashSet::new()
    };

    graph
        .variables
        .iter()
        .filter(|(id, node)| {
            (!has_factors || active_vars.contains(id))
                && match node.kind {
                    VariableKind::Ambiguity { .. } => true,
                    VariableKind::DdAmbiguity { constellation_id, .. } => constellation_id != 1,
                    _ => false,
                }
        })
        .map(|(id, _)| *id)
        .collect()
}

// ---- Tests ----------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_ambiguity_state_from_graph() {
        let mut graph = EstimationGraph::new();
        let amb1 = graph.add_variable(VariableKind::Ambiguity { constellation_id: 0, satellite: 1, frequency: 1, arc: 0 });
        let amb2 = graph.add_variable(VariableKind::Ambiguity { constellation_id: 0, satellite: 2, frequency: 1, arc: 0 });
        graph.add_variable(VariableKind::Pose { epoch: 0 });

        // Set ambiguity values
        graph.set_value(amb1, &[1.5]);
        graph.set_value(amb2, &[-0.3]);

        let amb_ids = collect_ambiguity_variables(&graph);
        assert_eq!(amb_ids.len(), 2);

        // Build a simple Hessian: identity (each variable independent, var=1)
        let values = VariableValues::build(&graph.variables);
        let hessian = DMatrix::identity(values.total_dim(), values.total_dim());

        let (float_amb, amb_cov) = extract_ambiguity_state(&graph, &hessian, &amb_ids);
        assert!((float_amb[0] - 1.5).abs() < 1e-12);
        assert!((float_amb[1] + 0.3).abs() < 1e-12);
        // With damped identity Hessian, marginal covariance is approximately identity
        assert!((amb_cov[(0, 0)] - 1.0).abs() < 1e-5);
    }

    #[test]
    fn ar_too_few_ambiguities_returns_float() {
        let float_amb = DVector::from_vec(vec![1.0, 2.0]);
        let amb_cov = DMatrix::identity(2, 2);
        let result = attempt_ar_fix(&float_amb, &amb_cov, 2.5);
        assert!(matches!(result, ArResult::Float { .. }));
    }

    #[test]
    fn inject_fixed_priors_adds_factors() {
        let mut graph = EstimationGraph::new();
        let amb1 = graph.add_variable(VariableKind::Ambiguity { constellation_id: 0, satellite: 1, frequency: 1, arc: 0 });
        let amb2 = graph.add_variable(VariableKind::Ambiguity { constellation_id: 0, satellite: 2, frequency: 1, arc: 0 });
        let ids = vec![amb1, amb2];
        let fixed = DVector::from_vec(vec![1.0, -1.0]);

        assert_eq!(graph.n_factors(), 0);
        inject_fixed_priors(&mut graph, &ids, &fixed);
        assert_eq!(graph.n_factors(), 2);
    }

    #[test]
    fn validate_fix_geometry_rejects_large_jump() {
        assert!(validate_fix_geometry(
            &EstimationGraph::new(),
            &[0.0, 0.0, 0.0],
            &[0.0, 0.0, 0.0],
            2.0,
        ));
        assert!(!validate_fix_geometry(
            &EstimationGraph::new(),
            &[0.0, 0.0, 0.0],
            &[10.0, 0.0, 0.0],
            2.0,
        ));
    }
}
