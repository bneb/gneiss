//! Normal equations accumulation and linear system solution.

use nalgebra::{DMatrix, DVector};
use crate::swfg::factor::Factor;
use crate::swfg::variables::VariableValues;

pub(crate) fn find_active_columns(
    factor: &dyn Factor,
    values: &VariableValues,
    m_rows: usize,
    j: &DMatrix<f64>,
) -> Vec<usize> {
    let mut cols = Vec::new();
    for &id in factor.variables() {
        if let Some((start, dim)) = values.index_of(id) {
            for d in 0..dim {
                cols.push(start + d);
            }
        }
    }
    if cols.is_empty() {
        for c in 0..values.total_dim() {
            if (0..m_rows).any(|m| j[(m, c)] != 0.0) {
                cols.push(c);
            }
        }
    }
    cols
}

pub(crate) fn accumulate_factor_normal_equations(
    factor: &dyn Factor,
    values: &VariableValues,
    jtj: &mut DMatrix<f64>,
    jtr: &mut DVector<f64>,
) -> f64 {
    let r = factor.residual(values);
    let j = factor.jacobian(values);
    let w = factor.information();
    let r_norm = r.norm();

    let weight = factor
        .robust_threshold()
        .map_or(1.0, |k| crate::swfg::factor::compute_robust_weight(r_norm, k, factor.use_cauchy()));
    let w_scaled = &w * weight;
    let m_rows = r.len();
    let cols = find_active_columns(factor, values, m_rows, &j);

    let mut w_j_act = Vec::with_capacity(cols.len());
    for &c in &cols {
        let mut col = vec![0.0; m_rows];
        for m in 0..m_rows {
            col[m] = (0..m_rows).map(|k| w_scaled[(m, k)] * j[(k, c)]).sum();
        }
        w_j_act.push(col);
    }

    for (i, &ci) in cols.iter().enumerate() {
        jtr[ci] += (0..m_rows).map(|m| w_j_act[i][m] * r[m]).sum::<f64>();
        for (j_idx, &cj) in cols.iter().enumerate() {
            jtj[(ci, cj)] += (0..m_rows).map(|m| j[(m, ci)] * w_j_act[j_idx][m]).sum::<f64>();
        }
    }

    let raw_s = (&r.transpose() * &w * &r)[(0, 0)];
    factor.robust_threshold().map_or(raw_s, |k| {
        if r_norm < 1e-12 {
            raw_s
        } else {
            raw_s * crate::swfg::factor::compute_robust_error(r_norm, k, factor.use_cauchy())
                / (r_norm * r_norm)
        }
    })
}

/// Sparse per-factor contribution: indices + values for JtJ/Jtr.
struct FactorContribution {
    entries: Vec<(usize, usize, f64)>, // (row, col, value) for JtJ
    grad: Vec<(usize, f64)>,           // (idx, value) for Jtr
    error: f64,
}

fn compute_factor_contribution(
    factor: &dyn Factor,
    values: &VariableValues,
) -> FactorContribution {
    let r = factor.residual(values);
    let j = factor.jacobian(values);
    let w = factor.information();
    let r_norm = r.norm();

    let weight = factor.robust_threshold().map_or(1.0, |k| {
        crate::swfg::factor::compute_robust_weight(r_norm, k, factor.use_cauchy())
    });
    let w_scaled = &w * weight;
    let m_rows = r.len();
    let cols = find_active_columns(factor, values, m_rows, &j);

    let mut w_j_act = Vec::with_capacity(cols.len());
    for &c in &cols {
        let mut col = vec![0.0; m_rows];
        for m in 0..m_rows {
            col[m] = (0..m_rows).map(|k| w_scaled[(m, k)] * j[(k, c)]).sum();
        }
        w_j_act.push(col);
    }

    let mut entries = Vec::with_capacity(cols.len() * cols.len());
    let mut grad = Vec::with_capacity(cols.len());
    for (i, &ci) in cols.iter().enumerate() {
        let g: f64 = (0..m_rows).map(|m| w_j_act[i][m] * r[m]).sum();
        grad.push((ci, g));
        for (j_idx, &cj) in cols.iter().enumerate() {
            let v: f64 = (0..m_rows).map(|m| j[(m, ci)] * w_j_act[j_idx][m]).sum();
            entries.push((ci, cj, v));
        }
    }

    let raw_s = (&r.transpose() * &w * &r)[(0, 0)];
    let error = factor.robust_threshold().map_or(raw_s, |k| {
        if r_norm < 1e-12 { raw_s }
        else { raw_s * crate::swfg::factor::compute_robust_error(r_norm, k, factor.use_cauchy()) / (r_norm * r_norm) }
    });

    FactorContribution { entries, grad, error }
}

/// Minimum factor count to justify rayon parallel overhead.
const PAR_THRESHOLD: usize = 8;

/// Accumulate normal equations from factors, using rayon when beneficial.
pub(crate) fn accumulate_factors_parallel(
    factors: &[Box<dyn Factor>],
    values: &VariableValues,
    jtj: &mut DMatrix<f64>,
    jtr: &mut DVector<f64>,
) -> f64 {
    if factors.len() < PAR_THRESHOLD {
        let mut total = 0.0;
        for f in factors {
            total += accumulate_factor_normal_equations(f.as_ref(), values, jtj, jtr);
        }
        return total;
    }

    use rayon::prelude::*;
    let contribs: Vec<FactorContribution> = factors
        .par_iter()
        .map(|f| compute_factor_contribution(f.as_ref(), values))
        .collect();

    let mut total_error = 0.0_f64;
    for c in contribs {
        total_error += c.error;
        for (idx, val) in c.grad {
            jtr[idx] += val;
        }
        for (r, col, val) in c.entries {
            jtj[(r, col)] += val;
        }
    }
    total_error
}


pub(crate) fn solve_linear_system(a: &DMatrix<f64>, b: &DVector<f64>) -> Option<DVector<f64>> {
    if let Some(chol) = a.clone().cholesky() {
        return Some(chol.solve(b));
    }
    if let Some(x) = a.clone().qr().solve(b) {
        return Some(x);
    }
    let lu = a.clone().lu();
    if let Some(x) = lu.solve(b) {
        return Some(x);
    }
    let svd = a.clone().svd(true, true);
    let mut x = DVector::zeros(b.len());
    for (i, &sigma) in svd.singular_values.iter().enumerate() {
        if sigma > 1e-12 {
            let u_col = svd.u.as_ref()?.column(i);
            let v_col = svd.v_t.as_ref()?.row(i).transpose();
            x += u_col.dot(b) / sigma * v_col;
        }
    }
    Some(x)
}
