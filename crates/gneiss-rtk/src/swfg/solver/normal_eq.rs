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

#[cfg(test)]
mod tests {
    //! The normal equations are the heart of the solver, so they are checked
    //! against the textbook definitions rather than against the solver's own
    //! output:
    //!   JtJ = J^T W J,   Jtr = J^T W r,   chi2 = r^T W r
    //! with every entry worked out by hand in the test that uses it.

    #![allow(clippy::unwrap_used)]

    use super::*;
    use crate::swfg::variables::{VariableId, VariableKind, VariableNode};
    use std::collections::BTreeMap;

    /// A factor with hand-picked residual, Jacobian and information matrix.
    #[derive(Debug)]
    struct Linear {
        vars: Vec<VariableId>,
        r: DVector<f64>,
        j: DMatrix<f64>,
        w: DMatrix<f64>,
        threshold: Option<f64>,
        cauchy: bool,
    }

    impl Linear {
        fn new(vars: &[VariableId], r: &[f64], j: &[&[f64]], w: &[f64]) -> Self {
            Self {
                vars: vars.to_vec(),
                r: DVector::from_row_slice(r),
                j: DMatrix::from_fn(j.len(), j[0].len(), |i, k| j[i][k]),
                w: DMatrix::from_diagonal(&DVector::from_row_slice(w)),
                threshold: None,
                cauchy: false,
            }
        }

        fn robust(mut self, k: f64, cauchy: bool) -> Self {
            self.threshold = Some(k);
            self.cauchy = cauchy;
            self
        }
    }

    impl Factor for Linear {
        fn variables(&self) -> &[VariableId] {
            &self.vars
        }
        fn residual(&self, _values: &VariableValues) -> DVector<f64> {
            self.r.clone()
        }
        fn jacobian(&self, _values: &VariableValues) -> DMatrix<f64> {
            self.j.clone()
        }
        fn information(&self) -> DMatrix<f64> {
            self.w.clone()
        }
        fn robust_threshold(&self) -> Option<f64> {
            self.threshold
        }
        fn use_cauchy(&self) -> bool {
            self.cauchy
        }
    }

    /// Two scalar variables with the given values, giving a 2-column state.
    fn two_scalars(a: f64, b: f64) -> (BTreeMap<VariableId, VariableNode>, VariableId, VariableId) {
        let mut vars = BTreeMap::new();
        for (i, v) in [a, b].iter().enumerate() {
            let id = VariableId::new(i as u64);
            let mut node = VariableNode::new(id, VariableKind::TropoZwd { epoch: 0 });
            node.value = DVector::from_element(1, *v);
            vars.insert(id, node);
        }
        (vars, VariableId::new(0), VariableId::new(1))
    }

    #[test]
    fn normal_equations_equal_j_transpose_w_j_and_j_transpose_w_r() {
        // J = [[2, 0], [0, 3]],  W = diag(1, 4),  r = (1, 2)
        //   W J   = [[2, 0], [0, 12]]
        //   Jt J  = [[4, 0], [0, 36]]
        //   W r   = (1, 8)
        //   Jt r  = (2, 24)
        //   r'Wr  = 1 + 2*8 = 17
        let (vars, x, y) = two_scalars(0.0, 0.0);
        let values = VariableValues::build(&vars);
        let f = Linear::new(&[x, y], &[1.0, 2.0], &[&[2.0, 0.0], &[0.0, 3.0]], &[1.0, 4.0]);

        let mut jtj = DMatrix::zeros(2, 2);
        let mut jtr = DVector::zeros(2);
        let chi2 = accumulate_factor_normal_equations(&f, &values, &mut jtj, &mut jtr);

        let want_jtj = DMatrix::from_diagonal(&DVector::from_vec(vec![4.0, 36.0]));
        assert!((&jtj - &want_jtj).norm() < 1e-12, "JtJ = {jtj:?}, want {want_jtj:?}");
        assert!((jtr[0] - 2.0).abs() < 1e-12, "Jtr[0] = {}", jtr[0]);
        assert!((jtr[1] - 24.0).abs() < 1e-12, "Jtr[1] = {}", jtr[1]);
        assert!((chi2 - 17.0).abs() < 1e-12, "chi2 = {chi2}");
    }

    #[test]
    fn huber_weighting_matches_the_closed_form() {
        // One variable, r = 3, J = 2, W = 1, k = 1.
        //   |r| = 3 > k      -> robust weight w = k/|r| = 1/3
        //   JtJ = J' W J = 2 * (1/3) * 2      = 4/3
        //   Jtr = J' W r = 2 * (1/3) * 3      = 2
        //   rho = |r|^2 = 9 inside the quadratic region would be 9, but the
        //        Huber loss gives 2 k |r| - k^2 = 5, and the reported error is
        //        raw * rho / |r|^2 = 9 * 5 / 9 = 5.
        let (vars, x, _y) = two_scalars(0.0, 0.0);
        let values = VariableValues::build(&vars);
        let f = Linear::new(&[x], &[3.0], &[&[2.0]], &[1.0]).robust(1.0, false);

        let mut jtj = DMatrix::zeros(2, 2);
        let mut jtr = DVector::zeros(2);
        let chi2 = accumulate_factor_normal_equations(&f, &values, &mut jtj, &mut jtr);
        assert!((jtj[(0, 0)] - 4.0 / 3.0).abs() < 1e-12, "JtJ = {}", jtj[(0, 0)]);
        assert!((jtr[0] - 2.0).abs() < 1e-12, "Jtr = {}", jtr[0]);
        assert!((chi2 - 5.0).abs() < 1e-12, "robust chi2 = {chi2}");
    }

    #[test]
    fn quadratic_loss_is_untouched_below_the_threshold() {
        // r = 0.5 with k = 1: weight 1 and the error stays r'Wr = 0.25.
        let (vars, x, _y) = two_scalars(0.0, 0.0);
        let values = VariableValues::build(&vars);
        let f = Linear::new(&[x], &[0.5], &[&[2.0]], &[1.0]).robust(1.0, false);
        let mut jtj = DMatrix::zeros(2, 2);
        let mut jtr = DVector::zeros(2);
        let chi2 = accumulate_factor_normal_equations(&f, &values, &mut jtj, &mut jtr);
        assert!((jtj[(0, 0)] - 4.0).abs() < 1e-12, "JtJ = {}", jtj[(0, 0)]);
        assert!((jtr[0] - 1.0).abs() < 1e-12, "Jtr = {}", jtr[0]);
        assert!((chi2 - 0.25).abs() < 1e-12, "chi2 = {chi2}");
    }

    #[test]
    fn cauchy_weighting_matches_the_closed_form() {
        // r = 3, k = 1: weight = 1/(1 + 9) = 1/10, so JtJ = 4/10 = 0.4 and
        // Jtr = 2 * 0.1 * 3 = 0.6.  The Cauchy loss is
        //   k^2 ln(1 + s/k^2) = ln(10)
        // and the reported error is raw * loss / s = 9 * ln(10) / 9 = ln(10).
        let (vars, x, _y) = two_scalars(0.0, 0.0);
        let values = VariableValues::build(&vars);
        let f = Linear::new(&[x], &[3.0], &[&[2.0]], &[1.0]).robust(1.0, true);
        let mut jtj = DMatrix::zeros(2, 2);
        let mut jtr = DVector::zeros(2);
        let chi2 = accumulate_factor_normal_equations(&f, &values, &mut jtj, &mut jtr);
        assert!((jtj[(0, 0)] - 0.4).abs() < 1e-12, "JtJ = {}", jtj[(0, 0)]);
        assert!((jtr[0] - 0.6).abs() < 1e-12, "Jtr = {}", jtr[0]);
        assert!((chi2 - 10.0_f64.ln()).abs() < 1e-12, "Cauchy chi2 = {chi2}");
    }

    #[test]
    fn active_columns_fall_back_to_scanning_the_jacobian() {
        // A factor that names a variable which is not in the state must still
        // contribute through the Jacobian scan, otherwise its information would
        // be silently dropped.
        let (vars, _x, _y) = two_scalars(0.0, 0.0);
        let values = VariableValues::build(&vars);
        let ghost = VariableId::new(77);
        let j = DMatrix::from_fn(2, 2, |i, k| if i == 0 && k == 0 { 5.0 } else { 0.0 });
        let f = Linear { vars: vec![ghost], r: DVector::zeros(2), j: j.clone(), w: DMatrix::identity(2, 2), threshold: None, cauchy: false };

        let cols = find_active_columns(&f, &values, 2, &j);
        assert_eq!(cols, vec![0], "only column 0 has a non-zero Jacobian entry");

        let mut jtj = DMatrix::zeros(2, 2);
        let mut jtr = DVector::zeros(2);
        let chi2 = accumulate_factor_normal_equations(&f, &values, &mut jtj, &mut jtr);
        assert!((jtj[(0, 0)] - 25.0).abs() < 1e-12, "fallback must still accumulate 5^2");
        assert!(chi2.abs() < 1e-12);
    }

    #[test]
    fn parallel_accumulation_matches_the_serial_sum() {
        // Ten factors is above the parallel threshold of 8; both paths must
        // produce the same normal equations and the same total error.
        let (vars, x, y) = two_scalars(0.0, 0.0);
        let values = VariableValues::build(&vars);
        let factors: Vec<Box<dyn Factor>> = (0..10)
            .map(|i| {
                let a = 1.0 + i as f64;
                Box::new(Linear::new(&[x, y], &[a, 2.0 * a], &[&[a, 0.0], &[0.0, a]], &[1.0, 4.0]))
                    as Box<dyn Factor>
            })
            .collect();

        let mut jtj_par = DMatrix::zeros(2, 2);
        let mut jtr_par = DVector::zeros(2);
        let err_par = accumulate_factors_parallel(&factors, &values, &mut jtj_par, &mut jtr_par);

        let mut jtj_ser = DMatrix::zeros(2, 2);
        let mut jtr_ser = DVector::zeros(2);
        let mut err_ser = 0.0;
        for f in &factors {
            err_ser += accumulate_factor_normal_equations(f.as_ref(), &values, &mut jtj_ser, &mut jtr_ser);
        }

        let djtj = (&jtj_par - &jtj_ser).norm();
        assert!(djtj < 1e-12, "JtJ differs by {djtj}");
        let djtr = (&jtr_par - &jtr_ser).norm();
        assert!(djtr < 1e-12, "Jtr differs by {djtr}");
        assert!((err_par - err_ser).abs() < 1e-9, "error differs by {}", (err_par - err_ser).abs());

        // Hand values over a = 1..10, with J = diag(a, a), W = diag(1, 4),
        // r = (a, 2a):
        //   J'WJ = diag(a^2, 4 a^2)      -> sum = diag(385, 1540)
        //   J'Wr = J'(a, 8a) = (a^2, 8a^2) -> sum = (385, 3080)
        //   chi2 = sum (a^2 + 4 (2a)^2) = 17 * 385 = 6545
        assert!((jtj_par[(0, 0)] - 385.0).abs() < 1e-9, "JtJ[0,0] = {}", jtj_par[(0, 0)]);
        assert!((jtj_par[(1, 1)] - 1540.0).abs() < 1e-9, "JtJ[1,1] = {}", jtj_par[(1, 1)]);
        assert!((jtr_par[0] - 385.0).abs() < 1e-9, "Jtr[0] = {}", jtr_par[0]);
        assert!((jtr_par[1] - 3080.0).abs() < 1e-9, "Jtr[1] = {}", jtr_par[1]);
        assert!((err_par - 6545.0).abs() < 1e-6, "chi2 = {err_par}");
    }

    #[test]
    fn serial_and_parallel_paths_agree_below_the_parallel_threshold() {
        let (vars, x, y) = two_scalars(0.0, 0.0);
        let values = VariableValues::build(&vars);
        let factors: Vec<Box<dyn Factor>> = (0..PAR_THRESHOLD - 1)
            .map(|i| {
                let a = 1.0 + i as f64;
                Box::new(Linear::new(&[x, y], &[a], &[&[a, 1.0]], &[2.0])) as Box<dyn Factor>
            })
            .collect();
        let mut jtj = DMatrix::zeros(2, 2);
        let mut jtr = DVector::zeros(2);
        let err = accumulate_factors_parallel(&factors, &values, &mut jtj, &mut jtr);
        // J'WJ for factor a is 2 * [[a^2, a], [a, 1]] and J'Wr = [2 a^2, 2 a],
        // so with sum(a^2) = 140 and sum(a) = 28 over a = 1..7:
        //   JtJ = [[280, 56], [56, 14]],  Jtr = (280, 56),  chi2 = 2*140 = 280
        assert!((jtj[(0, 0)] - 280.0).abs() < 1e-9, "JtJ[0,0] = {}", jtj[(0, 0)]);
        assert!((jtj[(0, 1)] - 56.0).abs() < 1e-9, "JtJ[0,1] = {}", jtj[(0, 1)]);
        assert!((jtj[(1, 1)] - 14.0).abs() < 1e-9, "JtJ[1,1] = {}", jtj[(1, 1)]);
        assert!((jtr[0] - 280.0).abs() < 1e-9, "Jtr[0] = {}", jtr[0]);
        assert!((jtr[1] - 56.0).abs() < 1e-9, "Jtr[1] = {}", jtr[1]);
        assert!((err - 280.0).abs() < 1e-9, "chi2 = {err}");
    }

    #[test]
    fn linear_solver_recovers_a_known_solution() {
        // [[2, 0], [0, 4]] x = (6, 8)  ->  x = (3, 2)
        let a = DMatrix::from_diagonal(&DVector::from_vec(vec![2.0, 4.0]));
        let x = solve_linear_system(&a, &DVector::from_vec(vec![6.0, 8.0])).expect("SPD system");
        assert!((x[0] - 3.0).abs() < 1e-12, "x0 = {}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-12, "x1 = {}", x[1]);
    }

    #[test]
    fn linear_solver_still_solves_a_singular_consistent_system() {
        // A = [[1, 1], [1, 1]] is singular (Cholesky fails), so the QR branch
        // must take over.  b = (1, 1) lies in the column space, so the system is
        // consistent and any correct solve must reproduce it exactly:
        //   A x = (x0 + x1, x0 + x1) = (1, 1)  =>  x0 + x1 = 1.
        // The minimum-norm solution would be (0.5, 0.5), but back substitution
        // on the rank-deficient factor legitimately picks a different point of
        // the solution line, so the invariant asserted here is that the
        // returned vector actually solves the system rather than a particular
        // representative of it.
        let a = DMatrix::from_element(2, 2, 1.0);
        let b = DVector::from_vec(vec![1.0, 1.0]);
        let x = solve_linear_system(&a, &b).expect("rank deficient system");
        assert!(x.iter().all(|v| v.is_finite()), "singular solve must stay finite: {x:?}");
        let residual = &a * &x - &b;
        assert!(residual.norm() < 1e-12, "A x - b = {residual:?}, the solve must satisfy the system");
        assert!((x[0] + x[1] - 1.0).abs() < 1e-12, "x = {x:?}");
    }
}
