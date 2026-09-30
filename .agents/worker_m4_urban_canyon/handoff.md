# Handoff Report: Worker M4 (Milestone 4: Prioritized PAR, DOP Guard, Positive-Definite Ambiguity Conditioning)

**Working Directory**: `/Users/kevin/projects/gneiss/.agents/worker_m4_urban_canyon`  
**Report File**: `/Users/kevin/projects/gneiss/.agents/worker_m4_urban_canyon/handoff.md`  
**Parent Conversation ID**: `5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0`  

---

## 1. Observation

1. **Prior Baseline Weaknesses Observed in Ambiguity Subset Ranking**:
   - `crates/gneiss-rtk/src/ambiguity/par.rs:24-25`:
     ```rust
     let cond_stdevs = compute_conditional_stdevs(q);
     let mut indexed: Vec<(usize, f64)> = (0..n).map(|i| (i, cond_stdevs[i])).collect();
     indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
     ```
     Subsets were ranked strictly by diagonal standard deviation, with no mechanism to ingest physical observation metrics (elevation, SNR, tracking lock duration, CMC multipath variance).
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_subsets.rs:80-84`:
     ```rust
     let score = |i: usize| {
         let frac = (a_float[i] - a_float[i].round()).abs();
         frac + q_amb[(i, i)].sqrt() * 0.5
     };
     c.sort_by(|&i, &j| score(i).total_cmp(&score(j)));
     ```
     In severe urban canyons, multipath-contaminated satellites with coincidental near-integer float values were prioritized ahead of clean, high-elevation, high-SNR satellites.

2. **DOP Geometry and Subset Size Deficiencies**:
   - In `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs:226-243`, `eval_par_subset` evaluated candidate subsets without auditing whether the selected subset of satellites maintained sufficient spatial geometry (PDOP $\le 10.0$) or minimum satellite count ($\ge 4$ unique satellites). Collinear street-canyon subsets risked passing the LAMBDA ratio test while degrading conditional 3D position estimates.

3. **Off-Diagonal Covariance Leakage in Fixed Ambiguity Conditioning**:
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs:349-353`:
     ```rust
     let mut new_cov = &state.cov - cov_red;
     new_cov = 0.5 * (&new_cov + &new_cov.transpose());
     for i in 0..state.dim() { new_cov[(i, i)] = new_cov[(i, i)].max(1e-6); }
     for &idx in &indices { new_cov[(idx, idx)] = 1e-4; }
     state.cov = new_cov;
     ```
     Setting `new_cov[(idx, idx)] = 1e-4` while leaving cross-covariances intact (`state.cov[(0, idx)]` around $0.01$) violated the Cauchy-Schwarz inequality $|P_{ij}| \le \sqrt{P_{ii} P_{jj}}$, producing negative eigenvalues in $P$.
   - In `crates/gneiss-rtk/src/composite/tc_ambiguity.rs:340-344`, while rows and columns were zeroed, post-fix $Q_{aa}$ lacked an explicit numerical symmetrization and spectral floor to ensure positive semi-definiteness ($Q_{aa} \succeq 0$).

4. **Strict Standards Enforcement (AGENTS.md)**:
   - Line count verification:
     - `crates/gneiss-rtk/src/ambiguity/par.rs`: 377 LOC (< 500 LOC)
     - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_subsets.rs`: 423 LOC (< 500 LOC)
     - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs`: 484 LOC (< 500 LOC)
     - `crates/gneiss-rtk/src/composite/tc_ambiguity.rs`: 492 LOC (< 500 LOC)
   - Function length: strictly $\le 32$ LOC across all functions.
   - Nesting depth: strictly $< 3$ levels across all functions.
   - Zero `unwrap()` calls in production code.
   - Zero compiler warnings and zero clippy warnings (`cargo clippy --workspace --all-targets -- -D warnings`).

---

## 2. Logic Chain

1. **CQM Formulation and Fallback Consistency**:
   - The Composite Quality Metric is formulated as:
     $$CQM_i = 0.25 \sin\theta_i + 0.25 \frac{S_i - 20}{30} + 0.20 \min\left(1, \frac{t_{lock}}{30}\right) - 0.15 \frac{\sigma_{cmc}}{2.0} - 0.15 \left(\sigma_a + 0.5 \text{frac}\right)$$
   - When metadata is unavailable, `compute_cqm(None, q_ii, float_val)` returns $-( \text{frac} + 0.5 \sigma_a )$. Sorting descending by CQM preserves identical ordering to the previous baseline ($frac + 0.5 \sigma_a$ ascending).
   - When metadata is present, clean satellites ($\theta = 60^\circ, S = 45\text{ dB-Hz}, \sigma_{cmc} = 0.15\text{m}$) achieve CQM $\approx +0.57$, while multipath-contaminated satellites ($\theta = 12^\circ, S = 22\text{ dB-Hz}, \sigma_{cmc} = 3.5\text{m}$) yield CQM $\approx -0.13$. Even when the corrupted satellite exhibits a coincidental float offset of $0.001$, the clean satellite with $0.25$ float offset is prioritized first for integer fixing.

2. **DOP Geometry Guard and Minimal Size Enforcement**:
   - For 3D positioning, at least 4 non-coplanar satellites ($\ge 3$ double-difference pairs) are required.
   - `validate_subset_geometry` computes PDOP using `gneiss_core::dop::compute_dop_from_positions`.
   - Subsets with $< 4$ satellites or collinear/degenerate geometry ($\text{PDOP} > 10.0$ or singular matrix) are rejected before integer fixing, preventing conditional position divergence.

3. **Definiteness and Elimination of Covariance Leakage**:
   - In `condition_state_on_integers`, when ambiguity index `idx` is fixed:
     `new_cov[(r, idx)] = 0.0; new_cov[(idx, r)] = 0.0;` for all $r \ne idx$, and `new_cov[(idx, idx)] = 1e-4`.
   - The decoupling decouples the fixed ambiguity mode.
   - Minimum eigenvalue across $P$ is verified using `nalgebra::linalg::SymmetricEigen::new(new_cov.clone()).eigenvalues.min()`. If $\lambda_{\min}(P) < 10^{-6}$, the integer conditioning is rejected to prevent filter instability.
   - In `tc_ambiguity.rs`, $Q_{aa}$ is symmetrized and floored at $10^{-6}$, guaranteeing $Q_{aa} \succeq 0$.

---

## 3. Caveats

- **No Caveats**: All changes strictly honor interface contracts, maintain backward compatibility, and operate within the existing Kalman filter and ESKF architectures without modifying non-owned files.

---

## 4. Conclusion

Milestone 4 (R4: Prioritized PAR, DOP Guard, Positive-Definite Ambiguity Conditioning) is fully implemented, verified, and passes all workspace tests and CI smoke guards:
1. **F10**: `AmbiguityMetadata` and `compute_cqm` added in `par.rs`; `compute_metadata_from_dd`, `select_par_candidates_with_metadata`, and `select_par_candidates_with_dd` added in `ar_subsets.rs`.
2. **F11**: `validate_subset_geometry` and `validate_dd_subset_geometry` added in `ar_subsets.rs`; DOP guard and size $\ge 4$ enforced in `eval_par_subset` in `ar.rs`.
3. **F12**: Off-diagonal covariance leakage zeroed out and $\lambda_{\min}(P) \ge 10^{-6}$ verified in `ar.rs:condition_state_on_integers`; $Q_{aa} \succeq 0$ symmetry and spectral floor guaranteed in `tc_ambiguity.rs`.
4. **Code Quality**: All 4 files strictly $< 500$ LOC, functions $\le 32$ LOC, nesting $< 3$, 0 `unwrap()`, 0 clippy warnings.

---

## 5. Verification Method

### 1. Verification Commands & Results

- **Clippy**:
  ```bash
  cargo clippy --workspace --all-targets -- -D warnings
  ```
  Result: `Finished dev profile [unoptimized + debuginfo] target(s)` — **0 warnings**.

- **Library Unit & Integration Tests**:
  ```bash
  cargo test -p gneiss-rtk --lib
  ```
  Result: **449 passed; 0 failed**.

- **Urban Canyon E2E Test Suite**:
  ```bash
  cargo test --test test_urban_canyon_e2e
  ```
  Result: **51 passed; 0 failed**.

- **Network RTK/PPK Smoke Benchmark**:
  ```bash
  python3 scripts/check_network_benchmark.py --smoke
  ```
  Result: **ALL CHECKS PASSED** (horizontal p50: 0.023m <= 0.04m, fix rates: 99.1%, 99.7%, 70.8%).

- **Multi-GNSS Smoke Benchmark**:
  ```bash
  python3 scripts/check_multignss_benchmark.py --smoke
  ```
  Result: **ALL CHECKS PASSED** (network fused fix rate: 99.90% >= 96.5%).

### 2. Files to Inspect

- `crates/gneiss-rtk/src/ambiguity/par.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_subsets.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs`
- `crates/gneiss-rtk/src/composite/tc_ambiguity.rs`

### 3. Invalidation Conditions

If a test case demonstrates that fixed ambiguities retain non-zero cross-covariances with position states, or that candidate ranking prioritizes a severe multipath satellite ahead of a clean satellite when metadata is provided, this conclusion would be invalidated.
