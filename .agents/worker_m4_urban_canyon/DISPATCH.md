# Task Assignment — Worker M4: C/N0- and Elevation-Prioritized Partial Ambiguity Resolution (PAR)

## Role & Mission
You are an implementation worker (`teamwork_preview_worker`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/worker_m4_urban_canyon/`.
Your parent orchestrator is: `5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0`.

## Authoritative Documents to Read Before Starting Work
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (authoritative user requirements, Section R4)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md` (project specifications, interface contracts, and code layout)
3. `/Users/kevin/projects/gneiss/.agents/survey_explorer_3/survey_r3_r4.md` (detailed R4 mathematical derivations, code inspection, and implementation specifications)
4. `/Users/kevin/projects/gneiss/.agents/survey_explorer_3/handoff.md` (summary of CQM ranking, DOP guard, and covariance leakage fix)
5. `/Users/kevin/projects/gneiss/AGENTS.md` (strict code quality standards)

## MANDATORY INTEGRITY WARNING
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## Code Quality Standards (AGENTS.md)
- File size strictly < 500 LOC (CRITICAL: `crates/gneiss-rtk/src/composite/tc_ambiguity.rs` is at 489 LOC, `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs` is at 463 LOC — do not exceed 500 LOC!)
- Function size strictly <= 32 LOC
- Nesting depth strictly < 3 levels
- Exactly 0 `unwrap()` calls in production code (`match`, `if let`, `ok_or()?`, or descriptive `.expect()` only; `unwrap()` in `#[cfg(test)]` is acceptable)
- Zero compiler warnings, zero clippy warnings (`cargo clippy --workspace --all-targets -- -D warnings`)
- All tests pass (`cargo test --workspace`)
- Both CI smoke guard scripts pass:
  - `python3 scripts/check_network_benchmark.py --smoke`
  - `python3 scripts/check_multignss_benchmark.py --smoke`

## Exclusive Write Ownership
You own ONLY:
- `crates/gneiss-rtk/src/ambiguity/par.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_subsets.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs`
- `crates/gneiss-rtk/src/composite/tc_ambiguity.rs`
Do NOT edit any files outside this set without coordinating.

## Scope & Implementation Objectives (Milestone 4: R4)

1. **Composite Quality Metric (CQM) Candidate Prioritization (`ar_subsets.rs` & `par.rs`)**:
   - In `ar_subsets.rs` / `par.rs`:
     - Currently `select_par_candidates` and `select_ils_subset` rank candidates purely by fractional float offset or variance $Q_{ii}$, ignoring elevation, SNR, lock duration, and CMC multipath variance.
     - Implement candidate ranking using a Composite Quality Metric (CQM):
       $$CQM_i = w_{el} \sin\theta_i + w_{snr} \frac{S_i - 20}{30} + w_{lock} \min\left(1, \frac{t_{lock}}{30}\right) - w_{cmc} \frac{\sigma_{cmc}}{2.0} - w_{var} \sigma_{a_i}$$
       (e.g., $w_{el}=0.25, w_{snr}=0.25, w_{lock}=0.20, w_{cmc}=0.15, w_{var}=0.15$ or calibrated weights).
     - Provide backward-compatible ranking when metadata is omitted (fallback to variance/fractional score).
     - Prioritize candidate subsets in descending order of CQM quality.

2. **DOP Geometry Guard & Minimum Subset Size Constraint (`ar_subsets.rs` & `tc_ambiguity.rs`)**:
   - Enforce minimum subset size $\ge 4$ satellites for any PAR candidate subset.
   - Add geometric Dilution of Precision (DOP) check (or reject degenerate subsets with collinear geometry) to ensure fixed subsets maintain stable 3D positioning.

3. **Strict Positive-Definite Ambiguity Conditioning (`ar.rs` & `tc_ambiguity.rs`)**:
   - In `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs:349-353` (`condition_state_on_integers`):
     - Fix off-diagonal covariance leakage: When setting fixed ambiguity diagonal variance `new_cov[(idx, idx)] = 1e-4`, explicitly zero out all cross-covariance row and column elements:
       `new_cov[(r, idx)] = 0.0;`
       `new_cov[(idx, r)] = 0.0;` for all $r \ne idx$.
     - Verify $\lambda_{\min}(P) \ge 10^{-6}$ and enforce positive definiteness ($P \succ 0$).
   - In `tc_ambiguity.rs`: Confirm integer conditioning preserves positive semi-definiteness ($Q_{aa} \succeq 0$) without off-diagonal leakage.

4. **Unit Tests & Invariants**:
   - Add unit tests verifying:
     - CQM ranks high-elevation, high-SNR, clean CMC satellites ahead of low-elevation or multipath-corrupted satellites even if the corrupted float offset is coincidental near an integer.
     - Subsets with $< 4$ satellites or bad DOP are rejected.
     - `condition_state_on_integers` produces strictly positive-definite covariance matrix with zero cross-covariances on fixed ambiguities.
   - Verify all files < 500 LOC, functions <= 32 LOC, nesting < 3, 0 unwrap in prod.

## Verification Commands
Before delivering your handoff, execute and verify:
```bash
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p gneiss-rtk --lib
cargo test --test test_urban_canyon_e2e
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke
```
## 2026-09-24T23:50:55Z
You are Worker M4 for Milestone 4: C/N0- and Elevation-Prioritized Partial Ambiguity Resolution (PAR).
Your working directory is: /Users/kevin/projects/gneiss/.agents/worker_m4_urban_canyon/
Read your task assignment in /Users/kevin/projects/gneiss/.agents/worker_m4_urban_canyon/DISPATCH.md, the authoritative user request in /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md, the project specifications in /Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md, and the deep survey findings in /Users/kevin/projects/gneiss/.agents/survey_explorer_3/survey_r3_r4.md and handoff.md.

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Strictly adhere to all AGENTS.md rules:
- File size < 500 LOC (CRITICAL: tc_ambiguity.rs is at 489 LOC, ar.rs at 463 LOC, par.rs at 284 LOC, ar_subsets.rs at 285 LOC — do not exceed 500 LOC!)
- Function size <= 32 LOC
- Nesting depth < 3 levels
- Zero unwrap() in production code
- Zero clippy/compiler warnings (cargo clippy --workspace --all-targets -- -D warnings)
- All workspace tests pass (cargo test -p gneiss-rtk --lib, cargo test --test test_urban_canyon_e2e)
- Both CI smoke guard scripts pass:
  python3 scripts/check_network_benchmark.py --smoke
  python3 scripts/check_multignss_benchmark.py --smoke

Implement F10 (CQM candidate ranking by elevation, C/N0, lock time, CMC variance in ar_subsets.rs / par.rs), F11 (DOP guard and minimum subset size >= 4), and F12 (fix off-diagonal covariance leakage in ar.rs and tc_ambiguity.rs to enforce P > 0 and Q_aa >= 0).

When all implementation and verifications pass, write /Users/kevin/projects/gneiss/.agents/worker_m4_urban_canyon/handoff.md and send a completion message to your parent.

## 2026-09-25T07:11:23Z
**Context**: Resuming Milestone 4 (R4: Prioritized PAR)
**Content**: The parent orchestrator has checked in. CQM ranking is implemented in `crates/gneiss-rtk/src/ambiguity/par.rs`.
Please finalize any remaining items for M4:
1. Ensure F10 (CQM ranking), F11 (DOP guard / min subset size >= 4 in `ar_subsets.rs` / `par.rs` / `tc_ambiguity.rs`), and F12 (zeroing off-diagonal cross-covariances on fixed ambiguities in `ar.rs:condition_state_on_integers` and `tc_ambiguity.rs`) are complete.
2. Verify all AGENTS.md rules (< 500 LOC per file, functions <= 32 LOC, nesting < 3, 0 unwrap in prod).
3. Run all verification commands:
   - `cargo clippy --workspace --all-targets -- -D warnings`
   - `cargo test -p gneiss-rtk --lib`
   - `cargo test --test test_urban_canyon_e2e`
   - `python3 scripts/check_network_benchmark.py --smoke`
   - `python3 scripts/check_multignss_benchmark.py --smoke`
4. Write your complete handoff report to `/Users/kevin/projects/gneiss/.agents/worker_m4_urban_canyon/handoff.md`.
**Action**: Execute the remaining verifications, write `handoff.md`, and reply with your completion report.
