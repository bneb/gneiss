# Dispatch: Survey Explorer 1 (R1 & Baseline Benchmarks)

## Objective
Map the full scope and existing codebase baseline for:
1. R1: Adaptive C/N0 (SNR) and elevation observation covariance weighting in `crates/gneiss-core/src/variance.rs`, `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`, and `formation.rs`.
2. Baseline benchmarks: Run and document baseline numbers for:
   - `python3 scripts/check_network_benchmark.py --smoke`
   - `python3 scripts/check_multignss_benchmark.py --smoke`
   - `cargo run --release --bin eval_odaiba_ins` (or relevant UrbanNav / F9P harness)
   - `cargo test --workspace` and `cargo clippy --workspace --all-targets -- -D warnings`

## Required Reading
- `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (specifically timestamp 2026-09-24T13:30:49Z)
- `/Users/kevin/projects/gneiss/AGENTS.md`

## Instructions
- Investigate how variance/covariance is currently formed for double differences.
- Check how SNR (C/N0) and elevation are currently represented and whether C/N0 is already passed into `formation_cov.rs` or if observation structs need enhancement.
- Derive the exact SIGMA-SNR mathematical formulation:
  $$\sigma^2(\theta, S) = \left( a^2 + \frac{b^2}{\sin^2 \theta} \right) \cdot f_{\text{SNR}}(S)$$
  with smooth derivatives and exponential noise scaling below 38-40 dB-Hz.
- Identify all files and structs that need modification for R1.
- Write your comprehensive findings to `/Users/kevin/projects/gneiss/.agents/survey_explorer_1/survey_r1_bench.md` and deliver `handoff.md`.
- Send a message to parent upon completion.

## 2026-09-24T13:33:14Z
You are Survey Explorer 1 for the Gneiss Urban Canyon Fix Rate Expansion and Multipath Mitigation task.
Your working directory is: /Users/kevin/projects/gneiss/.agents/survey_explorer_1
Your parent conversation ID is: c1309e2d-6c95-4b14-a86d-d26a13f2a150

Read your instructions in:
- /Users/kevin/projects/gneiss/.agents/survey_explorer_1/DISPATCH.md
- /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md (specifically 2026-09-24T13:30:49Z)
- /Users/kevin/projects/gneiss/AGENTS.md

Task:
1. Map the codebase baseline and implementation points for R1: Adaptive C/N0 (SNR) and elevation observation covariance weighting.
2. Formulate the exact SIGMA-SNR mathematical model and inspect variance.rs, formation_cov.rs, and formation.rs.
3. Check and run the baseline benchmarks:
   - python3 scripts/check_network_benchmark.py --smoke
   - python3 scripts/check_multignss_benchmark.py --smoke
   - cargo test --workspace
   - check any existing eval_f9p_rover or eval_odaiba_ins results
4. Document all findings in /Users/kevin/projects/gneiss/.agents/survey_explorer_1/survey_r1_bench.md.
5. Write your handoff.md and send a message back to parent when complete.

