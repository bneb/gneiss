# Dispatch: Challenger 1 (Milestone 1 Numerical Stress Testing)

## Objective
Empirically verify the correctness, numerical stability, and derivative continuity of Milestone 1 implementations.

## Target Modules
- `crates/gneiss-core/src/variance.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`

## Tasks
1. Stress test the SIGMA-SNR elevation and noise variance function across millions of random and grid sample evaluations:
   - Elevation $\theta \in [-10^\circ, 90^\circ]$
   - C/N0 $S \in [-20, 60]$ dB-Hz
2. Validate finite difference derivatives: confirm $\frac{\partial \sigma^2}{\partial \theta} \le 0$ and $\frac{\partial \sigma^2}{\partial S} \le 0$ everywhere.
3. Verify that variance is bounded in $(0, 1000.0]$ and never returns NaN or Inf.
4. Verify multi-satellite double-difference covariance matrices $R_{DD}$ are strictly positive definite ($\lambda_{\min}(R_{DD}) > 0$).
5. Deliver `handoff.md` with explicit verdict: `APPROVE` or `REJECT`.
6. Send completion message to parent.


## 2026-09-24T14:40:57Z
You are Challenger 1 for Milestone 1 of the Gneiss Urban Canyon project.
Your working directory is: /Users/kevin/projects/gneiss/.agents/challenger_1_m1
Your parent conversation ID is: c1309e2d-6c95-4b14-a86d-d26a13f2a150

Read your instructions in:
- /Users/kevin/projects/gneiss/.agents/challenger_1_m1/DISPATCH.md
- /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md
- /Users/kevin/projects/gneiss/AGENTS.md
- /Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon/PROJECT.md
- /Users/kevin/projects/gneiss/.agents/worker_m1/handoff.md

Empirically test and stress test variance.rs and formation_cov.rs across extreme elevation [-10, 90] and SNR [-20, 60].
Validate finite-difference derivatives, boundedness, positive-definiteness.
Deliver handoff.md with explicit verdict APPROVE or REJECT and send message back to parent.
