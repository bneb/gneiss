# Task Assignment — Sub-Orchestrator Milestone M1: 15-State ESKF/MEKF GNSS/INS

## Role & Mission
You are a sub-orchestrator (`teamwork_preview_orchestrator`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/sub_orch_m1/`.
Your parent orchestrator is: `1bd6ce81-03bf-4c40-b8b1-3b137333b5e7`.

## Authoritative Documents to Read Before Starting Work
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (authoritative user requirements)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md` (project architecture, interface contracts, and file layout)
3. `/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r1/survey_r1.md` (comprehensive mathematical derivations, state equations, and baseline metrics)
4. `/Users/kevin/projects/gneiss/AGENTS.md` (strict code quality standards)

## Scope & Milestone Objectives
Implement Frontier R1: 15-State Error-State Kalman Filter (ESKF/MEKF) for GNSS/INS:
1. State expansion to 15 states ($\delta \mathbf{p}^e, \delta \mathbf{v}^e, \delta \boldsymbol{\theta}, \delta \mathbf{b}_a, \delta \mathbf{b}_g$).
2. Transition matrix $\boldsymbol{\Phi}_{15\times 15}$ with strictly POSITIVE velocity-attitude coupling `vel_att = +f_e_skew * dt` (sign convention mandated by `AGENTS.md` and `predictor.rs:86-91`), and spectral process noise $\mathbf{Q}_{15\times 15}$.
3. Multiplicative error-quaternion feedback to nominal attitude ($\mathbf{q} \leftarrow \mathbf{q} \otimes \delta\mathbf{q}$) and state reset after measurement updates.
4. Closed-loop online accelerometer and gyroscope bias estimation driven by GNSS position and velocity innovations.
5. Dynamic vehicle Non-Holonomic Constraints (NHC) with attitude coupling Jacobian $\mathbf{E}_{23}(\mathbf{R}_b^e)^T[\mathbf{v}^e\times]$ and Zero-Velocity Updates (ZUPT) integrated into 15-state covariance.
6. Full 15-state backward Rauch-Tung-Striebel (RTS) smoother over forward filter history.
7. Benchmark validation: `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs` achieving $p_{50} < 2.5\text{ m}$ and $\text{RMS} < 5.2\text{ m}$ across the full 12,398-epoch trajectory (improving over the current 6-state baseline $p_{50}=2.907\text{ m}$, $\text{RMS}=5.508\text{ m}$).

## Exclusive Write Ownership
You and your dispatched workers own ONLY:
- `crates/gneiss-rtk/src/estimators/eskf/` (`mod.rs`, `types.rs`, `predict.rs`, `update.rs`, `constraints.rs`, `smoother.rs`)
- `crates/gneiss-rtk/src/estimators/mod.rs`
- `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs`
Do NOT write to files owned by other milestones.

## Orchestrator Procedure & Gate Verification
Execute via the standard iteration loop: Explorer -> Worker -> Reviewer -> Challenger -> Forensic Auditor -> Gate.
- Worker dispatch prompt MUST include the mandatory integrity warning verbatim.
- Forensic Auditor verdict is a non-negotiable binary veto.
- All code must satisfy AGENTS.md (< 500 LOC/file, < 32 LOC/function, < 3 nesting, 0 unwrap in prod, 0 warnings).
- `cargo check --bin eval_odaiba_ins` and `cargo test -p gneiss-rtk` must pass cleanly.
- `cargo run --release --bin eval_odaiba_ins` must demonstrate $p_{50} < 2.5\text{ m}$ and $\text{RMS} < 5.2\text{ m}$.

When complete, write `handoff.md` and notify parent orchestrator (`1bd6ce81-03bf-4c40-b8b1-3b137333b5e7`).
