# Task Assignment — Survey Frontier R1: 15-State ESKF/MEKF GNSS/INS

## Role & Mission
You are a read-only exploration agent (`teamwork_preview_explorer`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r1/`.
The authoritative user request is: `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md`.

## Objective
Survey the codebase for Frontier R1: 15-State Error-State Kalman Filter (ESKF/MEKF) for GNSS/INS:
1. State expansion to 15 states ($\delta \mathbf{p}^e, \delta \mathbf{v}^e, \delta \boldsymbol{\theta}, \delta \mathbf{b}_a, \delta \mathbf{b}_g$).
2. Error-quaternion feedback to nominal attitude ($\mathbf{q} \leftarrow \mathbf{q} \otimes \delta\mathbf{q}$).
3. Closed-loop online accelerometer and gyroscope bias estimation driven by GNSS position and velocity innovations.
4. Full 15-state backward Rauch-Tung-Striebel (RTS) smoother over forward filter history.
5. Dynamic vehicle Non-Holonomic Constraints (NHC) and Zero-Velocity Updates (ZUPT) integrated into 15-state covariance.
6. Benchmark: `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs` achieving $p_{50} < 2.5\text{ m}$ and RMS $< 5.2\text{ m}$ across full 12,398-epoch 10Hz trajectory.

## Scope of Investigation
- Locate and examine all existing INS/GNSS modules in `crates/gneiss-rtk` and `crates/gneiss-core` (e.g. `predictor.rs`, `ins/`, filter state representations, Jacobians, covariance propagation).
- Check `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs` and the Odaiba dataset: how are IMU data and GNSS observations ingested, what is the current baseline performance, what are the current state definitions and gaps.
- Map exact data structures, functions, interfaces, dependencies, and files that need to be created or modified.
- Verify adherence to AGENTS.md constraints: file size < 500 LOC, function size < 32 LOC, nesting < 3 levels, 0 unwrap() in prod, zero warnings.

## Deliverable
Write your complete survey report to:
`/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r1/survey_r1.md`
and write a standard `handoff.md` in your directory.
Send a message back to the orchestrator when complete.

## 2026-09-12T16:41:15Z
You are Survey Explorer R1.
Your working directory is: /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r1/
Please read your task instructions in /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r1/DISPATCH.md and the authoritative request in /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md.
Investigate the codebase for Frontier R1 (15-State ESKF/MEKF GNSS/INS, RTS smoother, NHC/ZUPT, eval_odaiba_ins benchmark).
Output your survey report to /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r1/survey_r1.md and your handoff.md. Send a completion message to the parent orchestrator when done.

