# Dispatch: Worker R1 — Frontier R1 (15-State ESKF Benchmark Verification)

Working Directory: /Users/kevin/projects/gneiss/.agents/worker_r1
Original Request: /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md
Master Plan: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md
Previous Progress: /Users/kevin/projects/gneiss/.agents/worker_m1/progress.md

## Objective
Verify and complete the benchmark for Frontier R1: 15-State Error-State Kalman Filter (ESKF/MEKF) for GNSS/INS in `crates/gneiss-rtk/src/estimators/eskf/` and `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs`.
Target metrics on Odaiba dataset (12,398 epochs):
- p50 horizontal error < 2.5 m
- RMS horizontal error < 5.2 m

## Exclusive File Ownership
- `crates/gneiss-rtk/src/estimators/eskf/**`
- `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs`

## Mandatory Rules & Integrity Warning
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Follow AGENTS.md standards:
- < 500 LOC per file
- < 32 LOC per function
- < 3 nesting depth
- 0 unwrap() in production code
- 0 clippy warnings (`cargo clippy -p gneiss-rtk --bin eval_odaiba_ins -- -D warnings`)
- All tests pass (`cargo test -p gneiss-rtk`)

## 2026-09-12T22:25:50Z
You are Worker R1 for the Gneiss positioning engine project.

Your assigned working directory is: /Users/kevin/projects/gneiss/.agents/worker_r1
The authoritative user request is in: /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md
The project plan is in: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md
Your dispatch instructions are in: /Users/kevin/projects/gneiss/.agents/worker_r1/DISPATCH.md
Previous progress on R1 is in: /Users/kevin/projects/gneiss/.agents/worker_m1/progress.md

MANDATORY FIRST STEPS:
1. Read ORIGINAL_REQUEST.md, PROJECT.md, and your DISPATCH.md.
2. Initialize BRIEFING.md and progress.md in your working directory.

OBJECTIVE:
Verify and complete the benchmark for Frontier R1: 15-State Error-State Kalman Filter (ESKF/MEKF) for GNSS/INS in `crates/gneiss-rtk/src/estimators/eskf/` and `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs`.
Target metrics on the Odaiba 12,398-epoch 10Hz trajectory:
- p50 horizontal error < 2.5 m
- RMS horizontal error < 5.2 m

EXCLUSIVE FILE OWNERSHIP:
- crates/gneiss-rtk/src/estimators/eskf/**
- crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

QUALITY STANDARDS (AGENTS.md):
- File size < 500 LOC
- Function size < 32 LOC
- Nesting depth < 3 levels
- 0 unwrap() in production code (unwrap is only allowed in #[cfg(test)])
- Zero compiler / clippy warnings (`cargo clippy -p gneiss-rtk --bin eval_odaiba_ins -- -D warnings`)
- Unit tests pass (`cargo test -p gneiss-rtk --lib estimators::eskf`)

EXECUTION & VERIFICATION:
1. Check existing implementation in crates/gneiss-rtk/src/estimators/eskf/ and eval_odaiba_ins.rs.
2. Run unit tests (`cargo test -p gneiss-rtk --lib estimators::eskf`).
3. Run the benchmark: `cargo run --release --bin eval_odaiba_ins`.
4. Inspect the output errors. If p50 < 2.5m and RMS < 5.2m are met, document the exact metrics. If not, inspect the filter tuning, lever-arm correction, NHC/ZUPT weighting, or backward RTS smoother integration and refine until the target metrics are achieved.
5. Ensure 0 clippy warnings and 0 unwrap() in production code.
6. Write a complete handoff report to `/Users/kevin/projects/gneiss/.agents/worker_r1/handoff.md` and report completion back via send_message.

