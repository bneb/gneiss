## 2026-09-13T02:06:28Z
You are Worker R1 tasked with achieving the Odaiba INS benchmark acceptance targets on Frontier R1.

Your working directory is: /Users/kevin/projects/gneiss/.agents/worker_r1_tune/

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Authoritative user request (MUST read first): /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md
Master project plan: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md
Explorer R1 Report: /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r1_status/report.md
Explorer R1 Handoff: /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r1_status/handoff.md

FILE OWNERSHIP:
You have exclusive write access to:
- crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs
- crates/gneiss-rtk/src/estimators/eskf/** (only if necessary)

GOALS & ACCEPTANCE CRITERIA:
1. Run `cargo run --release --bin eval_odaiba_ins` and verify that the Tokyo Odaiba benchmark achieves:
   - p50 < 2.5 m
   - RMS < 5.2 m
   across the full 12,398-epoch 10Hz trajectory (current performance: p50 = 2.761m, RMS = 5.472m).
2. Implement the concrete tuning recommended by Explorer R1:
   a. Adaptive measurement covariance inflation / Chi-square innovation gating in `update_gnss_innovation` to downweight/reject multipath spikes > 3-5m during elevated railway/highway occlusions.
   b. Doppler-gate or reduce reliance on noisy finite-differenced GNSS velocity in GNSS updates.
   c. Tune velocity process noise random walk $q_a$ (e.g. from 0.05 to ~0.005-0.01) to improve inertial coasting through urban canyon multipath.
   d. Ensure physical antenna lever arm is appropriately parameterized.
3. Verify all AGENTS.md standards:
   - All files < 500 LOC
   - All functions < 32 LOC
   - Nesting depth < 3
   - 0 unwrap() in production code
   - 0 clippy warnings (`cargo clippy -p gneiss-rtk --bin eval_odaiba_ins -- -D warnings`)
   - All unit tests pass (`cargo test -p gneiss-rtk --lib estimators::eskf`)
4. Document all results and write a complete handoff report to:
   /Users/kevin/projects/gneiss/.agents/worker_r1_tune/handoff.md
5. Send your completion message back to the orchestrator.
