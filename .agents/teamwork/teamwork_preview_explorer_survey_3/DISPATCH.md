## 2026-09-25T21:04:10Z

You are Survey Explorer 3 for the Gneiss Frame Safety & Epoch Alignment Refactoring.
Your working directory is: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_3
Authoritative request file: /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md

You MUST read /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md and /Users/kevin/projects/gneiss/AGENTS.md first.

Focus: R3 (Estimator State & Pipeline Refactoring) & R4 (Relational Coupling & Frame Safety Invariants) & Benchmarks.
1. Survey EskfState, SwfgEngine, PostProcessOptions, and eval_odaiba_ins in crates/gneiss-rtk.
2. Investigate how double-difference line-of-sight unit vectors and satellite elevations are computed for rover and base (see receiver_pcv.rs dd_correction_m citation in AGENTS.md, formation.rs, etc.).
3. Investigate relational coupling: how to derive rover and base geometry from a single shared ephemeris state / shared satellite position struct so disagreement is unrepresentable.
4. Survey existing test suite (789 tests), smoke benchmarks (scripts/check_network_benchmark.py --smoke, scripts/check_multignss_benchmark.py --smoke), and eval_odaiba_ins baseline performance.
5. Identify all modules/files that will need refactoring.
6. Save your full findings in /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_3/survey_r3_r4_estimators.md and write a complete self-contained handoff.md. Send a completion message when finished.
