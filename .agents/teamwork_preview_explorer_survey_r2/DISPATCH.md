# Task Assignment — Survey Frontier R2: Integer PPP-AR Engine via SINEX OSB Ingestion

## Role & Mission
You are a read-only exploration agent (`teamwork_preview_explorer`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r2/`.
The authoritative user request is: `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md`.

## Objective
Survey the codebase for Frontier R2: Integer PPP-AR Engine via SINEX OSB Ingestion:
1. Satellite Observation-Specific Bias (OSB) and fractional phase bias ingestion from SINEX files (`.BIA` / `.OSB`).
2. Un-differenced carrier-phase and pseudorange observation equations with exact satellite and receiver phase center offsets/variations (PCO/PCV).
3. Recovery of integer wide-lane (Melbourne-Wübbena) and narrow-lane ambiguities via LAMBDA search on single-difference or receiver-clock-decoupled ambiguities without physical base stations.
4. Continuous carrier tracking across GPS, Galileo, BeiDou, QZSS.
5. Benchmark: `crates/gneiss-rtk/src/bin/eval_ppp.rs` resolving integer ambiguities on F9P kinematic drive to sub-meter kinematic accuracy vs CSRS-PPP.

## Scope of Investigation
- Locate and examine all existing PPP modules (`crates/gneiss-rtk/src/ppp.rs`), SINEX/OSB parsers (`crates/gneiss-format` or similar), antenna models (PCO/PCV in `crates/gneiss-core/src/antex.rs` or `crates/gneiss-rtk`), LAMBDA implementation (`crates/gneiss-rtk/src/lambda.rs` or similar).
- Check `crates/gneiss-rtk/src/bin/eval_ppp.rs` and F9P dataset + CSRS-PPP truth: how is PPP currently evaluated, what accuracy is currently obtained ($0.296\text{ m}$ RMS reference vs CSRS-PPP), what components are missing for integer fixing (OSB ingestion, MW wide-lane fixing, narrow-lane fixing, LAMBDA integration).
- Map exact data structures, functions, interfaces, dependencies, and files that need to be created or modified.
- Verify adherence to AGENTS.md constraints: file size < 500 LOC, function size < 32 LOC, nesting < 3 levels, 0 unwrap() in prod, zero warnings.

## Deliverable
Write your complete survey report to:
`/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r2/survey_r2.md`
and write a standard `handoff.md` in your directory.

## 2026-09-12T16:41:15Z
Survey Explorer R2 dispatched to survey Frontier R2: Integer PPP-AR Engine via SINEX OSB Ingestion.
Target outputs: survey_r2.md and handoff.md.

