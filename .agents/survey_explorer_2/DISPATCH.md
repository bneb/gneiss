# Dispatch: Survey Explorer 2 (R2 CMC Multipath Mitigation)

## Objective
Map the full scope and existing codebase baseline for:
R2: Code-Minus-Carrier (CMC) Multipath Detection & Down-Weighting.

## Required Reading
- `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (specifically timestamp 2026-09-24T13:30:49Z)
- `/Users/kevin/projects/gneiss/AGENTS.md`

## Key Code Locations
- `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/update/robust.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_gate.rs`

## Instructions
- Investigate how double-difference carrier and pseudorange observations are formed and paired.
- Determine how continuous carrier tracking arcs are identified and tracked across epochs.
- Analyze how CMC ($CMC = P - \Phi - 2I$ or dual-frequency geometry-free combinations) can be maintained per satellite pair/arc.
- Investigate the robust IEKF update (`update/robust.rs`): how innovation screening, Huber/Tukey reweighting, or measurement variance inflation operates.
- Design the exact architecture for detecting code multipath jumps (5m–20m code steps) without flagging clean carrier phase.
- Ensure that down-weighting pseudorange preserves carrier phase weights and prevents corrupted measurements from entering the integer search.
- Write your comprehensive findings to `/Users/kevin/projects/gneiss/.agents/survey_explorer_2/survey_r2_cmc.md` and deliver `handoff.md`.
- Send a message to parent upon completion.

## 2026-09-24T13:33:14Z
Task:
1. Map the codebase baseline and implementation points for R2: Code-Minus-Carrier (CMC) Multipath Detection & Down-Weighting.
2. Inspect crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs, formation_cov.rs, update/robust.rs, mw.rs, and ar_gate.rs.
3. Detail how double differences are formed, how tracking arcs are identified, how CMC residuals are filtered/tracked, and how measurement down-weighting/variance inflation interacts with the robust IEKF update.
4. Formulate the algorithm for detecting 5m-20m pseudorange multipath steps without flagging clean carrier phase.
5. Document all findings in /Users/kevin/projects/gneiss/.agents/survey_explorer_2/survey_r2_cmc.md.
6. Write your handoff.md and send a message back to parent when complete.

