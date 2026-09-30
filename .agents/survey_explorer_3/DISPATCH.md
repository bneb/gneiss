# Dispatch: Survey Explorer 3 (R3 Doppler Slip Validation & R4 Prioritized PAR)

## Objective
Map the full scope and existing codebase baseline for:
1. R3: Doppler-assisted cycle slip detection & phase continuity validation.
2. R4: C/N0- and elevation-prioritized Partial Ambiguity Resolution (PAR).

## Required Reading
- `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (specifically timestamp 2026-09-24T13:30:49Z)
- `/Users/kevin/projects/gneiss/AGENTS.md`

## Key Code Locations
- `crates/gneiss-rtk/src/ambiguity/par.rs`
- `crates/gneiss-rtk/src/ambiguity/lambda.rs`
- `crates/gneiss-rtk/src/ambiguity/ffrt.rs`
- `crates/gneiss-rtk/src/composite/tc_ambiguity.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_gate.rs`
- `crates/gneiss-rtk/src/estimators/doppler/`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`

## Instructions
- For R3:
  - Check how Doppler measurements ($\dot{\Phi} = -\lambda f_D$) are currently parsed, stored, and filtered.
  - Determine how phase increments $\Delta\Phi_k = \Phi_k - \Phi_{k-1}$ are compared against integrated Doppler range rate $\int \dot{\rho} dt$ over epoch interval $\Delta t$.
  - Trace how cycle slips are flagged and how ambiguity states/covariances are re-seeded or inflated upon detection.
- For R4:
  - Examine `par.rs` and `tc_ambiguity.rs`. How is candidate subset selection currently implemented? Is it arbitrary truncation, or sorted by variance?
  - Detail how to order candidates by elevation, C/N0, lock duration, and CMC residual variance.
  - Check how DOP guards, minimum subset size ($\ge 4$), and positive semi-definiteness ($Q_{aa} \succeq 0, P \succ 0$) are preserved during subset fixing.
- Write your comprehensive findings to `/Users/kevin/projects/gneiss/.agents/survey_explorer_3/survey_r3_r4.md` and deliver `handoff.md`.
- Send a message to parent upon completion.

## 2026-09-24T13:33:14Z
Task:
1. Map the codebase baseline and implementation points for R3: Doppler-assisted cycle slip detection & phase continuity validation (inspect doppler, rtk_iekf, mw.rs, Rinex/obs parsing).
2. Map the codebase baseline and implementation points for R4: C/N0- and elevation-prioritized Partial Ambiguity Resolution (PAR) (inspect par.rs, lambda.rs, ffrt.rs, tc_ambiguity.rs, ar.rs).
3. Detail how Doppler range rate compares with phase increments, how slip flags are propagated and how ambiguity covariance is re-seeded.
4. Detail how PAR subset candidate ordering should incorporate elevation, C/N0, tracking lock duration, and CMC variance, while preserving DOP and covariance definiteness (Q_aa >= 0, P > 0).
5. Document all findings in /Users/kevin/projects/gneiss/.agents/survey_explorer_3/survey_r3_r4.md.
6. Write your handoff.md and send a message back to parent when complete.
