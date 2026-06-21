# BRIEFING — 2026-06-21T04:02:00Z

## Mission
Investigate Bug 1: Melbourne-Wübbena Dimensional Typo in `gneiss-rtk/src/engine/ppp_math.rs`

## 🔒 My Identity
- Archetype: Teamwork explorer
- Roles: Read-only investigation: analyze problems, synthesize findings, produce structured reports
- Working directory: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_mw
- Original parent: f16afb25-c177-42fe-985d-6840e173046f
- Milestone: sub_orch_milestone_1_tier_1_bugs

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- Propose a precise fix strategy without implementing it
- Output findings to handoff.md and notify the orchestrator via send_message

## Current Parent
- Conversation ID: f16afb25-c177-42fe-985d-6840e173046f
- Updated: not yet

## Investigation State
- **Explored paths**:
  - `crates/gneiss-rtk/src/engine/ppp_math.rs`
  - `crates/gneiss-rtk/src/engine/ppp.rs`
  - `crates/gneiss-rtk/src/engine/ambiguity.rs`
  - `crates/gneiss-rtk/src/engine/tcar.rs`
  - `crates/gneiss-rtk/src/measurements/combinations.rs`
- **Key findings**:
  - In `crates/gneiss-rtk/src/engine/ppp_math.rs`, the Melbourne-Wübbena combination in cycles was implemented using:
    `let mw = (cp1 - cp2) - (p1 / lam1 + p2 / lam2) * (lam1 * lam2) / (lam1 + lam2);`
  - The term `(lam1 * lam2) / (lam1 + lam2)` is $\lambda_{NL}$ (units of meters), creating a dimensional mismatch by multiplying a cycles term `(p1 / lam1 + p2 / lam2)` by meters.
  - The correct dimensionless ratio is $\lambda_{NL} / \lambda_{WL} = (\lambda_2 - \lambda_1) / (\lambda_1 + \lambda_2)$.
  - This typo was fixed in commit `da013e27a4be9319e98e5389afa792140f4b49f4` by changing the scaling factor to `(lam2 - lam1) / (lam1 + lam2)`.
- **Unexplored areas**: None

## Key Decisions Made
- Confirmed the dimensional mismatch mathematically and analyzed the git commit history to verify how it was resolved in the workspace.

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_mw/handoff.md — Handoff report containing observations, logic chain, caveats, conclusion, and verification method.
