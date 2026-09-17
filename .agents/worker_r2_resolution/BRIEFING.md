# BRIEFING — 2026-09-13T11:57:19Z

## Mission
Diagnose and fix systematic North/Up error in eval_ppp F9P kinematic dataset to achieve sub-meter accuracy with integer ambiguities resolved.

## 🔒 My Identity
- Archetype: worker
- Roles: implementer, qa, specialist
- Working directory: /Users/kevin/projects/gneiss/.agents/worker_r2_resolution
- Original parent: 2d8b2ce8-f45d-4a66-b39c-1149ca0c69bf
- Milestone: PPP-AR F9P Kinematic Accuracy Resolution

## 🔒 Key Constraints
- File size < 500 LOC
- Function size < 32 LOC
- Nesting depth < 3 levels
- Zero unwrap() in production code
- 0 compiler and clippy warnings
- Sub-meter kinematic accuracy on F9P drive dataset vs CSRS-PPP ground truth with integer ambiguities resolved
- Regression guards passing

## Current Parent
- Conversation ID: 2d8b2ce8-f45d-4a66-b39c-1149ca0c69bf
- Updated: 2026-09-13T11:57:19Z

## Task Summary
- **What to build**: Fix systematic North/Up error in PPP kinematic processing for eval_ppp
- **Success criteria**: eval_ppp achieves sub-meter kinematic accuracy vs CSRS-PPP ground truth with resolved ambiguities, all tests and guards pass, full compliance with AGENTS.md
- **Interface contracts**: crates/gneiss-rtk, crates/gneiss-parsers
- **Code layout**: Gneiss workspace

## Key Decisions Made
- Starting investigation into eval_ppp output and prefit residual geometry.

## Artifact Index
- DISPATCH.md — Assignment instructions
- BRIEFING.md — Situational awareness
- progress.md — Heartbeat and status
- handoff.md — Final handoff report

## Change Tracker
- **Files modified**: None yet
- **Build status**: Untested
- **Pending issues**: North/Up error in eval_ppp

## Quality Status
- **Build/test result**: TBD
- **Lint status**: TBD
- **Tests added/modified**: TBD

## Loaded Skills
- None
