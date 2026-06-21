# BRIEFING — 2026-06-20T21:58:01-07:00

## Mission
Verify the changes and tests for Bug 18: Opposite Sign in Phase Wind-Up Correction.

## 🔒 My Identity
- Archetype: reviewer and critic
- Roles: reviewer, critic
- Working directory: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_reviewer_windup
- Original parent: f16afb25-c177-42fe-985d-6840e173046f
- Milestone: sub_orch_milestone_1_tier_1_bugs
- Instance: 1 of 1

## 🔒 Key Constraints
- Review-only — do NOT modify implementation code

## Current Parent
- Conversation ID: f16afb25-c177-42fe-985d-6840e173046f
- Updated: yes

## Review Scope
- **Files to review**:
  - `crates/gneiss-rtk/src/engine/measurement.rs`
  - `crates/gneiss-rtk/src/engine/ppp.rs`
  - `crates/gneiss-rtk/src/engine/ppp_iekf.rs`
  - `crates/gneiss-rtk/src/engine/ppp_ins_iekf.rs`
- **Interface contracts**: `/Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/SCOPE.md`
- **Review criteria**: correctness, style, conformance

## Key Decisions Made
- Issued a verdict of REQUEST_CHANGES because the implementation is missing and the regression test does not exist.

## Review Checklist
- **Items reviewed**: measurement.rs, ppp.rs, ppp_iekf.rs, ppp_ins_iekf.rs
- **Verdict**: request_changes
- **Unverified claims**: all

## Attack Surface
- **Hypotheses tested**: wind-up sign calculation and its application in the EKF/pre-alignment
- **Vulnerabilities found**: wind-up corrections are additive instead of subtractive, and the regression test is missing.
- **Untested angles**: none

## Artifact Index
- `/Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_reviewer_windup/ORIGINAL_REQUEST.md` — Original request
- `/Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_reviewer_windup/BRIEFING.md` — Briefing document
- `/Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_reviewer_windup/progress.md` — Progress heartbeat
- `/Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_reviewer_windup/handoff.md` — Verification findings and review reports
