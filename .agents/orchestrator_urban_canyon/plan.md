# Execution Plan: Urban Canyon Fix Rate Expansion and Multipath Mitigation

## Architecture & Track Division

### 1. E2E Testing Track (Parallel, Opaque-Box)
- Owner: `sub_orch_e2e` / `teamwork_preview_test_writer`
- Scope:
  - Generate `TEST_INFRA.md` covering all 13 features across Tiers 1-4.
  - Implement Tier 1 (Feature Isolation, >= 5 tests per feature).
  - Implement Tier 2 (Boundary & Corner Cases, >= 5 tests per feature).
  - Implement Tier 3 (Cross-Feature Combinations, pairwise interactions).
  - Implement Tier 4 (Real-World Urban Canyon Workload Scenarios).
  - Publish `TEST_READY.md` upon completion.

### 2. Implementation Track (Sequential Milestones)
- **Milestone 1 (M1)**: Adaptive SNR & Elevation Observation Covariance Weighting
  - Scope: F1, F2, F3 (`obs.rs`, `variance.rs`, `formation_cov.rs`).
  - Loop: Explorer -> Worker -> Reviewers (2) -> Challengers (2) -> Auditor -> Gate.
- **Milestone 2 (M2)**: Code-Minus-Carrier (CMC) Multipath Screening & Down-Weighting
  - Scope: F4, F5, F6, F7 (`screen.rs`, `formation.rs`, `mw.rs`, `formation_cov.rs`).
  - Loop: Explorer -> Worker -> Reviewers (2) -> Challengers (2) -> Auditor -> Gate.
- **Milestone 3 (M3)**: Doppler-Assisted Cycle Slip Detection & Phase Continuity Validation
  - Scope: F8, F9 (`screening.rs`, `formation.rs`, `mw.rs`).
  - Loop: Explorer -> Worker -> Reviewers (2) -> Challengers (2) -> Auditor -> Gate.
- **Milestone 4 (M4)**: C/N0- and Elevation-Prioritized PAR Engine
  - Scope: F10, F11, F12 (`par.rs`, `ar_subsets.rs`, `ar.rs`, `tc_ambiguity.rs`).
  - Loop: Explorer -> Worker -> Reviewers (2) -> Challengers (2) -> Auditor -> Gate.
- **Milestone 5 (M5)**: Final Milestone — E2E Test Suite Pass (Tiers 1-4) & Adversarial Hardening (Tier 5)
  - Phase 1: Pass 100% of E2E tests created by E2E Testing Track.
  - Phase 2: Adversarial coverage hardening (Tier 5) via Challenger -> Worker -> Reviewer.
  - Verification on full UrbanNav datasets (Shinjuku, Whampoa, TST1, Odaiba INS).

---

## Execution Schedule
1. **Launch Phase 2 Parallel Tracks**:
   - Spawn E2E Testing Track sub-orchestrator (`sub_orch_e2e`) to build the 4-tier test suite.
   - Spawn Milestone 1 sub-orchestrator (`sub_orch_m1`) to implement Adaptive SNR & Elevation Observation Covariance.
2. **Sequential Flow for Remaining Milestones**:
   - When M1 passes its gate -> Launch Milestone 2 (`sub_orch_m2`) and Milestone 3 (`sub_orch_m3`). Note: M3 (Doppler cycle slips in `screening.rs`) can run concurrently or sequentially after M2 touching `formation.rs`.
   - When M2 & M3 pass their gates -> Launch Milestone 4 (`sub_orch_m4`).
   - When M4 passes and `TEST_READY.md` is published -> Launch Milestone 5 (`sub_orch_m5`).
