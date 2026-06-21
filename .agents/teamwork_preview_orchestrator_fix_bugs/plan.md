# Plan: Gneiss Navigation Engine Bug Fixes

This plan outlines the steps for executing the bug fixes across 3 main milestones, matching the tier and priority requirements.

## Execution Strategy

We will use a multi-agent orchestration pattern to delegate the implementation of these fixes. Because there are 25 bugs, we will spawn **sub-orchestrators** for each milestone to keep context sizes reasonable and respect spawn count thresholds (self-succession at 16 spawns).

### Milestone 1: Tier 1 Bugs (Ranks 1-8)
- **Objective**: Fix the 8 most critical, high-impact bugs.
- **Verification**: Ensure regression tests fail when the buggy code is restored and pass when fixed.
- **Delegation**: Spawn a sub-orchestrator `teamwork_preview_orchestrator` with scope `Milestone 1`.

### Milestone 2: Tier 2 & 3 Bugs (Ranks 9-18)
- **Objective**: Fix the 10 medium-to-high severity systematic biases in priority order.
- **Verification**: Crate-level unit tests and EKF validations.
- **Delegation**: Spawn a sub-orchestrator `teamwork_preview_orchestrator` with scope `Milestone 2`.

### Milestone 3: Tier 4 Features (Ranks 19-25)
- **Objective**: Implement the 7 geophysical and constellation enhancements. For complex items where full implementation is not feasible within scope constraints, add robust stubs and `#[test] #[ignore]` regression tests.
- **Verification**: Unit and regression tests.
- **Delegation**: Spawn a sub-orchestrator `teamwork_preview_orchestrator` with scope `Milestone 3`.

---

## Verifiability & CI
At the end of each milestone:
- Run `cargo build --workspace` to ensure no warnings or errors.
- Run `cargo test --workspace` to ensure all tests pass.
- Each fix must have a corresponding test matching naming rules (e.g. `test_...`).
