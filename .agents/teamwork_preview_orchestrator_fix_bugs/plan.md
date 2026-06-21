# Plan: Gneiss Navigation Engine Bug Fixes

This plan outlines the steps for executing the remaining bug fixes across the milestones.

## Execution Strategy

We will use a multi-agent orchestration pattern to delegate the implementation of these fixes. Because there are 25 bugs, we decompose them into sequential sub-orchestrations to manage context size, limit spawn counts per agent, and run specialized implementation and review steps.

### Milestone 1: Tier 1 Bugs (Ranks 1-8)
- **Objective**: Fix the 8 most critical, high-impact bugs.
- **Completed**:
  - Bug 17: GLONASS Time Scale Discrepancy
  - Bug 1: Melbourne-Wübbena Dimensional Typo
  - Bug 9: Sequential AR Covariance Mismatch
  - Bug 2: Velocity-Attitude Transition Sign Mismatch (with positive sign on disk)
- **In-Progress**:
  - Bug 18: Opposite Sign in Phase Wind-Up Correction
  - Bug 15: Incorrect Broadcast Clock TGD Correction
  - Bug 24: Outlier Tolerance in Precise Clock Gaps
  - Bug 6: GMF Legendre Unnormalized Polynomials
- **Delegation**: Spawn a sub-orchestrator (`self`) with scope `Milestone 1b (Bugs 18, 15, 24, 6)`.

### Milestone 2: Tier 2 & 3 Bugs (Ranks 9-18)
- **Objective**: Fix the 10 medium-to-high severity systematic biases in priority order.
- **Delegation**: Spawn a sub-orchestrator (`self`) with scope `Milestone 2`.

### Milestone 3: Tier 4 Features (Ranks 19-25)
- **Objective**: Implement the 7 geophysical and constellation enhancements. For complex items where full implementation is not feasible within scope constraints, add robust stubs and `#[test] #[ignore]` regression tests.
- **Delegation**: Spawn a sub-orchestrator (`self`) with scope `Milestone 3`.

---

## Verifiability & CI
At the end of each milestone:
- Run `cargo build --workspace` to ensure no warnings or errors.
- Run `cargo test --workspace` to ensure all tests pass.
- Each fix must have a corresponding regression test matching the naming rules that passes on the fix but would fail if the bug were restored.
