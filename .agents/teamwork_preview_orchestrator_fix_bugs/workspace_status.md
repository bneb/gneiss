# Workspace Status Report

Generated at: 2026-06-21T15:05:00Z
Cwd: `/Users/kevin/projects/gneiss`

## 1. Git Status and Diff Summary

### Git Status
```
On branch main
Your branch is ahead of 'origin/main' by 74 commits.

Changes not staged for commit:
	modified:   .agents/sub_orch_milestone_1_tier_1_bugs/BRIEFING.md
	modified:   .agents/sub_orch_milestone_1_tier_1_bugs/progress.md
	modified:   .agents/sub_orch_milestone_1b_tier_1_bugs/BRIEFING.md
	modified:   .agents/sub_orch_milestone_1b_tier_1_bugs/progress.md
	modified:   .agents/teamwork_preview_orchestrator_fix_bugs/progress.md
	modified:   crates/gneiss-rtk/src/engine/ppp_iekf.rs
	modified:   crates/gneiss-rtk/src/engine/processor/mod.rs

Untracked files:
	.agents/teamwork_preview_reviewer_2_bug_15_rep/
	.agents/teamwork_preview_worker_windup_fixed_retry/
```

### Git Diff Summary
The changes are grouped into two source code modifications and metadata updates:

1. **`crates/gneiss-rtk/src/engine/ppp_iekf.rs`**
   - Introduces a small regularization matrix `r_wl` of `1e-6` m² per pair to prevent `p_wl` (wide-lane constrained covariance) from collapsing to rank-singular when applying the Joseph covariance update with zero measurement noise.
   - Constrains the downstream narrow-lane (NL) LAMBDA search using the regularized wide-lane (WL) covariance matrix.

2. **`crates/gneiss-rtk/src/engine/processor/mod.rs`**
   - Adds `consecutive_rejections` count tracker to the `ProcessingEngine`.
   - Modifies the epoch processing error handler: resets EKF state (`self.current_state = None`) if `EngineError::InsufficientSatellites` occurs for `3` consecutive epochs. Non-satellite errors or successful epochs reset the rejection counter to `0`.

3. **Metadata Updates under `.agents/`**
   - Keeps milestone and sub-orchestrator progress/briefings updated for Milestone 1 / Milestone 1b.

---

## 2. Output and Results of Cargo Test

The test suite was run workspace-wide using `cargo test --workspace`.

### Summary
* **Total Crates tested**: `gneiss-rtk`, `gneiss-tests`, `gneiss-core`, `gneiss-fetch`, `gneiss-geodesy`, `gneiss-ntrip`, `gneiss-parsers`
* **Test execution outcome**: `SUCCESS`
* **Crate-specific counts**:
  - `gneiss_rtk`: **257 passed**; 0 failed; 1 ignored
  - `gneiss_tests`: **3 passed**; 0 failed; 0 ignored
  - All other crates (`gneiss_core`, `gneiss_fetch`, `gneiss_geodesy`, `gneiss_ntrip`, `gneiss_parsers`): 0 tests (doc-tests: 0 passed; 0 failed; 0 ignored)

### Detailed Cargo Test Output (Relevant portions)
```
running 258 tests
...
test result: ok. 257 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.21s

     Running unittests src/bin/test_chol.rs (target/debug/deps/test_chol-f9d7656c7adc879a)
running 0 tests
test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s

     Running unittests src/bin/test_nalgebra.rs (target/debug/deps/test_nalgebra-be1e1072b98b4e77)
running 0 tests
test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s

     Running unittests src/bin/test_svd.rs (target/debug/deps/test_svd-7c8dc4047811aa66)
running 0 tests
test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s

     Running unittests src/bin/train_gnn_raim.rs (target/debug/deps/train_gnn_raim-107f63ef9120582a)
running 0 tests
test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s

     Running unittests src/lib.rs (target/debug/deps/gneiss_tests-515d94ead4e464c7)
running 3 tests
test integration::test_cross_crate_fusion_initialization ... ok
test ppp_integration::test_ppp_skeleton ... ok
test urbannav_integration::test_urbannav_tst_replay_skeleton ... ok
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
```

---

## 3. Failing Tests
* **Failing tests**: **None** (0 failed). All 260 workspace tests passed successfully.
