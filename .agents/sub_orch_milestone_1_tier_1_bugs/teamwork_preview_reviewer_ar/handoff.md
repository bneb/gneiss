# Handoff Report — Bug 9: Sequential AR Covariance Mismatch

## 1. Observation
- **File reviewed**: `/Users/kevin/projects/gneiss/crates/gneiss-rtk/src/engine/ppp_iekf.rs`.
- **Method modified**: `resolve_cascade_ar` (lines 150-336), `resolve_widelane_ar` (lines 402-506), `resolve_narrowlane_ar` (lines 508-622).
- **Unit test verified**: `test_sequential_ar_mismatch_regression` (lines 1856-2023 of `crates/gneiss-rtk/src/engine/ppp_iekf.rs`).
- **Commands executed**:
  - `cargo test -p gneiss-rtk test_sequential_ar_mismatch_regression` -> `test engine::ppp_iekf::mutant_killer_tests::test_sequential_ar_mismatch_regression ... ok`
  - `cargo test -p gneiss-rtk` -> `254 passed; 0 failed`
  - `cargo test --workspace` -> `257 passed; 0 failed`
  - `cargo fmt --check` -> exit code 0 (no formatting errors).
  
We observed that `resolve_cascade_ar` (lines 150-336) now keeps track of state improvements during the per-constellation loop in local variable `x_current` and `p_current` (copied from `state.covariance`), rather than invoking `apply_state_vector` in-place inside the loop or in the fallback branch. The actual EKF state update is deferred until the final global validation check (lines 319-326) succeeds:

```rust
319:         // Global position validation against original float
320:         let float_pos = Vector3::new(x[0], x[1], x[2]);
321:         let fixed_pos = Vector3::new(x_current[0], x_current[1], x_current[2]);
322:         let jump = (fixed_pos - float_pos).norm();
323:         if jump > 20.0 {
324:             tracing::warn!("PPP-AR rejected: 3D position jump {:.2}m > 20m", jump);
325:             return Err("Position jump too large after AR fix");
326:         }
327: 
328:         tracing::info!(
329:             "PPP Cascade AR Fixed! N_Sats: {} jump={:.2}m",
330:             total_fixed_sats,
331:             jump
332:         );
333:         apply_state_vector(state, &x_current, p_current);
334:         state.is_fixed = true;
335:         Ok(())
```

## 2. Logic Chain
1. In the prior implementation of `resolve_cascade_ar`, each iteration of the constellation loop or the fallback update would call `apply_state_vector` directly on `state`, mutating the EKF's active state vector and covariance.
2. If subsequent validation (such as the final global 3D position jump check) failed, the method returned an `Err`.
3. Because the state had already been mutated in-place, the filter was left with a highly constrained covariance matrix (corresponding to resolved integer ambiguities) but remained in a float operating state (`state.is_fixed = false`).
4. This mismatch caused the filter to degrade or diverge in subsequent epochs.
5. In the current implementation, `state` is treated as read-only throughout the sequential loop and fallback steps. Instead, `resolve_widelane_ar` is updated to take the covariance matrix `p` as a parameter and avoid reading from `state.covariance` directly.
6. The state vector is only mutated at the very end of `resolve_cascade_ar` using `apply_state_vector(state, &x_current, p_current)` and setting `state.is_fixed = true` after the global position validation successfully passes.
7. If validation fails, `resolve_cascade_ar` returns an `Err` and `state` remains completely unmodified.

## 3. Caveats
- No caveats. The fix is clean and completely isolates the state mutations. The implementation of `resolve_cascade_ar` in `ppp_ins_iekf.rs` does not suffer from sequential updating as it is a single-step resolution, but we have verified it as well to ensure there are no similar gaps.

## 4. Conclusion
The implementation correctly eliminates the sequential AR covariance mismatch by deferring EKF state mutation until all sequential and fallback steps, including final global validation, are successfully completed.

## 5. Verification Method
To independently verify the fix and regression test, run:
```bash
cargo test -p gneiss-rtk test_sequential_ar_mismatch_regression
cargo test --workspace
```
These tests verify that:
- When a global validation failure occurs, the state vector, covariance, and fixed flag in the filter state remain unmodified.
- Code formatting and workspace layouts conform to project standards.

---

# Quality Review Report

## Review Summary
**Verdict**: APPROVE

## Findings
- **No findings of concern**. The implementation is correct, logically complete, and has high test coverage. No integrity violations (hardcoded tests, dummy functions, or fake validation output) were detected.

## Verified Claims
- Sequential AR updates do not modify state on validation failure -> verified via `test_sequential_ar_mismatch_regression` -> **PASS**
- Overall cargo workspace builds and passes all unit tests -> verified via `cargo test --workspace` -> **PASS**
- Workspace code formatting is compliant -> verified via `cargo fmt --check` -> **PASS**

## Coverage Gaps
- None. `ppp_ins_iekf.rs` was also inspected and does not have the sequential constellation loop.

## Unverified Items
- None.

---

# Adversarial Challenge Report

## Challenge Summary
**Overall risk assessment**: LOW

## Challenges

### [Low] Challenge 1: Local Copy Allocation Overhead
- **Assumption challenged**: Cloning the covariance matrix and state vector multiple times is memory/CPU efficient enough for a real-time system.
- **Attack scenario**: High rate EKF updates under resource pressure.
- **Blast radius**: Increased stack/heap allocations.
- **Mitigation**: The dimensions are extremely small (state dimension is typically <100), meaning the allocations represent negligible CPU/memory load and are fully acceptable for RTK performance.

## Stress Test Results
- Simulated sequential constellation fixes with sub-gate jumps adding up to a global validation failure -> expected to reject fix and preserve float state -> observed to reject and preserve state -> **PASS** (via `test_sequential_ar_mismatch_regression`).

## Unchallenged Areas
- None.
