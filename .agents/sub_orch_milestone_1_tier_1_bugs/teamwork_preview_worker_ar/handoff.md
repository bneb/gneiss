# Handoff Report — Bug 9: Sequential AR Covariance Mismatch

## 1. Observation
- File to modify: `crates/gneiss-rtk/src/engine/ppp_iekf.rs`.
- In the original implementation of `resolve_cascade_ar`:
  ```rust
  apply_state_vector(state, &x_current, p_current.clone());
  ```
  was called inside the `for (constellation, group_cands) in &const_groups` loop, mutating the `state` in-place before global validation checks were run.
- In `resolve_widelane_ar`:
  ```rust
  let mut d_wl_full = DMatrix::zeros(subset.len(), state.covariance.nrows());
  ```
  directly used `state.covariance` instead of the local updated covariance `p_current`.
- Ran baseline check and workspace tests using `cargo test --workspace` and they all passed:
  ```
  test result: ok. 253 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.13s
  ```

## 2. Logic Chain
- Mutating the `state` inside the per-constellation loop meant intermediate successes were written directly to `state`. If subsequent constellations failed or if the final global position validation check failed (`jump > 20.0` meters), `resolve_cascade_ar` returned `Err`, but the `state` vector and covariance were already partially or fully mutated.
- To prevent this, the updates should be accumulated in local variables `x_current` and `p_current`. `resolve_widelane_ar` needs to operate on `p_current` rather than `state.covariance`, so we pass it as a parameter `p`.
- At the end of `resolve_cascade_ar`, if and only if all validation checks pass, we call `apply_state_vector` once to commit the accumulated updates to `state`.
- The regression test `test_sequential_ar_mismatch_regression` mocks sequential AR succeeding and accumulating a 28.5m jump (above the 20.0m global limit), asserting that `resolve_cascade_ar` returns the correct jump error and leaves the input state/covariance completely unmodified.

## 3. Caveats
- No caveats. The fix strictly adheres to the requested strategy and maintains correct mathematical/logic execution.

## 4. Conclusion
- The sequential covariance mismatch bug is resolved. State/covariance updates during PPP-AR are now fully transactional: either the entire cascade fix succeeds and commits, or it fails with no side-effects.

## 5. Verification Method
### Compilation and test commands
- To build and run only the RTK engine tests:
  ```bash
  cargo test -p gneiss-rtk
  ```
  Output:
  ```
  test result: ok. 254 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.26s
  ```
- To run the entire workspace test suite:
  ```bash
  cargo test --workspace
  ```
  Output:
  ```
  test result: ok. 254 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.21s
  ```

### Files to Inspect
- `crates/gneiss-rtk/src/engine/ppp_iekf.rs`:
  - `resolve_cascade_ar` (lines 130–288)
  - `resolve_widelane_ar` (lines 352–429)
  - `resolve_narrowlane_ar` (lines 449–507)
  - `test_sequential_ar_mismatch_regression` (lines 1888–1996)
