# Plan — Fix Bug 9: Sequential AR Covariance Mismatch

## Goal
Ensure `resolve_cascade_ar` operates on local copies of state/covariance during sequential and fallback updates, and applies the state vector update only at the very end after all validations (including global position validation) pass.

## Detailed Steps
1. **Modify `resolve_widelane_ar` signature & body**:
   - Signature: `fn resolve_widelane_ar(&self, state: &RtkState, p: &DMatrix<f64>, subset: &..., x: &DVector<f64>) -> Result<...>`
   - Body: replace all usages of `state.covariance` with `p`.
2. **Modify `resolve_narrowlane_ar` body**:
   - Replace `state.covariance.nrows()` with `p_wl.nrows()`.
3. **Modify `resolve_cascade_ar` logic**:
   - Remove the `apply_state_vector(state, &x_current, p_current.clone());` call within the loop.
   - Pass `&p_current` into the loop call to `resolve_widelane_ar`.
   - Pass `&p_current` into the fallback call to `resolve_widelane_ar`.
   - Apply `apply_state_vector(state, &x_current, p_current)` only at the end after global validation passes.
4. **Build and verify existing tests**:
   - Run `cargo test -p gneiss-rtk` to verify it compiles and runs.
5. **Implement regression unit test**:
   - Write a unit test that mocks sequential AR succeeding but final global validation failing (position jump > 20m).
   - Verify `resolve_cascade_ar` returns `Err` and `state` remains unchanged (both state vector and covariance).
6. **Verify the new test and linting**:
   - Run `cargo test --workspace` and check code style with `cargo fmt`.
