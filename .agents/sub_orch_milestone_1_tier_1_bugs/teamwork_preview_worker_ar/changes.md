# Changes for Bug 9: Sequential AR Covariance Mismatch

## Fix Description
- Modified the signature of `resolve_widelane_ar` in `crates/gneiss-rtk/src/engine/ppp_iekf.rs` to accept `p: &DMatrix<f64>` as a parameter instead of directly using `state.covariance`.
- Replaced all references to `state.covariance` with the parameter `p` in `resolve_widelane_ar`.
- In `resolve_narrowlane_ar`, replaced `state.covariance.nrows()` with `p_wl.nrows()` since `p_wl` represents the correct covariance dimension.
- In `resolve_cascade_ar`:
  - Removed the immediate/intermediate `apply_state_vector` call from the sequential per-constellation AR loop (`for (constellation, group_cands) in &const_groups`).
  - Updated calls to `resolve_widelane_ar` inside the loop and fallback block to pass `&p_current`.
  - Applied the final state vector update via `apply_state_vector` only at the end of the method after all validations (including the global position validation) pass.
  - This ensures that if sequential AR succeeds but the global position validation fails (producing a jump > 20m), the state is not mutated and remains completely unchanged.

## Regression Test
- Added `test_sequential_ar_mismatch_regression` to the `mutant_killer_tests` module in `crates/gneiss-rtk/src/engine/ppp_iekf.rs`.
- Set up `ArMock` to mock sequential wide-lane and narrow-lane AR succeeding with a cumulative position jump of 28.5 meters (> 20 meters limit) across 3 constellations.
- Checked that `resolve_cascade_ar` returns `Err("Position jump too large after AR fix")` and verified that the input `state`'s state vector, covariance, and `is_fixed` flag remain completely unchanged.
