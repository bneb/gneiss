# Handoff Report — Bug 9: Sequential AR Covariance Mismatch

## Observation
In `/Users/kevin/projects/gneiss/crates/gneiss-rtk/src/engine/ppp_iekf.rs`, the method `resolve_cascade_ar` (lines 130–285) processes satellite ambiguity resolution per-constellation in a sequential loop:

```rust
164:         let mut x_current = x.clone();
165:         let mut p_current = state.covariance.clone();
...
169:         for (constellation, group_cands) in &const_groups {
...
191:             let wl_result = self.resolve_widelane_ar(state, &subset, &x_current);
...
205:             let nl_result = self.resolve_narrowlane_ar(state, &subset, &keep_indices, &x_wl, &p_wl);
...
230:             x_current = x_fixed;
231:             p_current = p_fixed;
232:             // Update state in-place so subsequent constellations use
233:             // the constrained covariance from this fix
234:             apply_state_vector(state, &x_current, p_current.clone());
...
```

At the end of the function, a global position validation check is performed against the original float position:

```rust
272:         // Global position validation against original float
273:         let float_pos = Vector3::new(x[0], x[1], x[2]);
274:         let fixed_pos = Vector3::new(x_current[0], x_current[1], x_current[2]);
275:         let jump = (fixed_pos - float_pos).norm();
276:         if jump > 20.0 {
277:             tracing::warn!("PPP-AR rejected: 3D position jump {:.2}m > 20m", jump);
278:             return Err("Position jump too large after AR fix");
279:         }
```

In `resolve_widelane_ar` (lines 342–427), the method signature takes `state: &RtkState` but does not take the covariance matrix `p` as an argument. Instead, it reads `state.covariance` directly:

```rust
359:         let q_wl_full = &d_wl_full * &state.covariance * d_wl_full.transpose();
...
398:         let mut q_wl = &d_wl * &state.covariance * d_wl.transpose();
...
415:         let k_wl = &state.covariance * d_wl.transpose() * s_inv;
416:         let dx_wl = &k_wl * (res_wl.best_integers - a_wl);
417:         Ok((
418:             x + dx_wl,
419:             crate::math::covariance::apply_joseph_covariance_update(
420:                 &state.covariance,
421:                 &k_wl,
422:                 &d_wl,
423:                 &DMatrix::zeros(keep_indices.len(), keep_indices.len()),
424:             ),
...
```

## Logic Chain
1. During the execution of `resolve_cascade_ar`, each constellation loop iteration executes `resolve_widelane_ar` and `resolve_narrowlane_ar`.
2. To allow subsequent iterations to use the updated (constrained) covariance matrix from previous successful iterations, `apply_state_vector` is called inside the loop (line 234). This immediately updates `state.covariance` in the mutable `state` reference.
3. If a subsequent constellation fails its validation, or if the final global position validation check (line 276) fails (e.g. 3D position jump > 20m), the function terminates early and returns `Err`.
4. Because `apply_state_vector` was already called inside the loop, the `state` object (both its internal state vector and `state.covariance`) remains modified with the constrained values from the successful constellation updates.
5. However, since the method returned `Err`, the caller (the `solve` loop at line 55) assumes ambiguity resolution failed and treats the state as a float solution (setting/retaining `state.is_fixed = false`).
6. This results in a critical mismatch: the filter covariance matrix is now highly constrained (having very small variances for resolved ambiguities), but the filter operates as a float solution, causing severe filter degradation or divergence in subsequent epochs.

## Caveats
1. No caveats. The problem is isolated to the in-place state mutation side-effects of `resolve_cascade_ar` prior to final acceptance and global validation of the cascade ambiguity resolution step.
2. The implementation of `resolve_cascade_ar` in `ppp_ins_iekf.rs` does not suffer from this issue because it does not run a per-constellation loop (it resolves all constellations in a single step and mutates state only at the end).

## Conclusion
The sequential AR covariance mismatch is caused by in-place mutation of the filter state (via `apply_state_vector`) during the sequential constellation loop in `resolve_cascade_ar`, prior to verifying global position validation and final acceptance of the fixed solution. If the cascade AR is ultimately rejected, these mutations leak, leaving the filter with a contaminated, over-constrained covariance matrix while remaining in float mode.

To resolve this issue without implementing it, we propose the following precise fix strategy:
1. Modify `resolve_widelane_ar` to accept the covariance matrix `p: &DMatrix<f64>` as a parameter:
   ```rust
   fn resolve_widelane_ar(
       &self,
       state: &RtkState,
       p: &DMatrix<f64>,
       subset: &[(
           (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
           (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
       )],
       x: &DVector<f64>,
   ) -> Result<(DVector<f64>, DMatrix<f64>, Vec<usize>), &'static str>
   ```
   Inside `resolve_widelane_ar`, replace all references to `state.covariance` with `p`.
2. In `resolve_narrowlane_ar`, replace `state.covariance.nrows()` with `p_wl.nrows()`. This eliminates the need for `state.covariance` inside `resolve_narrowlane_ar` entirely.
3. In `resolve_cascade_ar`:
   - Remove the in-place state updates (`apply_state_vector(state, &x_current, p_current.clone());`) inside the per-constellation loop.
   - Update the loop calls to pass `&p_current` to `resolve_widelane_ar`.
   - Update the fallback block to also pass `&p_current` to `resolve_widelane_ar` (rather than `state`).
   - Call `apply_state_vector` only at the very end of `resolve_cascade_ar` after all validations (including global position validation) have successfully passed:
     ```rust
     tracing::info!("PPP Cascade AR Fixed! N_Sats: {} jump={:.2}m", total_fixed_sats, jump);
     apply_state_vector(state, &x_current, p_current);
     state.is_fixed = true;
     ```

This ensures that the RTK state is treated as read-only throughout the cascade resolution and is only mutated once the entire solution is verified.

## Verification Method
1. Run `cargo test -p gneiss-rtk` to verify all tests compile and pass.
2. A new unit test should be added to `ppp_iekf.rs` that mocks a scenario where per-constellation AR succeeds (e.g. updating `x_current` and `p_current`), but final global validation fails. The test should verify that:
   - `resolve_cascade_ar` returns an `Err`.
   - The input `state`'s covariance and state vector remain completely unchanged (equal to their initial float values).
