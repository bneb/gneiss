# Change Log - Phase Wind-up Sign Fix (Bug 18)

## Modified Files

### `crates/gneiss-rtk/src/engine/measurement.rs`
- **Change**: In `apply_windup_to_obs` (lines 42-49), changed `+=` to `-=` to subtract the wind-up correction instead of adding it.
- **Change**: Added `test_phase_windup_correction_sign_rtk` unit test under `mod tests` block.
- **Rationale**: Phase wind-up correction rotates the effective phase by `wup` cycles. Correcting the carrier phase observation requires subtracting the wind-up correction. Adding it doubles the error.

### `crates/gneiss-rtk/src/engine/ppp_iekf.rs`
- **Change**: In `push_cp_measurement` (line 904), changed `+ windup` to `- windup`.
- **Change**: In UDUC carrier phase residual calculations (lines 1032 and 1049), changed `+ windup` to `- windup` to correctly subtract the phase wind-up correction.
- **Rationale**: Align prediction and residual models in the PPP iterated EKF to subtract the phase wind-up correction, matching the physical definition.

### `crates/gneiss-rtk/src/engine/ppp.rs`
- **Verification**: Verified that `wup` subtraction is already correctly implemented (`- wup` and `-=`) throughout the file. No changes were needed here.

## Verification & Build Status
- **Test execution**: Run `cargo test -p gneiss-rtk` and `cargo test --workspace`. All 258 tests passed successfully.
- **Regression verification**: Reverting the sign in `apply_windup_to_obs` to addition (`+=`) causes `test_phase_windup_correction_sign_rtk` to fail as expected, confirming the test is sensitive to the bug.
- **Code style**: Run `cargo fmt` and verified formatting with `cargo fmt --check` (exit status: 0).
