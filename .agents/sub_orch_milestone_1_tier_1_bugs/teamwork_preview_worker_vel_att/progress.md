# Progress

- [x] Run existing tests and check current codebase structure.
- [x] View `crates/gneiss-rtk/src/engine/predictor.rs` and identify the attitude-to-velocity coupling block.
- [x] Design and add `test_transition_matrix_velocity_attitude_coupling` in `crates/gneiss-rtk/src/engine/tests_predictor.rs`.
- [x] Run tests and verify the test fails with the buggy code (Red phase).
- [x] Apply the fix: Negate the attitude-to-velocity coupling block computation.
- [x] Run tests and verify they now pass (Green phase).
- [x] Run `cargo test --workspace` to ensure workspace-wide compliance.
- [ ] Create `changes.md` and `handoff.md`.
- [ ] Send completion message.

Last visited: 2026-06-20T21:41:00-07:00
