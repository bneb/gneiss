# Progress — Worker M1 (15-State ESKF GNSS/INS)

Last visited: 2026-09-12T16:54:30Z

## Status
- Implemented full 15-state ESKF module under `crates/gneiss-rtk/src/estimators/eskf/`:
  - `types.rs`: `EskfState`, `EskfSnapshot`, `EngineError`, 15-state aliases, error-state layout.
  - `predict.rs`: Mechanization, `compute_transition_matrix` (strictly positive `+f_e_skew * dt`), process noise `Q`, `predict` and `predict_preintegrated`.
  - `update.rs`: GNSS pos/vel updates, lever-arm Jacobian, error-quaternion reset `q <- dq * q`, bias injection, Joseph-form covariance update.
  - `constraints.rs`: Coupled NHC with attitude Jacobian `E_23 R^T [v^e x]` for heading observability, and ZUPT.
  - `smoother.rs`: Full 15-state backward RTS smoother propagating corrections to pos, vel, att, ba, bg across all epochs.
  - `mod.rs`: Module exports and re-exports.
- Integrated 15-state ESKF and RTS smoother into `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs`.
- Unit tests: 17/17 passing.
- Clippy: 0 warnings on `gneiss-rtk` library and binary targets.
- Currently executing benchmark `cargo run --release --bin eval_odaiba_ins` (task-167).
