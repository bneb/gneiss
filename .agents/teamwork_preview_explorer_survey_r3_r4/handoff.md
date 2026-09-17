# Handoff Report — Survey Frontier R3/R4

## 1. Observation
- **Regression Guard Scripts**:
  - `python3 scripts/check_network_benchmark.py --smoke`:
    ```
    [SMOKE MODE] Evaluating 1800 epochs
      ok     network fused horizontal p50 (m)             0.024 <= 0.04
      ok     network fused horizontal RMS (m)             0.041 <= 0.06
      ok     network fused vertical RMS (m)               0.043 <= 0.08
      ok     P181 smoothed fixed-only p50 (m)             0.022 <= 0.03
      ok     P222 smoothed fixed-only p50 (m)             0.060 <= 0.09
      ok     SLAC smoothed fixed-only p50 (m)             0.112 <= 0.12
      ok     OHLN smoothed fix rate (%)                  94.000 >= 79.0
      ok     P181 smoothed fix rate (%)                  92.100 >= 85.0
      ok     SLAC smoothed fix rate (%)                  70.600 >= 60.0
    ALL CHECKS PASSED
    ```
  - `python3 scripts/check_multignss_benchmark.py --smoke`:
    ```
    binary sha256[:12] = 5dbc8b82ddfe
    [SMOKE MODE] Evaluating 1800 epochs (~900.0 min)
      ok P181 fix rate (%)                    99.00 >= 97.5
      ok P181 h_p95 (mm)                     121.00 <= 145.0
      ok P181 v_p95 (mm)                     199.00 <= 290.0
      ok P225 fix rate (%)                    87.60 >= 71.0
      ok P225 h_p95 (mm)                      89.00 <= 245.0
      ok P225 v_p95 (mm)                     204.00 <= 370.0
      ok P222 fix rate (%)                    99.40 >= 86.0
      ok P222 h_p95 (mm)                     198.00 <= 295.0
      ok P222 v_p95 (mm)                      84.00 <= 135.0
      ok network fused fix rate (%)           99.30 >= 96.5
    ALL CHECKS PASSED
    ```
- **Frontier R3 (Network RTK / VRS)**:
  - `crates/gneiss-rtk/src/post_process/vrs.rs:78-98`: Currently only implements `fit_plane_gradient` (least-squares linear plane). Lacks Delaunay triangulation models for station troposphere and per-satellite ionospheric pierce points.
  - `crates/gneiss-rtk/src/bin/eval_network_ppk.rs`: Currently evaluates independent single-base PPK runs and fuses them with `network::fuse_network_solutions` (lines 693-703). It does not form inter-CORS multi-baseline network adjustments and does not feed synthesized VRS epochs into rover PPK.
  - `datasets/cors_sf_bay_network`: Contains 11 regional CORS base observation files (`cabl`, `capo`, `mhcb`, `ohln`, `p181`, `p222`, `p225`, `p261`, `p271`, `slac`, `tibb`) and rover `p224`.
- **Frontier R4 (Unified Composite Integration / ESKF)**:
  - `crates/gneiss-rtk/src/swfg/imu_preintegration/smoother.rs:14-15`: State is strictly 6-DOF: `pub type State6 = Vector6<f64>;` (pos 3, vel 3). Attitude is dead-reckoned open loop (`self.attitude *= preint.dq;`, line 85). Sensor biases $\mathbf{b}_a, \mathbf{b}_g$ are absent from the state vector and covariance.
  - `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs`: Uses loosely-coupled updates (`filter.update_gnss(pos, vel, ...)` at line 140) and lacks closed-loop error-state feedback.
  - `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`: `PppArSolver` provides wide-lane and narrow-lane integer fixing, but is not yet coupled to an inertial filter.
  - `crates/gneiss-parsers/src/sinex_bia.rs`: `SinexBias` parses and provides satellite phase biases, ready for uncombined carrier phase modeling.

## 2. Logic Chain
1. From Observation 1, the existing baseline passes both regression guards cleanly. Since both guard scripts strictly parse stdout for section headers like `=== Smoothed RTK [station] ===` and `=== NETWORK FUSED ===`, any additions to `eval_network_ppk.rs` must preserve these exact header formats to prevent breaking CI regressions.
2. From Observation 2, `vrs.rs` has the observable shifting logic (`shift_observables`), but planar fitting cannot accurately capture localized atmospheric variations over 15–50 km baselines. Replacing planar gradients with 2D Delaunay triangulation on CORS stations (troposphere) and ionospheric pierce points (ionosphere) provides localized, continuous barycentric interpolation without edge divergence.
3. From Observation 3, the current inertial smoother in `smoother.rs` has 402 lines. Expanding it to 15 states in that file would violate the 500 LOC ceiling in `AGENTS.md`. Therefore, the 15-state ESKF must be structured in a new, dedicated modular submodule (`crates/gneiss-rtk/src/estimators/eskf/`) split across cohesive files: `state.rs`, `predict.rs`, `constraints.rs`, `updates.rs`, and `smoother.rs`.
4. From Observation 4, composing the 15-state ESKF with `PppArSolver` (Tightly-Coupled PPP/INS) and with VRS synthesis (Tightly-Coupled Network RTK/INS) requires dedicated composite modules (`crates/gneiss-rtk/src/composite/tc_ppp.rs` and `tc_rtk.rs`) that ingest raw/double-difference carrier phase innovations and apply the calibrated antenna lever-arm cross-coupling.

## 3. Caveats
- No production source code was modified during this exploration turn (read-only investigation).
- While `eval_network_ppk` and `eval_odaiba_ins` datasets are locally available and verified, running the full 12,398-epoch Odaiba evaluation takes several minutes; smoke tests should be used during iterative development.
- Inter-station CORS ambiguity fixing on real-world datasets requires robust cycle-slip detection on base stations.

## 4. Conclusion
The codebase is in an ideal state to receive Frontier R3 and Frontier R4 implementations. The comprehensive survey report has been generated at:
`/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r3_r4/survey_r3_r4.md`.
The implementation should proceed in four sequential phases:
1. Spatial 2D Delaunay triangulation (`delaunay.rs`).
2. Multi-baseline inter-CORS network adjustment (`network_adj.rs`) and Delaunay VRS synthesis (`vrs.rs`).
3. 15-state ESKF/MEKF with closed-loop quaternion reset, online bias estimation, and 15-state RTS smoother (`estimators/eskf/`).
4. Unified Composite Integration (`composite/tc_ppp.rs` and `composite/tc_rtk.rs`).

## 5. Verification Method
- Code Quality:
  - `cargo clippy --workspace --all-targets -- -D warnings` (must be 0 warnings).
  - `cargo test --workspace` (all tests pass).
- Regression Guards:
  - `python3 scripts/check_network_benchmark.py --smoke`
  - `python3 scripts/check_multignss_benchmark.py --smoke`
- Performance Benchmarks:
  - `cargo run --release --bin eval_network_ppk` (confirm baseline ppm error reduction across P181, P225, P222).
  - `cargo run --release --bin eval_odaiba_ins` (confirm 15-state ESKF achieves $p_{50} < 2.5\text{ m}$ and RMS $< 5.2\text{ m}$).
