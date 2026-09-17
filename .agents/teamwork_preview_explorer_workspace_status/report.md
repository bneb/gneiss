# Comprehensive Workspace Health & AGENTS.md Compliance Report

**Date**: 2026-09-13T02:05:40Z  
**Investigator**: Explorer Workspace (`teamwork_preview_explorer_workspace_status`)  
**Repository**: `/Users/kevin/projects/gneiss`  
**Working Mode**: Read-Only Analysis  

---

## Executive Summary

| Category | Status | Details |
|---|:---:|---|
| **Workspace Test Suite (`cargo test --workspace`)** | **PASS** | **1,056 passed**, 0 failed, 1 ignored (annotated). All crate unit, integration, and doc tests pass. |
| **Workspace Clippy (`cargo clippy --workspace --all-targets -- -D warnings`)** | **FAIL (1 issue)** | **1 failure**: `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs:478:27` (`clippy::useless_vec`). Zero warnings in all other crates! |
| **Network RTK Benchmark (`check_network_benchmark.py --smoke`)** | **PASS** | All 9 performance and fix-rate regression metrics within thresholds. |
| **Multi-GNSS Benchmark (`check_multignss_benchmark.py --smoke`)** | **PASS** | All 10 constellation and baseline fix-rate regression metrics within thresholds. |
| **Frontiers E2E Suite (`test_frontiers_e2e`)** | **PASS** | **205 passed**, 0 failed, 0 ignored across Tiers 1–4. |
| **Zero `unwrap()` in Production Code** | **PASS** | **0 unwrap() calls** found across all 262 active workspace production files. |
| **Empty `#[ignore]` Tests** | **PASS** | 0 empty ignores. Sole ignored test (`sources::noaa::test_fetch_station_coordinate`) explicitly documented ("Requires live NOAA CORS network connection"). |
| **Frontier Files AGENTS.md Compliance** | **PARTIAL** | ESKF and Spatial modules 100% compliant (<32 LOC fns, <3 nesting depth, <500 LOC). `ppp_ar.rs` is 506 LOC and has 2 functions >= 32 LOC; minor nesting and function size debt in `tc_rtk.rs`, `tc_ppp.rs`, `network_adj.rs`, and `vrs.rs`. |

---

## 1. Full Workspace Tests (`cargo test --workspace`)

**Result**: **PASS** (Exit Code: 0)  
**Execution Time**: ~3 minutes (all crates + heavy convergence tests)  
**Breakdown by Target**:

- `gneiss_core`: **133 passed**, 0 failed, 0 ignored
- `gneiss_fetch`: **3 passed**, 0 failed, **1 ignored** (justified live network test)
- `gneiss_geodesy`: **20 passed**, 0 failed, 0 ignored
- `gneiss_parsers`: **267 passed**, 0 failed, 0 ignored
- `gneiss_parsers` (integration tests): **8 passed**, 0 failed, 0 ignored
- `gneiss_ntrip`: **1 passed**, 0 failed, 0 ignored
- `gneiss_rtk`: **385 passed**, 0 failed, 0 ignored
- `gneiss_rtk` (receiver_pcv_contract): **4 passed**, 0 failed, 0 ignored
- `gneiss_rtk` (benchmark_matrix): **15 passed**, 0 failed, 0 ignored
- `gneiss_cli`: **17 passed**, 0 failed, 0 ignored
- `gneiss_cli` (tests): **2 passed**, 0 failed, 0 ignored
- `gneiss_tests` (`test_frontiers_e2e`): **205 passed**, 0 failed, 0 ignored
- Doc tests across all crates: **0 failed**
- **Total**: **1,056 passed**, 0 failed, 1 ignored

### Ignored Test Audit:
- File: `crates/gneiss-fetch/src/sources/noaa.rs:218`
- Test: `test_fetch_station_coordinate`
- Attribute: `#[ignore = "Requires live NOAA CORS network connection"]`
- Status: **Compliant with AGENTS.md** (has explicit explanatory comment and requirement for external live NOAA connectivity).

---

## 2. Workspace Clippy (`cargo clippy --workspace --all-targets -- -D warnings`)

**Result**: **FAIL** (Exit Code: 101)  
**Issue**: Exactly 1 compiler/clippy error in the entire repository:

```text
error: useless use of `vec!`
   --> crates/gneiss-rtk/src/ambiguity/ppp_ar.rs:478:27
    |
478 |         let nl_integers = vec![12, -8, 15, 7];
    |                           ^^^^^^^^^^^^^^^^^^^ help: you can use an array directly: `[12, -8, 15, 7]`
    |
    = help: for further information visit https://rust-lang.github.io/rust-clippy/rust-1.97.0/index.html#useless_vec
    = note: `-D clippy::useless-vec` implied by `-D warnings`
```

### Proposed Fix:
In `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs:478`:
```rust
// Before:
let nl_integers = vec![12, -8, 15, 7];
// After:
let nl_integers = [12, -8, 15, 7];
```
*(All other workspace targets, including `gneiss-core`, `gneiss-parsers`, `gneiss-fetch`, `gneiss-geodesy`, `gneiss-ntrip`, `gneiss-cli`, and `gneiss-tests`, pass `cargo clippy` with zero warnings!)*

---

## 3. Regression Guard Scripts

### 3.1 Network RTK Benchmark (`python3 scripts/check_network_benchmark.py --smoke`)
**Result**: **PASS** (Exit Code: 0)  
**Output Summary**:
```text
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

### 3.2 Multi-GNSS Benchmark (`python3 scripts/check_multignss_benchmark.py --smoke`)
**Result**: **PASS** (Exit Code: 0)  
**Output Summary**:
```text
binary sha256[:12] = f1a65ce39c45
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

---

## 4. Frontiers E2E Test Suite (`cargo test -p gneiss-tests --test test_frontiers_e2e`)

**Result**: **PASS** (Exit Code: 0)  
**Test Results**: **205 passed; 0 failed; 0 ignored; finished in 0.01s**  
**Coverage across Tiers 1–4**:
- `tier1_eskf.rs`: 30 tests covering 15-state transition matrix positive `vel_att`, error quaternion feedback, online biases, RTS smoother, coupled NHC & ZUPT, Odaiba metrics.
- `tier1_ppp.rs`: 20 tests covering SINEX OSB $O(1)$ lookups, satellite nadir PCV, multi-constellation support, SD LAMBDA AR, kinematic PPP.
- `tier1_vrs.rs`: 25 tests covering CORS ingestion, DD network adjustment, 2D Delaunay IPP/ZWD interpolation, VRS epoch synthesis, Leica specs.
- `tier1_composite.rs`: 15 tests covering Tightly-Coupled PPP/INS and TC Network RTK/INS pipelines and mode switches.
- `tier2_eskf.rs`: 30 tests verifying boundary conditions, numerical stability, saturation limits, lever arms, outlier rejection.
- `tier2_ppp.rs`: 25 tests verifying empty tables, unknown sats, limb nadir angles, high-dimension LAMBDA, turn dynamics.
- `tier2_vrs.rs`: 25 tests verifying 3-station network, single base fallback, collinear geometry guards, storm front gradients.
- `tier2_composite.rs`: 20 tests verifying zero IMU fallback, 30s PPP outages, wheel slip, stationary ZUPT clamping.
- `tier3_pairwise.rs`: 10 tests verifying pairwise interactions between sub-engines.
- `tier4_scenarios.rs`: 5 tests verifying full end-to-end mission workflows (Odaiba 10Hz INS, F9P kinematic PPP-AR, regional CORS VRS, aerial survey, failover).

---

## 5. AGENTS.md Code Standards Audit

### 5.1 Zero `unwrap()` in Production Code
- **Status**: **100% PASS**
- Scanned 262 active workspace Rust files.
- **0 `unwrap()` calls** exist in production code paths. All production code employs `?`, `match`, `if let`, `unwrap_or`, or descriptive `.expect("invariant: ...")`. All occurrences of `.unwrap()` are isolated to unit/integration tests and test harnesses (`tests.rs`, `#[test]`, `#[cfg(test)]`).

### 5.2 File Size Standard (< 500 LOC)
- **Status**: **8 files in active workspace >= 500 LOC**

| File | LOC | Category | Note |
|---|---|---|---|
| `crates/gneiss-rtk/src/estimators/spp/tests.rs` | 1539 | Test | Pre-existing SPP unit test suite |
| `crates/gneiss-parsers/src/rinex/obs/tests.rs` | 884 | Test | Pre-existing RINEX test suite |
| `crates/gneiss-parsers/src/rtcm3/msm/tests.rs` | 861 | Test | Pre-existing RTCM3 MSM test suite |
| `crates/gneiss-rtk/src/bin/eval_network_ppk.rs` | 805 | Harness | Benchmark executable |
| `crates/gneiss-core/src/ephemeris/tests.rs` | 790 | Test | Pre-existing ephemeris test suite |
| `crates/gneiss-rtk/src/estimators/rtk_iekf/tests.rs` | 638 | Test | Pre-existing IEKF test suite |
| `crates/gneiss-parsers/src/rinex/nav/tests.rs` | 510 | Test | Pre-existing NAV test suite |
| **`crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`** | **506** | **Frontier Production** | **Exceeds 500 LOC by 6 lines due to embedded unit tests** |

**Frontier File LOC Breakdown**:
- `crates/gneiss-rtk/src/estimators/eskf/mod.rs`: 33 LOC (< 500)
- `crates/gneiss-rtk/src/estimators/eskf/types.rs`: 181 LOC (< 500)
- `crates/gneiss-rtk/src/estimators/eskf/predict.rs`: 240 LOC (< 500)
- `crates/gneiss-rtk/src/estimators/eskf/update.rs`: 130 LOC (< 500)
- `crates/gneiss-rtk/src/estimators/eskf/constraints.rs`: 160 LOC (< 500)
- `crates/gneiss-rtk/src/estimators/eskf/smoother.rs`: 188 LOC (< 500)
- `crates/gneiss-rtk/src/spatial/delaunay.rs`: 416 LOC (< 500)
- `crates/gneiss-rtk/src/post_process/network_adj.rs`: 408 LOC (< 500)
- `crates/gneiss-rtk/src/post_process/vrs.rs`: 417 LOC (< 500)
- `crates/gneiss-rtk/src/composite/tc_ppp.rs`: 496 LOC (< 500)
- `crates/gneiss-rtk/src/composite/tc_rtk.rs`: 492 LOC (< 500)
- **`crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`**: **506 LOC** (⚠️ > 500 LOC)
- `crates/gneiss-parsers/src/sinex_bia.rs`: 370 LOC (< 500)
- `tests/tests/test_frontiers_e2e/*.rs`: All 11 files between 112 and 304 LOC (< 500)

### 5.3 Function Size Standard (< 32 LOC)
- **Status in New Modules (ESKF, Spatial, SinexBia)**: **100% PASS** (Zero functions >= 32 LOC).
- **Status in Remaining Frontier Files**:

| File | Line | Function | LOC | Recommendation |
|---|---|---|---|---|
| `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs` | 232 | `resolve_sd_with_fixed_wl` | 59 | Extract SD matrix setup and ratio evaluation into helper functions. |
| `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs` | 93 | `fix_narrow_lane` | 33 | Split LAMBDA problem construction from post-fix validation. |
| `crates/gneiss-rtk/src/composite/tc_rtk.rs` | 206 | `formulate_dd_measurements` | 52 | Extract double-difference Jacobian formation into standalone function. |
| `crates/gneiss-rtk/src/composite/tc_rtk.rs` | 109 | `process_epoch` | 33 | Delegate IMU propagation and measurement dispatch to helpers. |
| `crates/gneiss-rtk/src/composite/tc_ppp.rs` | 203 | `formulate_measurements` | 35 | Decompose phase and code innovation equations. |
| `crates/gneiss-rtk/src/post_process/network_adj.rs` | 82 | `adjust_epoch` | 42 | Extract normal equations accumulation into a helper. |

*(Legacy code across the entire workspace contains 198 pre-existing functions >= 32 LOC, primarily in CLI evaluation harnesses and older solvers).*

### 5.4 Nesting Depth Standard (< 3 Levels)
- **Status in New Modules (ESKF, Spatial)**: **100% PASS** (Zero nesting depth >= 3).
- **Frontier Code Nesting Locations**:
  - `network_adj.rs`: lines 216, 219 (`find_common_sats`) and 300–306 (`extract_l1_l2`) have depth 3–4 loops/matches over station observables.
  - `vrs.rs`: lines 378–386 (`shift_observables`) have nested matches on observation types and carrier frequencies.
  - `tc_rtk.rs`: lines 249, 316–317 have nested iterations over base and rover satellites.
  - `ppp_ar.rs`: lines 262–277, 348–355 have nested loops over satellite candidates and conditional debug printouts.

---

## 6. Actionable Next Steps for Engineering Team

1. **Fix Clippy Failure**:
   - Change `vec![12, -8, 15, 7]` to `[12, -8, 15, 7]` in `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs:478`.
2. **Bring `ppp_ar.rs` under 500 LOC**:
   - Decompose `resolve_sd_with_fixed_wl` (59 LOC) and move test module helpers or extract tests into `ppp_ar/tests.rs` to reduce file size from 506 LOC to ~360 LOC.
3. **Refactor 6 Frontier Functions exceeding 32 LOC**:
   - Break down `resolve_sd_with_fixed_wl` (59 LOC), `formulate_dd_measurements` (52 LOC), `adjust_epoch` (42 LOC), `formulate_measurements` (35 LOC), `fix_narrow_lane` (33 LOC), `process_epoch` (33 LOC).
4. **Flatten Depth >= 3 Blocks**:
   - Refactor observable iterators in `network_adj.rs`, `vrs.rs`, and `tc_rtk.rs` using functional iterator chains (`find`, `filter_map`) rather than nested `for` / `if let` blocks.
