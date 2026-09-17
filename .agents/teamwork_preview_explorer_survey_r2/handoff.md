# Handoff Report — Frontier R2: Integer PPP-AR Engine via SINEX OSB Ingestion

**Agent**: Survey Explorer R2 (`teamwork_preview_explorer_survey_r2`)  
**Parent**: Orchestrator (`1bd6ce81-03bf-4c40-b8b1-3b137333b5e7`)  
**Date**: 2026-09-12  
**Deliverable Document**: `/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r2/survey_r2.md`  

---

## 1. Observation

1. **SINEX BIA Parser & File Presence**:
   - `crates/gneiss-parsers/src/sinex_bia.rs`: Lines 85–196 implement `SinexBias::parse` for `+BIAS/SOLUTION` blocks, storing `BiasRecord`s. Lines 198–210 implement `get_bias` via linear scan over `self.records`.
   - `datasets/rtkexplorer/sample_1/f9p_ppp_1224/com21374.bia` (686,779 bytes) exists and contains CODE MGEX IAR phase and code OSB records for day 359, 2020 (`OSB G063 G01 C1C ... ns -1.6145`, `OSB G063 G01 L1C ... ns -0.05369`, etc.).
   - In `crates/gneiss-rtk/src/bin/eval_ppp.rs` (line 420), `f9p_spec` has `bia_path: None`, leaving the F9P kinematic dataset uncorrected by satellite phase biases.

2. **Observation Equations & Antenna Models**:
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/satpos.rs`: Lines 171–206 implement `compute_phase_centre_3d` projecting ANTEX 3D satellite PCO into ECEF via nominal GNSS yaw-attitude modeling.
   - `crates/gneiss-parsers/src/antex.rs`: Lines 20–25 define `FrequencyPcv` with `pco: Vector3<f64>` and `noazi: Vec<f64>` (nadir-dependent PCV). No method currently interpolates satellite nadir PCV during observation modeling.
   - `crates/gneiss-parsers/src/receiver_pcv/mod.rs`: Implements 1D and 2D receiver PCV interpolation (`interpolate_az_zen`). However, `post_process/mod.rs:127` only defines `ReceiverPcvPair` (paired rover + base for RTK), with no standalone receiver PCV path for rover-only PPP.
   - `crates/gneiss-rtk/src/swfg/pipeline/factors/uduc.rs`: Full Undifferenced Uncombined (UDUC) pseudorange and carrier-phase factors with slant ionosphere states exist, but `swfg/engine/uduc_builder.rs:27` is never invoked by `SwfgEngine::process_impl` (`swfg/engine/mod.rs:313`).

3. **Ambiguity Resolution & LAMBDA**:
   - `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`: Lines 62–162 define `PppArSolver` with `fix_wide_lane`, `fix_sd_wide_lane`, and `fix_narrow_lane` (calling `resolve_lambda`).
   - `PppArSolver` is not called anywhere in the active execution pipeline outside its own unit tests.
   - `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs`: Lines 15–65 execute `execute_ar_step`, which extracts un-differenced float ambiguities in meters and calls `attempt_ar_fix` directly with LAMBDA, which fails to fix integer ambiguities for PPP.

4. **Multi-Constellation Gaps**:
   - `crates/gneiss-rtk/src/swfg/engine/epoch.rs`: Line 364 has `let is_supp = matches!(sat_obs.sat.constellation, Constellation::Gps | Constellation::Galileo); if !is_supp { return None; }`, discarding BeiDou and QZSS.
   - `satpos.rs:93` and `epoch.rs:211`: Omit `'J'` (QZSS) in `sys_char` matching.

5. **Benchmark Performance**:
   - Execution of `cargo run --release --bin eval_ppp` on the F9P kinematic drive yields:
     - CSRS-PPP vs RTK Truth: $p_{50} = \mathbf{0.287\text{ m}}$, RMS $= \mathbf{0.288\text{ m}}$ (33-epoch sample), $0.296\text{ m}$ RMS (full dataset).
     - Gneiss Float PPP vs CSRS-PPP: $p_{50} = \mathbf{10.459\text{ m}}$, RMS $= \mathbf{10.572\text{ m}}$.
     - Gneiss Float PPP vs RTK Truth: $p_{50} = \mathbf{10.077\text{ m}}$, RMS $= \mathbf{10.114\text{ m}}$.

---

## 2. Logic Chain

1. From Observation 1, `com21374.bia` provides the exact phase and code OSBs required for integer PPP-AR, but `eval_ppp.rs` omits it (`bia_path: None`). Consequently, the carrier-phase and pseudorange measurements enter the solver contaminated by satellite hardware and initial phase biases.
2. From Observation 2, satellite nadir PCV and standalone receiver antenna PCO/PCV are omitted in PPP mode, inducing systematic centimeter-to-decimeter modeling errors across varying satellite elevations.
3. From Observation 3, carrier phase ambiguities in un-differenced PPP absorb the receiver clock and initial phase offset, making them non-integer. Because `PppArSolver` is disconnected and `execute_ar_step` runs LAMBDA directly on un-differenced meter ambiguities, no integer fixing occurs. Single-differencing between satellites is mathematically required to eliminate the common receiver phase bias.
4. From Observation 4, BeiDou and QZSS are artificially dropped by `is_supp` in `epoch.rs`, reducing satellite geometry and satellite availability for ambiguity fixing.
5. Synthesizing Steps 1–4 explains Observation 5: without OSB biases, with truncated constellations, and without integer fixing, the Gneiss PPP filter operates as an uncalibrated float estimator dominated by pseudorange code noise, resulting in the observed $\approx 10\text{ m}$ horizontal error.
6. Implementing OSB ingestion, satellite nadir PCV, standalone receiver PCO/PCV, multi-constellation tracking (GPS/Galileo/BDS/QZSS), and single-differenced Wide-Lane/Narrow-Lane LAMBDA integer fixing directly addresses each root cause, enabling sub-meter kinematic accuracy matching CSRS-PPP ($0.296\text{ m}$ RMS).

---

## 3. Caveats

1. The F9P rover antenna in `rover.obs` is a vehicle patch antenna (`ANT # / TYPE` header is blank); its phase center offset relative to the vehicle roof monument requires using default patch antenna offsets (or known calibration) or estimating receiver clock offsets consistently.
2. Product alignment is essential: CODE OSB biases (`com21374.bia`) must be paired with CODE orbit/clock products (`COD0MGXFIN`) to maintain orbit/clock/bias consistency, rather than mixing with ESA products.

---

## 4. Conclusion

Frontier R2 is mathematically well-defined and architecturally ready for implementation. The core algorithms (`SinexBias`, `PppArSolver`, `resolve_lambda`, `UducCarrierPhaseFactor`) exist in the codebase, and the path to achieving sub-meter kinematic PPP-AR vs CSRS-PPP on `eval_ppp.rs` requires:
1. Indexing `SinexBias` and exposing wide-lane/narrow-lane satellite bias helpers.
2. Adding satellite nadir PCV and standalone receiver PCO/PCV to un-differenced observation processing.
3. Enabling BeiDou and QZSS in `epoch.rs`, `satpos.rs`, and band selection.
4. Hooking `PppArCascade` into `ar_handler::execute_ar_step` to execute single-differenced Wide-Lane (MW) rounding and Narrow-Lane LAMBDA integer search.
5. Updating `crates/gneiss-rtk/src/bin/eval_ppp.rs` to ingest `com21374.bia` and evaluate integer PPP-AR against CSRS-PPP.

---

## 5. Verification Method

1. **Compilation & Clippy**:
   ```bash
   cargo build --workspace
   cargo clippy --workspace --all-targets -- -D warnings
   ```
   *Expected*: 0 errors, 0 warnings.
2. **Unit Tests**:
   ```bash
   cargo test -p gneiss-rtk --lib ambiguity::ppp_ar
   cargo test -p gneiss-parsers --lib sinex_bia
   ```
   *Expected*: All tests pass with 0 failures.
3. **Kinematic PPP-AR Benchmark**:
   ```bash
   PPP_ONLY=f9p cargo run --release --bin eval_ppp
   ```
   *Expected*: Resolves integer ambiguities, achieving horizontal error $p_{50} < 1.0\text{ m}$ and RMS $< 1.0\text{ m}$, closing the discrepancy against CSRS-PPP ($0.296\text{ m}$ RMS).
4. **Invalidation Condition**:
   If `eval_ppp` achieves horizontal RMS $\ge 1.0\text{ m}$ or fails to fix integer ambiguities, the integer fixing cascade or product alignment is invalid.
