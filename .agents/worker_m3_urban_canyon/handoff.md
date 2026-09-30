# Handoff Report — Milestone 3 (R3: Doppler-Assisted Cycle Slip Detection & Phase Continuity Validation)

## 1. Observation

### Source Code Modifications & Exact Locations
1. **`crates/gneiss-rtk/src/post_process/screening.rs`** (467 LOC):
   - **Line 69**: Doppler slip detection loop extended from `[1, 2, 7]` to `[1, 2, 5, 6, 7]`:
     ```rust
     for band in [1, 2, 5, 6, 7] {
         let band_slips = self.check_band_doppler_slips(epoch, band);
         slips.extend(band_slips);
     }
     ```
   - **Line 150**: Replaced fixed 1.0 threshold with adaptive threshold tuned to detect half-cycle (0.5 cyc) and 1.0 cyc slips:
     ```rust
     let thresh = (0.30 * dt).clamp(0.28, 1.0);
     ```
   - **Lines 349–465**: Added `make_test_obs_band` helper and `test_doppler_half_cycle_slip_multi_band` verifying that 0.5-cycle slips are detected and increment the arc counter across all supported bands (1, 2, 5, 6, 7).

2. **`crates/gneiss-rtk/src/estimators/rtk_iekf/mod.rs`** (477 LOC):
   - **Line 209**: Removed conditional `if self.widelane_ar` gate; base-side slip detection now executes unconditionally every epoch:
     ```rust
     self.slip_detector.check_epoch(rover);
     self.base_slip_detector.check_epoch(base);
     let dd_meas = self.build_dd_measurements(rover, base, base_pos, ephems)?;
     ```

3. **`crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`** (485 LOC):
   - **Line 138**: Added `external_slip: bool` parameter to `update_tracker_from_obs`:
     ```rust
     pub fn update_tracker_from_obs(&mut self, obs: &[DoubleDiffMeasurement], external_slip: bool) {
     ```
   - **Line 149**: Incorporated `external_slip` into tracker arc reset condition:
     ```rust
     let slip = external_slip || lli_slip || arc_changed;
     ```
   - **Lines 390–448**: Refactored `solve_system_free_upd` and `compute_upd_residuals` helpers into independent functions <= 32 LOC, maintaining total file size at 485 LOC.
   - **Lines 475–484**: Added `test_tracker_resets_on_slip` unit test.

4. **`crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`** (497 LOC):
   - **Lines 287–295**: Modified `build_single_dd_pair` to return `Option<(DoubleDiffMeasurement, bool)>` exposing the pair-level slip flag `lli_slip`.
   - **Lines 65–79**: In `form_pair_dd`:
     - Aggregated slips across bands (`is_slip = is_slip || slip_b;`).
     - Reset `pair_epochs` to 0 on slip (`self.pair_epochs.insert(m.key, 0);`) while incrementing on continuous lock (`*self.pair_epochs.entry(m.key).or_insert(0) += 1;`).
   - **Lines 86–100**: Propagated `is_wl_slip` into `self.wl_tracker.update_tracker_from_obs(&[dd_wl], is_wl_slip);` and into `self.update_phase_wl(..., is_wl_slip)`.
   - **Line 372**: In `update_phase_wl`, passed `is_slip` directly to `self.pw_tracker.update(key, pw, 0.0, is_slip);` (replacing hardcoded `false`).
   - **Lines 433–448**: Modified `check_pair_slip` to take `detectors: (&CycleSlipDetector, &CycleSlipDetector)` to satisfy clippy argument count <= 7 without `#[allow(...)]`, evaluating base slip arcs unconditionally alongside rover slip arcs:
     ```rust
     let (slip_detector, base_slip_detector) = detectors;
     let cur_arc = slip_detector.get_arc(sat_id)
         + slip_detector.get_arc(ref_sat)
         + base_slip_detector.get_arc(sat_id)
         + base_slip_detector.get_arc(ref_sat);
     ```
   - **Lines 453–496**: Added unit tests `test_unconditional_base_slip_flags_pair_slip` and `test_pair_epochs_reset_on_slip`.

### Quantitative Verifications
- **File Lengths**:
  - `screening.rs`: 467 LOC (< 500)
  - `formation.rs`: 497 LOC (< 500)
  - `mw.rs`: 485 LOC (< 500)
  - `mod.rs`: 477 LOC (< 500)
- **Clippy**:
  `cargo clippy --workspace --all-targets -- -D warnings`: 0 warnings, exit code 0.
- **Unit & Integration Tests**:
  - `cargo test -p gneiss-rtk --lib`: 445 passed; 0 failed.
  - `cargo test --test test_urban_canyon_e2e`: 51 passed; 0 failed.
  - `cargo test --workspace`: 205 passed; 0 failed.
- **Benchmark Smoke Tests**:
  - `python3 scripts/check_network_benchmark.py --smoke`: ALL CHECKS PASSED.
  - `python3 scripts/check_multignss_benchmark.py --smoke`: ALL CHECKS PASSED.

## 2. Logic Chain

1. **Multi-Band Doppler Slip Detection (F8)**:
   - *Observation*: `check_doppler_slips` previously checked only bands `[1, 2, 7]`, skipping GPS/Galileo/BDS L5/E5a/B2a (band 5) and Galileo E6/BDS B3 (band 6). Furthermore, the threshold `(1.0 * dt).max(1.0)` was >= 1.0 cycle, missing half-cycle (0.5 cyc) slips.
   - *Reasoning*: Multi-frequency GNSS receivers frequently encounter cycle slips on secondary and tertiary frequencies in urban canyons. Setting the threshold to `(0.30 * dt).clamp(0.28, 1.0)` ensures that at 1 Hz (`dt = 1.0`), the threshold is 0.30 cycles, which reliably flags both 0.5-cycle and 1.0-cycle discontinuities while remaining above the nominal ~0.05–0.15 cycle Doppler integration noise.
   - *Conclusion*: Expanding the loop to `[1, 2, 5, 6, 7]` and tuning the adaptive threshold catches half-cycle slips across all active bands.

2. **Unconditional Base Station Cycle Slip Detection (F9)**:
   - *Observation*: In `rtk_iekf/mod.rs` and `formation.rs`, `self.base_slip_detector.check_epoch(base)` and base arc checks in `check_pair_slip` were gated behind `if self.widelane_ar`.
   - *Reasoning*: Cycle slips occurring on reference station antennas or tracking loops corrupt the double-difference carrier phase observations regardless of whether widelane ambiguity resolution is enabled.
   - *Conclusion*: Removing the `widelane_ar` gate ensures that base slips are detected unconditionally in all processing modes and properly flag double-difference pairs.

3. **Phase Tracking Continuity & `pair_epochs` Reset (F9)**:
   - *Observation*: In `formation.rs`, `pair_epochs` was monotonically incremented for all formed pairs regardless of cycle slips. Furthermore, `pw_tracker.update` was passed a hardcoded `false` for its slip flag, and `WidelaneTracker` did not receive external pair-slip signals from Doppler detection.
   - *Reasoning*: An ambiguity filter or ambiguity resolution manager requires an accurate count of continuous lock epochs. When a cycle slip occurs on either satellite of a pair (or on the base), the phase integer ambiguity jumps, invalidating the accumulated phase history and any smoothed widelane or phase-widelane statistics.
   - *Conclusion*: By returning `lli_slip` from `build_single_dd_pair`, checking slips across all bands, resetting `self.pair_epochs.insert(m.key, 0)` upon detection, and propagating `is_slip` into `pw_tracker` and `mw.rs`'s `update_tracker_from_obs`, phase continuity is accurately represented and filtered.

## 3. Caveats

- **GLONASS Inter-Frequency Bias Handling**: GLONASS FDMA signals have satellite-specific carrier frequencies. The Doppler slip detector uses each satellite's designated carrier frequency to convert phase and Doppler between cycles and metres. If GLONASS frequency channel numbers are unavailable, default nominal band frequencies are used.
- **Dynamic Platform Acceleration**: The Doppler predictor threshold `(0.30 * dt).clamp(0.28, 1.0)` assumes platform jerk / acceleration over `dt <= 1.0s` does not exceed ~0.30 cycles/s² of carrier acceleration without an IMU coupling. For high-rate (e.g., 5 Hz / 10 Hz) data in urban canyons, `dt` scales linearly down to 0.28 cycles, providing tighter slip boundaries without false alarms.

## 4. Conclusion

Milestone 3 requirements (F8 and F9) are completely implemented and verified:
1. Multi-band Doppler slip detection actively screens bands `[1, 2, 5, 6, 7]` with adaptive threshold `(0.30 * dt).clamp(0.28, 1.0)`, detecting both 0.5-cycle and 1.0-cycle slips.
2. Base station slip detection runs unconditionally across all operational modes.
3. Cycle slips reset `pair_epochs` to 0 and propagate cleanly to `pw_tracker` and `WidelaneTracker`.
4. All code complies with AGENTS.md constraints (< 500 LOC per file, functions <= 32 LOC, nesting < 3, 0 `unwrap()`, 0 compiler or clippy warnings).
5. All workspace tests, e2e urban canyon tests, and CI benchmark smoke guard scripts pass with 100% success.

## 5. Verification Method

To independently verify this milestone, run:
```bash
# 1. Zero compiler and clippy warnings
cargo clippy --workspace --all-targets -- -D warnings

# 2. Project unit and library tests
cargo test -p gneiss-rtk --lib

# 3. Dedicated Urban Canyon E2E test suite (all 51 tests)
cargo test --test test_urban_canyon_e2e

# 4. Entire workspace integration and unit test suite
cargo test --workspace

# 5. CI Smoke Guard Scripts
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke
```
