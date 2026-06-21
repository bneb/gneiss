# Handoff Report — Bug 15: Incorrect Broadcast Clock TGD Correction Fix

## 1. Observation
- **File**: `crates/gneiss-core/src/ephemeris.rs`
  - The variant-specific `position` methods (lines 477, 509, 541, 579) call `calc_keplerian`, passing constellation-specific group delay parameters (e.g., `self.tgd`, `self.bgd_e1_e5a`, `self.tgd1`).
  - At line 465, `calc_keplerian` subtracts this delay: `let clk_err = af0 + af1 * tc + af2 * tc * tc + F * e * sqrt_a * libm::sin(ek) - tgd;`.
- **File**: `crates/gneiss-rtk/src/engine/ppp.rs`
  - In `compute_sat_state`, the code previously called `eph.position(t_nom)` (line 427) and manually added back `eph.tgd()` (line 434) to undo the single-frequency correction: `dt_s = brdc_clk + eph.tgd();`.
- **File**: `crates/gneiss-rtk/src/estimators/spp.rs`
  - In `compute_sat_state`, the code called `eph.position(t_tx_sat_gps)` (line 185) and `eph.position(t_tx_true_gps)` (line 190) directly without checking if the measurements were iono-free (`m.is_iono_free`), leading to incorrect clock correction for dual-frequency iono-free SPP observations.

## 2. Logic Chain
1. For dual-frequency iono-free combinations, the group delay (TGD/BGD) cancels out in the linear combination, meaning the clock error should not have TGD/BGD subtracted.
2. In order to avoid manual and constellation-dependent TGD addition hacks in estimators, we implement a `position_iono_free` API on the `Ephemeris` enum and variants.
3. This API delegates to `calc_keplerian` with a `0.0` value for the group delay parameter. For GLONASS (which has no group delay), it delegates directly to `position`.
4. In `ppp.rs`, we update `compute_sat_state` to use `position_iono_free` directly and remove the manual `tgd()` addition.
5. In `spp.rs`, we update `compute_sat_state` to conditionally use `position_iono_free` when `m.is_iono_free` is true, and `position` when it is false.
6. A unit test `test_broadcast_clock_tgd_correct` is added to check that `position_iono_free` and `position` differ by exactly `tgd()` / group delay for GPS, Galileo, BeiDou, and QZSS.

## 3. Caveats
- No caveats.

## 4. Conclusion
The broadcast clock TGD correction bug has been fully resolved by implementing `position_iono_free` on the `Ephemeris` enum and variants, and updating the PPP engine and SPP estimator to utilize it correctly. All changes align with the minimal change principle and preserve the original codebase architecture.

## 5. Verification Method
- **Verification Commands**:
  - Run the specific test for TGD correction:
    `cargo test ephemeris::tests::test_broadcast_clock_tgd_correct`
  - Run the full workspace test suite:
    `cargo test`
- **Files to Inspect**:
  - `crates/gneiss-core/src/ephemeris.rs`
  - `crates/gneiss-rtk/src/engine/ppp.rs`
  - `crates/gneiss-rtk/src/estimators/spp.rs`
