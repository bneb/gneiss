# Handoff Report — Bug 15: Incorrect Broadcast Clock TGD Correction

## 1. Observation
- In `crates/gneiss-core/src/ephemeris.rs`:
  - `calc_keplerian` (line 465) subtracts `tgd` directly from the clock error:
    ```rust
    let clk_err = af0 + af1 * tc + af2 * tc * tc + F * e * sqrt_a * libm::sin(ek) - tgd;
    ```
  - The `position` implementations for `GpsEphemeris`, `GalileoEphemeris`, `BeidouEphemeris`, and `QzssEphemeris` pass their respective group delays to `calc_keplerian` (lines 500, 532, 570, 602).
- In `crates/gneiss-rtk/src/engine/ppp.rs`:
  - `compute_sat_state` (lines 431-434) adds back the group delay to undo the subtraction:
    ```rust
    dt_s = brdc_clk + eph.tgd();
    ```
- In `crates/gneiss-rtk/src/estimators/spp.rs`:
  - `compute_sat_state` (lines 185-191) queries `m.eph.position(...)` but does not add back the group delay when `m.is_iono_free` is `true`.

## 2. Logic Chain
- **Step 1**: The broadcast clock correction parameters ($a_{f0}, a_{f1}, a_{f2}$) are referenced to the dual-frequency (ionosphere-free) combination of signals. Therefore, the raw clock correction represents the dual-frequency clock error directly.
- **Step 2**: For single-frequency users, the group delay (`tgd`) must be subtracted from the raw clock correction to get the single-frequency clock correction.
- **Step 3**: For dual-frequency or ionosphere-free users, the group delay cancels in the ionosphere-free combination and must not be subtracted.
- **Step 4**: Because `position()` in `ephemeris.rs` always subtracts `tgd`, it represents a single-frequency corrected clock.
- **Step 5**: The PPP engine in `ppp.rs` correctly cancels the subtraction by adding back `tgd` (observed at `ppp.rs:434`). However, this is an API round-trip that relies on manually adding the value back in the estimator.
- **Step 6**: The SPP engine in `spp.rs` implements a dual-frequency ionosphere-free mode (`is_iono_free = true`), but fails to undo the TGD subtraction (observed at `spp.rs:185`), causing incorrect clock corrections for dual-frequency SPP users.

## 3. Caveats
- Single-frequency L2 or other frequency bands were not investigated in-depth as SPP currently only supports L1/E1 single frequency or dual-frequency iono-free. However, the proposed fix allows future implementation of other single-frequency corrections by exposing `position_iono_free()`.

## 4. Conclusion
The codebase needs a clean way to obtain the raw/ionosphere-free broadcast clock bias without subtracting group delay. Propose introducing `position_iono_free()` on `Ephemeris` that calls `calc_keplerian` with `tgd = 0.0`. This avoids the manual round-trip in PPP and fixes the incorrect clock corrections in dual-frequency SPP.

## 5. Verification Method
- **Verification Command**:
  Run the test suite using `cargo test --workspace` to ensure no existing tests are broken.
- **Regression Test**:
  Verify the implementation of a new unit test named `test_broadcast_clock_tgd_correct` in `crates/gneiss-core/src/ephemeris.rs` that asserts:
  - `eph.position(t).2` is equal to the raw clock minus the group delay.
  - `eph.position_iono_free(t).2` is equal to the raw clock.
  - The difference between the two is exactly the group delay value.
