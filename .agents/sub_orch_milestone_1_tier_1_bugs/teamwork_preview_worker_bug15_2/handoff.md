# Handoff Report: Bug 15 — Incorrect Broadcast Clock TGD Correction

## 1. Observation
We observed that the broadcast clock timing group delay (TGD) corrections for the Beidou constellation were incorrectly implemented across several files in the workspace:
- In `crates/gneiss-parsers/src/rinex.rs` (lines 573-574), the RINEX parser mapped `vals[23]` (representing `TGD2`) directly to `aodc` of `BeidouEphemeris`, and completely discarded the actual `AODC` value from `vals[25]`. It also lacked support for a `tgd2` field:
  ```rust
  tgd1: vals[22],
  aodc: vals[23] as u32,
  ```
- In `crates/gneiss-core/src/ephemeris.rs` (lines 155-180), `BeidouEphemeris` did not have a `tgd2` field.
- In `crates/gneiss-core/src/ephemeris.rs` (lines 707), `BeidouEphemeris::position_iono_free` passed a hardcoded `0.0` as the timing group delay `tgd` to `calc_keplerian`:
  ```rust
  calc_keplerian(
      ...
      0.0,
      ...
  )
  ```
- In `crates/gneiss-core/src/ephemeris.rs` (line 1524), `test_broadcast_clock_tgd_correct` asserted:
  ```rust
  assert!((clk_if - clk_pos - bds_eph.tgd1).abs() < 1e-15);
  ```
  which assumed `clk_if` did not have any TGD correction applied.
- The `cargo test --workspace` run was initially blocked by file locks from concurrent processes executing target build/test operations, and failed with directory creation errors after compile interruptions.

## 2. Logic Chain
- **Step 1**: To correctly align the dual-frequency (B1I/B2I) ionosphere-free broadcast clock corrections, the combined group delay $T_{GD\_IF}$ must be subtracted:
  $$T_{GD\_IF} = \frac{f_{B1I}^2 T_{GD1} - f_{B2I}^2 T_{GD2}}{f_{B1I}^2 - f_{B2I}^2}$$
  Therefore, the parser must correctly map both `tgd1` and `tgd2` fields, and `BeidouEphemeris` must hold both values.
- **Step 2**: We added the `tgd2` field to `BeidouEphemeris` in `crates/gneiss-core/src/ephemeris.rs` and updated `build_beidou_ephemeris` in `crates/gneiss-parsers/src/rinex.rs` to map `tgd2` from `vals[23]` and `aodc` from `vals[25] as u32`.
- **Step 3**: We updated `Ephemeris::bgd_e5b` to return `tgd2` for Beidou, and implemented the combined group delay calculation in `BeidouEphemeris::position_iono_free` using frequencies `FREQ_BDS_B1I` and `FREQ_GAL_E5B`.
- **Step 4**: We updated all instantiations of `BeidouEphemeris` across the workspace tests and SPP estimator tests to initialize `tgd2`. We corrected `test_broadcast_clock_tgd_correct` to assert:
  `clk_if - clk_pos = tgd1 - tgd_if`
- **Step 5**: To prevent file lock and directory cleanup conflicts with other concurrent workspace compilation tasks, we ran our checks and tests using a custom target directory: `--target-dir target_bug15`.

## 3. Caveats
- No caveats. The implementation covers all tasks, updates all occurrences of `BeidouEphemeris`, and passes all test assertions successfully.

## 4. Conclusion
The broadcast clock TGD correction for Beidou has been correctly fixed. All dual-frequency ionosphere-free calculations for Beidou now use the combined group delay $T_{GD\_IF}$, and all tests build and pass successfully.

## 5. Verification Method
To verify the implementation independently, run the following commands in the workspace root `/Users/kevin/projects/gneiss`:
```bash
cargo test --workspace --target-dir target_bug15
```
Inspect the files:
- `crates/gneiss-core/src/ephemeris.rs`
- `crates/gneiss-parsers/src/rinex.rs`
- `crates/gneiss-rtk/src/estimators/spp.rs`
Invalidation condition: If any test fails, or the Beidou ionosphere-free clock correction does not match `tgd1 - tgd_if`.
