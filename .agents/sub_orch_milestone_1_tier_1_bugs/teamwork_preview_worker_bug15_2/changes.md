# Changes: Fix Bug 15 — Incorrect Broadcast Clock TGD Correction

## Modifications

### 1. `crates/gneiss-core/src/ephemeris.rs`
- Added `pub tgd2: f64` field to the `BeidouEphemeris` struct.
- Modified `Ephemeris::bgd_e5b` to return `e.tgd2` when matching `Ephemeris::Beidou(e)` instead of using `other.tgd()` (which defaults to `tgd1`).
- Modified `BeidouEphemeris::position_iono_free` to calculate the combined ionosphere-free group delay $T_{GD\_IF}$:
  $$T_{GD\_IF} = \frac{f_{B1I}^2 T_{GD1} - f_{B2I}^2 T_{GD2}}{f_{B1I}^2 - f_{B2I}^2}$$
  using frequencies `FREQ_BDS_B1I` and `FREQ_GAL_E5B`. Passed this value to `calc_keplerian` instead of `0.0`.
- Initialized `tgd2: 0.0` in `bds_geo_eph`, `bds_igso_eph`, and `test_bds_boundaries` test cases.
- In `test_broadcast_clock_tgd_correct`, set `tgd1: 4.0e-9` and `tgd2: 2.0e-9` in `bds_eph`, calculated `tgd_if` using the same formula, and corrected the assertion:
  `assert!((clk_if - clk_pos - (bds_eph.tgd1 - tgd_if)).abs() < 1e-15);`

### 2. `crates/gneiss-parsers/src/rinex.rs`
- In `build_beidou_ephemeris`, mapped `tgd2` from `vals[23]` and corrected `aodc` mapping to use `vals[25] as u32` (instead of using `vals[23]` for `aodc` and discarding the real `aodc` value).

### 3. `crates/gneiss-rtk/src/estimators/spp.rs`
- Initialized `tgd2: 0.0` in the test case `eph_bds` instantiation.

## Verification Outcomes
- Build run via `cargo check --workspace --target-dir target_bug15` completed successfully.
- Test run via `cargo test --workspace --target-dir target_bug15` successfully passed 683/683 tests.
