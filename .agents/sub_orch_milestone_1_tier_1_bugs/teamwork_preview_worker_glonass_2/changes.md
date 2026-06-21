# Changes - Bug 17: GLONASS Time Scale Discrepancy Fix

## Modified Files
### `crates/gneiss-parsers/src/rinex.rs`
- **GLONASS Epoch Conversion**:
  Updated the GLONASS epoch time parsing logic to properly subtract the 3-hour Moscow Time offset (`10800.0` seconds) when converting the constellation's time to GPST. This is achieved by:
  ```rust
  let mut toc_gpst = parse_rinex_nav_epoch_time(&line, is_rinex_3);
  match current_constellation {
      Constellation::Glonass => toc_gpst = toc_gpst + (18.0 - 10800.0),
      Constellation::Beidou => toc_gpst = toc_gpst + 14.0,
      _ => {}
  }
  current_toc = toc_gpst;
  ```
  This ensures that addition of the negative offset utilizes `GpsTime`'s `Add<f64>` trait implementation, which automatically normalizes the `tow` and shifts the continuous GPS `week` accordingly.
  
- **Unit/Regression Test**:
  Updated `test_parse_rinex_3_nav_date` to assert the corrected TOW value of `411318.0` (shifted back by 3 hours from the uncorrected value of `422118.0`).

## Verification Results
- `cargo test -p gneiss-parsers` successfully passes all tests.
- `cargo test --workspace` successfully passes all 256 tests in the workspace with no regressions.
