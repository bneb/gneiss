# Handoff Report: GLONASS Time Scale Discrepancy Fix

## 1. Observation
- In `crates/gneiss-parsers/src/rinex.rs`, the original code converted the GLONASS epoch time to GPST without accounting for the 3-hour Moscow Time offset:
  ```rust
  let mut toc_gpst = parse_rinex_nav_epoch_time(&line, is_rinex_3);
  match current_constellation {
      Constellation::Glonass => toc_gpst.tow += 18.0,
      Constellation::Beidou => toc_gpst.tow += 14.0,
      _ => {}
  }
  ```
- The unit test `test_parse_rinex_3_nav_date` asserted the uncorrected TOW of `422118.0`:
  ```rust
  // Year 2020, Month 12, Day 24, Hour 21, Min 15
  assert_eq!(eph.toe().week, 2137);
  // Thursday 21:15 UTC + 18s leap seconds
  assert!(
      (eph.toe().tow - 422118.0).abs() < 1e-4,
      "Expected TOW near 422118.0, got {}",
      eph.toe().tow
  );
  ```
- Command `cargo test -p gneiss-parsers` and `cargo test --workspace` initially compiled and passed, but the GLONASS conversion logic was mathematically incorrect by 3 hours (10,800 seconds).

## 2. Logic Chain
- GLONASST is linked to Moscow Time, which is UTC+3 (10,800 seconds ahead of UTC).
- Converting GLONASST to GPST requires subtracting the 3-hour (10,800 seconds) offset and adding the leap second offset (18.0 seconds).
- Using direct manipulation on `toc_gpst.tow` (e.g. `toc_gpst.tow += 18.0 - 10800.0`) would leave the `GpsTime` struct in an unnormalized state with negative `tow`.
- `GpsTime` implements the `Add<f64>` trait, which internally normalizes the week and TOW:
  ```rust
  impl Add<f64> for GpsTime {
      type Output = Self;
      fn add(self, seconds: f64) -> Self::Output {
          GpsTime::new(self.week, self.tow + seconds)
      }
  }
  ```
- Therefore, replacing the logic with `toc_gpst = toc_gpst + (18.0 - 10800.0)` guarantees that the new offset is added correctly and the resulting `GpsTime` is automatically normalized.
- The unit test `test_parse_rinex_3_nav_date` must assert the corrected TOW `411318.0` (which is `422118.0 - 10800.0`).

## 3. Caveats
- No caveats. The fix aligns perfectly with GNSS time conversion specifications.

## 4. Conclusion
- The GLONASS time scale discrepancy bug has been successfully resolved. GLONASS ephemeris epochs are now correctly converted to GPST by subtracting the 3-hour Moscow Time offset and properly normalising the `GpsTime` struct. All tests pass successfully.

## 5. Verification Method
- Execute the unit tests specifically for `gneiss-parsers`:
  ```bash
  cargo test -p gneiss-parsers
  ```
- Verify workspace tests:
  ```bash
  cargo test --workspace
  ```
- Verify that `test_parse_rinex_3_nav_date` asserts `411318.0` in `crates/gneiss-parsers/src/rinex.rs`.
