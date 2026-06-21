# Handoff Report: GLONASS Time Scale Discrepancy (Bug 17)

## 1. Observation

In the file `crates/gneiss-parsers/src/rinex.rs` (lines 686-692), during navigation file parsing, the time of clock/ephemeris (`toc_gpst`) is adjusted based on the constellation:
```rust
            let mut toc_gpst = parse_rinex_nav_epoch_time(&line, is_rinex_3);
            match current_constellation {
                Constellation::Glonass => toc_gpst.tow += 18.0,
                Constellation::Beidou => toc_gpst.tow += 14.0,
                _ => {}
            }
            current_toc = toc_gpst;
```

Additionally, in `crates/gneiss-parsers/src/rinex.rs` (lines 945-971), the unit test `test_parse_rinex_3_nav_date` asserts this behavior for a GLONASS navigation record:
```rust
    #[test]
    fn test_parse_rinex_3_nav_date() {
        let file_contents =
            "     3.03           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
R 6 2020 12 24 21 15  0  .189751386642E-03  .000000000000E+00  .422910000000E+06
     -.740158740234E+04 -.212037086487E+00  .000000000000E+00  .000000000000E+00
     -.206682856445E+05 -.176755714417E+01  .931322574615E-09 -.400000000000E+01
      .129489067383E+05 -.294115734100E+01 -.186264514923E-08  .000000000000E+00
";
        use std::io::BufReader;
        let mut reader = BufReader::new(file_contents.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();

        assert_eq!(ephemerides.len(), 1);
        let eph = &ephemerides[0];

        // Year 2020, Month 12, Day 24, Hour 21, Min 15
        assert_eq!(eph.toe().week, 2137);
        // Thursday 21:15 UTC + 18s leap seconds
        assert!(
            (eph.toe().tow - 422118.0).abs() < 1e-4,
            "Expected TOW near 422118.0, got {}",
            eph.toe().tow
        );
    }
```
The test comment states `// Thursday 21:15 UTC + 18s leap seconds`, and verifies that TOW is `422118.0` (which is `422100.0 + 18.0`).

## 2. Logic Chain

1. **GLONASS Time Reference:** The GLONASS time scale (GLONASST) is reference-linked to Russian UTC(SU) with a constant offset of +3 hours (Moscow Time), i.e., $GLONASST = UTC + 3\text{ hours}$ (10,800 seconds).
2. **RINEX Ephemeris Time Scale:** Per the RINEX specification (e.g., RINEX 3.03, Section 5.3), the epoch time fields (`toc` / `toe`) in a GLONASS navigation file are written in GLONASS Time (GLONASST).
3. **GPS Time Conversion:** 
   - GPS Time (GPST) is a continuous time scale that is ahead of UTC by leap seconds (18.0 seconds for the relevant epoch since 2017).
   - Thus, $GPST = UTC + 18.0\text{ seconds}$.
   - Substituting UTC: $GPST = (GLONASST - 3\text{ hours}) + 18.0\text{ seconds} = GLONASST - 10,800.0\text{ s} + 18.0\text{ s} = GLONASST - 10,782.0\text{ s}$.
4. **Code Discrepancy:**
   - The parser `parse_rinex_nav_epoch_time` extracts the date-time numbers directly from the file (which are in GLONASST) and creates a `GpsTime` object from them.
   - For `Constellation::Glonass`, the code then does `toc_gpst.tow += 18.0;`.
   - This completely misses the subtraction of 3 hours (10,800 seconds), causing the resulting `toe` / `toc` for GLONASS ephemerides to be exactly 3 hours too late in the GPS Time scale.
5. **Propagation Impact:**
   - In `GlonassEphemeris::position`, the orbit is integrated numerically using `let dt = t - self.toe;`.
   - An error of $+10,800$ seconds in `toe` translates to a $-10,800$ seconds error in the propagation time `dt`.
   - Since GLONASS satellites travel at approximately $2.5\text{ km/s}$, integrating the orbit with a $3\text{ hour}$ timing discrepancy results in a position error of roughly $10,800\text{ s} \times 2.5\text{ km/s} \approx 27,000\text{ kilometers}$.
   - This massive coordinate discrepancy causes GLONASS observations to be completely discarded (or fail to converge) in solvers like SPP or PPP.

## 3. Caveats

- We assume that the receiver observations are correctly aligned to GPST and do not suffer from similar timezone offsets (which has been verified: standard receiver logs and UBX measurements output raw `tow` aligned to GPS/GNSS system time, not local Moscow Time).
- The leap second value of `18.0` is hardcoded in the parser and tests. While leap seconds can change, the codebase currently standardizes on this constant offset, so the fix strategy preserves this standard for consistency with the rest of the engine.

## 4. Conclusion

The GLONASS navigation parsing code has a 3-hour timezone discrepancy because it does not subtract the Moscow Time offset ($3\text{ hours} = 10,800\text{ seconds}$) when converting the GLONASS epoch time (GLONASST) to GPS Time (GPST).

### Proposed Fix Strategy

1. **Modify Time Conversion in `crates/gneiss-parsers/src/rinex.rs`:**
   Instead of just adding `18.0` to the TOW for GLONASS, subtract `10782.0` (which is $18.0 - 10800.0$) using the `+` operator to ensure the value is correctly normalized (as `GpsTime` implements `Add<f64>` which automatically normalizes the week and time-of-week).
   
   *Before:*
   ```rust
             let mut toc_gpst = parse_rinex_nav_epoch_time(&line, is_rinex_3);
             match current_constellation {
                 Constellation::Glonass => toc_gpst.tow += 18.0,
                 Constellation::Beidou => toc_gpst.tow += 14.0,
                 _ => {}
             }
   ```
   *After:*
   ```rust
             let mut toc_gpst = parse_rinex_nav_epoch_time(&line, is_rinex_3);
             match current_constellation {
                 Constellation::Glonass => toc_gpst = toc_gpst + (18.0 - 10800.0),
                 Constellation::Beidou => toc_gpst = toc_gpst + 14.0,
                 _ => {}
             }
   ```

2. **Update the Unit Test in `crates/gneiss-parsers/src/rinex.rs`:**
   Update `test_parse_rinex_3_nav_date` to assert the correct GPS TOW after subtracting the 3-hour Moscow Time offset.
   
   - The GLONASST epoch in the test is `2020-12-24 21:15:00` (Thursday).
   - Calendar TOW for `21:15:00` is $4\text{ days} \times 86400\text{ s} + 21 \times 3600\text{ s} + 15 \times 60\text{ s} = 422,100.0\text{ s}$.
   - Correct GPST TOW is $422,100.0 - 10,800.0\text{ (offset)} + 18.0\text{ (leap)} = 411,318.0\text{ s}$.
   
   *Before:*
   ```rust
         // Thursday 21:15 UTC + 18s leap seconds
         assert!(
             (eph.toe().tow - 422118.0).abs() < 1e-4,
             "Expected TOW near 422118.0, got {}",
             eph.toe().tow
         );
   ```
   *After:*
   ```rust
         // Thursday 21:15 Moscow Time (GLONASST) = 18:15:00 UTC + 18s leap seconds = 18:15:18 GPST
         assert!(
             (eph.toe().tow - 411318.0).abs() < 1e-4,
             "Expected TOW near 411318.0, got {}",
             eph.toe().tow
         );
   ```

## 5. Verification Method

To independently verify the fix strategy:
1. Inspect `crates/gneiss-parsers/src/rinex.rs` at the lines indicated.
2. Run the existing tests using `cargo test -p gneiss-parsers`.
3. Apply the changes proposed in the Fix Strategy.
4. Run `cargo test -p gneiss-parsers` again to confirm that the updated test `test_parse_rinex_3_nav_date` passes with the corrected expected value.
