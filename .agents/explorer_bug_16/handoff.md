# Handoff Report: Bug 16 - Mismatched Galileo BGD Correction

## Summary
The codebase contains a mismatched Galileo Group Delay (BGD) correction where `bgd_e1_e5a` is applied to Galileo E5b (band 7) observations instead of `bgd_e1_e5b`. This investigation identifies the locations where this occurs and proposes a clean diff patch to correct it.

---

## 1. Observation

Direct code observations from the repository:

1. **`crates/gneiss-core/src/ephemeris.rs:614-642`**:
   `GalileoEphemeris` defines a function `position_e5b` which correctly uses `bgd_e1_e5b` to compute Keplerian parameters and corrects the clock error:
   ```rust
   pub fn position_e5b(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
       calc_keplerian(
           t,
           self.toe,
           self.toc,
           self.af0,
           self.af1,
           self.af2,
           self.crs,
           self.crc,
           self.cuc,
           self.cus,
           self.cic,
           self.cis,
           self.m0,
           self.e,
           self.sqrt_a,
           self.delta_n,
           self.omega0,
           self.omega_dot,
           self.i0,
           self.idot,
           self.omega,
           self.bgd_e1_e5b, // Bug 16: use E5b BGD, not E5a
           MU_GAL,
           OMEGA_E_GAL,
           false,
       )
   }
   ```

2. **`crates/gneiss-core/src/ephemeris.rs:37-45`**:
   The general enum `Ephemeris` delegates `position` to the constellation-specific structs. For Galileo, it delegates to `e.position(t)`, which always corrects using `bgd_e1_e5a`:
   ```rust
   impl Ephemeris {
       pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
           match self {
               Ephemeris::Gps(e) => e.position(t),
               Ephemeris::Galileo(e) => e.position(t),
               Ephemeris::Beidou(e) => e.position(t),
               Ephemeris::Qzss(e) => e.position(t),
               Ephemeris::Glonass(e) => e.position(t),
           }
       }
   }
   ```
   No `position_e5b` delegate exists on the `Ephemeris` enum.

3. **`crates/gneiss-rtk/src/estimators/spp.rs:185-198`** (and similarly in `crates/gneiss-rtk/src/engine/spp_tight.rs:214-220`):
   When computing single-frequency satellite clock and coordinates, SPP calls `position`:
   ```rust
   let (_, _, sat_clk_err_rough, _) = if m.is_iono_free {
       m.eph.position_iono_free(t_tx_sat_gps)
   } else {
       m.eph.position(t_tx_sat_gps)
   };
   ```

4. **`crates/gneiss-core/src/signal.rs:27-53`**:
   The helper function `get_frequency` lacks support for band 7 (Galileo E5b/Beidou B2b), causing it to default to `FREQ_GPS_L1` via the wildcard match `_ => FREQ_GPS_L1`:
   ```rust
   pub fn get_frequency(sat: SatelliteId, freq_band: u8, freq_num: i8) -> f64 {
       match freq_band {
           1 => ...
           2 => ...
           5 => ...
           _ => FREQ_GPS_L1,
       }
   }
   ```

---

## 2. Logic Chain

1. **Galileo Single-Frequency Processing**: Galileo satellites broadcast two separate Broadcast Group Delays (BGDs):
   - `bgd_e1_e5a` for Galileo E5a (band 5) single-frequency users.
   - `bgd_e1_e5b` for Galileo E5b (band 7) single-frequency users.
2. **Standard delegation always uses E5a BGD**: The general `Ephemeris::position()` method calls `GalileoEphemeris::position()`, which always subtracts `bgd_e1_e5a`.
3. **Mismatched correction in estimators**: Since the single-frequency estimators (`spp.rs`, `spp_tight.rs`) only have access to `Ephemeris::position()`, any Galileo E5b (band 7) single-frequency observations will have their satellite clocks corrected using `bgd_e1_e5a` instead of `bgd_e1_e5b`.
4. **Missing frequency mapping**: Without matching band 7 in `get_frequency`, any code path querying frequency/wavelength for band 7 observations will get the incorrect L1 frequency, introducing modeling errors.

---

## 3. Caveats

- **Double-differenced (DD) EKF RTK**: Satellite clock errors cancel out completely in double-difference combinations. Thus, this mismatch has negligible impact on standard RTK positioning, but directly affects Single Point Positioning (SPP) and tight coupling estimators.
- **Iono-free dual-frequency PPP**: The iono-free clock combination removes group delay effects and uses `position_iono_free()` (where `tgd = 0.0`), which is already correct.

---

## 4. Conclusion

To use `bgd_e1_e5b` for band 7 Galileo E5b observations, we must:
1. Expose `position_e5b` on the `Ephemeris` enum.
2. Update `get_frequency` to support band 7.
3. Add a `freq_band` field to `SppMeasurement` so SPP estimators can distinguish frequency bands.
4. Update `build_single_measurement` in `spp.rs` to track the frequency band and support single-frequency E5b.
5. Update `compute_sat_state` in `spp.rs` and `process_measurement` in `spp_tight.rs` to invoke `position_e5b` when `m.freq_band == 7`.

A complete, machine-applicable patch file has been created at:
`/Users/kevin/projects/gneiss/.agents/explorer_bug_16/fix_bug_16.patch`

---

## 5. Verification Method

To verify the proposed fix independently:
1. Run standard project tests to ensure no regressions:
   ```bash
   cargo test --workspace
   ```
2. Write a regression unit test in `crates/gneiss-core/src/ephemeris.rs` to verify that `position_e5b` produces the expected clock correction difference compared to `position` (matching `bgd_e1_e5b - bgd_e1_e5a` in meters).
3. Verify that `get_frequency(sat_galileo, 7, 0)` returns `FREQ_GAL_E5B` (1207.14 MHz).
