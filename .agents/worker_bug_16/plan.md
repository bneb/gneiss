# Implementation Plan - Bug 16 Fix

This plan details the steps to implement the fix for Bug 16: Mismatched Galileo BGD Correction.

## Steps

### Step 1: Define `position_e5b` on `Ephemeris` and `GalileoEphemeris`
- **File**: `crates/gneiss-core/src/ephemeris.rs`
- **Implementation**:
  - In `GalileoEphemeris` impl, add `pub fn position_e5b(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64)` which calculates keplerian orbit projection using `bgd_e1_e5b` instead of `bgd_e1_e5a`.
  - In `Ephemeris` impl, add `pub fn position_e5b(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64)`:
    ```rust
    pub fn position_e5b(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        match self {
            Ephemeris::Galileo(e) => e.position_e5b(t),
            other => other.position(t),
        }
    }
    ```
- **Verification**: Run `cargo check -p gneiss-core` to verify it compiles.

### Step 2: Handle Band 7 in `get_frequency`
- **File**: `crates/gneiss-core/src/signal.rs`
- **Implementation**:
  - Add `7` arm to `get_frequency`:
    ```rust
        7 => match sat.constellation {
            Constellation::Galileo | Constellation::Beidou => FREQ_GAL_E5B,
            _ => FREQ_GPS_L2, // Default fallback
        },
    ```
- **Verification**: Run `cargo check -p gneiss-core` to verify it compiles.

### Step 3: Add `freq_band` field to `SppMeasurement`
- **File**: `crates/gneiss-rtk/src/estimators/spp.rs`
- **Implementation**:
  - Add `pub freq_band: u8` field to the `SppMeasurement` struct definition.
- **Verification**: Run `cargo check -p gneiss-rtk` to verify compilation.

### Step 4: Set `freq_band` in `build_single_measurement`
- **File**: `crates/gneiss-rtk/src/estimators/spp.rs`
- **Implementation**:
  - Update `build_single_measurement` to set `freq_band` based on:
    - If `p2_opt` is used (meaning either iono-free or single-freq p2_opt), determine what band it corresponds to.
      Wait! Let's check how `freq_band` is determined:
      - "1 if p1_opt, 7 if Galileo/Beidou and p2_opt, 2 for other p2_opt, etc."
      Wait, is it referring to the band used by the pseudorange observation?
      Wait, if we use iono-free, the combination is constructed from both bands. But wait, since we are doing iono-free, which band's BGD should we apply?
      Actually, the request says: "1 if p1_opt, 7 if Galileo/Beidou and p2_opt, 2 for other p2_opt, etc."
      Wait! If `p1_opt` is some value and `p2_opt` is `Some`, we use iono-free combination. Let's look at `build_single_measurement` again:
      ```rust
      let (raw_pr, is_iono_free) = if let (Some(p1), Some(p2)) = (p1_opt, p2_opt) {
          ...
          ((f1_sq * p1 - f2_sq * p2) / (f1_sq - f2_sq), true)
      } else {
          (p1_opt?, false)
      };
      ```
      Wait, if `is_iono_free` is true, we are combining p1 and p2. But which BGD correction does this combination need?
      Wait, in the dual-frequency iono-free combination, the Galileo broadcast clock is already corrected for E1/E5a iono-free. If we use the E1/E5b iono-free combination (band 7), we must apply the E1-E5b BGD correction!
      Wait, what if `p2_opt` is not present? Then we are single frequency, and it uses `p1_opt`. In that case, `freq_band` should be `1` (or whatever the band of `p1` is, usually `1`).
      So the logic for `freq_band` should be:
      - If `p2_opt.is_some()`:
        - If constellation is Galileo or Beidou, `freq_band = 7`.
        - Else `freq_band = 2` (or potentially 5 or other if another constellation uses it).
      - Else (only `p1_opt` is present):
        - `freq_band = 1` (or 2 for Beidou since Beidou `p1_opt` gets observable 2? Let's check).
        Let's read finding 4 again:
        "4. `build_single_measurement` needs to set `freq_band` based on observations (1 if p1_opt, 7 if Galileo/Beidou and p2_opt, 2 for other p2_opt, etc.)."
        Let's look at this carefully:
        - "1 if p1_opt" -> wait, does it mean "if we fall back to p1_opt, we set it to 1"?
          Let's re-read: "1 if p1_opt, 7 if Galileo/Beidou and p2_opt, 2 for other p2_opt, etc."
          Yes, if we use `p1_opt` (which means `p2_opt` is `None`), we set `freq_band = 1`.
          Wait, what if `p2_opt` is present? Then:
          - If Galileo or Beidou: `freq_band = 7`.
          - Else: `freq_band = 2`.
          Let's check if there are other cases, e.g. "etc."
          Wait, is there any other band? What about GPS L5? In SPP, is L5 ever used as `p2`?
          Let's look at `p2_opt` definition:
          ```rust
          let p2_opt = match sat_obs.sat.constellation {
              gneiss_core::sat::Constellation::Galileo => {
                  sat_obs.get_observable(7).or(sat_obs.get_observable(5))
              }
              ...
          ```
          Ah! For Galileo, `p2_opt` is `sat_obs.get_observable(7).or(sat_obs.get_observable(5))`.
          Wait! If `sat_obs.get_observable(7)` is not found, but `sat_obs.get_observable(5)` is found, then the observable used is band 5.
          In that case, should `freq_band` be 5 instead of 7?
          Yes! "7 if Galileo/Beidou and p2_opt" -> wait, if `sat_obs.get_observable(7)` is selected, then `freq_band = 7`. If `sat_obs.get_observable(5)` is selected, then `freq_band = 5`.
          Wait! How can we know which one was selected?
          Let's look at:
          ```rust
          let p2_opt = match sat_obs.sat.constellation {
              gneiss_core::sat::Constellation::Galileo => {
                  sat_obs.get_observable(7).or(sat_obs.get_observable(5))
              }
              ...
          ```
          Wait, we can rewrite this or track the selected band:
          ```rust
          let mut freq_band = 1;
          let p1_opt = match sat_obs.sat.constellation {
              gneiss_core::sat::Constellation::Beidou => sat_obs.get_observable(2),
              _ => sat_obs.get_observable(1),
          };
          // Actually, if we only use p1_opt, freq_band is 1. Wait, for Beidou, p1_opt is observable 2. But we can still use 1 or 2 as the band?
          // Finding 4 says: "1 if p1_opt, 7 if Galileo/Beidou and p2_opt, 2 for other p2_opt, etc."
          // Let's implement it exactly like:
          let p2_opt = match sat_obs.sat.constellation {
              gneiss_core::sat::Constellation::Galileo => {
                  if let Some(obs) = sat_obs.get_observable(7) {
                      freq_band = 7;
                      Some(obs)
                  } else if let Some(obs) = sat_obs.get_observable(5) {
                      freq_band = 5;
                      Some(obs)
                  } else {
                      None
                  }
              }
              gneiss_core::sat::Constellation::Beidou => {
                  if let Some(obs) = sat_obs.get_observable(7) {
                      freq_band = 7;
                      Some(obs)
                  } else if let Some(obs) = sat_obs.get_observable(6) {
                      freq_band = 6;
                      Some(obs)
                  } else {
                      None
                  }
              }
              _ => {
                  if let Some(obs) = sat_obs.get_observable(2) {
                      freq_band = 2;
                      Some(obs)
                  } else {
                      None
                  }
              }
          };

          let (raw_pr, is_iono_free) = if let (Some(p1), Some(p2)) = (p1_opt, p2_opt) {
              let f1_sq = f1 * f1;
              let f2_sq = f2 * f2;
              ((f1_sq * p1 - f2_sq * p2) / (f1_sq - f2_sq), true)
          } else {
              freq_band = 1;
              (p1_opt?, false)
          };
          ```
          Wait, is this correct?
          If `p2_opt` is `Some` and `p1_opt` is `Some`, we use iono-free combination. In that case, `freq_band` is set to the band of `p2_opt` (which is `7` or `5` or `6` or `2`).
          Wait! What if `p1_opt` is `None`? Then `p1_opt?` will return `None`, and the function returns `None`.
          If `p2_opt` is `None`, the `else` branch is taken: `freq_band = 1; (p1_opt?, false)`.
          This is extremely clean and matches exactly "1 if p1_opt, 7 if Galileo/Beidou and p2_opt, 2 for other p2_opt, etc."!
- **Verification**: Run `cargo check -p gneiss-rtk` to verify it compiles.

### Step 5: Update `compute_sat_state` and `process_measurement`
- **Files**:
  - `crates/gneiss-rtk/src/estimators/spp.rs`
  - `crates/gneiss-rtk/src/engine/spp_tight.rs`
- **Implementation**:
  - In `compute_sat_state`, call `position_e5b` if `m.freq_band == 7`, else `position`.
  - In `process_measurement`, call `position_e5b` if `m.freq_band == 7`, else `position`.
- **Verification**: Run `cargo check` and `cargo test`.

### Step 6: Add Regression Unit Test
- **File**: `crates/gneiss-core/src/ephemeris.rs` or similar test file.
- **Implementation**:
  - Add a unit test verifying that `Ephemeris::position_e5b` is called when `freq_band == 7` and corrects using `bgd_e1_e5b`.
  - The test should verify that for a Galileo ephemeris, `position_e5b` yields the position with `bgd_e1_e5b` clock correction, whereas `position` yields the position with `bgd_e1_e5a` correction.
- **Verification**: Run `cargo test` and ensure the regression test passes.
