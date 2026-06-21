# Analysis: Bug 18 — Opposite Sign in Phase Wind-Up Correction in `gneiss-rtk`

## Executive Summary
Phase wind-up is an electromagnetic propagation effect that shifts the measured carrier phase of circularly polarized GPS signals due to the relative rotation between the satellite and receiver antennas. The correct physical equation requires **subtracting** the computed phase wind-up correction (in cycles) from the raw carrier phase observations to obtain the corrected carrier phase matching the geometric propagation model. 

However, the current implementation in `gneiss-rtk` consistently **adds** the wind-up correction to the carrier phase observations in both RTK (double-differenced) and PPP (single-differenced/undifferenced) engines. This results in an opposite sign correction, leading to degraded positioning accuracy and filter convergence issues when phase wind-up is non-negligible.

---

## 1. Problem Investigation & Evidence

### Physical Principle of Phase Wind-Up Correction
The measured carrier phase $\phi_{\text{meas}}$ (in cycles) is modeled as:
$$\phi_{\text{meas}} = \phi_{\text{geom}} + N + \phi_{\text{windup}} + \epsilon$$

Where:
- $\phi_{\text{geom}}$ is the geometric range in cycles ($\rho / \lambda$).
- $N$ is the integer ambiguity.
- $\phi_{\text{windup}}$ is the phase wind-up correction in cycles.

To compare the observation directly to the geometric model, the corrected carrier phase observation $\phi_{\text{corrected}}$ must be:
$$\phi_{\text{corrected}} = \phi_{\text{meas}} - \phi_{\text{windup}}$$

### Codebase Investigation & Discrepancies
The codebase computes phase wind-up using `gneiss_core::windup::phase_windup`, which correctly implements the standard wind-up model (Wu et al. 1992). However, three areas in `gneiss-rtk` incorrectly apply the correction with a positive sign (addition instead of subtraction).

#### 1. RTK (Double-Differenced) Measurement Engine
In `crates/gneiss-rtk/src/engine/measurement.rs`:
```rust
fn apply_windup_to_obs(obs: &mut DdObservation, windup: f64) {
    if let Some(cp) = &mut obs.cp_l1 {
        *cp += windup; // <-- Bug: adds instead of subtracting
    }
    if let Some(cp2) = &mut obs.cp_l2 {
        *cp2 += windup; // <-- Bug: adds instead of subtracting
    }
}
```

#### 2. PPP Phase Ambiguity Resolution & Pre-alignment
In `crates/gneiss-rtk/src/engine/ppp.rs`:
- Inside `update_phase_ambiguities` (lines 495-504, 545-546):
```rust
        let l_meas = if sat.is_iono_free && sat.cp2.is_some() {
            crate::engine::ppp_math::compute_iono_free(
                (cp1 + wup) * sat.lam1,             // <-- Bug: cp1 + wup
                (sat.cp2.unwrap() + wup) * sat.lam2, // <-- Bug: cp2 + wup
                sat.f1,
                sat.f2,
            )
        } else {
            (cp1 + wup) * sat.lam1 // <-- Bug: cp1 + wup
        };
        // ...
        if !sat.is_iono_free && sat.cp2.is_some() && sat.p2.is_some() {
            let l1_m = (cp1 + wup) * sat.lam1;             // <-- Bug: cp1 + wup
            let l2_m = (sat.cp2.unwrap() + wup) * sat.lam2; // <-- Bug: cp2 + wup
```
- Inside `add_uduc_ambiguities` (lines 588-589):
```rust
    let l1_meas = (cp1 + wup) * sat.lam1;             // <-- Bug: cp1 + wup
    let l2_meas = (sat.cp2.unwrap() + wup) * sat.lam2; // <-- Bug: cp2 + wup
```

#### 3. PPP Iterated EKF (IEKF) Solver
In `crates/gneiss-rtk/src/engine/ppp_iekf.rs`:
- Inside carrier phase prediction (lines 897-898):
```rust
            let windup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
            let l_meas = (cp1 + windup) * sat.lam1; // <-- Bug: cp1 + windup
```
- Inside UDUC phase measurement residual construction (lines 1026, 1043):
```rust
        let res_l1 = (sat.cp1.unwrap() + windup) * sat.lam1 - (expected_base - i1 + n1); // <-- Bug: cp1 + windup
        // ...
        let res_l2 = (sat.cp2.unwrap() + windup) * sat.lam2 - (expected_base - gamma * i1 + n2); // <-- Bug: cp2 + windup
```

---

## 2. Recommended Fix Strategy

All additions of the phase wind-up correction to the carrier phase observations must be changed to subtractions. 

### Step 1: Fix RTK double-differenced observations
In `crates/gneiss-rtk/src/engine/measurement.rs`, change the operators to `-=`:
```rust
fn apply_windup_to_obs(obs: &mut DdObservation, windup: f64) {
    if let Some(cp) = &mut obs.cp_l1 {
        *cp -= windup;
    }
    if let Some(cp2) = &mut obs.cp_l2 {
        *cp2 -= windup;
    }
}
```

### Step 2: Fix PPP pre-alignment and ambiguity resolution calculations
In `crates/gneiss-rtk/src/engine/ppp.rs`:
1. In `update_phase_ambiguities`, modify the calculation of `l_meas`:
```rust
        let l_meas = if sat.is_iono_free && sat.cp2.is_some() {
            crate::engine::ppp_math::compute_iono_free(
                (cp1 - wup) * sat.lam1,
                (sat.cp2.unwrap() - wup) * sat.lam2,
                sat.f1,
                sat.f2,
            )
        } else {
            (cp1 - wup) * sat.lam1
        };
```
2. In the same function, modify the calculation of `l1_m` and `l2_m`:
```rust
        if !sat.is_iono_free && sat.cp2.is_some() && sat.p2.is_some() {
            let l1_m = (cp1 - wup) * sat.lam1;
            let l2_m = (sat.cp2.unwrap() - wup) * sat.lam2;
```
3. In `add_uduc_ambiguities`, modify the calculation of `l1_meas` and `l2_meas`:
```rust
    let l1_meas = (cp1 - wup) * sat.lam1;
    let l2_meas = (sat.cp2.unwrap() - wup) * sat.lam2;
```

### Step 3: Fix PPP IEKF solver phase residuals
In `crates/gneiss-rtk/src/engine/ppp_iekf.rs`:
1. In `predict_carrier_phase` or similar, modify `l_meas`:
```rust
            let windup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
            let l_meas = (cp1 - windup) * sat.lam1;
```
2. In the residual calculations inside UDUC phase measurement construction:
```rust
        let res_l1 = (sat.cp1.unwrap() - windup) * sat.lam1 - (expected_base - i1 + n1);
        // ...
        let res_l2 = (sat.cp2.unwrap() - windup) * sat.lam2 - (expected_base - gamma * i1 + n2);
```

---

## 3. Unit Regression Test Design

To prevent future regression, we must rewrite/supplement unit tests to explicitly assert the direction of the wind-up correction (that a positive wind-up value results in a decreased corrected carrier phase).

### 1. RTK Engine Regression Test
In `crates/gneiss-rtk/src/engine/measurement.rs`, supplement the unit test suite with:
```rust
    #[test]
    fn test_phase_windup_correction_sign_rtk() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        
        let mut obs = DdObservation {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            pr_l1: 0.0,
            pr_l2: None,
            cp_l1: Some(10.0),
            cp_l2: Some(20.0),
            doppler: 0.0,
            snr: 45.0,
            locktime: None,
        };

        let windup = 0.25; // 0.25 cycles of positive wind-up
        apply_windup_to_obs(&mut obs, windup);

        // Corrected carrier phase = raw_cp - windup
        assert_eq!(obs.cp_l1.unwrap(), 9.75);
        assert_eq!(obs.cp_l2.unwrap(), 19.75);
    }
```

### 2. PPP Engine Regression Test
We can also verify the PPP pre-alignment logic. In `crates/gneiss-rtk/src/engine/ppp.rs`:
```rust
    #[test]
    fn test_phase_windup_correction_sign_ppp() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use crate::filter::RtkState;
        use gneiss_core::coords::{Coordinate, Datum, Frame};
        use gneiss_core::time::GpsTime;

        let time = GpsTime::new(2137, 422922.0);
        let mut state = RtkState::new(
            time,
            Coordinate::new(
                Vector3::new(1000.0, 2000.0, 3000.0),
                Datum::WGS84,
                Frame::ECEF,
                time,
            ),
            10.0,
        );

        // Inject a known positive wind-up correction into the state
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let windup_val = 0.25;
        state.windup.insert(sat, windup_val);

        // Create a dummy ProcessedSat
        // When update_phase_ambiguities runs, it will read `windup_val` (or compute it, 
        // but we can test that the correction subtraction logic works for generating l_meas).
        // Alternatively, test the mathematical subtraction itself.
    }
```
