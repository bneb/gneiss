# Project: Gneiss Navigation Engine Bug Fixes

## Architecture
The Gneiss navigation engine is structured as a multi-crate Rust workspace:
- `gneiss-core`: General GNSS algorithms, troposphere/ionosphere atmospheric models, tides, wind-up corrections, broadcast ephemeris orbit/clock calculations.
- `gneiss-rtk`: High-precision PPP and RTK EKF estimation engine, measurement models, predictions, sequential ambiguity resolution (AR), and Melbourne-Wübbena (MW) combination calculations.
- `gneiss-parsers`: RINEX observation and precise clock/orbit (SP3/CLK) parser logic.

## Milestones

| # | Milestone Name | Scope (Bugs to Fix) | Dependencies | Status |
|---|----------------|---------------------|--------------|--------|
| 1 | **Milestone 1: Tier 1 Bugs** | Fix ranks 1 to 8 (Bugs 17, 1, 9, 2, 18, 15, 24, 6) | None | `IN_PROGRESS` (f16afb25-c177-42fe-985d-6840e173046f) |
| 2 | **Milestone 2: Tier 2 & 3 Bugs** | Fix ranks 9 to 18 (Bugs 10, 25, 12, 23, 16, 5, 11, 8, 22, 3) | Milestone 1 | `PLANNED` |
| 3 | **Milestone 3: Tier 4 Features** | Implement/Stub ranks 19 to 25 (Bugs 13, 19, 20, 21, 14, 7, 4) | Milestone 2 | `PLANNED` |

## Detailed Bug Mapping

### Milestone 1: Tier 1 Bugs (Ranks 1–8)
1. **GLONASS Time Scale Discrepancy (Bug 17)**: Ensure GLONASS orbit calculation handlesGPST vs GLONASST time differences consistently.
2. **Melbourne-Wübbena Typo (Bug 1)**: Correct MW combination scaling factor to `(lam2 - lam1) / (lam1 + lam2)`.
3. **Sequential AR Covariance (Bug 9)**: Compute Narrowlane covariance from post-Widelane update state covariance.
4. **Velocity-Attitude Jacobian (Bug 2)**: Fix EKF transition matrix attitude skew sign to positive (`+f_e_skew`).
5. **Opposite Wind-Up Sign (Bug 18)**: Subtract wind-up correction `wup` instead of adding it.
6. **Broadcast Clock TGD Correction (Bug 15)**: Stop subtracting TGD from broadcast clocks for dual-frequency/iono-free measurements.
7. **Outlier Precise Clock Gap (Bug 24)**: Return `None` (not stale bias) when gap in precise clock records exceeds 900s limit.
8. **GMF Legendre Polynomials (Bug 6)**: Implement fully normalized Associated Legendre Functions for GMF tropospheric mapping.

### Milestone 2: Tier 2 & 3 Bugs (Ranks 9–18)
9. **TOF Sat Position Clock Bias (Bug 10)**: Account for receiver clock bias `cdt_r` in Time-of-Flight (TOF) calculation.
10. **Covariance Re-Init on Slip Omission (Bug 25)**: Inflate receiver position/velocity covariances when cycle slip occurs and ambiguity is reset.
11. **Receiver Antenna PCV Omission (Bug 12)**: Implement receiver PCV correction as function of azimuth/elevation.
12. **Klobuchar Evaluated at Receiver (Bug 23)**: Evaluate Klobuchar model at Ionospheric Pierce Point (IPP) instead of receiver position.
13. **Mismatched Galileo BGD (Bug 16)**: Use Galileo E1/E5b group delay correction instead of E1/E5a for E5b observations.
14. **GMF Troposphere Longitude Omission (Bug 5)**: Include longitude parameter in GMF spherical harmonics evaluator.
15. **Sat PCV Zenith-Dependent Correction Omission (Bug 11)**: Project and apply zenith-dependent PCV from ANTEX on line-of-sight.
16. **INS State AR Update Cutoff (Bug 8)**: Stop zeroing out the Kalman gain for INS states in the Narrowlane constraint update.
17. **Saastamoinen Dry Delay Pressure (Bug 22)**: Use actual/standard surface pressure instead of a constant sea-level pressure scaled by height.
18. **L2C Phase Shift Bias (Bug 3)**: Verify L2C phase shift behavior against receiver model/standard, and avoid hardcoding static 0.25 shift.

### Milestone 3: Tier 4 Features (Ranks 19–25)
19. **Lack of Ocean Tide Loading (Bug 13)**: Parse BLQ files and implement 11-constituent Ocean Tide Loading (OTL) corrections.
20. **Missing Yaw Steering Models (Bug 19)**: Implement constellation-specific yaw steering laws for BeiDou, Galileo, and GLONASS.
21. **Beidou B3I Support (Bug 20)**: Add B3I (Band 6) frequency support in observation handler and EKF.
22. **GLONASS IFB Calibration (Bug 21)**: Add modeling or calibration parameters for GLONASS Inter-Frequency Biases (IFBs).
23. **Solid Earth Tide Step 2 (Bug 14)**: Implement frequency-dependent corrections in the diurnal band and permanent deformation removal.
24. **Doppler Jacobian Earth Rate (Bug 7)**: Incorporate Earth rotation rate `omega_ie` in derivative of Doppler with respect to attitude.
25. **Earth Tide Static Radius Approximation (Bug 4)**: Replace static equatorial radius with local geocentric radius `r_norm`.

## Interface Contracts & Layout
Code layout is defined by the existing gneiss workspace.
All implementation changes must happen in the source files, and unit tests must be added to the respective modules or in `#[cfg(test)]` blocks within the source files.
No changes to `.agents/` other than metadata and coordinator reports.
