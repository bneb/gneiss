# Project: Gneiss Navigation Engine Bug Fixes

## Architecture
The Gneiss navigation engine is structured as a multi-crate Rust workspace:
- `gneiss-core`: General GNSS algorithms, troposphere/ionosphere atmospheric models, tides, wind-up corrections, broadcast ephemeris orbit/clock calculations.
- `gneiss-rtk`: High-precision PPP and RTK EKF estimation engine, measurement models, predictions, sequential ambiguity resolution (AR), and Melbourne-Wübbena (MW) combination calculations.
- `gneiss-parsers`: RINEX observation and precise clock/orbit (SP3/CLK) parser logic.

## Milestones

| # | Milestone Name | Scope (Bugs to Fix) | Dependencies | Status |
|---|----------------|---------------------|--------------|--------|
| 1 | **Milestone 1: Tier 1 Bugs** | Fix ranks 1 to 8 (Bugs 17, 1, 9, 2, 18, 15, 24, 6) | None | `IN_PROGRESS` (Bugs 17, 1, 9, 2, 18, 6 are DONE; Bug 15, 24 remaining) |
| 2 | **Milestone 2: Tier 2 & 3 Bugs** | Fix ranks 9 to 18 (Bugs 10, 25, 12, 23, 16, 5, 11, 8, 22, 3) | Milestone 1 | `PLANNED` (Bugs 25, 12 are DONE) |
| 3 | **Milestone 3: Tier 4 Features** | Implement/Stub ranks 19 to 25 (Bugs 13, 19, 20, 21, 14, 7, 4) | Milestone 2 | `PLANNED` |

## Detailed Bug Mapping

### Milestone 1: Tier 1 Bugs (Ranks 1–8)
1. **GLONASS Time Scale Discrepancy (Bug 17)**: `DONE` (Fixed, verified, and audited CLEAN)
2. **Melbourne-Wübbena Typo (Bug 1)**: `DONE` (Fixed, verified, and audited CLEAN)
3. **Sequential AR Covariance (Bug 9)**: `DONE` (Fixed, verified, and audited CLEAN)
4. **Velocity-Attitude Jacobian (Bug 2)**: `DONE` (Fixed, verified, with positive sign `+f_e_skew` preserved on disk as requested)
5. **Opposite Wind-Up Sign (Bug 18)**: `DONE` (Verified already fixed in working tree)
6. **Broadcast Clock TGD Correction (Bug 15)**: `PLANNED` (Stop subtracting TGD from broadcast clocks for dual-frequency/iono-free measurements)
7. **Outlier Precise Clock Gap (Bug 24)**: `PLANNED` (Return `None` when gap in precise clock records exceeds 900s limit)
8. **GMF Legendre Polynomials (Bug 6)**: `DONE` (Verified already fixed in working tree)

### Milestone 2: Tier 2 & 3 Bugs (Ranks 9–18)
9. **TOF Sat Position Clock Bias (Bug 10)**: `PLANNED` (Account for receiver clock bias in TOF calculation)
10. **Covariance Re-Init on Slip Omission (Bug 25)**: `DONE` (Verified already fixed in working tree)
11. **Receiver Antenna PCV Omission (Bug 12)**: `DONE` (Verified already fixed in working tree)
12. **Klobuchar Evaluated at Receiver (Bug 23)**: `PLANNED` (Evaluate Klobuchar model at IPP instead of receiver)
13. **Mismatched Galileo BGD (Bug 16)**: `PLANNED` (Use Galileo E1/E5b group delay correction instead of E1/E5a for E5b observations)
14. **GMF Troposphere Longitude Omission (Bug 5)**: `PLANNED` (Include longitude parameter in GMF spherical harmonics evaluator)
15. **Sat PCV Zenith-Dependent Correction Omission (Bug 11)**: `PLANNED` (Project and apply zenith-dependent PCV from ANTEX on line-of-sight)
16. **INS State AR Update Cutoff (Bug 8)**: `PLANNED` (Stop zeroing out Kalman gain for INS states in Narrowlane constraint update)
17. **Saastamoinen Dry Delay Pressure (Bug 22)**: `PLANNED` (Use actual/standard surface pressure instead of a constant sea-level pressure scaled by height)
18. **L2C Phase Shift Bias (Bug 3)**: `PLANNED` (Verify L2C phase shift behavior and avoid hardcoding static 0.25 shift)

### Milestone 3: Tier 4 Features (Ranks 19–25)
19. **Lack of Ocean Tide Loading (Bug 13)**: `PLANNED` (Parse BLQ files and implement OTL corrections)
20. **Missing Yaw Steering Models (Bug 19)**: `PLANNED` (Implement constellation-specific yaw steering laws)
21. **Beidou B3I Support (Bug 20)**: `PLANNED` (Add B3I/Band 6 frequency support)
22. **GLONASS IFB Calibration (Bug 21)**: `PLANNED` (Add modeling or calibration parameters for GLONASS IFBs)
23. **Solid Earth Tide Step 2 (Bug 14)**: `PLANNED` (Implement Step 2 Solid Earth Tide corrections)
24. **Doppler Jacobian Earth Rate (Bug 7)**: `PLANNED` (Incorporate Earth rotation rate in derivative of Doppler w.r.t attitude)
25. **Earth Tide Static Radius Approximation (Bug 4)**: `PLANNED` (Replace static equatorial radius with local geocentric radius)

## Interface Contracts & Layout
Code layout is defined by the existing gneiss workspace.
All implementation changes must happen in the source files, and unit tests must be added to the respective modules or in `#[cfg(test)]` blocks within the source files.
No changes to `.agents/` other than metadata and coordinator reports.
