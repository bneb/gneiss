# Scope: Milestone 2 — Tier 2 & 3 Bugs

## Architecture
Gneiss is a GNSS/INS navigation and RTK positioning engine written in Rust. It is divided into several crates:
- `gneiss-core`: Contains core definitions for orbits, observations, signals, coordinates, time, atmospheric models (Klobuchar, troposphere models), and basic math/geodesy.
- `gneiss-geodesy`: Reference frames, antenna calibrations (ANTEX), Helmert transformations.
- `gneiss-parsers`: Rinex, SP3, ANTEX, UBX parser implementations.
- `gneiss-rtk`: The main positioning engine, including SPP, PPP, RTK estimators, Kalman filters (IEKF), and measurements.

## Milestones
| # | Name | Scope | Dependencies | Status |
|---|------|-------|-------------|--------|
| 1 | Bug 16: Mismatched Galileo BGD | Correct Galileo BGD correction selection logic for E5b observations (band 7) | None | IN_PROGRESS |
| 2 | Bug 23: Klobuchar Receiver Eval | Evaluate Klobuchar ionosphere model at Ionospheric Pierce Point (IPP) instead of receiver position | M1 | PLANNED |
| 3 | Bug 5: GMF Troposphere Longitude | Include receiver longitude in spherical harmonics annual evaluator for GMF | M2 | PLANNED |
| 4 | Bug 10: TOF Sat Position Clock Bias | Account for receiver clock bias `cdt_r` in signal Time-of-Flight calculation | M3 | PLANNED |
| 5 | Bug 11: Sat PCV Zenith-Dependent | Project and apply zenith-dependent PCV from ANTEX along the line-of-sight | M4 | PLANNED |
| 6 | Bug 8: INS State AR Update Cutoff | Stop zeroing out the Kalman gain for INS states in the Narrowlane constraint update | M5 | PLANNED |
| 7 | Bug 22: Saastamoinen Dry Delay | Use actual/standard surface pressure instead of constant sea-level pressure scaled by height | M6 | PLANNED |
| 8 | Bug 3: L2C Phase Shift Bias | Avoid hardcoding static 0.25 shift if receiver tracks in-phase, or handle properly | M7 | PLANNED |

## Interface Contracts
- No new cross-module interfaces are required; changes will conform to existing function signatures and structs in `gneiss-core` and `gneiss-rtk`.
