# Scope: Gneiss Tier 1 Bug Fixes

## Architecture
The gneiss navigation engine consists of several packages:
- `gneiss-core`: ephemeris calculation, atmospheric delay models (GMF/Legendre), broadcast clock correction.
- `gneiss-rtk`: RTK navigation engine including PPP math, predictor, IEKF state estimator, and phase wind-up.
- `gneiss-parsers`: parsers for GNSS metadata and observations, including RINEX clock files.

## Milestones
These bugs are fixed sequentially to avoid conflicts and verify correctness step-by-step.

| # | Name | Scope | Dependencies | Status |
|---|------|-------|-------------|--------|
| 1 | Bug 17: GLONASS Time Scale Discrepancy | Correct time scale conversion in gneiss-core ephemeris | None | DONE |
| 2 | Bug 1: Melbourne-Wübbena Dimensional Typo | Correct dimensional/unit typo in gneiss-rtk ppp_math | M1 | DONE |
| 3 | Bug 9: Sequential AR Covariance Mismatch | Fix covariance mismatch in gneiss-rtk ppp_iekf | M2 | DONE |
| 4 | Bug 2: Velocity-Attitude Transition Sign Mismatch | Fix transition matrix sign in gneiss-rtk predictor | M3 | DONE |
| 5 | Bug 18: Opposite Sign in Phase Wind-Up Correction | Correct wind-up sign in gneiss-rtk ppp | M4 | IN_PROGRESS |
| 6 | Bug 15: Incorrect Broadcast Clock TGD Correction | Fix TGD correction logic in gneiss-core ephemeris | M5 | PLANNED |
| 7 | Bug 24: Outlier Tolerance in Precise Clock Gaps | Fix gap handling / outlier tolerance in gneiss-parsers rinex_clk | M6 | PLANNED |
| 8 | Bug 6: GMF Legendre Unnormalized Polynomials | Correct Legendre normalization in gneiss-core atmosphere | M7 | PLANNED |

## Interface Contracts
No cross-module interface changes are expected for these bug fixes. The internal behavior will be corrected while maintaining existing public signatures.
