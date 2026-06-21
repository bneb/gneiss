# Scope: Milestone 1b — Resolve Remaining Tier 1 Bugs

## Architecture
This milestone targets four distinct components of the `gneiss` navigation engine:
1. `gneiss-rtk` Engine: Phase Wind-up correction in `crates/gneiss-rtk/src/engine/ppp.rs`. (ALREADY FIXED)
2. `gneiss-core` Ephemeris: Broadcast clock TGD correction in `crates/gneiss-core/src/ephemeris.rs`. (TO IMPLEMENT)
3. `gneiss-parsers` Rinex Clock: Outlier tolerance in precise clock gaps in `crates/gneiss-parsers/src/rinex_clk.rs`. (TO IMPLEMENT)
4. `gneiss-core` Atmosphere: GMF Legendre normalization in `crates/gneiss-core/src/atmosphere.rs`. (ALREADY FIXED)

There are no direct code dependencies between these four bug fixes, but they must be resolved in priority order as listed.

## Milestones
| # | Name | Scope | Dependencies | Status |
|---|---|---|---|---|
| 1 | Bug 18: Opposite Sign in Phase Wind-Up Correction | Fix `wup` sign correction in `ppp.rs` and add regression test `test_windup_sign_correct`. | None | DONE |
| 2 | Bug 15: Incorrect Broadcast Clock TGD Correction | Fix dual-frequency/iono-free clock TGD subtraction logic in `ephemeris.rs` and add regression test `test_broadcast_clock_tgd_correct`. | None | PLANNED |
| 3 | Bug 24: Outlier Tolerance in Precise Clock Gaps | Fix outlier tolerance by returning `None` instead of `Some(r1.bias)` in `rinex_clk.rs` when gap exceeds 900s, and add regression test `test_precise_clock_gap_tolerance`. | Bug 15 | PLANNED |
| 4 | Bug 6: GMF Legendre Unnormalized Polynomials | Implement fully normalized Associated Legendre Functions in `atmosphere.rs` and add regression test `test_gmf_legendre_normalization`. | Bug 24 | DONE |

## Interface Contracts
- No new interfaces are introduced. Existing function signatures and structures must remain compatible.
- All fixes must compile cleanly with `cargo build`.
- All tests must pass cleanly with `cargo test`.
