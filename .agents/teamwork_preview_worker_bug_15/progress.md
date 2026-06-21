## Current Status
Last visited: 2026-06-21T03:05:00-07:00
- [x] Initialized BRIEFING.md and task setup
- [x] Investigate files: `crates/gneiss-core/src/ephemeris.rs`, `crates/gneiss-rtk/src/engine/ppp.rs`, `crates/gneiss-rtk/src/estimators/spp.rs`
- [x] Implement `position_iono_free` and variant methods
- [x] Add unit test `test_broadcast_clock_tgd_correct`
- [x] Update PPP engine to use `position_iono_free` and remove manual `tgd` addition
- [x] Update SPP estimator to use `position_iono_free` conditionally
- [x] Run build and test suite
- [x] Write handoff report and notify orchestrator
