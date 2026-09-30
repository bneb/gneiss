# Progress: Survey 2 — Temporal Frame & Epoch Alignment Safety (R2)

Last visited: 2026-09-25T21:12:00Z

## Status
- [x] Initialized DISPATCH.md and BRIEFING.md
- [x] Survey time representations across crates/gneiss-core and crates/gneiss-rtk
- [x] Identify ad-hoc float rounding or integer truncation for epoch matching (16+ sites cataloged)
- [x] Investigate leap-second handling (GPS vs UTC 18s offset, GLONASS leap seconds, BeiDou 14s ephemeris discrepancy)
- [x] Analyze IMU / GNSS synchronization (eval_odaiba_ins, backward.rs, streaming.rs)
- [x] Propose strictly typed epoch types & conversions (Epoch<Scale>, TimeDelta, EpochKey)
- [x] Identify all modules/files requiring refactoring (across core, parsers, fetch, rtk)
- [x] Write survey_r2_temporal.md and handoff.md
- [x] Verify all 789 workspace tests pass
- [ ] Send completion message to parent
