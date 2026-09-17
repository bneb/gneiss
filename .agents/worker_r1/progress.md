# Progress — Worker R1 (15-State ESKF Benchmark Verification)

Last visited: 2026-09-12T22:48:30Z

## Status
- Unit tests: 17/17 passing in `estimators::eskf` (0.01s).
- Clippy passes with 0 warnings on library and benchmark binary.
- All files in `crates/gneiss-rtk/src/estimators/eskf/` verified strictly < 500 LOC, functions < 32 LOC, nesting < 3, 0 unwraps.
- Fixed stationary ZUPT detection by wiring wheel-speed CAN signal, reducing Forward Filter p50 from 6.08m to 3.96m.
- Currently generating AR-enabled GNSS fixes (`target/gnss_fixes_odaiba_ar.csv`, task-690) to evaluate ambiguity-fixed RTK performance.
