# gneiss-tests

`tests` contains the workspace integration test suite, benchmark regression harnesses, and synthetic simulation scenarios.

## 4-Tier Testing Architecture

Integration tests in `tests/tests/test_frontiers_e2e/` are structured into four progressive validation tiers:

### Tier 1: Component Contract & Boundary Tests
- Verifies mathematical contracts and input validation for core subsystems:
  - `tier1_vrs.rs`: Delaunay triangulation, IPP geometry, and spatial interpolation bounds.
  - `tier1_eskf.rs`: 15-state error quaternion mechanization, covariance propagation, and matrix symmetries.
  - `tier1_ppp.rs`: Melbourne-Wübbena combination, ionosphere-free residual formation, and OSB ingestion.
  - `tier1_composite.rs`: Pipeline lifecycle transitions and time synchronization.

### Tier 2: Differential Invariant & Property Tests
- Enforces physical conservation laws and estimator invariants:
  - `tier2_vrs.rs`: Baseline collapse invariants and atmospheric gradient continuity.
  - `tier2_eskf.rs`: Stationary zero-velocity updates (ZUPT) and accelerometer bias observability.
  - `tier2_ppp.rs`: Integer ambiguity resolution discrimination ratios ($R \ge 3.0$).
  - `tier2_composite.rs`: State continuity across measurement drops.

### Tier 3: Pairwise Subsystem Integration Tests
- Validates cross-crate interactions:
  - `tier3_pairwise.rs`: Coupled Network RTK / INS and PPP / INS state updates under simulated sensor drift.

### Tier 4: Mission Scenario & Stress Tests
- End-to-end mission verification under challenging operational environments:
  - `tier4_scenarios.rs`: Multi-minute urban canyon drives with bridge overpass outages, severe ionospheric cycle slips, and antenna lever-arm dynamics.

## Benchmark Matrix Harnesses

- `src/benchmark_matrix.rs`: Validates solution percentiles and fix rates across NOAA CORS baselines.
- `src/f9p_rover_benchmark.rs`: Low-cost receiver multipath rejection and convergence testing.
- `src/post_process_simulation.rs`: Forward/backward RTS smoothing and multi-pass initialization verification.
- `src/inertial_outage_simulation.rs`: Evaluates dead-reckoning drift during GNSS outages.

## Running Tests

Execute the full workspace test suite:
```bash
cargo test --workspace
```

Execute integration tests specifically:
```bash
cargo test -p tests
```
