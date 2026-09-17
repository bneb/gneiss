# Test Infrastructure: Tier 1–4 Frontiers E2E Test Suite

## 1. Overview & Architecture
This document defines the test architecture, test tier hierarchy, execution commands, and verification criteria for the Gneiss GNSS/INS Tier-1 Commercial Frontiers test suite.

The test suite is built strictly as a requirement-driven, opaque-box validation framework covering all 19 features in `PROJECT.md § Feature Inventory` across four hierarchical tiers:
- **Tier 1 (Feature Coverage)**: $\ge 5$ test cases per feature in isolation, testing fundamental mechanics, kinematics, and algorithmic correctness.
- **Tier 2 (Boundary & Corner Cases)**: $\ge 5$ edge/boundary cases per feature testing numerical limits, zero/extreme inputs, missing resources, and degenerate geometry.
- **Tier 3 (Cross-Feature Combinations)**: 10 pairwise and multi-feature interaction tests validating integration contracts and state coupling.
- **Tier 4 (Real-World Application Scenarios)**: 5 comprehensive mission workflows validating end-to-end execution on realistic GNSS/INS datasets and synthetic flight/drive profiles.

## 2. Test Directory & File Layout
In strict compliance with `AGENTS.md` (< 500 LOC per file, < 32 LOC per function, < 3 nesting depth, 0 clippy warnings), the test suite is partitioned into modular units under `tests/tests/test_frontiers_e2e/`:

```
tests/
├── Cargo.toml
└── tests/
    ├── test_frontiers_e2e.rs                 # Integration test entry point (< 50 LOC)
    └── test_frontiers_e2e/
        ├── common.rs                         # Common oracles, geodesy math, fixtures (< 350 LOC)
        ├── tier1_eskf.rs                     # Tier 1: Features 1–6 (30 tests, < 400 LOC)
        ├── tier1_ppp.rs                      # Tier 1: Features 7–11 (25 tests, < 400 LOC)
        ├── tier1_vrs.rs                      # Tier 1: Features 12–16 (25 tests, < 400 LOC)
        ├── tier1_composite.rs                # Tier 1: Features 17–19 (15 tests, < 300 LOC)
        ├── tier2_eskf.rs                     # Tier 2: Features 1–6 boundaries (30 tests, < 400 LOC)
        ├── tier2_ppp.rs                      # Tier 2: Features 7–11 boundaries (25 tests, < 400 LOC)
        ├── tier2_vrs.rs                      # Tier 2: Features 12–16 boundaries (25 tests, < 400 LOC)
        ├── tier2_composite.rs                # Tier 2: Features 17–19 boundaries (15 tests, < 300 LOC)
        ├── tier3_pairwise.rs                 # Tier 3: Cross-feature pairwise interactions (10 tests, < 350 LOC)
        └── tier4_scenarios.rs                # Tier 4: End-to-end mission workflows (5 tests, < 350 LOC)
```

## 3. Runner Commands & Verification
The primary verification command for this test suite is:
```bash
cargo test -p gneiss-tests --test test_frontiers_e2e
```

To run individual tiers:
```bash
cargo test -p gneiss-tests --test test_frontiers_e2e tier1_
cargo test -p gneiss-tests --test test_frontiers_e2e tier2_
cargo test -p gneiss-tests --test test_frontiers_e2e tier3_
cargo test -p gneiss-tests --test test_frontiers_e2e tier4_
```

Lint and code quality verification:
```bash
cargo clippy -p gneiss-tests --test test_frontiers_e2e -- -D warnings
```

## 4. Coverage Criteria & Thresholds
1. **Total Test Cases**: $\ge 205$ unique requirement-driven test cases (95 Tier 1 + 95 Tier 2 + 10 Tier 3 + 5 Tier 4).
2. **Feature Coverage**: 100% of the 19 features in `PROJECT.md` covered with at least 5 isolated tests and 5 boundary tests each.
3. **Pass Rate**: 100% pass rate with zero failures and zero ignored tests.
4. **Code Quality**:
   - Every file strictly $< 500$ LOC.
   - Every function strictly $< 32$ LOC.
   - Maximum nesting depth strictly $< 3$ levels.
   - Zero clippy warnings.

## 5. Feature Inventory Mapping
| # | Feature Name | Tier 1 Module | Tier 2 Module | Tier 3 Interaction |
|---|--------------|---------------|---------------|-------------------|
| 1 | 15-State ESKF Formulation | `tier1_eskf.rs` | `tier2_eskf.rs` | `tier3_pairwise.rs` (ESKF + NHC) |
| 2 | Error-Quaternion Feedback | `tier1_eskf.rs` | `tier2_eskf.rs` | `tier3_pairwise.rs` (ESKF + RTS) |
| 3 | Online Bias Estimation | `tier1_eskf.rs` | `tier2_eskf.rs` | `tier3_pairwise.rs` (ESKF + NHC) |
| 4 | 15-State RTS Smoother | `tier1_eskf.rs` | `tier2_eskf.rs` | `tier3_pairwise.rs` (ESKF + RTS) |
| 5 | Coupled NHC & ZUPT | `tier1_eskf.rs` | `tier2_eskf.rs` | `tier3_pairwise.rs` (ESKF + NHC) |
| 6 | Odaiba Benchmark Target | `tier1_eskf.rs` | `tier2_eskf.rs` | `tier4_scenarios.rs` (Scenario 1) |
| 7 | Fast SINEX OSB Ingestion | `tier1_ppp.rs` | `tier2_ppp.rs` | `tier3_pairwise.rs` (OSB + AR) |
| 8 | Antenna PCO/PCV Corrections | `tier1_ppp.rs` | `tier2_ppp.rs` | `tier3_pairwise.rs` (PCV + Constellation) |
| 9 | Multi-Constellation PPP | `tier1_ppp.rs` | `tier2_ppp.rs` | `tier3_pairwise.rs` (PCV + Constellation) |
| 10 | Single-Differenced LAMBDA AR | `tier1_ppp.rs` | `tier2_ppp.rs` | `tier3_pairwise.rs` (OSB + AR) |
| 11 | Kinematic PPP-AR Benchmark | `tier1_ppp.rs` | `tier2_ppp.rs` | `tier4_scenarios.rs` (Scenario 2) |
| 12 | Multi-Station CORS Ingestion | `tier1_vrs.rs` | `tier2_vrs.rs` | `tier3_pairwise.rs` (CORS + DD Net) |
| 13 | Network Baseline Adjustment | `tier1_vrs.rs` | `tier2_vrs.rs` | `tier3_pairwise.rs` (CORS + DD Net) |
| 14 | Delaunay Atmospheric Models | `tier1_vrs.rs` | `tier2_vrs.rs` | `tier3_pairwise.rs` (Delaunay + VRS) |
| 15 | Localized VRS Synthesis | `tier1_vrs.rs` | `tier2_vrs.rs` | `tier3_pairwise.rs` (Delaunay + VRS) |
| 16 | Network RTK Benchmark | `tier1_vrs.rs` | `tier2_vrs.rs` | `tier4_scenarios.rs` (Scenario 3) |
| 17 | Tightly-Coupled PPP/INS | `tier1_composite.rs` | `tier2_composite.rs` | `tier3_pairwise.rs` (TC-PPP + ESKF) |
| 18 | Tightly-Coupled Network RTK/INS | `tier1_composite.rs` | `tier2_composite.rs` | `tier3_pairwise.rs` (TC-RTK + ESKF) |
| 19 | Composite Execution Tests | `tier1_composite.rs` | `tier2_composite.rs` | `tier3_pairwise.rs` (Multi-Mode Switch) |
