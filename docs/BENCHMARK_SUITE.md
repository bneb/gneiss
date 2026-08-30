# Gneiss Benchmark Suite

This document outlines the standard real-world benchmark datasets assembled in `datasets/` to ensure Gneiss's performance across diverse challenging GNSS scenarios. 

## Benchmark Matrix

### 1. Profile A: Kinematic UAV / Drone Aerial Survey
*   **Location**: `datasets/profile_a_uav/`
*   **Description**: Aerial flight trajectory with high roll/pitch dynamics, camera shutter events, and lever arms.
*   **Objective**: Validates sensor fusion (IMU + GNSS), accurate timestamp interpolation for camera shutter events, and strict lever-arm compensation during high dynamics.
*   **Guard Script**: `scripts/check_kinematic_uav_benchmark.py`

### 2. Profile B: Space Weather & Solar Storm Day
*   **Location**: `datasets/profile_b_storm/`
*   **Description**: Real severe geomagnetic storm dataset (e.g., May 2024 solar storm) with high TEC (Total Electron Content) gradients, rapid phase scintillation, and severe cycle slips.
*   **Objective**: Validates cycle-slip repair algorithms, atmospheric gradient modeling, and dual-frequency ionosphere-free (IF) / geometry-free (GF) robust combinations.
*   **Guard Script**: `scripts/check_storm_benchmark.py`

### 3. Profile C: Global MGEX PPP / PPP-AR Network
*   **Location**: `datasets/profile_c_mgex/`
*   **Description**: Multi-continent distributed stations (equatorial, mid-latitude, polar) for worldwide base-free Precise Point Positioning (PPP) testing.
*   **Objective**: Benchmarks convergence time, PPP Ambiguity Resolution (PPP-AR) success rate, and global product ingestion (SP3/CLK/BIA).
*   **Guard Script**: `scripts/check_mgex_benchmark.py`

### 4. Profile D: Low-Cost Multi-Band Hardware
*   **Location**: `datasets/profile_d_f9p/`
*   **Description**: u-blox ZED-F9P + helical patch antenna data recorded under challenging environments like foliage and urban multipath.
*   **Objective**: Validates multipath mitigation, SNR-based weighting, and AR reliability on low-cost single/dual-band hardware.
*   **Guard Script**: `scripts/check_f9p_benchmark.py`

## Usage Instructions

To verify the engine against the benchmark suite, execute the guard scripts:

```bash
# Run all benchmark guards
./scripts/check_kinematic_uav_benchmark.py
./scripts/check_storm_benchmark.py
./scripts/check_mgex_benchmark.py
./scripts/check_f9p_benchmark.py
```

These scripts compile the core engine via `cargo build --release` (or mock the verification) and enforce strict pass/fail criteria on budgets for fix rates and error percentiles.

## Dataset Ingestion

Because GNSS files are large, they are represented via cropped samples in the `datasets/` folders. If you need to re-initialize the baseline test files, run:

```bash
python3 scripts/ingest_benchmark_datasets.py
```

This will link or mock the required `.ubx`, `.rtcm3`, and `.rnx` files needed by the testing pipeline while keeping dataset sizes reasonable (cropped to ~15-30 min or 24h at 30s intervals).

## Quality Invariants

*   **Zero Warnings**: All benchmark checks must execute without producing any internal compiler warnings or unhandled exceptions.
*   **Speed**: End-to-end evaluation must stay fast. Datasets should strictly be cropped to represent the challenging phenomena (e.g. max 30-min during peak solar storm).
