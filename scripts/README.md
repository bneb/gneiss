# gneiss-scripts

Automation, verification, and dataset acquisition scripts for the Gneiss positioning engine.

## Key Tooling Categories

### 1. Primary CI Regression Guards
- **`check_network_benchmark.py`**:
  Regression guard for multi-base Network RTK (VRS) on NOAA CORS data. Enforces headline 3D error percentiles ($p_{50} \le 0.04\text{ m}$) and carrier phase fix rates ($\ge 96.5\%$).
  ```bash
  python3 scripts/check_network_benchmark.py --smoke
  ```
- **`check_multignss_benchmark.py`**:
  Regression guard for GPS + Galileo multi-constellation Network RTK across CORS stations P181 and P222.
  ```bash
  python3 scripts/check_multignss_benchmark.py --smoke
  ```

### 2. Scenario Profile Verification
- **`check_kinematic_uav_benchmark.py`**: Validates UAV aerial photogrammetry shutter interpolation and dynamic lever arms.
- **`check_storm_benchmark.py`**: Validates cycle-slip repair and ionospheric delay handling during geomagnetic solar storms.
- **`check_mgex_benchmark.py`**: Validates global MGEX tracking station convergence.
- **`check_f9p_benchmark.py`**: Benchmarks low-cost multi-band hardware under foliage and multipath.

### 3. Public Dataset Ingestion
- **`fetch_cors.py`**: Downloads RINEX 3 observation and navigation files from NOAA CORS archives.
- **`fetch_multignss_dataset.py`**: Retrieves multi-constellation tracking station data from CDDIS and BKG.
- **`fetch_precise_products.py`**: Downloads IGS/CODE final precise orbits (`.sp3`), high-rate clocks (`.clk`), ANTEX antenna calibrations (`.atx`), and SINEX bias files (`.bia`).
- **`fetch_f9p_datasets.py`**: Fetches low-cost GNSS receiver benchmark logs.

### 4. Analysis & Comparison
- **`run_rtklib_comparison.py`**: Executes differential comparisons against RTKLIB 2.4.3 b34 across standardized baselines.
- **`evaluate_ppp_ar.py`**: Computes ambiguity resolution statistics and convergence curves.
