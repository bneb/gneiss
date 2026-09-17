# Gneiss Benchmark Suite

This document defines the automated benchmark and verification suite used to evaluate positioning performance, fix reliability, and regression invariance across GNSS processing engines.

## 1. Primary CI Regression Guards

The CI pipeline enforces strict performance gates on accuracy and ambiguity resolution (AR) fix rates using two dual-mode guard scripts:

### `scripts/check_network_benchmark.py`
- **Dataset**: NOAA CORS network (P181, P222, P225, SLAC, OHLN) over short to medium regional baselines.
- **Engine**: Network RTK / PPK VRS solver.
- **Pass Criteria**:
  - Smoke mode (`--smoke`): Median 3D position error $p_{50} \le 0.04\text{ m}$, carrier phase fix rate $\ge 96.5\%$.
  - Full evaluation: Complete multi-baseline VRS network convergence across all epochs.

### `scripts/check_multignss_benchmark.py`
- **Dataset**: 2025 DOY 160 multi-GNSS (GPS + Galileo) CORS network observations.
- **Engine**: Multi-constellation Network RTK / PPK.
- **Pass Criteria**:
  - Smoke mode (`--smoke`):
    - Single baseline P181 fix rate $\ge 97.5\%$.
    - Single baseline P222 fix rate $\ge 86.0\%$.
    - Network-fused VRS fix rate $\ge 96.5\%$.
  - Full evaluation: Multi-frequency ambiguity resolution stability across constellation boundaries.

---

## 2. Production Evaluation Harnesses

Five dedicated binaries located in `crates/gneiss-rtk/src/bin/` provide end-to-end evaluation against authoritative commercial baselines and ground truth:

| Binary | Scenario & Dataset | Reference Truth | Primary Metrics |
|---|---|---|---|
| `eval_network_ppk` | NOAA CORS regional network (P181, P225, P222, SLAC, OHLN) | CORS monument coordinates (NAD83(2011)) | Horizontal/3D error ($p_{50}, p_{95}$), VRS fix % |
| `eval_ppp` | Multi-GNSS kinematic & static PPP (RTK Explorer F9P, WTZR, ALIC) | Canada Geodetic Service CSRS-PPP & RTK PPK | Median 3D error vs CSRS-PPP ($0.262\text{ m}$), Calibrated RTK error ($1.7\text{ cm}$) |
| `eval_odaiba_ins` | Urban canyon Tokyo Odaiba (12,398 epochs, 10 Hz GNSS + 50 Hz IMU) | NovAtel SPAN-CPT tactical-grade ground truth | 15-state ESKF horizontal error ($p_{50} = 2.31\text{ m}, \text{RMS} = 4.64\text{ m}$) |
| `eval_qinertia_ppk` | NGS short baseline (112.5 m) and RTK Explorer kinematic drive | NGS published coordinates and RTK PPK | Baseline length agreement, integer closure |
| `eval_f9p_rover` | Low-cost u-blox ZED-F9P across 6 urban/suburban environments | High-grade multi-frequency reference station | Ambiguity resolution ratio, multipath rejection |

Run an evaluation binary via Cargo:
```bash
cargo run --release --bin eval_ppp
cargo run --release --bin eval_network_ppk
cargo run --release --bin eval_odaiba_ins
```

---

## 3. Integration Profile Guards

Four supplementary profile scripts validate specific environmental stresses:

- **Profile A: Kinematic UAV / Drone Survey** (`scripts/check_kinematic_uav_benchmark.py`):
  Validates attitude dynamics, camera shutter time-tag interpolation, and antenna lever-arm compensation.
- **Profile B: Space Weather & Solar Storm** (`scripts/check_storm_benchmark.py`):
  Evaluates cycle-slip repair, ionospheric delay gradients, and geometry-free (GF) / ionosphere-free (IF) phase combinations during geomagnetic disturbances (e.g. May 2024 storm).
- **Profile C: Global MGEX PPP** (`scripts/check_mgex_benchmark.py`):
  Tests global tracking station convergence (WTZR, ALIC) using precise IGS orbits (SP3), satellite clocks (CLK), and SINEX code/phase biases (BIA).
- **Profile D: Low-Cost Multi-Band Hardware** (`scripts/check_f9p_benchmark.py`):
  Validates low-cost patch antenna performance under foliage and multipath.

---

## 4. Benchmark Dataset Directory Layout

```
datasets/
├── noaa_cors/         # NOAA CORS RINEX 3 observation and broadcast navigation files
├── rtkexplorer/       # u-blox ZED-F9P kinematic rover and base station RINEX logs
├── urbannav/          # Tokyo Odaiba UrbanNav dataset (u-blox UBX, IMU CSV, NovAtel SPAN truth)
├── igs/               # Precise ephemeris products: IGS/COD SP3, CLK, BIA, ATX, Bernese DCB
├── profile_a_uav/     # UAV aerial survey logs and event markers
├── profile_b_storm/   # Space weather storm observation snippets
├── profile_c_mgex/    # Multi-GNSS global station observation files
└── profile_d_f9p/     # Challenging environment low-cost GNSS logs
```

---

## 5. Verification Commands

Run the primary CI guards:
```bash
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke
```

Run profile benchmarks:
```bash
python3 scripts/check_kinematic_uav_benchmark.py
python3 scripts/check_storm_benchmark.py
python3 scripts/check_mgex_benchmark.py
python3 scripts/check_f9p_benchmark.py
```
