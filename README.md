<p align="center">
  <img src="assets/logo.svg" alt="Gneiss Logo" width="200"/>
</p>

# Gneiss Navigation Engine

Gneiss is a tightly-coupled GNSS and Inertial Navigation System (INS) engine written in Rust. It fuses raw satellite observations (pseudorange, carrier phase, and Doppler) with high-rate inertial sensors (accelerometer and gyroscope) using an Extended Kalman Filter (EKF).

The engine is designed for robust operation in multi-path environments, utilizing empirical statistical methods for outlier rejection and variance scaling.

## Features

- **Tightly-Coupled Integration**: Direct fusion of raw GNSS observations and IMU data within the primary state vector.
- **Adaptive Estimation**: Implements Innovation-based Adaptive Estimation (IAE) and Median Absolute Deviation (MAD) RAIM to scale observation variances dynamically.
- **Ambiguity Resolution**: Uses the LAMBDA (Least-squares AMBiguity Decorrelation Adjustment) algorithm for carrier-phase integer ambiguity resolution across GPS, Galileo, BeiDou, and GLONASS.
- **Precise Point Positioning (PPP-AR)**: Supports PPP utilizing RTCM SSR streams with solid earth tides and phase wind-up physical modeling.
- **Protection Levels & ARAIM**: Employs Solution Separation to rigorously calculate Horizontal and Vertical Protection Levels (HPL/VPL) bounding faults to target integrity risks.
- **Hardware-Agnostic Calibrations**: Supports injection of custom Temperature-Calibrated IMU misalignments and Antenna Phase Center (APC) models.
- **Post-Processing (PPK)**: Supports forward-backward Rauch-Tung-Striebel (RTS) smoothing to produce continuous trajectories from static files.
- **Real-Time Streaming**: `#![no_std]` compatible core engine, with a `tokio`-based `live` subcommand for streaming data via local UART and NTRIP casters.

## Project Status & Roadmap

The engine's architecture provides a mathematical scaffold for multiple GNSS processing modes. Currently, the project is heavily focused on optimizing **RTK + INS** workflows for urban canyon and challenging multipath environments.

**Supported Modes:**
  - Single Point Positioning (SPP)
  - Real-Time Kinematic (RTK)
  - RTK + INS (Tightly-Coupled)
  - Precise Point Positioning (PPP)
  - PPP + INS

**Features Included:**
  - Ionosphere-Free Linear Combinations for multi-frequency correction.
  - Zenith Wet Delay (ZWD) tropospheric estimation.
  - Geophysical corrections (Solid Earth Tides, Satellite Phase Wind-Up).
  - CDDIS SP3 and precise clock (.clk) parsing via 10th-order Lagrange interpolation.
  - Clock Jump State Preservation algorithms to prevent EKF divergence during TCXO adjustments.

## Architecture

The engine operates on a causal, recursive filtering architecture:

```mermaid
graph TD
    subgraph Inputs
    A[Raw Satellite Data]
    C[Raw Inertial Data]
    S[RTCM SSR Stream]
    end

    subgraph "Pass 1: Auto-Calibration"
    P1[GNSS Intrinsics Profiler]
    P2[6-DOF Extrinsics Optimizer]
    P1 -->|Clock Variances| P2
    P2 -->|Lever Arm & Angles| B
    end

    A --> P1
    C --> P2

    subgraph "Pass 2: Execution Engine"
    B(Gneiss Engine)
    B --> D{Extended Kalman Filter}
    D -->|Float State| H[LAMBDA Ambiguity Resolution]
    H -->|Fixed Ambiguities| I{FFRT Validation}
    I -->|Pass| J[Apply Fix & Hold]
    I -->|Fail| D
    D --> K[ARAIM Solution Separation]
    end

    A --> B
    C --> B
    S --> B

    subgraph Outputs
    D --> E[Position Trajectory]
    D --> F[Calibrated Sensor Biases]
    D --> G[Attitude & Heading]
    K --> L[HPL & VPL Bounds]
    end
```

## Quickstart

Gneiss provides a command-line interface for testing and processing datasets. To run a tightly-coupled fusion pipeline with backward smoothing on a static dataset:

```bash
cargo run --release -p gneiss-cli -- process \
    --rover datasets/my_data/rover.obs \
    --base datasets/my_data/base.obs \
    --output trajectory.pos \
    --enable-imu-fusion \
    --enable-backward-smoothing \
    --lever-arm "0.1,0.0,-0.2"
```

To run the engine in real-time streaming mode using a serial port and an NTRIP base station:

```bash
cargo run --release -p gneiss-cli -- live \
    --port /dev/ttyACM0 \
    --baud 460800 \
    --ntrip-url rtk2go.com \
    --ntrip-mount MOUNT \
    --ntrip-user USER \
    --ntrip-pass PASS
```

### CLI Configuration (`process` subcommand)

- `--rover <PATH>`: **(Required)** Path to the rover's raw GNSS observations (RINEX `.obs`).
- `--base <PATH>`: Path to the base station's raw GNSS observations (RINEX `.obs`). Required for RTK/PPK. If omitted, the engine defaults to SPP (Single Point Positioning).
- `--output <PATH>`: **(Required)** Path for the resulting `.pos` trajectory file.
- `--config <PATH>`: Path to a JSON configuration file to override default EKF process noise and tuning parameters.
- `--enable-backward-smoothing`: Enables the Rauch-Tung-Striebel (RTS) backward smoother.

- `--lambda-ratio <FLOAT>`: Minimum ratio for the LAMBDA Partial Ambiguity Resolution (PAR) test (default: `3.0`).
- `--lambda-subset <INT>`: Minimum number of satellites required for PAR (default: `7`). 
- `--lever-arm <X,Y,Z>`: Translation vector (in meters) from the IMU center of navigation to the GNSS antenna phase center in the vehicle body frame.
- `--calibrate`: **(New)** Enables the automated 2-pass calibration pipeline. In Pass 1, the engine dynamically estimates GNSS receiver clock variance profiles and uses a Nelder-Mead simplex optimizer to auto-detect the 6-DOF IMU mounting angles and lever arms. In Pass 2, it executes the tight-coupling with the detected parameters. This is highly recommended for cold-starts without a-priori lever arm measurements.
- `--calibrate-imu`: Enables dynamic state-estimation of IMU mounting rotations relative to the vehicle frame.
- `--raim-outlier-m <FLOAT>`: SPP Receiver Autonomous Integrity Monitoring (RAIM) threshold in meters.
- `--chi-square-pr <FLOAT>`: EKF Chi-Square threshold for pseudorange measurement rejection.
- `--chi-square-cp <FLOAT>`: EKF Chi-Square threshold for carrier phase measurement rejection.

## Accuracy Benchmarks
In deep urban canyons (e.g., Tokyo Shinjuku), the tight coupling engine produces sub-3m accuracy from an entirely blind cold-start using the `--calibrate` mode. If ground-truth lever arm parameters are precisely known, the engine can achieve sub-1.5m accuracy.

| Configuration | Median 3D Error | 95% 3D Error |
|---------------|-----------------|--------------|
| Commercial Baseline (NovAtel) | 3.60 m | >10 m |
| `gneiss` (Auto-Calibrated) | 2.93 m | 9.29 m |
| `gneiss` (Manually Tuned) | 1.34 m | 3.77 m |

## Documentation

For technical implementation details, see the following documents:
- [Architecture Details](./ARCHITECTURE.md)
- [Benchmark Methodology](./BENCHMARKS.md)
- [Precise Point Positioning (PPP-AR) Explained](./docs/PPP_AR_EXPLAINED.md)

## Workspace Structure

| Crate | Purpose |
| :--- | :--- |
| [`gneiss-core`](./crates/gneiss-core) | Core data structures, physical constants, and geometric models. |
| [`gneiss-geodesy`](./crates/gneiss-geodesy) | Earth reference frames, datum transformations, and gravity models. |
| [`gneiss-parsers`](./crates/gneiss-parsers) | Decoders for standard positioning formats (RINEX, UBX, RTCM3). |
| [`gneiss-rtk`](./crates/gneiss-rtk) | The Extended Kalman Filter, mechanization, and ambiguity resolution logic. |
| [`gneiss-ntrip`](./crates/gneiss-ntrip) | Asynchronous networking client for RTK corrections. |
| [`gneiss-cli`](./bin/gneiss-cli) | Command-line interface for dataset processing and real-time execution. |
