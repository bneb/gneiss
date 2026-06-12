# Gneiss Benchmarking Guide

This document outlines the standard operating procedure for evaluating the Gneiss navigation engine. We rely on a centralized benchmarking infrastructure to ensure empirical, regression-proof development across various challenging environments (urban canyons, static datasets, varying hardware).

## 1. Acquiring Datasets

Gneiss tests against several major public datasets. Due to their massive size (raw `.ubx`, `.obs`, and `.rtcm3` files), these are **not** committed to the repository.

Instead, use the fetching scripts to initialize the directory structures and download what is available automatically:

```bash
./scripts/fetch_datasets.sh
```

**Note on Massive Datasets:** Datasets like **UrbanLoco**, **TEX-CUP**, **WHU-Smartphone**, and **smartLoc** are hosted on Google Drive or university servers. You must download the raw sequences manually into their respective `datasets/` subdirectories.

## 2. The Unified Orchestrator

All evaluation runs through a single Python entry point: `scripts/benchmark.py`. 

This orchestrator compiles the `gneiss-cli`, runs the specified datasets through the engine, evaluates the resulting `.pos` files against the dataset's ground truth, and generates Markdown reports.

### Running Benchmarks

**1. Gneiss Baseline Suite**
Runs all GNSS and INS processing modes (SPP, RTK, PPP, Loosely/Tightly coupled) for Gneiss and generates `BENCHMARKS.md`.
```bash
./scripts/benchmark.py --suite gneiss
```

**2. RTKLIB Comparison Suite**
Runs Gneiss side-by-side against RTKLIB (demo5) across baseline modes (SPP, RTK Forward, RTK Smoothed) to compare position error metrics. Generates `COMPARISON.md`.
```bash
./scripts/benchmark.py --suite rtklib
```

**3. Full 18-Grid Matrix**
Runs the comprehensive matrix combining Base Modes $\times$ INS Coupling $\times$ Filter Direction against RTKLIB baselines. Generates `BENCHMARKS_MATRIX.md`.
```bash
./scripts/benchmark.py --suite matrix
```

### Useful Flags

- `--dataset <NAME>`: Filter the run to a specific dataset (e.g., `--dataset UrbanLoco`). Useful for rapid iteration on a single sequence.
- `--dry-run`: Prints the underlying `gneiss-cli` and `rnx2rtkp` commands without executing them.
- `--eval-only`: Skips the heavy processing step and strictly regenerates the Markdown reports from existing `.pos` files.

## 3. Adding New Datasets

To integrate a new dataset into the automated pipeline:

1. Download the `.obs`, `.nav`, and ground truth `.csv/.pos` files into a subdirectory within `datasets/`.
2. Open `scripts/benchmark.py`.
3. Locate the `DATASETS` dictionary at the top of the file.
4. Add your new dataset as a dictionary block:

```python
    "My New Dataset": {
        "rover": "datasets/my_new_dataset/rover.obs",
        "base": "datasets/my_new_dataset/base.obs",
        "nav": "datasets/my_new_dataset/brdc.nav",
        "truth": "datasets/my_new_dataset/ground_truth.csv",
        "gneiss_config": "datasets/my_new_dataset/config.json", # Optional EKF overrides
        "is_static": False, # Set to True if the rover is stationary
    },
```

The orchestrator will execute it in the next run.
