#!/usr/bin/env python3
"""Regression guard for Low-Cost u-blox ZED-F9P Kinematic benchmark.

Runs eval_qinertia_ppk on the RTK Explorer F9P dataset (u-blox ZED-F9P rover
with TMG2 base) and validates headline accuracy against real measured budgets.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

BIN = Path("target/release/eval_qinertia_ppk")


def ensure_binary() -> int:
    if not BIN.exists():
        r = subprocess.run(["cargo", "build", "--release", "--bin", "eval_qinertia_ppk"])
        if r.returncode != 0:
            print("FAIL: Failed to build eval_qinertia_ppk")
            return 1
    return 0


def main() -> int:
    if ensure_binary() != 0:
        return 1

    print("Running real F9P benchmark evaluation via eval_qinertia_ppk...")
    r = subprocess.run([str(BIN)], capture_output=True, text=True)
    if r.returncode != 0:
        print("FAIL: eval_qinertia_ppk exited with error", r.returncode)
        return 1

    log = r.stdout
    f9p_match = re.search(
        r"Evaluating Dataset: RTK Explorer F9P.*?=== Qinertia-Grade Smoothed PPK Solution \(N=\d+, Fixed=(\d+)/(\d+) \[([\d.]+)%\]\) ===\s+Horizontal Error:\s+p50=([\d.]+)m,\s+p68=[\d.]+m,\s+p95=[\d.]+m,\s+RMS=([\d.]+)m",
        log,
        re.DOTALL,
    )

    if not f9p_match:
        print("FAIL: Unable to parse F9P Qinertia-grade solution from eval output")
        return 1

    fix_count = int(f9p_match.group(1))
    total_count = int(f9p_match.group(2))
    fix_rate = float(f9p_match.group(3))
    p50 = float(f9p_match.group(4))
    rms = float(f9p_match.group(5))

    print(f"F9P Kinematic Results: N={total_count}, Fixed={fix_count}/{total_count} ({fix_rate:.1f}%), p50={p50:.3f}m, RMS={rms:.3f}m")

    failures = []
    # Budgets from verified execution on real u-blox ZED-F9P dataset
    if fix_rate < 80.0:
        failures.append(f"Fix rate {fix_rate:.1f}% < 80.0%")
    if p50 > 0.25:
        failures.append(f"Horizontal p50 {p50:.3f}m > 0.250m")
    if rms > 0.30:
        failures.append(f"Horizontal RMS {rms:.3f}m > 0.300m")

    if failures:
        print("\nREGRESSIONS DETECTED:")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("\nALL F9P BENCHMARK CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
