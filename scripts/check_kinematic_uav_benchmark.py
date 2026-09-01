#!/usr/bin/env python3
"""Regression guard for High-Dynamic Kinematic & UAV Trajectory benchmarks.

Executes integration benchmark suite covering high-angular-rate circular vehicle
motion and open-sky kinematic PPK trajectories against exact ground truth.
"""

from __future__ import annotations

import re
import subprocess
import sys


def main() -> int:
    print("Running Kinematic & UAV benchmark integration suite...")
    cmd = [
        "cargo",
        "test",
        "-p",
        "gneiss-tests",
        "--",
        "test_simulation_high_dynamic_circular_kinematic_sub_centimeter",
        "test_post_process_open_sky_kinematic_sub_centimeter",
        "--nocapture",
    ]

    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("FAIL: Kinematic benchmark execution failed\n", r.stdout, r.stderr)
        return 1

    out = r.stdout
    circ_match = re.search(r"Simulation High-Dynamic Circular Kinematic:\s+RMS=([\d.]+)m", out)
    sky_match = re.search(r"Open-sky Kinematic Post-Processed RTK:\s+RMS=([\d.]+)m,\s+Fixed=(\d+)/(\d+)", out)

    if not circ_match or not sky_match:
        print("FAIL: Could not parse kinematic benchmark output metrics\n", out)
        return 1

    circ_rms = float(circ_match.group(1))
    sky_rms = float(sky_match.group(1))
    sky_fixed = int(sky_match.group(2))
    sky_total = int(sky_match.group(3))

    print(f"High-Dynamic Circular Kinematic: RMS = {circ_rms:.4f}m")
    print(f"Open-Sky Kinematic PPK: RMS = {sky_rms:.4f}m, Fixed = {sky_fixed}/{sky_total}")

    failures = []
    if circ_rms > 0.010:
        failures.append(f"High-dynamic circular RMS {circ_rms:.4f}m > 0.010m")
    if sky_rms > 0.010:
        failures.append(f"Open-sky kinematic RMS {sky_rms:.4f}m > 0.010m")
    if sky_fixed < 28:
        failures.append(f"Open-sky fixed epochs {sky_fixed} < 28")

    if failures:
        print("\nREGRESSIONS DETECTED:")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("\nALL KINEMATIC & UAV BENCHMARK CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
