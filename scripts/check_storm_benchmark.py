#!/usr/bin/env python3
"""Regression guard for Severe Solar Storm / Cycle Slip & Outage Recovery benchmark.

Executes integration benchmark suite covering cycle slip recovery under severe
ionospheric disturbance and positioning continuity across satellite constellation outages.
"""

from __future__ import annotations

import re
import subprocess
import sys


def main() -> int:
    print("Running Storm / Cycle Slip & Outage Recovery benchmark suite...")
    cmd = [
        "cargo",
        "test",
        "-p",
        "gneiss-tests",
        "--",
        "test_post_process_cycle_slip_recovery",
        "test_post_process_outage_continuity",
        "--nocapture",
    ]

    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("FAIL: Storm/recovery benchmark execution failed\n", r.stdout, r.stderr)
        return 1

    out = r.stdout
    slip_match = re.search(r"Cycle Slip Post-Processed RTK:\s+RMS=([\d.]+)m", out)
    outage_match = re.search(r"Outage Post-Processed RTK:\s+RMS=([\d.]+)m", out)

    if not slip_match or not outage_match:
        print("FAIL: Could not parse storm/recovery benchmark output metrics\n", out)
        return 1

    slip_rms = float(slip_match.group(1))
    outage_rms = float(outage_match.group(1))

    print(f"Cycle Slip Recovery PPK: RMS = {slip_rms:.4f}m")
    print(f"Outage Continuity PPK:   RMS = {outage_rms:.4f}m")

    failures = []
    if slip_rms > 0.010:
        failures.append(f"Cycle slip recovery RMS {slip_rms:.4f}m > 0.010m")
    if outage_rms > 0.015:
        failures.append(f"Outage continuity RMS {outage_rms:.4f}m > 0.015m")

    if failures:
        print("\nREGRESSIONS DETECTED:")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("\nALL STORM & RECOVERY BENCHMARK CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
