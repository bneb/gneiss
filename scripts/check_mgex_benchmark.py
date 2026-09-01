#!/usr/bin/env python3
"""Regression guard for Global MGEX / IGS Standalone Float PPP Tracking benchmark.

Executes integration benchmark suite on real IGS tracking data (Wettzell WTZR
and Alice Springs ALIC) using CODE precise orbit (SP3) and clock (CLK) products.
"""

from __future__ import annotations

import re
import subprocess
import sys


def main() -> int:
    print("Running Global MGEX / IGS Standalone PPP benchmark suite...")
    cmd = [
        "cargo",
        "test",
        "-p",
        "gneiss-tests",
        "--",
        "test_real_rinex_wtzr_float_ppp_convergence",
        "test_real_rinex_alic_float_ppp_convergence",
        "--nocapture",
    ]

    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("FAIL: MGEX/IGS PPP benchmark execution failed\n", r.stdout, r.stderr)
        return 1

    out = r.stdout
    wtzr_match = re.search(r"WTZR Real RINEX Float PPP:\s+p50=([\d.]+)m,\s+min=[\d.]+m,\s+last_horiz=[\d.]+m,\s+last_vert=([\d.]+)m", out)
    alic_match = re.search(r"ALIC Real RINEX Float PPP:\s+p50=([\d.]+)m", out)

    if not wtzr_match or not alic_match:
        print("FAIL: Could not parse MGEX/IGS PPP benchmark output metrics\n", out)
        return 1

    wtzr_p50 = float(wtzr_match.group(1))
    wtzr_last_v = float(wtzr_match.group(2))
    alic_p50 = float(alic_match.group(1))

    print(f"WTZR Standalone Float PPP: p50 = {wtzr_p50:.4f}m, final vert dU = {wtzr_last_v:.4f}m")
    print(f"ALIC Standalone Float PPP: p50 = {alic_p50:.4f}m (from 3m perturbed seed)")

    failures = []
    if wtzr_p50 > 0.60:
        failures.append(f"WTZR horizontal p50 {wtzr_p50:.4f}m > 0.60m")
    if wtzr_last_v > 0.10:
        failures.append(f"WTZR final vertical error {wtzr_last_v:.4f}m > 0.10m")
    if alic_p50 > 1.50:
        failures.append(f"ALIC horizontal p50 {alic_p50:.4f}m > 1.50m")

    if failures:
        print("\nREGRESSIONS DETECTED:")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("\nALL GLOBAL MGEX / IGS PPP BENCHMARK CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
