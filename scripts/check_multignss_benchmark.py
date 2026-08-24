#!/usr/bin/env python3
"""Regression guard for the multi-GNSS benchmark (2025 DOY160).

Runs the release binary against datasets/multignss_2025d160 with GPS+
Galileo systems and checks headline metrics against budgets from the
verified state. Complements check_network_benchmark.py (dataset A).
Exits non-zero on regression or unparseable output.
"""

import os
import re
from pathlib import Path
import subprocess
import sys

BIN = Path("target/release/eval_network_ppk")

# Budgets from verified GE run @ c8a85f0+ (see docs/NETWORK_RTK_NEXT_STEPS.md)
BUDGETS = {
    "P181": {"fix_min": 96.0, "h_p95_max": 150.0, "v_p95_max": 320.0},
    "P225": {"fix_min": 68.0, "h_p95_max": 260.0, "v_p95_max": 400.0},
    "P222": {"fix_min": 84.0, "h_p95_max": 300.0, "v_p95_max": 140.0},
}
NETWORK_FIX_MIN = 94.0


def main() -> int:
    env = dict(os.environ,
               GNEISS_DATASET="multi2025",
               GNEISS_SYSTEMS="GE")
    r = subprocess.run([str(BIN)], env=env, capture_output=True, text=True)
    if r.returncode != 0:
        print("FAIL: eval exited", r.returncode)
        return 1
    log = r.stdout

    failures = []

    def check(name, value, op, bound):
        if value is None:
            failures.append(f"{name}: unparseable")
            print(f"  FAIL {name}: unparseable")
            return
        ok = (value >= bound) if op == ">=" else (value <= bound)
        sym = ">=" if op == ">=" else "<="
        print(f"  {'ok' if ok else 'REGRESS'} {name:<34} {value:7.2f} {sym} {bound}")
        if not ok:
            failures.append(f"{name}: {value} {op} {bound}")

    for block in re.split(r"^=== ", log, flags=re.M):
        m = re.match(r"Smoothed RTK \[(\w+)\]", block)
        if not m:
            continue
        b = m.group(1)
        if b not in BUDGETS:
            continue
        bud = BUDGETS[b]
        fm = re.search(r"Fixed=\d+/\d+ \[([\d.]+)%\]", block)
        h = re.search(r"Horizontal Error:\s+p50=[\d.]+m,\s+p68=[\d.]+m,\s+p95=([\d.]+)m", block)
        v = re.search(r"Vertical Error:\s+p50=[+-][\d.]+m,\s+RMS=[\d.]+m,\s+p95=([\d.]+)m", block)
        if fm:
            check(f"{b} fix rate (%)", float(fm.group(1)), ">=", bud["fix_min"])
        if h:
            check(f"{b} h_p95 (mm)", float(h.group(1)) * 1000, "<=", bud["h_p95_max"])
        if v:
            check(f"{b} v_p95 (mm)", float(v.group(1)) * 1000, "<=", bud["v_p95_max"])

    m = re.search(r"=== NETWORK FUSED.*?Fixed=\d+/\d+ \[([\d.]+)%\]", log)
    if m:
        check("network fused fix rate (%)", float(m.group(1)), ">=", NETWORK_FIX_MIN)
    else:
        failures.append("network fused: missing")

    print()
    if failures:
        print("REGRESSIONS:", failures)
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
