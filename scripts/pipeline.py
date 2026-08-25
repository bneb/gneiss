#!/usr/bin/env python3
"""Automated improvement pipeline for Gneiss PPK engine.

Implements the RUNBOOK.md procedure as a single command:
  1. Verify baseline green (guards + tests)
  2. Report current metrics
  3. Exit non-zero if baseline broken

Usage:
  python3 scripts/pipeline.py          # verify + report
  python3 scripts/pipeline.py --fix    # verify + report + auto-fix trivial issues
"""

import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).parent.parent


def run(cmd, timeout=300):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout, cwd=ROOT)
    return r.returncode, r.stdout, r.stderr


def step_verify_baseline():
    """Verify all tests pass and both guards are green."""
    print("=== STEP 1: VERIFY BASELINE ===")
    
    rc, out, _ = run("cargo build --workspace 2>&1")
    errors = out.count("error[") if rc != 0 else 0
    print(f"  build: {'OK' if rc == 0 else f'FAILED ({errors} errors)'}")
    if rc != 0:
        print(out[-500:] if len(out) > 500 else out)
        return False
    
    rc, out, _ = run("cargo test --workspace 2>&1 | grep 'test result'", timeout=600)
    total_pass = sum(int(line.split("passed")[0].strip().split()[-1]) 
                     for line in out.split('\n') if 'passed' in line and '0 failed' in line)
    any_fail = 'failed' in out and out.count('failed') > out.count('0 failed')
    print(f"  tests: {total_pass} passed, {'FAILURES DETECTED' if any_fail else 'all green'}")
    if any_fail:
        return False
    
    for guard, name in [("scripts/check_network_benchmark.py", "dataset A"),
                        ("scripts/check_multignss_benchmark.py", "dataset B")]:
        rc, out, _ = run(f"python3 {guard}", timeout=600)
        ok = "PASSED" in out
        print(f"  guard {name}: {'GREEN' if ok else 'RED'}")
        if not ok:
            print(out[-300:] if len(out) > 300 else out)
            return False
    
    return True


def step_report_metrics():
    """Report current accuracy metrics."""
    print("\n=== STEP 2: CURRENT METRICS ===")
    rc, out, _ = run(
        "GNEISS_DATASET=multi2025 GNEISS_SYSTEMS=GE ./target/release/eval_network_ppk",
        timeout=600
    )
    if rc == 0:
        for line in out.split('\n'):
            if 'Smoothed RTK' in line or 'NETWORK FUSED' in line:
                print(f"  {line.strip()}")
            elif 'Horizontal Error' in line or 'Vertical Error' in line:
                print(f"    {line.strip()}")
    
    rc, out, _ = run("./target/release/eval_network_ppk", timeout=600)
    if rc == 0:
        # dataset A metrics
        for m in __import__('re').finditer(
            r'=== Smoothed RTK \[(\w+)\].*?Horizontal Error:\s+p50=([\d.]+)m.*?RMS=([\d.]+)m', out
        ):
            b, p50, rms = m.group(1), float(m.group(2))*1000, float(m.group(3))*1000
            print(f"    A:{b}: p50={p50:.0f}mm rms={rms:.0f}mm")


def main():
    print("=" * 60)
    print("Gneiss PPK Engine — Automated Pipeline")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    if not step_verify_baseline():
        print("\nBASELINE BROKEN — fix before proceeding.")
        return 1
    
    step_report_metrics()
    
    print("\n=== READY FOR IMPROVEMENT ===")
    print("Pick top item from docs/PROJECT_STATUS.md roadmap.")
    print("Implement with TDD. Re-run this script after.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
