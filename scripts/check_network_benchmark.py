#!/usr/bin/env python3
"""Regression guard for the multi-base network RTK benchmark.

Runs the full-day eval (release binary - build it first) and checks
headline metrics against budgets verified on the
network-rtk-long-baseline branch. Exits non-zero on any regression.

Usage:
    cargo build --release --bin eval_network_ppk
    scripts/check_network_benchmark.py

Requires datasets/cors_short_baseline/ (see scripts/fetch_datasets.sh).
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

BIN = Path("target/release/eval_network_ppk")


def run_eval() -> str:
    if not BIN.exists():
        sys.exit(
            f"FAIL: {BIN} not found - "
            "run 'cargo build --release --bin eval_network_ppk' first"
        )
    result = subprocess.run(
        [str(BIN)], capture_output=True, text=True, check=True
    )
    return result.stdout


def iter_blocks(lines: list[str]):
    """Yield (header, body-lines) per '===' delimited section."""
    header: str | None = None
    body: list[str] = []
    for line in lines:
        if line.startswith("=== "):
            if header is not None:
                yield header, body
            header, body = line, []
        elif header is not None:
            body.append(line)
    if header is not None:
        yield header, body


def find_block(log: str, product: str) -> list[str] | None:
    for header, body in iter_blocks(log.splitlines()):
        if product in header:
            return body
    return None


def line_value(block: list[str], prefix: str, index: int = 0) -> float | None:
    """Nth float on the line starting with <prefix>, else None."""
    for line in block:
        if line.startswith(prefix):
            nums = re.findall(r"-?\d+\.\d+", line[len(prefix):])
            return float(nums[index]) if index < len(nums) else None
    return None


def fixed_only_p50(block: list[str]) -> float | None:
    """p50 from the '[...] fixed-only horizontal p50: X' summary line."""
    for line in block:
        if "fixed-only" in line:
            m = re.search(r"p50:\s*(\d+\.\d+)", line)
            return float(m.group(1)) if m else None
    return None


def fix_rate(block_header: str) -> float | None:
    m = re.search(r"\[([\d.]+)%\]", block_header)
    return float(m.group(1)) if m else None


def main() -> int:
    log = run_eval()

    def net_block() -> list[str]:
        b = find_block(log, "NETWORK FUSED")
        if b is None:
            sys.exit("FAIL: NETWORK FUSED block missing from eval output")
        return b

    failures: list[str] = []

    def check(name: str, value: float | None, op: str, bound: float) -> None:
        if value is None:
            failures.append(f"{name}: value not parsed")
            print(f"  FAIL {name:<44} <no value>")
            return
        ok = value <= bound if op == "<=" else value >= bound
        sym = "<=" if op == "<=" else ">="
        print(f"  {'ok' if ok else 'REGRESS':<6} {name:<42} {value:7.3f} {sym} {bound}")
        if not ok:
            failures.append(f"{name}: {value} {op} {bound}")

    net = net_block()
    check("network fused horizontal p50 (m)",
          line_value(net, "Horizontal Error:", 0), "<=", 0.04)
    check("network fused horizontal RMS (m)",
          line_value(net, "Horizontal Error:", 3), "<=", 0.12)
    check("network fused vertical RMS (m)",
          line_value(net, "Vertical Error:", 1), "<=", 0.15)

    for base, budget in [("P181", 0.05), ("P222", 0.12), ("SLAC", 0.15)]:
        blk = find_block(log, f"Smoothed RTK [{base}]")
        check(f"{base} smoothed fixed-only p50 (m)",
              None if blk is None else fixed_only_p50(blk), "<=", budget)

    # Fix-rate floors reflect HONEST quality accounting on the
    # long-baseline branch: disputed epochs are flagged float rather than
    # counted as fixed, so rates sit below the legacy (dishonest) numbers.
    # Fix-rate floors reflect honest accounting with ZWD estimation:
    # fewer false fixes = lower rate but higher quality.
    for base, floor in [("OHLN", 85.0), ("P181", 88.0), ("SLAC", 40.0)]:
        hdr = next((h for h in log.splitlines()
                    if h.startswith("=== ") and f"[{base}]" in h and "Smoothed" in h), "")
        check(f"{base} smoothed fix rate (%)",
              fix_rate(hdr) if hdr else None, ">=", floor)

    if failures:
        print("\nREGRESSIONS DETECTED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("\nALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
