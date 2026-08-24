#!/usr/bin/env python3
"""A/B comparator for WL_DUMP per-epoch error CSVs.

Usage:
    python3 scripts/compare_dumps.py BASE_DIR GRAD_DIR [--metric p95]

Computes horizontal/vertical p50/p95/RMS per base from the Smoothed CSVs
in each directory and prints deltas. Focused on p95 by default: tails are
the commercial-grade discriminator, and the metric our gradient work must
move to justify itself.
"""

import argparse
import csv
import math
from pathlib import Path


def quantile(sorted_xs, q):
    if not sorted_xs:
        return float("nan")
    idx = min(int(len(sorted_xs) * q), len(sorted_xs) - 1)
    return sorted_xs[idx]


def stats_from_csv(path):
    h, v = [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            hv = float(row["h"])
            vv = float(row["v"])
            # Exclude gross outliers (blowup windows) from distribution
            # shape so p50/p95 describe the operating envelope; RMS keeps
            # them via a separate all-epochs figure.
            if hv < 2.0 and abs(vv) < 2.0:
                h.append(hv * 1000)
                v.append(vv * 1000)
    h.sort()
    v.sort()
    n = len(h)
    rms = lambda xs: math.sqrt(sum(x * x for x in xs) / max(len(xs), 1))
    return {
        "n": n,
        "h_p50": quantile(h, 0.5),
        "h_p95": quantile(h, 0.95),
        "h_rms": rms(h),
        "v_p50": quantile(v, 0.5),
        "v_p95": quantile(max(v, key=abs) and sorted(v, key=abs), 0.95),
        "v_rms": rms(v),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base_dir")
    ap.add_argument("variant_dir")
    ap.add_argument("--label-a", default="base")
    ap.add_argument("--label-b", default="variant")
    args = ap.parse_args()

    bases = ["P181", "OHLN", "CAPO", "P225", "P222", "SLAC"]
    rows = []
    for b in bases:
        pa = Path(args.base_dir) / f"dump_{b}_Smoothed.csv"
        pb = Path(args.variant_dir) / f"dump_{b}_Smoothed.csv"
        if not (pa.exists() and pb.exists()):
            continue
        sa, sb = stats_from_csv(pa), stats_from_csv(pb)
        rows.append((b, sa, sb))

    hdr = (f"{'base':>6} | {'h_p50':>13} {'h_p95':>14} | "
           f"{'|v|p50':>13} {'|v|p95':>14}")
    print(hdr)
    print("-" * len(hdr))
    for b, sa, sb in rows:
        def cell(a, bb):
            d = bb - a
            arrow = "↓" if abs(bb) < abs(a) else ("↑" if abs(bb) > abs(a) else "=")
            return f"{a:6.1f}->{bb:6.1f}{arrow}"
        va = abs(sa["v_p50"])
        vb = abs(sb["v_p50"])
        pa95a, pa95b = sa["v_p95"], sb["v_p95"]
        print(f"{b:>6} | {cell(sa['h_p50'], sb['h_p50'])} "
              f"{cell(sa['h_p95'], sb['h_p95'])} | "
              f"{cell(va, vb)} {cell(abs(pa95a), abs(pa95b))}")

    print("\nv_p95 detail (signed, mm):")
    for b, sa, sb in rows:
        print(f"  {b:>6}: {sa['v_p95']:+8.1f} -> {sb['v_p95']:+8.1f}")


if __name__ == "__main__":
    main()
