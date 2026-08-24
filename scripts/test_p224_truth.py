#!/usr/bin/env python3
"""Tests for p224_truth_2025 truth extraction utilities.

Assertions are grounded in independently verifiable facts:
- WGS84 on-equator prime meridian: llh(0,0,h) -> (a+h, 0, 0) exactly.
- Date tokens follow UNR's DDMMMYY convention (e.g. 09JUN25 = 2025-06-09).
"""

import math
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from p224_truth_2025 import llh_to_ecef, extract_rows, parse_ddmmmyy

FAILURES = []


def check(name, cond, detail=""):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" ({detail})" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def approx(a, b, tol):
    return abs(a - b) <= tol


def test_parse_ddmmmyy():
    check("09JUN25 -> 2025-06-09", parse_ddmmmyy("09JUN25") == (2025, 6, 9))
    check("01JAN99 -> 1999-01-01", parse_ddmmmyy("01JAN99") == (1999, 1, 1))
    check("31DEC26 -> 2026-12-31", parse_ddmmmyy("31DEC26") == (2026, 12, 31))
    try:
        parse_ddmmmyy("32JAN25")
        check("invalid day rejected", False)
    except ValueError:
        check("invalid day rejected", True)


def test_llh_to_ecef_prime_meridian_equator():
    a = 6378137.0
    x, y, z = llh_to_ecef(0.0, 0.0, 0.0)
    check("llh(0,0,0) -> (a,0,0)", approx(x, a, 1e-6) and approx(y, 0, 1e-6) and approx(z, 0, 1e-6),
          f"{x:.4f},{y:.4f},{z:.4f}")
    x, y, z = llh_to_ecef(0.0, 0.0, 100.0)
    check("equator height adds along X", approx(x, a + 100.0, 1e-6))


def test_llh_to_ecef_pole():
    # At the pole: |r| = b = a(1-f) at height 0, z-axis aligned
    f = 1 / 298.257223563
    b = 6378137.0 * (1 - f)
    x, y, z = llh_to_ecef(90.0, 0.0, 0.0)
    check("pole: x,y ~ 0", approx(x, 0, 1e-3) and approx(y, 0, 1e-3))
    check("pole: z == b", approx(z, b, 1e-3), f"{z:.4f} vs {b:.4f}")


def test_extract_rows_filters_and_medians(tmpdir):
    # Synthetic tenv3 rows: columns 21..23 hold lat/lon/h (0-based split).
    def row(ddmonyy, lat, lon, h, site="P224"):
        # real tenv3 tail: 23 tokens, llh at indices 20..22
        head = ["P224" if i == 0 else "0" for i in range(20)]
        head[1] = ddmonyy
        return f"{site} {' '.join(head[1:20])} {lat:.9f} {lon:.9f} {h:.4f}"

    lines = [
        row("15MAY25", 37.8638959, -122.2190591, 407.3500),
        row("09JUN25", 37.8638974, -122.2190658, 407.3551),
        row("10JUN25", 37.8638975, -122.2190659, 407.3549),
        row("11JUN25", 37.8638976, -122.2190660, 407.3553),
        row("05JUL25", 37.8638980, -122.2190670, 407.3560),
        # foreign station must be ignored
        row("09JUN25", 37.0, -122.0, 100.0, site="P181"),
    ]
    p = Path(tmpdir) / "mini.tenv3"
    p.write_text("\n".join(lines) + "\n")

    rows = extract_rows(p, "P224", (2025, 6))
    check("only June rows extracted", len(rows) == 3, str(len(rows)))
    check("rows sorted by day", [r[0] for r in rows] == [9, 10, 11])
    med_lat = sum(r[1] for r in rows) / len(rows)
    check("median latitude matches fixture",
          approx(med_lat, (37.8638974 + 37.8638975 + 37.8638976) / 3, 1e-12))


def main():
    with tempfile.TemporaryDirectory() as td:
        test_parse_ddmmmyy()
        test_llh_to_ecef_prime_meridian_equator()
        test_llh_to_ecef_pole()
        test_extract_rows_filters_and_medians(td)
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        return 1
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
