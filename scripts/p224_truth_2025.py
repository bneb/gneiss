#!/usr/bin/env python3
"""Extract truth coordinates from UNR IGS20 tenv3 time series.

P224 (the multi-GNSS benchmark rover) has no single published coordinate
file for 2025; UNR's daily IGS20 solution gives per-day lat/lon/h with
real scatter. Median over the target month -> truth ECEF for benchmark
scoring. Works identically for any CORS in the series.

Usage:
    python3 scripts/p224_truth_2025.py [tenv3-path] [site] [YYYY-MM]
Defaults: /tmp/ds2025/P224_IGS20.tenv3, P224, 2025-06
"""

import math
import statistics
import sys
from pathlib import Path

MONTHS = {m: i + 1 for i, m in enumerate(
    ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
     "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"])}

WGS84_A = 6378137.0
WGS84_E2 = 6.69437999014e-3


def parse_ddmmmyy(tok: str):
    """UNR date token DDMMMYY -> (year, month, day)."""
    dd, mon, yy = tok[:2], tok[2:5], tok[5:]
    if mon not in MONTHS or not dd.isdigit() or not yy.isdigit():
        raise ValueError(f"bad date token {tok!r}")
    day = int(dd)
    if not 1 <= day <= 31:
        raise ValueError(f"day out of range: {tok!r}")
    yy_i = int(yy)
    # RINEX/UNR two-digit-year pivot: 80-99 -> 19xx, 00-79 -> 20xx
    year = 1900 + yy_i if yy_i >= 80 else 2000 + yy_i
    return year, MONTHS[mon], day


def llh_to_ecef(lat_deg: float, lon_deg: float, h_m: float):
    """Geodetic (WGS84) -> ECEF metres."""
    la, lo = math.radians(lat_deg), math.radians(lon_deg)
    N = WGS84_A / math.sqrt(1 - WGS84_E2 * math.sin(la) ** 2)
    return (
        (N + h_m) * math.cos(la) * math.cos(lo),
        (N + h_m) * math.cos(la) * math.sin(lo),
        (N * (1 - WGS84_E2) + h_m) * math.sin(la),
    )


def extract_rows(path: Path, site: str, year_month: tuple):
    """Daily (day, lat, lon, h) rows for `site` within `year_month`, sorted."""
    y_t, m_t = year_month
    rows = []
    for line in path.open():
        p = line.split()
        if len(p) < 23 or p[0] != site:
            continue
        try:
            y, m, d = parse_ddmmmyy(p[1])
        except ValueError:
            continue
        if (y, m) == (y_t, m_t):
            # tenv3 row tail: ... corr_nu, latitude, longitude, height
            rows.append((d, float(p[20]), float(p[21]), float(p[22])))
    rows.sort()
    return rows


def median_llh(rows):
    return tuple(statistics.median(r[i] for r in rows) for i in (1, 2, 3))


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        Path("/tmp/ds2025/P224_IGS20.tenv3")
    site = sys.argv[2] if len(sys.argv) > 2 else "P224"
    y, m = (int(x) for x in (sys.argv[3].split("-") if len(sys.argv) > 3
                             else ["2025", "06"]))

    rows = extract_rows(path, site, (y, m))
    if not rows:
        print(f"no {site} epochs in {y}-{m:02d}")
        return 1
    lat, lon, h = median_llh(rows)
    sig_h = statistics.pstdev(r[3] for r in rows) * 1000
    x, y_, z = llh_to_ecef(lat, lon, h)
    print(f"{site}: epochs={len(rows)} "
          f"(day {rows[0][0]}..{rows[-1][0]} of {y}-{m:02d})")
    print(f"IGS20 median: lat={lat:.9f} lon={lon:.9f} h={h:.4f} "
          f"(sigma_h={sig_h:.1f} mm)")
    print(f"ECEF: [{x:.4f}, {y_:.4f}, {z:.4f}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
