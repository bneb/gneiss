#!/usr/bin/env python3
"""Generate eval-ready constant truth + station coords for the 2025 DOY160
multi-GNSS benchmark from UNR IGS20 medians.

Outputs (datasets/multignss_2025d160/):
  p224_truth.pos   2880 rows @30 s, constant ECEF (existing Truth::new parser)
  station_coords.json  {site: [x,y,z]} for rover+bases

Self-checks: row count, first/last timestamps, coordinate invariance.
"""

import json
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from p224_truth_2025 import extract_rows, median_llh, llh_to_ecef  # noqa: E402

SERIES_DIR = Path("/tmp/ds2025")
OUT_DIR = Path("datasets/multignss_2025d160")
SITES = ["P224", "P181", "P222", "P225"]
EPOCHS = 2880
DT = 30.0


def fail(msg):
    print(f"  [FAIL] {msg}")
    sys.exit(1)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    coords = {}
    for site in SITES:
        src = SERIES_DIR / f"{site}_IGS20.tenv3"
        rows = extract_rows(src, site, (2025, 6))
        if len(rows) < 15:
            fail(f"{site}: only {len(rows)} June epochs")
        lat, lon, h = median_llh(rows)
        x, y, z = llh_to_ecef(lat, lon, h)
        coords[site] = [round(x, 4), round(y, 4), round(z, 4)]
        sig_h = statistics.pstdev(r[3] for r in rows) * 1000
        print(f"  [{site}] sigma_h={sig_h:.1f} mm")

    # truth file for rover P224
    x, y, z = coords["P224"]
    out = OUT_DIR / "p224_truth.pos"
    lines = []
    # GPS week/time-of-week: 2025-06-09 00:00 == week 2322, tow 345600? Verify:
    # 2025-01-05 was start of GPS week 2322? Simpler: emit UTC-style date rows;
    # existing parser reads columns by position (date x y z Q ns).
    import datetime as dt
    t0 = dt.datetime(2025, 6, 9, 0, 0, 0)
    for k in range(EPOCHS):
        t = t0 + dt.timedelta(seconds=k * DT)
        lines.append(f"{t.strftime('%Y/%m/%d %H:%M:%S.000000')} "
                     f"{x:.4f} {y:.4f} {z:.4f} 1 20")
    out.write_text("\n".join(lines) + "\n")

    # --- self-checks -----------------------------------------------------
    got = out.read_text().strip().split("\n")
    if len(got) != EPOCHS:
        fail(f"truth rows {len(got)} != {EPOCHS}")
    if not got[0].startswith("2025/06/09 00:00:00"):
        fail(f"first timestamp wrong: {got[0][:30]}")
    if not got[-1].startswith("2025/06/09 23:59:30"):
        fail(f"last timestamp wrong: {got[-1][:34]}")
    f0 = got[0].split()
    if abs(float(f0[2]) - x) > 1e-3 or abs(float(f0[4]) - z) > 1e-3:
        fail("coordinates drifted between generation and file")
    (OUT_DIR / "station_coords.json").write_text(json.dumps(coords, indent=2))

    print(f"  wrote {out} ({len(got)} rows)")
    print(f"  wrote {OUT_DIR / 'station_coords.json'}")
    print("ALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
