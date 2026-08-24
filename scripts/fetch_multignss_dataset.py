#!/usr/bin/env python3
"""Fetch the modern multi-GNSS CORS benchmark dataset (2025-06-09, DOY 160).

Sources (all open, no auth):
  obs : NGS corsdata  rinex/YYYY/DDD/<site>/<site>DDD0.YYo.gz   (RINEX 2.11
        MIXED, 20 observables, GPS+GLONASS+Galileo for upgraded stations)
  nav : BKG IGS mirror root_ftp/IGS/obs/YYYY/DDD/<STATION>_R_..._MN.rnx.gz
        (any MGEX station works — broadcast ephemerides are globally
        identical; we take the first listed)                    RINEX 3.04

Stations (P224 rover + three cross/different-length bases, all upgraded
to multi-GNSS receivers, all NGS mm-grade coordinates available):
  P224 rover TRM59800.00 SCIT   P181 15.0 km NW  TRM59800.80 SCIT
  P225      TRM29659.00 SCIT    P222 38.0 km SSE TRM59800.80 SCIT
"""

import gzip
import re
import subprocess
import sys
from pathlib import Path

YEAR, DOY = 2025, 160
SITES = ["p224", "p181", "p225", "p222"]
ADJACENT_DAYS = [157, 158, 159, 161, 162, 163]  # sidereal/multipath repeats
OUT = Path("/tmp/ds2025")


def fetch(url: str, dest: Path, min_size: int = 1000) -> bool:
    if dest.exists() and dest.stat().st_size >= min_size:
        print(f"  cached {dest.name}")
        return True
    subprocess.run(["curl", "-s", "--max-time", "900", "-o", str(dest), url])
    ok = dest.exists() and dest.stat().st_size >= min_size
    print(f"  {'ok ' if ok else 'FAIL'} {dest.name} ({dest.stat().st_size if dest.exists() else 0} B)")
    return ok


def main() -> int:
    OUT.mkdir(exist_ok=True)
    ok = True
    days = [DOY] + ADJACENT_DAYS
    for doy in days:
        for s in SITES:
            ok &= fetch(
                f"https://geodesy.noaa.gov/corsdata/rinex/{YEAR}/{doy:03d}/{s}/"
                f"{s}{doy:03d}0.{YEAR%100}o.gz",
                OUT / f"{s}{doy:03d}0.{YEAR%100}o.gz",
            )
    # mixed broadcast nav from BKG mirror (first station listing)
    listing = subprocess.run(
        ["curl", "-s", "--max-time", "30",
         f"https://igs.bkg.bund.de/root_ftp/IGS/obs/{YEAR}/{DOY}/"],
        capture_output=True, text=True).stdout
    mn = re.findall(r'href="([A-Z]{4}\d{2}[A-Z]{3}_R_'
                    f"{YEAR}{DOY:03d}_01D_MN\\.rnx\\.gz)\"", listing)
    if not mn:
        print("FAIL no MN nav found in BKG listing")
        return 1
    ok &= fetch(f"https://igs.bkg.bund.de/root_ftp/IGS/obs/{YEAR}/{DOY}/{mn[0]}",
                OUT / "station_mixed_nav.rnx.gz", min_size=100_000)
    # validate nav contains Galileo
    try:
        c = set()
        with gzip.open(OUT / "station_mixed_nav.rnx.gz", "rt", errors="replace") as f:
            for line in f:
                if len(line) > 3 and line[0] == "E" and line[1].isdigit():
                    c.add("E"); break
        print("  nav has Galileo:", bool(c))
        ok &= bool(c)
    except Exception as e:
        print("  nav validation error:", e)
        ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
