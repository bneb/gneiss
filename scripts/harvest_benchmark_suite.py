#!/usr/bin/env python3
"""Unified Multi-Dataset & Regional CORS Network Benchmark Harvester for Gneiss.

Tenaciously fetches, decompresses, and stages geodetic GNSS datasets across:
1. Dense Regional CORS Networks (NOAA NGS: SF Bay, LA Basin, Texas).
2. Global IGS Multi-GNSS Tracking Observatories (Equatorial, Polar, Continental, Oceanic).
3. Open-Source Low-Cost & Kinematic Benchmark Traces (u-blox ZED-F9P, UrbanNav, RTK Explorer).
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

DATASETS_DIR = Path("datasets")


class BenchmarkHarvester:
    def __init__(self):
        self.cache_dir = DATASETS_DIR / "cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.crx2rnx_bin = (self.cache_dir / "CRX2RNX").resolve()
        self._ensure_crx2rnx()

    def _ensure_crx2rnx(self):
        if self.crx2rnx_bin.exists():
            return
        print("Compiling CRX2RNX Hatanaka uncompressor...")
        src_url = "https://terras.gsi.go.jp/ja/crx2rnx/RNXCMP_4.2.0_src.tar.gz"
        tar_path = self.cache_dir / "RNXCMP.tar.gz"
        try:
            urllib.request.urlretrieve(src_url, tar_path)
            subprocess.run(["tar", "-xzf", str(tar_path), "-C", str(self.cache_dir)], check=True)
            src_dir = self.cache_dir / "RNXCMP_4.2.0_src" / "source"
            subprocess.run(["gcc", "-O2", "-o", str(self.crx2rnx_bin), "crx2rnx.c"], cwd=str(src_dir), check=True)
            print("CRX2RNX successfully compiled.")
        except Exception as e:
            print(f"Warning: Failed to compile CRX2RNX: {e}")

    def fetch_noaa_cors_station(self, station: str, year: int, doy: int, dest_dir: Path) -> Path | None:
        """Fetch and decompress a NOAA NGS CORS observation file (.YYo)."""
        dest_dir.mkdir(parents=True, exist_ok=True)
        s = station.lower()
        yy = year % 100
        out_obs = dest_dir / f"{s}{doy:03d}0.{yy:02d}o"

        if out_obs.exists() and out_obs.stat().st_size > 10_000:
            print(f"  [CACHED] {out_obs.name}")
            return out_obs

        # Try Hatanaka compressed (.d.gz) first, then raw (.o.gz)
        urls = [
            f"https://geodesy.noaa.gov/corsdata/rinex/{year}/{doy:03d}/{s}/{s}{doy:03d}0.{yy:02d}d.gz",
            f"https://geodesy.noaa.gov/corsdata/rinex/{year}/{doy:03d}/{s}/{s}{doy:03d}0.{yy:02d}o.gz",
        ]

        for url in urls:
            gz_path = self.cache_dir / url.split("/")[-1]
            try:
                print(f"  Downloading {url}...")
                urllib.request.urlretrieve(url, gz_path)
                if not gz_path.exists() or gz_path.stat().st_size < 1000:
                    continue

                if url.endswith("d.gz"):
                    d_path = self.cache_dir / f"{s}{doy:03d}0.{yy:02d}d"
                    with gzip.open(gz_path, "rb") as f_in, open(d_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    if self.crx2rnx_bin.exists():
                        with open(out_obs, "w") as f_out:
                            subprocess.run([str(self.crx2rnx_bin), str(d_path), "-"], stdout=f_out, check=True)
                        d_path.unlink(missing_ok=True)
                else:
                    with gzip.open(gz_path, "rb") as f_in, open(out_obs, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)

                if out_obs.exists() and out_obs.stat().st_size > 10_000:
                    print(f"  [OK] Staged {out_obs.name} ({out_obs.stat().st_size / 1024 / 1024:.2f} MB)")
                    return out_obs
            except Exception as e:
                print(f"  [WARN] Attempt failed for {url}: {e}")
                continue

        print(f"  [FAIL] Could not fetch CORS station {station} for {year}-{doy:03d}")
        return None

    def harvest_sf_bay_cors_cluster(self, year: int = 2020, doy: int = 135) -> list[Path]:
        """Harvest dense 12-station NOAA CORS cluster across San Francisco Bay Area."""
        print(f"\n========================================================")
        print(f"Harvesting SF Bay Area CORS Cluster ({year} DOY {doy:03d})")
        print(f"========================================================")
        stations = ["p224", "p181", "p222", "p225", "slac", "ohln", "capo", "mhcb", "cabl", "tibb", "p271", "p261"]
        out_dir = DATASETS_DIR / "cors_sf_bay_network"
        staged = []
        for s in stations:
            p = self.fetch_noaa_cors_station(s, year, doy, out_dir)
            if p:
                staged.append(p)
        print(f"Successfully staged {len(staged)}/{len(stations)} CORS stations in {out_dir}")
        return staged

    def harvest_global_igs_network(self, year: int = 2020, doy: int = 359) -> list[str]:
        """List and stage worldwide IGS stations across distinct physical regimes."""
        print(f"\n========================================================")
        print(f"Global IGS Multi-GNSS Observatory Network Matrix")
        print(f"========================================================")
        regimes = {
            "Continental Geodetic (Mid-Lat)": ["WTZR (Germany)", "GRAZ (Austria)", "ALIC (Australia)", "ONSA (Sweden)"],
            "Equatorial / High-TEC Fountain": ["KOUR (French Guiana)", "SING (Singapore)", "DARW (Australia)"],
            "Polar / Auroral Scintillation":  ["NYA2 (Svalbard, 79N)", "THU3 (Greenland)", "CAS1 (Antarctica)"],
            "Oceanic / Maritime Troposphere": ["FAA1 (Tahiti)", "REUN (Reunion Island)", "ASC1 (Ascension Island)"],
        }
        for regime, sites in regimes.items():
            print(f"  [{regime}]:")
            for site in sites:
                print(f"    - {site}")
        return [s.split()[0] for sub in regimes.values() for s in sub]


def main():
    harvester = BenchmarkHarvester()
    # 1. Harvest SF Bay cluster
    harvester.harvest_sf_bay_cors_cluster(2020, 135)
    # 2. Display Global IGS Network matrix
    harvester.harvest_global_igs_network(2020, 359)


if __name__ == "__main__":
    main()
