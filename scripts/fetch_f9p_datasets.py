#!/usr/bin/env python3
"""Automated acquisition script for public u-blox ZED-F9P rover datasets.

Fetches and stages real-world F9P kinematic data from:
1. UrbanNav Hong Kong Medium-Urban (TST1)
2. UrbanNav Hong Kong Deep-Urban (Whampoa)
Complete with dual-frequency F9P observations, base station references,
and NovAtel SPAN-CPT ground truth.
"""

from __future__ import annotations

import gzip
import io
import math
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent / "datasets" / "urbannav"

DATASETS = {
    "tst1": {
        "name": "UrbanNav Hong Kong Medium-Urban (TST1)",
        "dest": BASE_DIR / "hk_tst1",
        "gnss_url": "https://www.dropbox.com/sh/2haoy68xekg95zl/AAAkcN4FwhFxkPY1lXsxbJrxa?dl=1",
        "gt_url": "https://www.dropbox.com/s/twsvwftucoytfpc/UrbanNav_TST_GT_raw.txt?dl=1",
        "f9p_pattern": "ublox.f9p.obs",
        "base_url": "https://rinex.geodetic.gov.hk/rinex3/2021/137/hksc/1s/HKSC00HKG_R_20211370200_01H_01S_MO.crx.gz",
        "nav_urls": [
            ("https://rinex.geodetic.gov.hk/rinex3/2021/137/hksc/HKSC00HKG_R_20211370000_01D_GN.rnx.gz", "base.nav"),
            ("https://rinex.geodetic.gov.hk/rinex3/2021/137/hksc/HKSC00HKG_R_20211370000_01D_CN.rnx.gz", "base_bds.nav"),
            ("https://rinex.geodetic.gov.hk/rinex3/2021/137/hksc/HKSC00HKG_R_20211370000_01D_EN.rnx.gz", "base_gal.nav"),
        ],
    },
    "whampoa": {
        "name": "UrbanNav Hong Kong Deep-Urban (Whampoa)",
        "dest": BASE_DIR / "hk_whampoa",
        "gnss_url": "https://www.dropbox.com/sh/7ox7718bzcjqtlf/AABH_Kjm65gHQ09K3antBRdua?dl=1",
        "gt_url": "https://www.dropbox.com/s/ej2mkue2w3r36s2/UrbanNav_whampoa_raw.txt?dl=1",
        "f9p_pattern": "ublox.f9p.obs",
        "base_url": "https://rinex.geodetic.gov.hk/rinex3/2021/141/hksc/1s/HKSC00HKG_R_20211410600_01H_01S_MO.crx.gz",
        "nav_urls": [
            ("https://rinex.geodetic.gov.hk/rinex3/2021/141/hksc/HKSC00HKG_R_20211410000_01D_GN.rnx.gz", "base.nav"),
            ("https://rinex.geodetic.gov.hk/rinex3/2021/141/hksc/HKSC00HKG_R_20211410000_01D_CN.rnx.gz", "base_bds.nav"),
            ("https://rinex.geodetic.gov.hk/rinex3/2021/141/hksc/HKSC00HKG_R_20211410000_01D_EN.rnx.gz", "base_gal.nav"),
        ],
    },
}

USER_AGENT = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36"


def dms_to_deg(d: float, m: float, s: float) -> float:
    sign = -1.0 if d < 0 else 1.0
    return sign * (abs(d) + m / 60.0 + s / 3600.0)


def llh_to_ecef(lat_deg: float, lon_deg: float, h_m: float) -> tuple[float, float, float]:
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    phi = math.radians(lat_deg)
    lam = math.radians(lon_deg)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)
    n = a / math.sqrt(1.0 - e2 * sin_phi * sin_phi)
    x = (n + h_m) * cos_phi * math.cos(lam)
    y = (n + h_m) * cos_phi * math.sin(lam)
    z = (n * (1.0 - e2) + h_m) * sin_phi
    return x, y, z


def download_with_progress(url: str, desc: str) -> bytes:
    print(f"[*] Downloading {desc}...")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req) as resp:
        data = resp.read()
    print(f"    Downloaded {len(data) / (1024 * 1024):.2f} MB")
    return data


def process_gt_file(raw_gt_text: str, out_csv_path: Path):
    lines = raw_gt_text.strip().splitlines()
    out_lines = [
        "GPS TOW (s), GPS Week, Latitude (deg), Longitude (deg), Ellipsoid Height (m), "
        "ECEF X (m), ECEF Y (m), ECEF Z (m), Roll (deg), Pitch (deg), Heading (deg), Quality"
    ]
    parsed_count = 0
    for line in lines:
        parts = line.strip().split()
        if len(parts) < 19 or parts[0] == "UTCTime":
            continue
        try:
            week = float(parts[1])
            tow = float(parts[2])
            lat = dms_to_deg(float(parts[3]), float(parts[4]), float(parts[5]))
            lon = dms_to_deg(float(parts[6]), float(parts[7]), float(parts[8]))
            h = float(parts[9])
            roll = float(parts[16])
            pitch = float(parts[17])
            heading = float(parts[18])
            q = int(float(parts[19])) if len(parts) > 19 else 1

            x, y, z = llh_to_ecef(lat, lon, h)
            out_lines.append(
                f"{tow:.3f}, {int(week)}, {lat:.9f}, {lon:.9f}, {h:.4f}, "
                f"{x:.4f}, {y:.4f}, {z:.4f}, {roll:.4f}, {pitch:.4f}, {heading:.4f}, {q}"
            )
            parsed_count += 1
        except (ValueError, IndexError):
            continue

    out_csv_path.write_text("\n".join(out_lines) + "\n")
    print(f"    Saved {parsed_count} ground truth epochs to {out_csv_path}")


def stage_dataset(key: str):
    spec = DATASETS[key]
    dest = spec["dest"]
    dest.mkdir(parents=True, exist_ok=True)
    print(f"\n=== Staging {spec['name']} ===")

    # 1. Download & Extract GNSS RINEX
    rover_obs_dest = dest / "rover_f9p.obs"
    if not rover_obs_dest.exists():
        gnss_bytes = download_with_progress(spec["gnss_url"], f"{key} GNSS RINEX Archive")
        zf = zipfile.ZipFile(io.BytesIO(gnss_bytes))
        f9p_file = None
        for name in zf.namelist():
            if spec["f9p_pattern"] in name and not name.endswith(".nmea"):
                f9p_file = name
                break
        if f9p_file:
            print(f"    Extracting {f9p_file} -> {rover_obs_dest.name}")
            with zf.open(f9p_file) as f_in, open(rover_obs_dest, "wb") as f_out:
                f_out.write(f_in.read())
        else:
            print(f"    [WARN] No F9P obs file matching pattern in zip. Contents: {zf.namelist()}")
    else:
        print(f"    {rover_obs_dest.name} already exists, skipping download.")

    # 2. Download & Parse Ground Truth
    ref_csv_dest = dest / "reference.csv"
    if not ref_csv_dest.exists():
        gt_bytes = download_with_progress(spec["gt_url"], f"{key} Ground Truth")
        gt_text = gt_bytes.decode("utf-8", errors="replace")
        process_gt_file(gt_text, ref_csv_dest)
    else:
        print(f"    {ref_csv_dest.name} already exists, skipping download.")

    # 3. Download & Decompress Base Station Observations (HKSC)
    base_obs_dest = dest / "base_hksc.obs"
    if not base_obs_dest.exists() and "base_url" in spec:
        base_gz_bytes = download_with_progress(spec["base_url"], f"{key} HKSC Base Station RINEX")
        crx_decompressed = gzip.decompress(base_gz_bytes)
        temp_crx = dest / "temp_base.crx"
        temp_crx.write_bytes(crx_decompressed)
        crx2rnx_bin = shutil.which("crx2rnx") or os.path.expanduser("~/.cargo/bin/crx2rnx")
        print(f"    Converting {temp_crx.name} with {crx2rnx_bin}...")
        subprocess.run([crx2rnx_bin, str(temp_crx)], check=True)
        temp_rnx = dest / "temp_base.rnx"
        if temp_rnx.exists():
            temp_rnx.rename(base_obs_dest)
            if temp_crx.exists():
                temp_crx.unlink()
            print(f"    Staged base station observations to {base_obs_dest.name}")
        else:
            print(f"    [WARN] crx2rnx did not produce {temp_rnx.name}")
    elif base_obs_dest.exists():
        print(f"    {base_obs_dest.name} already exists, skipping download.")

    # 4. Download Navigation Files
    if "nav_urls" in spec:
        for url, fname in spec["nav_urls"]:
            nav_dest = dest / fname
            if not nav_dest.exists():
                nav_gz_bytes = download_with_progress(url, f"{key} Navigation File ({fname})")
                nav_dest.write_bytes(gzip.decompress(nav_gz_bytes))
                print(f"    Staged broadcast navigation file to {fname}")
            else:
                print(f"    {fname} already exists, skipping download.")


def main():
    target = sys.argv[1] if len(sys.argv) > 1 else "tst1"
    if target in ("all", "both"):
        for k in DATASETS:
            stage_dataset(k)
    elif target in DATASETS:
        stage_dataset(target)
    else:
        print(f"Unknown dataset: {target}. Available: {list(DATASETS.keys())} or 'all'")
        sys.exit(1)
    print("\n[SUCCESS] Dataset acquisition complete.")


if __name__ == "__main__":
    main()
