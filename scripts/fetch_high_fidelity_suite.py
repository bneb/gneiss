#!/usr/bin/env python3
import os
import sys
import urllib.request
import gzip
import shutil
import subprocess
from datetime import datetime

class HighFidelityDatasetFetcher:
    def __init__(self, target_dir: str):
        self.target_dir = os.path.abspath(target_dir)
        os.makedirs(self.target_dir, exist_ok=True)
        self.crx2rnx_bin = os.path.join(self.target_dir, "CRX2RNX")

    def _ensure_crx2rnx(self):
        if os.path.exists(self.crx2rnx_bin):
            return
        
        print("Compiling CRX2RNX from source...")
        src_url = "https://terras.gsi.go.jp/ja/crx2rnx/RNXCMP_4.2.0_src.tar.gz"
        tar_path = os.path.join(self.target_dir, "RNXCMP.tar.gz")
        urllib.request.urlretrieve(src_url, tar_path)
        
        subprocess.run(["tar", "-xzf", tar_path, "-C", self.target_dir], check=True)
        src_dir = os.path.join(self.target_dir, "RNXCMP_4.2.0_src", "source")
        
        compile_cmd = ["gcc", "-O2", "-o", self.crx2rnx_bin, "crx2rnx.c"]
        subprocess.run(compile_cmd, cwd=src_dir, check=True)

    def fetch_station(self, station: str, year: int, doy: int) -> str:
        yy = year % 100
        url = f"https://geodesy.noaa.gov/corsdata/rinex/{year}/{doy:03d}/{station.lower()}/{station.lower()}{doy:03d}0.{yy:02d}d.gz"
        filename = f"{station.lower()}{doy:03d}0.{yy:02d}d.gz"
        compressed_path = os.path.join(self.target_dir, filename)
        d_file = compressed_path[:-3]
        o_file = d_file[:-1] + "o"

        if not os.path.exists(o_file):
            print(f"Downloading {url}...")
            urllib.request.urlretrieve(url, compressed_path)
            
            print(f"Decompressing {filename}...")
            with gzip.open(compressed_path, 'rb') as f_in:
                with open(d_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            print(f"Converting {d_file} from Compact RINEX (Hatanaka) to RINEX...")
            self._ensure_crx2rnx()
            with open(o_file, 'w') as f_out:
                subprocess.run([self.crx2rnx_bin, d_file, "-"], stdout=f_out, check=True)
            
            if os.path.exists(d_file):
                os.remove(d_file)
            if os.path.exists(compressed_path):
                os.remove(compressed_path)

        return o_file

    def fetch_broadcast_nav(self, year: int, doy: int) -> str:
        yy = year % 100
        nav_file = os.path.join(self.target_dir, f"brdc{doy:03d}0.{yy:02d}n")
        if os.path.exists(nav_file):
            return nav_file

        urls = [
            f"https://geodesy.noaa.gov/corsdata/rinex/{year}/{doy:03d}/brdc{doy:03d}0.{yy:02d}n.gz",
            f"https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{year}/{doy:03d}/brdc{doy:03d}0.{yy:02d}n.gz",
        ]

        for url in urls:
            try:
                print(f"Attempting to download broadcast ephemeris from {url}...")
                gz_path = nav_file + ".gz"
                urllib.request.urlretrieve(url, gz_path)
                with gzip.open(gz_path, 'rb') as f_in:
                    with open(nav_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                if os.path.exists(gz_path):
                    os.remove(gz_path)
                return nav_file
            except Exception as e:
                print(f"Failed to download from {url}: {e}")

        # Fallback: check if rover.nav exists in gsdc or rtkexplorer
        fallback = os.path.join(os.path.dirname(self.target_dir), "gsdc", "rover.nav")
        if os.path.exists(fallback):
            print(f"Using local fallback nav file: {fallback}")
            shutil.copyfile(fallback, nav_file)
            return nav_file

        raise RuntimeError(f"Could not retrieve broadcast navigation for Year {year}, DOY {doy}")

def main():
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    out_dir = os.path.join(root, "datasets", "cors_short_baseline")
    fetcher = HighFidelityDatasetFetcher(out_dir)

    print("==================================================")
    print("Fetching Tier 1 Short-Baseline CORS Pair (P222 & P224)")
    print("==================================================")

    year = 2020
    doy = 135 # May 14, 2020

    base_obs = fetcher.fetch_station("p222", year, doy)
    rover_obs = fetcher.fetch_station("p224", year, doy)
    nav_file = fetcher.fetch_broadcast_nav(year, doy)

    print(f"\n[SUCCESS] Datasets staged successfully:")
    print(f"Base Station (P222):  {base_obs}")
    print(f"Rover Station (P224): {rover_obs}")
    print(f"Broadcast Nav:        {nav_file}")

if __name__ == "__main__":
    main()
