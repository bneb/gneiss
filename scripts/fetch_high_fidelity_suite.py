#!/usr/bin/env python3
import os
import sys
import urllib.request
import gzip
import shutil
import ssl
import subprocess
from datetime import datetime

try:
    import certifi
    _SSL_CTX = ssl.create_default_context(cafile=certifi.where())
    urllib.request.install_opener(
        urllib.request.build_opener(urllib.request.HTTPSHandler(context=_SSL_CTX))
    )
except ImportError:
    pass  # fall back to system CA bundle

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

    def fetch_broadcast_nav(self, year: int, doy: int, suffix: str = "n") -> str:
        yy = year % 100
        nav_file = os.path.join(self.target_dir, f"brdc{doy:03d}0.{yy:02d}{suffix}")
        if os.path.exists(nav_file):
            return nav_file

        urls = [
            f"https://geodesy.noaa.gov/corsdata/rinex/{year}/{doy:03d}/brdc{doy:03d}0.{yy:02d}{suffix}.gz",
            f"https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{year}/{doy:03d}/brdc{doy:03d}0.{yy:02d}{suffix}.gz",
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

        # Fallback: check if rover.nav exists in gsdc or rtkexplorer (GPS only)
        if suffix == "n":
            fallback = os.path.join(os.path.dirname(self.target_dir), "gsdc", "rover.nav")
            if os.path.exists(fallback):
                print(f"Using local fallback nav file: {fallback}")
                shutil.copyfile(fallback, nav_file)
                return nav_file

        raise RuntimeError(f"Could not retrieve broadcast navigation for Year {year}, DOY {doy} ({suffix})")

    def fetch_station_coordinate(self, station: str) -> str:
        coord_file = os.path.join(self.target_dir, f"{station.lower()}_14.coord.txt")
        if os.path.exists(coord_file):
            return coord_file
        url = f"https://geodesy.noaa.gov/corsdata/coord/coord_14/{station.lower()}_14.coord.txt"
        print(f"Downloading surveyed coordinates from {url}...")
        urllib.request.urlretrieve(url, coord_file)
        return coord_file

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
    glo_nav_file = fetcher.fetch_broadcast_nav(year, doy, suffix="g")

    print("==================================================")
    print("Fetching Network RTK Base Stations (multi-base)")
    print("==================================================")
    # Bay Area bases around rover P224 for multi-base network RTK (baselines 15-50 km)
    network_bases = ["p222", "slac", "p181", "ohln", "capo", "p225"]
    network_files = []
    for st in network_bases:
        try:
            obs = fetcher.fetch_station(st, year, doy)
            coord = fetcher.fetch_station_coordinate(st)
            network_files.append((st, obs, coord))
        except Exception as e:
            print(f"[WARN] Failed to fetch station {st}: {e}")

    print(f"\n[SUCCESS] Datasets staged successfully:")
    print(f"Base Station (P222):  {base_obs}")
    print(f"Rover Station (P224): {rover_obs}")
    print(f"Broadcast Nav (GPS):  {nav_file}")
    print(f"Broadcast Nav (GLO):  {glo_nav_file}")
    print(f"Network Bases ({len(network_files)}):")
    for st, obs, coord in network_files:
        print(f"  {st.upper():6s} obs={os.path.basename(obs)} coord={os.path.basename(coord)}")

if __name__ == "__main__":
    main()
