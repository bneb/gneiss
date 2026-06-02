import os
import urllib.request
import gzip
import shutil
import subprocess
from datetime import datetime

class CorsFetcher:
    def __init__(self, output_dir: str):
        self.output_dir = os.path.abspath(output_dir)
        os.makedirs(self.output_dir, exist_ok=True)
        self.crx2rnx_bin = os.path.join(self.output_dir, "CRX2RNX")

    def _get_doy(self, dt: datetime) -> int:
        return dt.timetuple().tm_yday

    def _build_ngs_url(self, station: str, dt: datetime) -> str:
        doy = self._get_doy(dt)
        year = dt.year
        yy = year % 100
        # https://geodesy.noaa.gov/corsdata/rinex/2020/135/p222/p2221350.20d.gz
        return f"https://geodesy.noaa.gov/corsdata/rinex/{year}/{doy:03d}/{station.lower()}/{station.lower()}{doy:03d}0.{yy:02d}d.gz"

    def _ensure_crx2rnx(self):
        if os.path.exists(self.crx2rnx_bin):
            return
        
        print("Compiling CRX2RNX from source...")
        src_url = "https://terras.gsi.go.jp/ja/crx2rnx/RNXCMP_4.2.0_src.tar.gz"
        tar_path = os.path.join(self.output_dir, "RNXCMP.tar.gz")
        urllib.request.urlretrieve(src_url, tar_path)
        
        subprocess.run(["tar", "-xzf", tar_path, "-C", self.output_dir], check=True)
        src_dir = os.path.join(self.output_dir, "RNXCMP_4.2.0_src", "source")
        
        # Compile CRX2RNX
        compile_cmd = ["gcc", "-O2", "-o", self.crx2rnx_bin, "crx2rnx.c"]
        subprocess.run(compile_cmd, cwd=src_dir, check=True)

    def fetch(self, station: str, dt: datetime) -> str:
        url = self._build_ngs_url(station, dt)
        filename = url.split("/")[-1]
        compressed_path = os.path.join(self.output_dir, filename)
        
        d_file = compressed_path[:-3] # remove .gz
        o_file = d_file[:-1] + "o" # change .20d to .20o
        
        if not os.path.exists(o_file):
            print(f"Downloading {url}...")
            urllib.request.urlretrieve(url, compressed_path)
            
            print(f"Decompressing {filename}...")
            with gzip.open(compressed_path, 'rb') as f_in:
                with open(d_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            print("Converting Hatanaka to RINEX...")
            self._ensure_crx2rnx()
            with open(o_file, 'w') as f_out:
                subprocess.run([self.crx2rnx_bin, d_file, "-"], stdout=f_out, check=True)
            
        return o_file

