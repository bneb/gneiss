import urllib.request
import os
import gzip
import shutil

# EarthData Credentials are in .env, let's assume netrc works
files = [
    "GRG0MGXFIN_20203590000_01D_05M_ORB.SP3.gz",
    "GRG0MGXFIN_20203590000_01D_30S_CLK.CLK.gz",
    "GRG0MGXFIN_20203590000_01D_01D_OSB.BIA.gz"
]

base_url = "https://cddis.nasa.gov/archive/gnss/products/mgex/2137/"

for f in files:
    url = base_url + f
    dest = "datasets/rtkexplorer/sample_1/f9p_ppp_1224/" + f
    print(f"Downloading {url}...")
    try:
        os.system(f"curl -n -L -c .urs_cookies -b .urs_cookies {url} -o {dest}")
        print("Done")
        
        # unzip
        if os.path.exists(dest):
            with gzip.open(dest, 'rb') as f_in:
                with open(dest[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(dest)
    except Exception as e:
        print(f"Failed: {e}")
