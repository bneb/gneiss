import urllib.request
import os
import gzip
import shutil

os.makedirs("datasets/wtzr_ppp_1224", exist_ok=True)

# 1. Download WTZR from CDDIS
url_obs = "https://cddis.nasa.gov/archive/gnss/data/daily/2020/359/20d/WTZR00DEU_R_20203590000_01D_30S_MO.crx.gz"
dest_obs = "datasets/wtzr_ppp_1224/WTZR00DEU_R_20203590000_01D_30S_MO.crx.gz"

token = os.environ.get("EARTHDATA_TOKEN", "eyJ0...")

print("Downloading WTZR observation...")
req = urllib.request.Request(url_obs)
req.add_header("Authorization", f"Bearer {token}")
try:
    with urllib.request.urlopen(req) as response:
        with open(dest_obs, 'wb') as f_out:
            shutil.copyfileobj(response, f_out)
    print("Done")
except Exception as e:
    print(f"Failed to download OBS: {e}")

# 2. Copy the BRDC from sample_1
shutil.copy("datasets/rtkexplorer/sample_1/f9p_ppp_1224/BRDC00IGS_R_20203590000_01D_MN.rnx", "datasets/wtzr_ppp_1224/")

# 3. Download CODE products from BKG (com21374.eph.Z, com21374.clk.Z)
base_url = "https://igs.bkg.bund.de/root_ftp/IGS/products/mgex/2137/"
files = ["com21374.eph.Z", "com21374.clk.Z", "com21374.bia.Z"]
for f in files:
    url = base_url + f
    dest = "datasets/wtzr_ppp_1224/" + f
    print(f"Downloading {url}...")
    try:
        urllib.request.urlretrieve(url, dest)
        os.system(f"gunzip -f {dest}")
        print("Done")
    except Exception as e:
        print(f"Failed: {e}")

