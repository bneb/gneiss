import urllib.request
import os
import gzip
import shutil

os.makedirs("datasets/wtzr_ppp_1224", exist_ok=True)

# 1. Download WTZR from BKG (avoids CDDIS Earthdata auth)
url_obs = "https://igs.bkg.bund.de/root_ftp/IGS/obs/2020/359/WTZR00DEU_R_20203590000_01D_30S_MO.crx.gz"
dest_obs = "datasets/wtzr_ppp_1224/WTZR00DEU_R_20203590000_01D_30S_MO.crx.gz"

print("Downloading WTZR observation...")
try:
    urllib.request.urlretrieve(url_obs, dest_obs)
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

