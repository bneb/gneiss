import urllib.request
import os

files = [
    "com21374.dcb.Z",
    "com21374.bia.Z",
    "grg21374.clk.Z",
    "grg21374.sp3.Z"
]

base_url = "https://igs.bkg.bund.de/root_ftp/IGS/products/mgex/2137/"

for f in files:
    url = base_url + f
    dest = "datasets/rtkexplorer/sample_1/f9p_ppp_1224/" + f
    print(f"Downloading {url}...")
    try:
        urllib.request.urlretrieve(url, dest)
        print("Done")
        os.system(f"gunzip -f {dest}")
    except Exception as e:
        print(f"Failed: {e}")
