import urllib.request
import re

url = "https://igs.bkg.bund.de/root_ftp/IGS/products/mgex/2137/"
req = urllib.request.Request(url)
try:
    with urllib.request.urlopen(req) as response:
        html = response.read().decode('utf-8')
        links = re.findall(r'href=[\'"]?([^\'" >]+)', html)
        for link in links:
            if "CAS0MGXRAP_2020359" in link:
                print(link)
except Exception as e:
    print(f"Failed to fetch directory: {e}")
