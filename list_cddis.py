import os
import requests

token = os.environ.get("EARTHDATA_TOKEN")

session = requests.Session()
class BearerAuth(requests.auth.AuthBase):
    def __init__(self, token):
        self.token = token
    def __call__(self, r):
        if "urs.earthdata.nasa.gov" in r.url:
            r.headers["Authorization"] = f"Bearer {self.token}"
        return r

url = "https://cddis.nasa.gov/archive/gnss/products/2137/"
print(f"Fetching {url}...")
response = session.get(url, auth=BearerAuth(token), allow_redirects=True)
print("Status Code:", response.status_code)
html = response.text
lines = html.split('\n')
for line in lines:
    if "com" in line.lower() or "cod" in line.lower():
        print(line.strip())

