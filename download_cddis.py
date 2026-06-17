import os
import requests

with open(".env") as f:
    for line in f:
        if line.startswith("EARTHDATA_TOKEN="):
            token = line.split("=", 1)[1].strip()

if not token:
    print("NO TOKEN")
    exit(1)

url = "https://cddis.nasa.gov/archive/gnss/products/2137/com21374.sp3.Z"

session = requests.Session()
def add_auth_header(request, *args, **kwargs):
    if "urs.earthdata.nasa.gov" in request.url:
        request.headers["Authorization"] = f"Bearer {token}"
    return request

session.hooks['response'] = lambda r, *args, **kwargs: r

print(f"Downloading {url}...")
# Follow redirects, but inject auth header on URS domain
# Requests allows auth tuple or a custom auth object
class BearerAuth(requests.auth.AuthBase):
    def __init__(self, token):
        self.token = token
    def __call__(self, r):
        if "urs.earthdata.nasa.gov" in r.url:
            r.headers["Authorization"] = f"Bearer {self.token}"
        return r

response = session.get(url, auth=BearerAuth(token), allow_redirects=True)
print("Final URL:", response.url)
print("Status Code:", response.status_code)

if response.status_code == 200:
    with open("com21374.sp3.Z", "wb") as f:
        f.write(response.content)
    print(f"Saved {len(response.content)} bytes.")
else:
    print(response.text[:500])

