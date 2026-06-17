import earthaccess
import os

token = os.environ.get("EARTHDATA_TOKEN")

# We can search without login
results = earthaccess.search_data(
    short_name="GNSS_IGS_MGEX_ORBIT_PRODUCTS",
    temporal=("2020-12-24", "2020-12-25")
)

for r in results:
    if "COD0" in str(r) or "com" in str(r):
        print(r)
print(f"Total results: {len(results)}")
