import os
import gzip
import shutil
from urllib import request

token = os.environ.get("EARTHDATA_TOKEN", "eyJ0eXAiOiJKV1QiLCJvcmlnaW4iOiJFYXJ0aGRhdGEgTG9naW4iLCJzaWciOiJlZGxqd3RwdWJrZXlfb3BzIiwiYWxnIjoiUlMyNTYifQ.eyJ0eXBlIjoiVXNlciIsInVpZCI6ImJuZWIiLCJleHAiOjE3ODU1OTE2NTMsImlhdCI6MTc4MDQwNzY1MywiaXNzIjoiaHR0cHM6Ly91cnMuZWFydGhkYXRhLm5hc2EuZ292IiwiaWRlbnRpdHlfcHJvdmlkZXIiOiJlZGxfb3BzIiwiYWNyIjoiZWRsIiwiYXNzdXJhbmNlX2xldmVsIjozfQ.YILVJGvjQVPDJTKonkmVh5tNPsvlP5tv6rJOd3NqumZVL6_uk2gTzi7RRU830-DMMnWqnkGPYumffaIIRzioL8kP_zz0OrBsWg-qPZACM20fRb5NnGe3SbZ6nH7yeXygSCu_8CTGT1qdz4_GvAlj0Jl-aWly78ftjOqkE0xTvRM-oTMDonJYocl-2JexIuYEaiWBev-tmN1pbYwMV66p_6-DlbzWnnItKvz3-YcHGW-G_45S-kt9WlcOzEHeQoAYJuEVDdHz_LDKl5vLL5OfbLrrMG11d7ewOTv9P9MIaR_r-jrX6SUu0WHV8wSSQUxPatYCscVL91yUmPlJq7UTWg")

files = [
    "CAS0MGXRAP_20203590000_01D_01D_DCB.BIA.gz"
]

base_url = "https://cddis.nasa.gov/archive/gnss/products/mgex/2137/"

for f in files:
    url = base_url + f
    dest = "datasets/rtkexplorer/sample_1/f9p_ppp_1224/" + f
    print(f"Downloading {url}...")
    try:
        req = request.Request(url)
        req.add_header('Authorization', f'Bearer {token}')
        with request.urlopen(req) as response:
            with open(dest, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
        print("Done")
        
        # unzip
        if os.path.exists(dest):
            with gzip.open(dest, 'rb') as f_in:
                with open(dest[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(dest)
    except Exception as e:
        print(f"Failed: {e}")
