import ftplib
import gzip
import shutil
import os

files = [
    "CAS0MGXRAP_20203590000_01D_01D_DCB.BIA.gz"
]

try:
    ftp = ftplib.FTP("igs.ign.fr")
    ftp.login("anonymous", "anonymous@")
    ftp.cwd("/pub/igs/products/mgex/2137")
    
    for f in files:
        dest = "datasets/rtkexplorer/sample_1/f9p_ppp_1224/" + f
        print(f"Downloading {f}...")
        with open(dest, 'wb') as local_file:
            ftp.retrbinary(f"RETR {f}", local_file.write)
        print("Done")
        
        # unzip
        if os.path.exists(dest):
            with gzip.open(dest, 'rb') as f_in:
                with open(dest[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(dest)
            
    ftp.quit()
except Exception as e:
    print(e)
