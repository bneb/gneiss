import ftplib
try:
    ftp = ftplib.FTP("garner.ucsd.edu")
    ftp.login("anonymous", "anonymous@")
    ftp.cwd("/pub/products/mgex/2137")
    files = ftp.nlst()
    for f in files:
        if "GRG" in f.upper() or "grg" in f.lower():
            print(f)
    ftp.quit()
except Exception as e:
    print(e)
