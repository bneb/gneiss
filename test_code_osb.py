def check_code_osb():
    lines = open("datasets/wtzr_ppp_1224/com21374.bia").readlines()
    for l in lines:
        if " C1W " in l and "G01" in l: print(l.strip())
        if " C2W " in l and "G01" in l: print(l.strip())

check_code_osb()
