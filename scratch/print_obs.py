import re
with open("datasets/gsdc/Pixel4_GnssLog.20o") as f:
    for i, line in enumerate(f):
        if i > 25 and line.startswith("G06"):
            print(line.strip())
