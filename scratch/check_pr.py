import math
import sys

# read some lines from the rinex file
with open('datasets/gsdc/Pixel4_GnssLog.20o') as f:
    for _ in range(50):
        print(f.readline().strip())
