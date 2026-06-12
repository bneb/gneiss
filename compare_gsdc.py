import sys
import math
import csv

def read_pos(path):
    pos = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('%'): continue
            parts = line.split()
            if len(parts) >= 6 and len(parts[0]) == 4:
                try:
                    week = int(parts[0])
                    tow = float(parts[1])
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    pos.append((week, tow, x, y, z))
                except ValueError:
                    pass
    return pos

gn = read_pos("benchmarks/rtklib_comparison/gneiss_GSDC_rtk_forward.pos")

import subprocess

def get_eval():
    result = subprocess.run(["target/release/gneiss-cli", "eval", "--solution", "benchmarks/rtklib_comparison/gneiss_GSDC_rtk_forward.pos", "--truth", "datasets/gsdc/reference.csv", "--json"], capture_epoch_errors=True)

# Actually I'll just write a quick Rust snippet to print largest errors.

