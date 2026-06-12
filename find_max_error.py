import csv
import math

ref_data = {}
with open('datasets/gsdc/reference.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        ref_data[int(float(row[0]))] = (float(row[3]), float(row[4]), float(row[5]))

max_err = 0.0
max_err_tow = 0

with open('/tmp/test_ppp_l1.pos', 'r') as f:
    for line in f:
        if line.startswith('%'): continue
        parts = line.split(',')
        if len(parts) < 6: continue
        tow = int(float(parts[1]))
        x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
        if tow in ref_data:
            rx, ry, rz = ref_data[tow]
            err = math.sqrt((x-rx)**2 + (y-ry)**2 + (z-rz)**2)
            if err > max_err:
                max_err = err
                max_err_tow = tow

print(f"Max error: {max_err} at TOW: {max_err_tow}")

with open('/tmp/test_ppp_l1.pos', 'r') as f:
    for line in f:
        if line.startswith('%'): continue
        parts = line.split(',')
        if len(parts) < 6: continue
        if int(float(parts[1])) == max_err_tow:
            print(f"File line: {line.strip()}")
