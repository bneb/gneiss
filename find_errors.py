import sys
import math

ref = {}
with open('datasets/gsdc/reference.csv') as f:
    for line in f:
        if line.startswith('GPS') or line.startswith('%'): continue
        parts = line.strip().split(',')
        if len(parts) > 7:
            # Lat, Lon are at indices 2, 3
            ref[float(parts[0])] = (float(parts[5]), float(parts[6]), float(parts[7]), float(parts[2]), float(parts[3]))

for filename in ['benchmarks/rtklib_comparison/gneiss_GSDC_rtk_forward.pos', 'benchmarks/rtklib_comparison/gneiss_GSDC_rtk_smoothed.pos']:
    print(f"\nEvaluating {filename}")
    with open(filename) as f:
        for line in f:
            if line.startswith('%'): continue
            parts = line.split()
            if len(parts) < 5: continue
            tow = float(parts[1])
            x = float(parts[2])
            y = float(parts[3])
            z = float(parts[4])
            
            if tow in ref:
                rx, ry, rz, rlat, rlon = ref[tow]
                dx = x - rx
                dy = y - ry
                dz = z - rz
                
                # Approximate horizontal error:
                # Up vector
                lat_r = math.radians(rlat)
                lon_r = math.radians(rlon)
                up_x = math.cos(lat_r) * math.cos(lon_r)
                up_y = math.cos(lat_r) * math.sin(lon_r)
                up_z = math.sin(lat_r)
                
                v_err = dx*up_x + dy*up_y + dz*up_z
                h_err = math.sqrt(dx*dx + dy*dy + dz*dz - v_err*v_err)
                
                if h_err > 15.0:
                    print(f"tow {tow}: h_err {h_err:.2f}m")
