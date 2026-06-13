import sys

tows = [425536.442, 425537.442, 425538.442]

for filename in ['benchmarks/rtklib_comparison/gneiss_GSDC_rtk_forward.pos', 'benchmarks/rtklib_comparison/gneiss_GSDC_rtk_smoothed.pos']:
    print(f"\n{filename}")
    with open(filename) as f:
        for line in f:
            if line.startswith('%'): continue
            parts = line.split()
            if len(parts) < 5: continue
            tow = float(parts[1])
            if tow in tows:
                print(f"TOW {tow}: X={parts[2]} Y={parts[3]} Z={parts[4]} Q={parts[5]}")
