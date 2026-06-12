import sys
import math

def parse_rtk(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('%'): continue
            parts = line.split()
            if len(parts) >= 15:
                # 2020/05/14 22:11:03.442
                m_str = parts[1].split(':')[1]
                s_str = parts[1].split(':')[2]
                t = int(m_str) * 60 + float(s_str)
                data.append((t, float(parts[2]), float(parts[3]), float(parts[4])))
    return data

def parse_gn(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('%'): continue
            parts = line.split()
            if len(parts) >= 5 and len(parts[0]) == 4:
                # 2105 425463.442
                tow = float(parts[1])
                data.append((tow, float(parts[2]), float(parts[3]), float(parts[4])))
    return data

rtk = parse_rtk("benchmarks/rtklib_comparison/rtklib_GSDC_Pixel_4_RTK_Kinematic.pos")
gn = parse_gn("benchmarks/rtklib_comparison/gneiss_GSDC_rtk_forward.pos")

print(f"Loaded {len(rtk)} RTK and {len(gn)} Gneiss")

gn_dict = {round(g[0] * 10): g for g in gn} # TOW * 10
base_tow = gn[0][0] - rtk[0][0] # offset between TOW and RTKLIB minute+second

diffs = []
for r in rtk:
    rs = round((r[0] + base_tow) * 10)
    if rs in gn_dict:
        g = gn_dict[rs]
        d = math.sqrt((r[1]-g[1])**2 + (r[2]-g[2])**2 + (r[3]-g[3])**2)
        diffs.append((d, r[0], r, g))

diffs.sort(reverse=True)
print("Top 20 diffs:")
for d, t, r, g in diffs[:20]:
    print(f"T={t:.1f}, Diff={d:.2f} m")

