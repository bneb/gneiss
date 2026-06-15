import pandas as pd
import numpy as np

def read_pos(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('%')]
    data = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 5:
            # Assuming format: week tow x y z ...
            try:
                tow = float(parts[1])
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                data.append({'tow': tow, 'x': x, 'y': y, 'z': z})
            except:
                pass
    return pd.DataFrame(data)

def read_ref(path):
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={'GPS TOW (s)': 'tow', 'ECEF X (m)': 'x_ref', 'ECEF Y (m)': 'y_ref', 'ECEF Z (m)': 'z_ref'})
    return df

res = read_pos('benchmarks/shinjuku_600.pos')
ref = read_ref('datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv')

merged = pd.merge_asof(res.sort_values('tow'), ref.sort_values('tow'), on='tow', direction='nearest', tolerance=0.1)

if 'x' in merged.columns and 'x_ref' in merged.columns:
    error_x = merged['x'] - merged['x_ref']
    error_y = merged['y'] - merged['y_ref']
    error_z = merged['z'] - merged['z_ref']
    merged['error'] = np.sqrt(error_x**2 + error_y**2 + error_z**2)
    
    # Print the first time the error exceeds 10m
    diverge = merged[merged['error'] > 10.0]
    if len(diverge) > 0:
        print("First epoch > 10m error:")
        print(diverge.iloc[0][['tow', 'error']])
    
    diverge100 = merged[merged['error'] > 100.0]
    if len(diverge100) > 0:
        print("\nFirst epoch > 100m error:")
        print(diverge100.iloc[0][['tow', 'error']])
