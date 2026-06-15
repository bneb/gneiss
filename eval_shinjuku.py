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
    # Strip whitespace from column names just in case
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={'GPS TOW (s)': 'tow', 'ECEF X (m)': 'x_ref', 'ECEF Y (m)': 'y_ref', 'ECEF Z (m)': 'z_ref'})
    return df

res = read_pos('benchmarks/shinjuku_1000_fixed.pos')
ref = read_ref('datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv')

merged = pd.merge_asof(res.sort_values('tow'), ref.sort_values('tow'), on='tow', direction='nearest', tolerance=0.1)

if 'x' in merged.columns and 'x_ref' in merged.columns:
    error_x = merged['x'] - merged['x_ref']
    error_y = merged['y'] - merged['y_ref']
    error_z = merged['z'] - merged['z_ref']
    error = np.sqrt(error_x**2 + error_y**2 + error_z**2)
    print(f"Mean Error: {error.mean():.3f} m")
    print(f"Median Error: {error.median():.3f} m")
    print(f"95th Percentile: {np.percentile(error, 95):.3f} m")
    print(f"Max Error: {error.max():.3f} m")
    
    # Dump worst epochs
    merged['error'] = error
    print("\nWorst Epochs:")
    print(merged.sort_values('error', ascending=False)[['tow', 'error']].head(20))
