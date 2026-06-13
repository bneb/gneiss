import numpy as np
def read_ecef(path):
    times, positions = [], []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('%') or not line.strip(): continue
            parts = line.split()
            try:
                times.append(float(parts[1]))
                positions.append([float(parts[2]), float(parts[3]), float(parts[4])])
            except: pass
    return np.array(times), np.array(positions)

t, p = read_ecef('scratch/test.pos')
idx = (t >= 427000) & (t <= 427100)
if np.sum(idx) > 0:
    t_win = t[idx]
    p_win = p[idx]
    dp = np.linalg.norm(p_win[1:] - p_win[:-1], axis=1)
    print(f"Max jump in 427000-427100: {np.max(dp):.2f}m at {t_win[np.argmax(dp)]}")
    
