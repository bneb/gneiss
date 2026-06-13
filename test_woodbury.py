import numpy as np

P_k_full = np.array([
    [7.869002, 50.0, 0.040804, 0.038019], # pos
    [50.0, 10000.0,  99.9,     99.9],   # clock
    [0.040804, 99.9, 2.84,     2.839],   # amb1
    [0.038019, 99.9, 2.839,    2.84]     # amb2
])

Q_full = np.array([
    [5.0, 0.0,   0.0, 0.0],
    [0.0, 100.0, 0.0, 0.0],
    [0.0, 0.0,   0.0, 0.0],
    [0.0, 0.0,   0.0, 0.0]
])

P_pred_full = P_k_full + Q_full
P_pred_inv_full = np.linalg.inv(P_pred_full)

print("P_pred_inv full:")
print(P_pred_inv_full)

