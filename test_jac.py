import numpy as np

sat_pos = np.array([20000000.0, 15000000.0, 10000000.0])
nominal_r = np.array([1000.0, 2000.0, 3000.0])
delta = np.array([5.0, -3.0, 2.0])
rx = nominal_r[0] + delta[0]
ry = nominal_r[1] + delta[1]
rz = nominal_r[2] + delta[2]

dx = sat_pos[0] - rx
dy = sat_pos[1] - ry
dz = sat_pos[2] - rz
dist = np.sqrt(dx**2 + dy**2 + dz**2)

anal_jac = np.array([dx/dist, dy/dist, dz/dist, -1.0, -1.0]) # zwd map_wet=1.0

def res(d):
    rx = nominal_r[0] + d[0]
    ry = nominal_r[1] + d[1]
    rz = nominal_r[2] + d[2]
    dx = sat_pos[0] - rx
    dy = sat_pos[1] - ry
    dz = sat_pos[2] - rz
    dist = np.sqrt(dx**2 + dy**2 + dz**2)
    return -dist - d[3] - d[4] * 1.0

eps = 1e-6
num_jac = []
for i in range(5):
    d1 = delta.copy()
    d2 = delta.copy()
    d1 = np.append(d1, [0.1, 0.05])
    d2 = np.append(d2, [0.1, 0.05])
    d1[i] += eps
    d2[i] -= eps
    r1 = res(d1)
    r2 = res(d2)
    num_jac.append((r1 - r2) / (2 * eps))

print("Anal:", anal_jac)
print("Num :", num_jac)
print("Diff:", np.array(anal_jac) - np.array(num_jac))
