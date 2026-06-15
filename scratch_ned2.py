import numpy as np

def ecef_to_ned_matrix(lat, lon):
    slat = np.sin(lat)
    clat = np.cos(lat)
    slon = np.sin(lon)
    clon = np.cos(lon)
    return np.array([
        [-slat * clon, -slat * slon, clat],
        [-slon,        clon,         0],
        [-clat * clon, -clat * slon, -slat]
    ])

# Epoch 85.0
lat, lon = np.radians(35.68700593), np.radians(139.69211781)
R = ecef_to_ned_matrix(lat, lon)
dx = -3955039.8671 - (-3955041.6035)
dy = 3355056.7136 - 3355054.7213
dz = 3700086.9867 - 3700086.8807

v_ecef = np.array([dx, dy, dz])
d_ned = R @ v_ecef

print(f"v_ecef: {v_ecef}")
print(f"d_ned: {d_ned}")
