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
dx, dy, dz = 1.5718, 1.2814, -0.8778

d_ned = R @ np.array([dx, dy, dz])

print(f"d_ned: {d_ned}")
