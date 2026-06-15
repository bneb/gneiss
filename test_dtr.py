import numpy as np

def relativistic_correction(x, y, z, vx, vy, vz):
    c = 299792458.0
    dtr = -2.0 * (x * vx + y * vy + z * vz) / (c * c)
    return dtr

# typical MEO orbit
r = 26500000.0
v = 3800.0
print(relativistic_correction(r, 0, 0, v, 0, 0))
