import sympy as sp

# We want to derive the Jacobians for the IMU preintegration residuals.
# Since sympy doesn't handle SO(3) elegantly, we'll derive the translational and velocity parts,
# and write down the SO(3) parts manually based on Forster 2015.

# States at i
p_i = sp.Matrix(sp.symbols('p_ix p_iy p_iz'))
v_i = sp.Matrix(sp.symbols('v_ix v_iy v_iz'))
# Rotation matrix R_i is treated as a 3x3 matrix symbol
R_i = sp.MatrixSymbol('R_i', 3, 3)

# States at j
p_j = sp.Matrix(sp.symbols('p_jx p_jy p_jz'))
v_j = sp.Matrix(sp.symbols('v_jx v_iy_j v_iz_j'))

# Preintegrated measurements
dp = sp.Matrix(sp.symbols('dp_x dp_y dp_z'))
dv = sp.Matrix(sp.symbols('dv_x dv_y dv_z'))

# Gravity and time
g = sp.Matrix(sp.symbols('g_x g_y g_z'))
dt = sp.Symbol('dt')

# Residuals
r_p = R_i.T * (p_j - p_i - v_i * dt - 0.5 * g * dt**2) - dp
r_v = R_i.T * (v_j - v_i - g * dt) - dv

print("r_p w.r.t p_i:")
print(sp.diff(r_p, p_i))

print("r_p w.r.t p_j:")
print(sp.diff(r_p, p_j))

print("r_p w.r.t v_i:")
print(sp.diff(r_p, v_i))

