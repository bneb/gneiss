import math

x = 1000
y = 1000
theta = 0.01

print("GNSS-RTK:")
print("x_new =", x * math.cos(theta) + y * math.sin(theta))
print("y_new =", -x * math.sin(theta) + y * math.cos(theta))

# Inverse rotation matrix (since signal travels from sat to rx, Earth rotates forward)
print("Rotation +theta (Earth rotates):")
print("x_new =", x * math.cos(theta) - y * math.sin(theta))
print("y_new =", x * math.sin(theta) + y * math.cos(theta))
