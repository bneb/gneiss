import math

def ecef_to_llh(x, y, z):
    a = 6378137.0
    b = 6356752.3142
    f = (a - b) / a
    e_sq = f * (2.0 - f)
    
    p = math.sqrt(x*x + y*y)
    theta = math.atan2(z*a, p*b)
    
    lon = math.atan2(y, x)
    lat = math.atan2(z + math.pow(1.0 - f, -2.0) * e_sq * b * math.pow(math.sin(theta), 3.0),
                     p - e_sq * a * math.pow(math.cos(theta), 3.0))
    N = a / math.sqrt(1.0 - e_sq * math.pow(math.sin(lat), 2.0))
    h = p / math.cos(lat) - N
    
    return math.degrees(lat), math.degrees(lon), h

print(ecef_to_llh(-2694569.5800, -4296509.3320, 3854822.9075))
