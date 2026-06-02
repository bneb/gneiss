import sys
import pandas as pd
import numpy as np
import pyproj

def llh_to_ecef(lat, lon, alt):
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    x, y, z = pyproj.transform(lla, ecef, lon, lat, alt, radians=False)
    return x, y, z

def millis_to_tow(millis):
    # GPS epoch is 1980-01-06
    # 1 week = 604800 seconds
    seconds = millis / 1000.0
    week = int(seconds // 604800)
    tow = seconds % 604800
    return tow

def convert(input_csv, output_csv):
    df = pd.read_csv(input_csv)
    
    with open(output_csv, 'w') as f:
        f.write("GPS TOW (s), GPS Week, Latitude (deg), Longitude (deg), Ellipsoid Height (m), ECEF X (m), ECEF Y (m), ECEF Z (m), Roll (deg), Pitch (deg), Heading (deg), Velocity X (m/s), Velocity Y (m/s), Velocity Z (m/s), Acceleration X (m/s^2), Acceleration Y (m/s^2), Acceleration Z (m/s^2), Angular rate X (rad/s), Angular rate Y (rad/s), Angular rate Z (rad/s)\n")
        
        for _, row in df.iterrows():
            millis = row['millisSinceGpsEpoch']
            tow = millis_to_tow(millis)
            
            lat = row['latDeg']
            lon = row['lngDeg']
            alt = row['heightAboveWgs84EllipsoidM']
            
            x, y, z = llh_to_ecef(lat, lon, alt)
            
            f.write(f"{tow:.3f}, 0, {lat:.8f}, {lon:.8f}, {alt:.4f}, {x:.4f}, {y:.4f}, {z:.4f}, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0\n")

if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2])
