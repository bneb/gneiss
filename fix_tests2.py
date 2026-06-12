import re

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    content = f.read()

# Fix Coordinate
content = content.replace(
    "Coordinate::new(1000.0, 1000.0, 1000.0)", 
    "Coordinate::new(Vector3::new(1000.0, 1000.0, 1000.0), gneiss_core::coords::Datum::Wgs84, gneiss_core::coords::Frame::Ecef, GpsTime::new(2000, 100000.0))"
)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(content)
print("Tests fixed again")
