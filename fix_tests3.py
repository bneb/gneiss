import re

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    content = f.read()

content = content.replace("gneiss_core::coords::Datum::Wgs84", "gneiss_core::coords::Datum::WGS84")
content = content.replace("gneiss_core::coords::Frame::Ecef", "gneiss_core::coords::Frame::ECEF")

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(content)
print("Fixed Datum and Frame")
