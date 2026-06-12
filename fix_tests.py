import re

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    content = f.read()

# Fix dummy state initialization
content = content.replace("RtkState::new(Coordinate::new(1000.0, 1000.0, 1000.0), GpsTime::new(2000, 100000.0))",
                          "RtkState::new(GpsTime::new(2000, 100000.0), Coordinate::new(1000.0, 1000.0, 1000.0), 10.0)")

# Fix Constellation::GPS to Constellation::Gps
content = content.replace("Constellation::GPS", "Constellation::Gps")

# Fix SatelliteId::new
content = content.replace("SatelliteId::new(Constellation::Gps, 1).unwrap()", 
                          "SatelliteId { constellation: Constellation::Gps, prn: 1 }")
content = content.replace("SatelliteId::new(Constellation::Gps, 2).unwrap()", 
                          "SatelliteId { constellation: Constellation::Gps, prn: 2 }")


with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(content)
print("Tests fixed")
