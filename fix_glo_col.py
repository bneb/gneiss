import sys

with open("crates/gneiss-rtk/src/spp.rs", "r") as f:
    content = f.read()

content = content.replace(
"""    let mut gps_col = None;
    let mut gal_col = None;
    let mut bds_col = None;

    if measurements.iter().any(|m| m.constellation == gneiss_core::sat::Constellation::Gps) {
        gps_col = Some(cols); cols += 1;
    }
    if measurements.iter().any(|m| m.constellation == gneiss_core::sat::Constellation::Galileo) {
        gal_col = Some(cols); cols += 1;
    }
    if measurements.iter().any(|m| m.constellation == gneiss_core::sat::Constellation::Beidou) {
        bds_col = Some(cols); cols += 1;
    }""",
"""    let mut gps_col = None;
    let mut gal_col = None;
    let mut bds_col = None;
    let mut glo_col = None;

    if measurements.iter().any(|m| m.constellation == gneiss_core::sat::Constellation::Gps) {
        gps_col = Some(cols); cols += 1;
    }
    if measurements.iter().any(|m| m.constellation == gneiss_core::sat::Constellation::Galileo) {
        gal_col = Some(cols); cols += 1;
    }
    if measurements.iter().any(|m| m.constellation == gneiss_core::sat::Constellation::Beidou) {
        bds_col = Some(cols); cols += 1;
    }
    if measurements.iter().any(|m| m.constellation == gneiss_core::sat::Constellation::Glonass) {
        glo_col = Some(cols); cols += 1;
    }"""
)

with open("crates/gneiss-rtk/src/spp.rs", "w") as f:
    f.write(content)
