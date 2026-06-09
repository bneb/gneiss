import sys

with open("crates/gneiss-rtk/src/spp.rs", "r") as f:
    content = f.read()

content = content.replace(
"""        match m.constellation {
            gneiss_core::sat::Constellation::Gps => { if let Some(c) = gps_col { h_matrix[(row, c)] = 1.0; } },
            gneiss_core::sat::Constellation::Galileo => { if let Some(c) = gal_col { h_matrix[(row, c)] = 1.0; } },
            gneiss_core::sat::Constellation::Beidou => { if let Some(c) = bds_col { h_matrix[(row, c)] = 1.0; } },
            _ => {},
        }""",
"""        match m.constellation {
            gneiss_core::sat::Constellation::Gps => { if let Some(c) = gps_col { h_matrix[(row, c)] = 1.0; } },
            gneiss_core::sat::Constellation::Galileo => { if let Some(c) = gal_col { h_matrix[(row, c)] = 1.0; } },
            gneiss_core::sat::Constellation::Beidou => { if let Some(c) = bds_col { h_matrix[(row, c)] = 1.0; } },
            gneiss_core::sat::Constellation::Glonass => { if let Some(c) = glo_col { h_matrix[(row, c)] = 1.0; } },
            _ => {},
        }"""
)

with open("crates/gneiss-rtk/src/spp.rs", "w") as f:
    f.write(content)
