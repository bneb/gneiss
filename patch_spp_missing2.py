import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

content = content.replace("let seed_cdt_bds = state.cdt_bds;", "let seed_cdt_bds = state.cdt_bds;\n    let seed_cdt_glo = state.cdt_glo;")
content = content.replace("let next_cdt_bds = if let Some(c) = bds_col { state.cdt_bds + dx_vec[c] } else { state.cdt_bds };", "let next_cdt_bds = if let Some(c) = bds_col { state.cdt_bds + dx_vec[c] } else { state.cdt_bds };\n    let next_cdt_glo = if let Some(c) = glo_col { state.cdt_glo + dx_vec[c] } else { state.cdt_glo };")

# glo_col might have failed to insert because my `target` string in previous Python script was slightly wrong
content = content.replace("    let mut bds_col = None;\n\n    // We ALWAYS want a master clock column", "    let mut bds_col = None;\n    let mut glo_col = None;\n\n    // We ALWAYS want a master clock column")

with open(filename, "w") as f:
    f.write(content)
