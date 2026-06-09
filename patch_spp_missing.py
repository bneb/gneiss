import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

content = content.replace("seed_cdt\n    )", "seed_cdt,\n        seed_cdt\n    )")
content = content.replace("seed_cdt_bds\n                );", "seed_cdt_bds,\n                    seed_cdt_glo\n                );")
content = content.replace("next_cdt_bds\n    ))", "next_cdt_bds,\n        next_cdt_glo\n    ))")

content = content.replace("""    let mut bds_col = None;

    // We ALWAYS want a master clock column, even if GPS isn't present,""", """    let mut bds_col = None;
    let mut glo_col = None;

    // We ALWAYS want a master clock column, even if GPS isn't present,""")

with open(filename, "w") as f:
    f.write(content)
