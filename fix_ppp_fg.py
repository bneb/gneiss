import re

content = open('crates/gneiss-rtk/src/engine/ppp_fg.rs').read()

content = content.replace("nominal_dt_gal: state.cdt_gal, nominal_dt_bds: state.cdt_bds, nominal_dt_glo: state.cdt_glo,", "nominal_dt_gal: 0.0, nominal_dt_bds: 0.0, nominal_dt_glo: 0.0,")
content = content.replace("state.cdt_gal += opt_delta[18];", "// state.cdt_gal += opt_delta[18];")
content = content.replace("state.cdt_bds += opt_delta[19];", "// state.cdt_bds += opt_delta[19];")
content = content.replace("state.cdt_glo += opt_delta[20];", "// state.cdt_glo += opt_delta[20];")

content = content.replace("state.reset_ambiguity_variance(s.sat_obs.sat, 0, amb_est, 100.0);", "// state.reset_ambiguity_variance(s.sat_obs.sat, 0, amb_est, 100.0);")

with open('crates/gneiss-rtk/src/engine/ppp_fg.rs', 'w') as f:
    f.write(content)
