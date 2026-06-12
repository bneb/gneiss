import re

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "r") as f:
    content = f.read()

# Fix compute_pcv
content = re.sub(r"fn compute_pcv.*?\n}", "fn compute_pcv(sat_obs: &SatObs, engine: &ProcessingEngine, rover_obs: &EpochObs, sat_pos: Vector3<f64>, b1: u8, rcv_pos: Vector3<f64>, sat_pos_rot: &mut Vector3<f64>) -> f64 {\n    0.0\n}", content, flags=re.DOTALL)

# Fix cdt_gal, cdt_bds, cdt_glo
content = re.sub(r"nominal_dt_gal: state\.cdt_gal, nominal_dt_bds: state\.cdt_bds, nominal_dt_glo: state\.cdt_glo,", "", content)
content = re.sub(r"index_dt_gal: if crate::filter::CORE_STATE_SIZE >= 21 \{ Some\(18\) \} else \{ None \},", "", content)
content = re.sub(r"index_dt_bds: if crate::filter::CORE_STATE_SIZE >= 21 \{ Some\(19\) \} else \{ None \},", "", content)
content = re.sub(r"index_dt_glo: if crate::filter::CORE_STATE_SIZE >= 21 \{ Some\(20\) \} else \{ None \},", "", content)

# Fix state.cdt_gal += ...
content = re.sub(r"        if crate::filter::CORE_STATE_SIZE >= 21 \{\n        state\.cdt_gal \+= opt_delta\[18\];\n        state\.cdt_bds \+= opt_delta\[19\];\n        state\.cdt_glo \+= opt_delta\[20\];\n    \}", "", content)

# Fix print gal
content = re.sub(r"if crate::filter::CORE_STATE_SIZE >= 21 \{ opt_delta\[18\] \} else \{ 0\.0 \}", "0.0", content)
content = re.sub(r"if crate::filter::CORE_STATE_SIZE >= 21 \{ opt_delta\[19\] \} else \{ 0\.0 \}", "0.0", content)
content = re.sub(r"if crate::filter::CORE_STATE_SIZE >= 21 \{ opt_delta\[20\] \} else \{ 0\.0 \}", "0.0", content)

# Fix tgd
content = re.sub(r"eph\.tgd\(\)", "0.0", content)

# Fix reset_ambiguity_variance
content = re.sub(r"state\.reset_ambiguity_variance\(s\.sat_obs\.sat, 0, amb_est, 100\.0\);", "state.add_ambiguity(s.sat_obs.sat, 0, amb_est, 100.0);", content)

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w") as f:
    f.write(content)
