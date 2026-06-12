import re

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "r") as f:
    code = f.read()

code = code.replace("let idx_i = 18 + i;", "let idx_i = crate::filter::CORE_STATE_SIZE + i;")
code = code.replace("let idx_j = 18 + j;", "let idx_j = crate::filter::CORE_STATE_SIZE + j;")
code = code.replace("index_amb: 18 + amb_idx", "index_amb: crate::filter::CORE_STATE_SIZE + amb_idx")
code = code.replace("18 + state.ambiguities.len()", "crate::filter::CORE_STATE_SIZE + state.ambiguities.len()")
code = code.replace("cov[(18,18)]", "cov[(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE)]")
code = code.replace("opt_delta[18 + i]", "opt_delta[crate::filter::CORE_STATE_SIZE + i]")

# Insert cdt update
tgt = "for i in 0..state.ambiguities.len() { state.ambiguities[i] += opt_delta[crate::filter::CORE_STATE_SIZE + i]; }"
cdt_update = """    if crate::filter::CORE_STATE_SIZE >= 21 {
        state.cdt_gal += opt_delta[18];
        state.cdt_bds += opt_delta[19];
        state.cdt_glo += opt_delta[20];
    }
    """ + tgt
code = code.replace(tgt, cdt_update)

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w") as f:
    f.write(code)

