import sys
import re

with open("crates/gneiss-rtk/src/engine/tight_fg.rs", "r") as f:
    text = f.read()

text = text.replace("let initial_delta = DVector::zeros(18 + state.ambiguities.len());", "let initial_delta = DVector::zeros(crate::filter::CORE_STATE_SIZE + state.ambiguities.len());")
text = text.replace("for i in 0..state.ambiguities.len() { state.ambiguities[i] += opt_delta[18 + i]; }", "for i in 0..state.ambiguities.len() { state.ambiguities[i] += opt_delta[crate::filter::CORE_STATE_SIZE + i]; }")

with open("crates/gneiss-rtk/src/engine/tight_fg.rs", "w") as f:
    f.write(text)
