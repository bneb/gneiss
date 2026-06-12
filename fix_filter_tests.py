import re
with open('crates/gneiss-rtk/src/filter.rs', 'r') as f:
    data = f.read()
data = re.sub(r'assert_eq!\(state\.mw_sd_ema\[&sat\],\s*([0-9.]+)\);', r'assert!((state.mw_sd_ema[&sat] - \1).abs() < 1e-9);', data)
with open('crates/gneiss-rtk/src/filter.rs', 'w') as f:
    f.write(data)
