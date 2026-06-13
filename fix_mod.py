with open("crates/gneiss-rtk/src/engine/mod.rs", "r") as f:
    lines = f.readlines()

new_lines = []
skip = False
for line in lines:
    if line.startswith("    pub fn run_combined_ppk(&mut self) -> Result<Vec<RtkState>, EngineError> {"):
        new_lines.append(line)
        new_lines.append("        smoother::run_combined_ppk(self)\n")
        new_lines.append("    }\n")
        skip = True
        continue
        
    if skip:
        if line.startswith("    }") and lines[lines.index(line) + 1].startswith("}"):
            skip = False
        continue
        
    if line.startswith("pub mod ppp_fg;"):
        new_lines.append("pub mod smoother;\n")
        
    new_lines.append(line)

with open("crates/gneiss-rtk/src/engine/mod.rs", "w") as f:
    f.writelines(new_lines)
