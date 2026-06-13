import re

with open("crates/gneiss-rtk/src/engine/mod.rs", "r") as f:
    content = f.read()

# Replace the run_combined_ppk method with a call to smoother::run_combined_ppk
pattern = re.compile(r'pub fn run_combined_ppk\(&mut self\) -> Result<Vec<RtkState>, EngineError> \{.*?\n    \}\n\}', re.DOTALL)
replacement = '''pub fn run_combined_ppk(&mut self) -> Result<Vec<RtkState>, EngineError> {
        smoother::run_combined_ppk(self)
    }
}'''
new_content = pattern.sub(replacement, content, count=1)

with open("crates/gneiss-rtk/src/engine/mod.rs", "w") as f:
    f.write(new_content)
