import re

with open("crates/gneiss-rtk/src/engine/mod.rs", "r") as f:
    content = f.read()

# Replace block 1
block1_regex = r"        if self\.config\.enable_nhc \{\n            if let Some\(state\) = self\.current_state\.as_mut\(\) \{\n(?:.*\n){1,70}?                    let _ = crate::nhc::apply_nhc\(state.*?tuning\);\n                \}\n            \}\n        \}"

new_block1 = """        if let Some(state) = self.current_state.as_mut() {
            apply_imu_constraints(state, &self.config, &self.imu_history);
        }"""

content, n1 = re.subn(block1_regex, new_block1, content, count=1, flags=re.MULTILINE)

# Replace block 2
block2_regex = r"        if self\.config\.enable_nhc && is_ins \{\n            if let Some\(state\) = self\.current_state\.as_mut\(\) \{\n(?:.*\n){1,70}?                    let _ = crate::nhc::apply_nhc\(state.*?tuning\);\n                \}\n            \}\n        \}"

new_block2 = """        if is_ins {
            if let Some(state) = self.current_state.as_mut() {
                apply_imu_constraints(state, &self.config, &self.imu_history);
            }
        }"""

content, n2 = re.subn(block2_regex, new_block2, content, count=1, flags=re.MULTILINE)

print(f"Replaced block 1: {n1}, block 2: {n2}")

with open("crates/gneiss-rtk/src/engine/mod.rs", "w") as f:
    f.write(content)
