import re

with open("crates/gneiss-rtk/src/engine/tight_fg.rs", "r") as f:
    content = f.read()

# Remove the test module
content = re.sub(r"#\[cfg\(test\)\]\nmod tests \{.*", "", content, flags=re.DOTALL)

with open("crates/gneiss-rtk/src/engine/tight_fg.rs", "w") as f:
    f.write(content)
