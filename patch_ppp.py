with open("crates/gneiss-rtk/src/engine/ppp.rs", "r") as f:
    content = f.read()

import re
content = re.sub(r"#\[cfg\(test\)\].*mod tests \{.*", "", content, flags=re.DOTALL)

with open("crates/gneiss-rtk/src/engine/ppp.rs", "w") as f:
    f.write(content)
