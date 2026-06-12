with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "r") as f:
    content = f.read()

# I will just erase the test_compute_dd I just added, and append a better one.
import re
content = re.sub(r'#\[test\]\s+fn test_compute_dd.*?\}\n\}', '}', content, flags=re.DOTALL)

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "w") as f:
    f.write(content)
