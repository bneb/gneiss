import re
import os

for root, _, files in os.walk("crates/gneiss-rtk/src"):
    for file in files:
        if file.endswith(".rs"):
            with open(os.path.join(root, file), "r") as f:
                text = f.read()

            functions = re.findall(r'((?:pub )?fn [^{]+)\{([^}]*)\}', text)
            for sig, body in functions:
                loc = len(body.strip().split('\n'))
                if loc > 30:
                    name = sig.strip().split('(')[0]
                    if not "test_" in name:
                        print(f"{os.path.join(root, file)}: {name} is {loc} LOC")

