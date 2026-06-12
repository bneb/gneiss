import json
import re

with open("ppp_fg_view.txt", "r") as f:
    for line in f:
        try:
            data = json.loads(line)
            if data.get("step_index") == 13 and data.get("type") == "VIEW_FILE":
                content = data.get("content", "")
                
                lines = content.split('\n')
                print(f"Total lines in content: {len(lines)}")
                
                matches = 0
                for i, l in enumerate(lines[:100]):
                    m = re.match(r'^(\d+):\s?(.*)', l)
                    if m:
                        matches += 1
                    else:
                        print(f"Failed to match line {i}: {repr(l[:50])}")
                print(f"Matches in first 100: {matches}")
                break
        except Exception as e:
            pass
