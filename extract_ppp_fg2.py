import json
import re

with open("ppp_fg_view.txt", "r") as f:
    for line in f:
        try:
            data = json.loads(line)
            if data.get("step_index") == 13 and data.get("type") == "VIEW_FILE":
                content = data.get("content", "")
                
                lines = content.split('\n')
                out_lines = []
                for l in lines:
                    m = re.match(r'^(\d+):\s?(.*)', l)
                    if m:
                        out_lines.append(m.group(2))
                
                with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w") as out_f:
                    out_f.write('\n'.join(out_lines))
                print(f"Extracted {len(out_lines)} lines.")
                break
        except Exception as e:
            pass
