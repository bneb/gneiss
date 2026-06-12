import json
import os
import glob

# Search all transcripts for "tight_fg.rs"
transcripts = glob.glob("/Users/kevin/.gemini/antigravity/brain/*/.system_generated/logs/transcript.jsonl")

best_content = ""
max_lines = 0

for ts in transcripts:
    try:
        with open(ts, 'r') as f:
            for line in f:
                try:
                    data = json.loads(line)
                    # Check tool calls
                    if "tool_calls" in data:
                        for tc in data["tool_calls"]:
                            if tc.get("function", {}).get("name") == "default_api:write_to_file":
                                args = tc.get("function", {}).get("arguments", "{}")
                                try:
                                    args_json = json.loads(args)
                                    if "tight_fg.rs" in args_json.get("TargetFile", ""):
                                        content = args_json.get("CodeContent", "")
                                        if len(content.split('\n')) > max_lines:
                                            max_lines = len(content.split('\n'))
                                            best_content = content
                                except:
                                    pass
                    
                    # Check tool responses
                    if "tool_responses" in data:
                        for tr in data["tool_responses"]:
                            if "tight_fg.rs" in tr.get("response", "") and "The following code has been modified to include a line number" in tr.get("response", ""):
                                # This is a view_file response
                                text = tr.get("response", "")
                                # Extract lines
                                lines = text.split('\n')
                                extracted = []
                                for l in lines:
                                    # Match "1: " format
                                    if ": " in l:
                                        parts = l.split(": ", 1)
                                        if parts[0].isdigit():
                                            extracted.append(parts[1])
                                if len(extracted) > max_lines:
                                    max_lines = len(extracted)
                                    best_content = "\n".join(extracted)
                except Exception as e:
                    pass
    except:
        pass

if best_content:
    with open("/Users/kevin/projects/gneiss/crates/gneiss-rtk/src/engine/tight_fg.rs", "w") as f:
        f.write(best_content)
    print(f"Recovered {max_lines} lines to tight_fg.rs")
else:
    print("Could not find full content.")
