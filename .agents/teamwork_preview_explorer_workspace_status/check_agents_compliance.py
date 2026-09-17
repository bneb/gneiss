#!/usr/bin/env python3
"""
Audit Rust source files against AGENTS.md metrics:
- File size < 500 LOC
- Function size < 32 LOC
- Nesting depth < 3 levels
- Zero unwrap() in production code
"""
import os
import re
import sys
from pathlib import Path

WORKSPACE_ROOT = Path("/Users/kevin/projects/gneiss")
CRATES_DIR = WORKSPACE_ROOT / "crates"

def analyze_file(filepath: Path):
    rel_path = filepath.relative_to(WORKSPACE_ROOT)
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    total_loc = len(lines)
    file_too_long = total_loc >= 500

    is_test_file = "tests.rs" in filepath.name or "/tests/" in str(filepath) or filepath.name.startswith("test_")
    in_test_module = is_test_file
    brace_depth = 0
    test_module_brace_depth = 0 if is_test_file else -1

    # Track function definitions
    fn_start_re = re.compile(r'^\s*(pub(\([^)]+\))?\s+)?(async\s+)?(unsafe\s+)?(extern\s+("[^"]+"\s+)?)?fn\s+([a-zA-Z0-9_]+)')

    functions = []
    current_fn = None
    current_fn_start_depth = 0

    nesting_violations = []
    unwrap_violations = []
    next_is_test = False

    for idx, raw_line in enumerate(lines):
        line_num = idx + 1
        line = raw_line.strip()

        # Check for test module or test attribute
        if "#[cfg(test)]" in line or "#[test]" in line:
            next_is_test = True

        # Check comment stripping
        cleaned_line = re.sub(r'"(\\.|[^"\\])*"', '""', raw_line)
        cleaned_line = re.sub(r'//.*$', '', cleaned_line)

        # Check for function start
        fn_match = fn_start_re.search(cleaned_line)
        if fn_match and current_fn is None:
            fn_name = fn_match.group(7)
            is_test_fn = in_test_module or next_is_test or fn_name.startswith("test_")
            current_fn = {
                "name": fn_name,
                "start": line_num,
                "is_test": is_test_fn,
                "brace_opened": False
            }

        # Check for unwrap() in production code
        if not in_test_module and (current_fn is None or not current_fn["is_test"]):
            if re.search(r'\.unwrap\(\)', cleaned_line):
                unwrap_violations.append((line_num, raw_line.strip()))

        # Count open/close braces
        for char in cleaned_line:
            if char == '{':
                brace_depth += 1
                if "#[cfg(test)]" in lines[max(0, idx-2):idx+1] or "mod test" in raw_line:
                    if not in_test_module:
                        in_test_module = True
                        test_module_brace_depth = brace_depth

                if current_fn is not None:
                    if not current_fn["brace_opened"]:
                        current_fn["brace_opened"] = True
                        current_fn_start_depth = brace_depth
                    else:
                        # Inside function, check nesting depth relative to function body
                        # Function body is at relative depth 0 (or 1)
                        # Statements inside if/for/match/loop are relative depth 1, nested block is 2, etc.
                        rel_nesting = brace_depth - current_fn_start_depth
                        if rel_nesting >= 3:
                            # Nesting depth >= 3 levels inside function
                            nesting_violations.append((line_num, rel_nesting, raw_line.strip(), current_fn["is_test"]))

            elif char == '}':
                if current_fn is not None and current_fn["brace_opened"]:
                    if brace_depth == current_fn_start_depth:
                        # Function ended
                        fn_loc = line_num - current_fn["start"] + 1
                        functions.append((
                            current_fn["name"],
                            current_fn["start"],
                            line_num,
                            fn_loc,
                            current_fn["is_test"]
                        ))
                        current_fn = None
                if in_test_module and brace_depth == test_module_brace_depth:
                    in_test_module = False
                    test_module_brace_depth = -1
                brace_depth -= 1

    fn_too_long = [f for f in functions if f[3] >= 32]

    return {
        "file": rel_path,
        "loc": total_loc,
        "file_too_long": file_too_long,
        "functions": functions,
        "fn_too_long": fn_too_long,
        "nesting_violations": nesting_violations,
        "unwrap_violations": unwrap_violations,
    }

def main():
    rs_files = sorted(WORKSPACE_ROOT.glob("**/*.rs"))
    # Exclude target/, .agents/, .worktrees/, scratch/
    rs_files = [
        f for f in rs_files
        if not any(part.startswith(".") for part in f.parts[:-1])
        and "target" not in f.parts
        and "scratch" not in f.parts
    ]
    print(f"Auditing {len(rs_files)} active workspace Rust files...")

    results = [analyze_file(f) for f in rs_files]

    long_files = [r for r in results if r["file_too_long"]]
    long_fns = [(r["file"], fn) for r in results for fn in r["fn_too_long"]]
    unwraps = [(r["file"], u) for r in results for u in r["unwrap_violations"]]
    nestings = [(r["file"], n) for r in results for n in r["nesting_violations"] if not n[3]]

    print(f"\n==================================================")
    print(f"--- File Size Violations (>= 500 LOC) ---: {len(long_files)}")
    print(f"==================================================")
    for r in sorted(long_files, key=lambda x: -x['loc']):
        print(f"  {r['file']}: {r['loc']} LOC")

    print(f"\n==================================================")
    print(f"--- unwrap() in Production Code ---: {len(unwraps)}")
    print(f"==================================================")
    for file, u in unwraps:
        print(f"  {file}:{u[0]}: {u[1]}")

    print(f"\n==================================================")
    prod_long_fns = [f for f in long_fns if not f[1][4]]
    print(f"--- Production functions >= 32 LOC ---: {len(prod_long_fns)}")
    print(f"==================================================")
    for file, fn in sorted(prod_long_fns, key=lambda x: -x[1][3]):
        print(f"  {file}:{fn[1]} '{fn[0]}' ({fn[3]} LOC)")

    print(f"\n==================================================")
    print(f"--- Nesting Depth Violations (>= 3 in prod) ---: {len(nestings)}")
    print(f"==================================================")
    print(f"Total nesting depth >= 3 lines in production: {len(nestings)}")
    # Group by file
    files_with_nesting = {}
    for file, n in nestings:
        files_with_nesting[file] = files_with_nesting.get(file, 0) + 1
    for file, count in sorted(files_with_nesting.items(), key=lambda x: -x[1]):
        print(f"  {file}: {count} occurrences")

if __name__ == "__main__":
    main()
