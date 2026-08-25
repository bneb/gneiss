#!/usr/bin/env python3
"""Round heartbeat: writes state/heartbeat.json with session metrics.

Called at end of each round. Makes monitoring a pure file-read —
no session access needed for external oversight.
"""

import json
import subprocess
import time
from pathlib import Path

STATE_DIR = Path("state")
HEARTBEAT = STATE_DIR / "heartbeat.json"
PREV = HEARTBEAT  # same file; previous beat is overwritten but unix_ts preserved


def _run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return r.stdout.strip()


def count_tests():
    out = _run("cargo test --workspace 2>&1 | grep 'test result'")
    passed = failed = 0
    for line in out.split('\n'):
        p = line.find('passed')
        f = line.find('failed')
        if p >= 0:
            try: passed += int(line[:p].strip().split()[-1])
            except: pass
        if f >= 0:
            try: failed += int(line[:f].strip().split()[-1])
            except: pass
    return passed, failed


def guard_status():
    net = _run("python3 scripts/check_network_benchmark.py 2>&1 | tail -1")
    mgx = _run("python3 scripts/check_multignss_benchmark.py 2>&1 | tail -1")
    return {
        "dataset_a": "PASS" if "PASSED" in net else "FAIL",
        "dataset_b": "PASS" if "PASSED" in mgx else "FAIL",
    }


def main():
    STATE_DIR.mkdir(exist_ok=True)
    now = time.time()
    
    tests_passed, tests_failed = count_tests()
    guards = guard_status()
    head = _run("git log --oneline -1")
    dirty = int(bool(_run("git status --porcelain")))

    # Read previous beat for delta
    prev_interval_min = None
    if PREV.exists():
        try:
            prev = json.loads(PREV.read_text())
            prev_interval_min = (now - prev["unix_ts"]) / 60
        except Exception:
            pass

    beat = {
        "unix_ts": now,
        "human_time": time.strftime("%Y-%m-%d %H:%M:%S"),
        "head": head,
        "dirty_files": dirty,
        "tests_passed": tests_passed,
        "tests_failed": tests_failed,
        "guards": guards,
        "prev_interval_min": round(prev_interval_min, 1) if prev_interval_min else None,
    }
    HEARTBEAT.write_text(json.dumps(beat, indent=2))
    print(json.dumps(beat, indent=2))


if __name__ == "__main__":
    main()
