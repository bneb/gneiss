# Gneiss PPK Engine — Round Runbook

## Every round follows this procedure. No exceptions.

### 1. VERIFY (2 min)
```bash
cargo test --workspace 2>&1 | grep "test result"
./target/release/eval_network_ppk >/dev/null 2>&1
python3 scripts/check_network_benchmark.py 2>&1 | tail -1
python3 scripts/check_multignss_benchmark.py 2>&1 | tail -1
```
- If ANY fail: fix before proceeding. Do not build on broken state.
- If ALL pass: continue to step 2.

### 2. IMPLEMENT (one change per round)
- Pick ONE item from the prioritized roadmap in PROJECT_STATUS.md
- Write tests FIRST (coordinator-written, not implementer-written)
- Implement against tests
- Run tests → must be green before committing

### 3. MEASURE (5 min)
```bash
GNEISS_DATASET=multi2025 GNEISS_SYSTEMS=GE ./target/release/eval_network_ppk >/dev/null 2>&1
# Compare against previous run's metrics
```
- If improved: keep, update guard budgets if warranted
- If neutral: keep only if it prevents a bug class or enables future work
- If regressed: revert immediately, document why

### 4. COMMIT + HEARTBEAT (2 min)
```bash
git add -A && git commit -m "descriptive message with measured impact"
python3 scripts/heartbeat.py
```

## Rules

1. **Never commit without running both guards**
2. **One change per round** — don't bundle unrelated improvements
3. **Revert on regression** — no exceptions, no "I'll fix it next round"
4. **Document negative results** — they prevent wasted future effort
5. **Frame safety first** — use types from gnss_time.rs, frames.rs, frequencies.rs for new code
6. **Sequential agents on shared crates** — parallel agents on mod.rs/update.rs create conflicts

## Known Dead Ends (do not re-attempt)

See docs/PROJECT_STATUS.md and docs/archive/NETWORK_RTK_NEXT_STEPS.md "Validated Negative Results" section.

## Key Files

| file | purpose |
|---|---|
| docs/PROJECT_STATUS.md | comprehensive status + sprint logs |
| docs/TIER1_ROADMAP.md | Tier-1 PPK roadmap |
| docs/FRAME_SAFETY_PLAN.md | frame-safety architecture |
| scripts/check_network_benchmark.py | dataset A guard |
| scripts/check_multignss_benchmark.py | dataset B guard |
| scripts/heartbeat.py | round telemetry |
| docs/archive/ | archived historical logs & superseded docs |
