# Gneiss — Claude Code Instructions

## Code Quality Standards

All code MUST meet these standards. Existing code that violates them is technical debt to be resolved.

| Metric | Standard | Rationale |
|--------|----------|-----------|
| **File size** | < 500 LOC | Files over 500 lines are hard to review and reason about |
| **Function size** | < 32 LOC | Functions over 32 lines do too many things |
| **Nesting depth** | < 3 levels | Deep nesting is a bug magnet |
| **Test coverage** | > 95% line coverage | Untested code is broken code |
| **Mutation testing** | 0 survivors | Every mutant must be caught by tests |
| **Compiler warnings** | 0 | Warnings are future errors |

## Hard Rules

- **No `unwrap()` in production code.** Use `match`, `if let`, `ok_or()?`, or `.expect()` with a descriptive message explaining why the invariant holds. `unwrap()` in `#[cfg(test)]` is acceptable.
- **No `.orig` or `.rej` files committed.** These are merge conflict artifacts. Delete them.
- **No empty `#[ignore]` tests.** Either implement the test or remove it. Ignored tests must have a comment explaining why and a condition for re-enabling.
- **No dead code.** Functions, structs, and modules that are never called must be removed or gated behind `#[cfg(test)]`.
- **No duplicate function definitions.** Shared utilities belong in a canonical location imported by all callers.

## Conventions

- **Sign conventions** are documented in predictor.rs:86-91 (attitude error state d_theta = -ψ, left-multiplied global-frame). All measurement Jacobians must be consistent.
- **Phase windup** is always **subtracted** from carrier phase: `cp - windup`, never added.
- **Error handling** uses `Result<T, EngineError>` with `?` propagation. Functions that can fail must not panic.
- **Tests** go in the same file as the code they test, inside `#[cfg(test)] mod tests { ... }`.

## Build & Test Commands

```bash
cargo build --workspace          # Must pass with 0 warnings
cargo test --workspace           # Must pass all tests
cargo clippy --workspace         # Must pass with 0 warnings
```

## Project Architecture

See [ARCHITECTURE.md](./ARCHITECTURE.md) for the full architecture document (note: describes an
earlier engine generation; current engine is network RTK/PPK, see docs/PROJECT_STATUS.md).
See [docs/archive/SPRINT_PLAN_V2.md](./docs/archive/SPRINT_PLAN_V2.md) for the historical PPP accuracy
investigation that motivated the SWFG rewrite (POST_MORTEM.md never existed; this was the closest
equivalent and is archived, not current).
See [docs/PROJECT_STATUS.md](./docs/PROJECT_STATUS.md) and [docs/NETWORK_RTK_NEXT_STEPS.md](./docs/NETWORK_RTK_NEXT_STEPS.md)
for the current sprint plan and roadmap.
See [docs/TIER1_ROADMAP.md](./docs/TIER1_ROADMAP.md) for the cross-cutting roadmap to tier-1 PPK
parity (performance, UX, features, code quality; accuracy is covered in PROJECT_STATUS.md's Sprint 16).
