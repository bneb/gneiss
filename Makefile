.PHONY: coverage mutants lint check test

# Run tarpaulin with HTML output
coverage:
	cargo tarpaulin --out Html --output-dir target/tarpaulin $(ARGS)

# Run cargo-mutants on all crates (skip noise from tests helper crate)
mutants:
	cargo mutants --workspace --exclude tests $(ARGS)

# Run clippy with strict settings across the full workspace
lint:
	cargo clippy --workspace --all-targets -- -D warnings -D clippy::all -W clippy::pedantic $(ARGS)

# Run tests, then lint, then coverage, then mutants
check: test lint coverage mutants

# Run workspace tests
test:
	cargo test --workspace $(ARGS)
