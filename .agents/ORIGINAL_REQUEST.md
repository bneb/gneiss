# Original User Request

## Initial Request — 2026-06-15T18:14:31Z

# Teamwork Project Prompt — Draft

> Status: Launched
> Goal: Craft prompt → get user approval → delegate to teamwork_preview

Federate a search through the `gneiss` workspace's test suite to find mistakes in test assertions that might be hiding actual bugs.

Working directory: /Users/kevin/projects/gneiss
Integrity mode: demo

## Requirements

### R1. Broad Workspace Scan
Scan the test suite across the entire `gneiss` workspace. Identify test assertions that contain logical errors, trivially pass when they shouldn't, shadow underlying bugs, or misuse mathematical/domain constraints.

### R2. Non-Destructive Reporting
Do not modify the test or production codebase. Focus entirely on analysis and generate a comprehensive markdown report of the uncovered assertion mistakes. 

## Acceptance Criteria

### Reporting Quality
- [ ] The team outputs a report named `suspicious_tests_report.md`.
- [ ] Every flagged test in the report explicitly cites the absolute file path and line number of the suspicious assertion.
- [ ] Every flagged test includes a brief, concrete explanation of why the assertion is logically flawed and what bug it might be obscuring.
