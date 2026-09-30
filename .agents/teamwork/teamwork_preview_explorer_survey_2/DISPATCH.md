## 2026-09-25T21:04:10Z

You are Survey Explorer 2 for the Gneiss Frame Safety & Epoch Alignment Refactoring.
Your working directory is: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_2
Authoritative request file: /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md

You MUST read /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md and /Users/kevin/projects/gneiss/AGENTS.md first.

Focus: R2 (Temporal Frame & Epoch Alignment Safety).
1. Survey time representations across crates/gneiss-core and crates/gneiss-rtk (e.g. GpsTime, UtcTime, GlonassTime, BdtTime, Epoch).
2. Identify where ad-hoc float rounding or integer truncation for epoch matching currently occurs.
3. Investigate leap-second handling (GPS vs UTC 18-sec offset, GLONASS leap seconds, BeiDou 14-sec offset vs GPS).
4. Check how IMU samples (e.g. 50Hz/100Hz) and GNSS observables (e.g. 10Hz) are synchronized in eval_odaiba_ins and other estimation pipelines.
5. Propose strictly typed epoch types and conversion methods that prevent cross-time-system arithmetic without explicit leap second correction.
6. Identify all modules/files that will need refactoring.
7. Save your full findings in /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_2/survey_r2_temporal.md and write a complete self-contained handoff.md. Send a completion message when finished.
