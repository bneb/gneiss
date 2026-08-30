# Gneiss Geodetic & Data Integrity Audit Report

## 1. Coordinate Frame Tags and Velocity Propagation
- Verified tags for ITRF2014, IGS20, and WGS84.
- Ensured correct velocity propagation to the data epoch.

## 2. ARP vs APC Offsets and ANTEX Calibrations
- Checked phase center offset logic (APC vs ARP).
- Verified `igs20.atx` calibrations to prevent double counting.

## 3. Data Leakage
- Verified that ground truth coordinates are isolated to evaluation scripts (`eval` command) and NEVER leaked into filter initialization, seeding, or ambiguity resolution states.
- Audited `run_comprehensive.py`.

## 4. Edge-Case Stress Tests
- Tested scintillation cycle slips handling.
- Evaluated missing broadcast ephemeris fallback.
- Validated high PDOP performance constraints.

## Code Standards
- Confirmed files are < 500 LOC.
- Confirmed functions are < 32 LOC.
- Confirmed 0 warnings.
