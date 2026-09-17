# Gneiss Geodetic & Data Integrity Audit Report

This document records the geodetic integrity, coordinate frame models, and physical correction audits implemented in Gneiss.

## 1. Reference Frame Inventory & Realizations
- **Supported Realizations**: `ITRF2014`, `ITRF2020`, `IGS14`, `IGS20`, `WGS84` (G1762/G2139), `Nad83_2011`, `Etrs89`, and `Gda2020`.
- **14-Parameter Time-Dependent Helmert Transformations**:
  - Full coordinate transformation including 3D translations ($T_x, T_y, T_z$), rotations ($R_x, R_y, R_z$), scale factor ($s$), and their secular rates of change ($\dot{T}_x, \dot{T}_y, \dot{T}_z, \dot{R}_x, \dot{R}_y, \dot{R}_z, \dot{s}$) evaluated at the epoch of observation.
  - Aligned with authoritative specifications: IOGP EPSG:8970 (Method 1056) and NOAA NGS HTDP (Horizontal Time-Dependent Positioning).
  - Validated against textbook benchmark point (Station SALT AIR at epochs 2010.0 and 2020.0).

## 2. Geodynamic & Physical Corrections
- **Solid Earth Tides (SET)**:
  - Formulated in accordance with IERS Conventions (2010), Chapter 7.
  - Evaluates degree-2 and degree-3 Love and Shida numbers driven by low-precision analytical solar and lunar ephemerides ($d_{J2000}$ astronomical century calculation).
- **Ocean Tide Loading (OTL)**:
  - Models 11 diurnal and semi-diurnal ocean tide constituents ($M_2, S_2, N_2, K_2, K_1, O_1, P_1, Q_1, M_f, M_m, S_{sa}$) from ocean loading grids (e.g. FES2014b / GOT4.10c).
- **Phase Windup Correction**:
  - Carrier phase windup angle computed from receiver and satellite dipole antenna unit vectors, strictly subtracted from carrier phase measurements (`cp - windup`).
- **Antenna Phase Center Offsets (PCO) & Variations (PCV)**:
  - Strict separation of Antenna Reference Point (ARP) and Antenna Phase Center (APC).
  - Frequency-dependent satellite and receiver PCO/PCV corrections from IGS ANTEX (`.atx`) tables.

## 3. Local Datum Ties & Site Calibration
- **`LocalDatumTie` Estimator**:
  - Computes rigid 3D spatial translations between local CORS monument coordinates (e.g. NAD83(2011)) and global satellite orbit frames (e.g. ITRF2014).
  - Protects against artificial residual inflation when comparing global PPP solutions to local RTK ground truth.

## 4. Benchmark Isolation & Invariant Checks
- **Data Leakage Safeguards**:
  - Ground truth coordinates in `tests/src/benchmark_matrix.rs` and evaluation binaries are strictly isolated for post-hoc validation; filter states are seeded from perturbed or unconstrained initial positions.
- **Standards Compliance**:
  - Files $< 500$ LOC, functions $< 32$ LOC, nesting depth $< 3$ levels, 0 compiler/clippy warnings, and zero `unwrap()` in production code.
