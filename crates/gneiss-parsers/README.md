# gneiss-parsers

`gneiss-parsers` provides zero-allocation, stream-oriented, and batch decoders for industry-standard GNSS observation, broadcast ephemeris, precise product, and receiver protocol formats.

## Supported Formats

### 1. Observations & Navigation Messages
- **RINEX (v2.xx, v3.xx)**:
  - Observation files (`.obs`, `.*o`): Multi-constellation, multi-frequency observation records with signal quality metrics.
  - Navigation files (`.nav`, `.*n`, `.*p`): GPS, GLONASS, Galileo, and BeiDou broadcast ephemerides.
- **RINEX Precise Clocks (`.clk`)**: High-rate satellite and receiver clock corrections.
- **RTCM Standard 10403.x (RTCM3)**:
  - Multiple Signal Messages (MSM1 through MSM7) for GPS, GLONASS, Galileo, and BeiDou.
  - Ephemeris messages (1019, 1020, 1045, 1046, 1042).
  - Station antenna description and reference point coordinates (1005, 1006, 1007, 1008).
  - State Space Representation (SSR) orbit and clock corrections.

### 2. Precise Geodetic Products
- **Precise Orbits (`.sp3`)**: SP3-c and SP3-d format precise satellite ephemerides with Lagrange polynomial interpolation.
- **Biases & Calibration**:
  - **SINEX OSB/BIA (`.bia`, `.bsx`)**: Observation-Specific Biases for code and carrier phase ambiguity resolution.
  - **Bernese DCB (`.dcb`)**: CODE/Bernese Differential Code Biases resolving inter-frequency biases (e.g. GLONASS $P_2-C_2$).
  - **ANTEX (`.atx`)**: Antenna Exchange format for satellite and receiver Phase Center Offsets (PCO) and Variations (PCV).
  - **BLQ Ocean Loading (`.blq`)**: Station-specific 11-constituent ocean tide loading amplitudes and phase angles.
  - **IONEX (`.ionex`)**: Global Ionospheric Map (GIM) total electron content grids.

### 3. Receiver Hardware Protocols & Reference Trajectories
- **u-blox UBX**: `RXM-RAWX` (carrier phase and pseudorange measurements), `RXM-SFRBX` (subframe navigation bits), and `ESF-RAW` (inertial measurements).
- **Septentrio SBF**: Raw measurement and tracking blocks.
- **Applanix SBET**: Binary trajectory exchange format for POSPac reference verification.
- **CSRS-PPP (`.pos`)**: Canada Geodetic Service Precise Point Positioning output files for commercial benchmark comparison.
