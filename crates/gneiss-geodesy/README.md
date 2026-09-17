# gneiss-geodesy

`gneiss-geodesy` implements reference frame transformations, geodynamic tidal physics, phase windup, and local site calibration for high-precision GNSS positioning.

## Capabilities

### 1. 14-Parameter Time-Dependent Helmert Transformations
- Full time-dependent transformation between global and continental reference frames (e.g. ITRF2014, ITRF2020, NAD83(2011), ETRS89, GDA2020).
- Incorporates 3D translations ($T_x, T_y, T_z$), rotations ($R_x, R_y, R_z$), scale factor ($s$), and their secular rates ($\dot{T}, \dot{R}, \dot{s}$) evaluated at the epoch of observation:
  $$\mathbf{x}(t) = \mathbf{T}(t) + [1 + s(t)] \mathbf{R}(t) \mathbf{x}_0$$
- Mathematically conforms to IOGP EPSG:8970 (Method 1056) and NOAA NGS HTDP standards.

### 2. Geodynamic Tidal Physics
- **Solid Earth Tides (SET)**: Formulated in accordance with IERS Conventions (2010), Chapter 7, utilizing degree-2 and degree-3 Love and Shida numbers driven by analytical solar and lunar ephemerides.
- **Ocean Tide Loading (OTL)**: Models 11 principal diurnal and semi-diurnal harmonic tidal constituents ($M_2, S_2, N_2, K_2, K_1, O_1, P_1, Q_1, M_f, M_m, S_{sa}$) derived from ocean loading grids.
- **Pole Tides**: Compensates for rotational elastic deformation caused by Earth's Chandler wobble and annual polar motion.

### 3. Carrier Phase Windup
- Computes phase rotation induced by relative orientation of satellite and receiver dipole antennas.
- Enforces strict subtraction invariant: $\Phi_{\text{corr}} = \Phi_{\text{meas}} - \delta\phi_{\text{windup}}$.

### 4. Local Datum Ties & Site Calibration (`site_calibration.rs`)
- `LocalDatumTie`: Estimates rigid 3D spatial translation between local CORS monument coordinates (NAD83(2011)) and global satellite orbit frames (ITRF2014).
- Eliminates apparent geodetic datum offsets (e.g. $0.29\text{ m}$ between NAD83 and ITRF2014), collapsing horizontal residuals against local RTK ground truth to the centimeter level.
