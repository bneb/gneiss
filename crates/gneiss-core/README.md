# gneiss-core

`gneiss-core` provides the foundational data structures, reference frame realizations, time systems, and physical models for the Gneiss positioning engine.

## Key Subsystems

### 1. Reference Frames & Realizations (`frames/`)
- Type-safe reference frame tags and 14-parameter time-dependent Helmert transformations with plate tectonic velocity models.
- Realizations include `ITRF2014`, `ITRF2020`, `IGS14`, `IGS20`, `WGS84`, `Nad83_2011`, `Etrs89`, and `Gda2020`.
- Mathematical alignment with IOGP EPSG:8970 (Method 1056) and NOAA NGS HTDP standards.

### 2. Coordinate Types & Geometric Transformations (`coordinates/`)
- Strongly-typed coordinate containers: `EcefPos`, `GeodeticPos`, `NedPos`, and `EpochPosition`.
- WGS84, GRS80, and IERS reference ellipsoid transformations with zero heap allocation.
- Local tangent plane projections (North-East-Down, East-North-Up).

### 3. Time Systems & Ephemerides (`time/`, `ephemeris/`)
- High-precision time representations (`GpsTime`, UTC conversions, leap-second offsets).
- Broadcast navigation message decoders and satellite orbit/clock propagation routines.
- Transmit-time iteration with Earth rotation (Sagnac effect) compensation.

### 4. Physical & Atmospheric Models
- Tropospheric delay modeling: Saastamoinen zenith hydrostatic and wet delays with Vienna/Niell mapping functions.
- Ionospheric delay estimation: Dual-frequency ionosphere-free linear combinations and Klobuchar single-frequency broadcast model.
- IERS 2010 Solid Earth Tide models driven by solar and lunar ephemerides.
