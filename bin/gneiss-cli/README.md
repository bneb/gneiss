# gneiss-cli

`gneiss-cli` is the command-line interface for the Gneiss positioning engine, providing tools for batch post-processing, real-time RTK streaming, accuracy evaluation, camera event interpolation, and local datum site calibration.

## Subcommands

### 1. `process`
Execute single-baseline RTK/PPK, Network RTK (VRS), or PPP post-processing:
```bash
# Post-process kinematic rover against base station with 2 calibration passes
gneiss process \
    --rover rover.obs \
    --base base.obs \
    --output trajectory.pos \
    --calibrate-passes 2

# Precise Point Positioning (PPP) with precise orbits and clocks
gneiss process \
    --rover rover.obs \
    --mode ppp \
    --sp3 cod23200.sp3 \
    --clk cod23200.clk \
    --antex igs14.atx \
    --output ppp_solution.pos
```

### 2. `calibrate`
Compute local site calibration parameters (7-parameter Helmert or 3D rigid translation) from paired GNSS and ground control point coordinates:
```bash
gneiss calibrate \
    --input control_points.csv \
    --output calibration.json
```

### 3. `batch`
Process an entire directory of rover sessions against common base stations:
```bash
gneiss batch \
    --input-dir ./rover_sessions/ \
    --base base.obs \
    --output-dir ./solutions/ \
    --format pos
```

### 4. `eval`
Compute empirical Cumulative Distribution Functions (CDF), percentiles ($p_{50}, p_{68}, p_{95}$), and RMS against ground truth trajectories:
```bash
gneiss eval \
    --solution trajectory.pos \
    --truth ground_truth.csv
```

### 5. `events`
Interpolate exact camera shutter event positions from raw trajectory solutions for UAV photogrammetry:
```bash
gneiss events \
    --trajectory trajectory.pos \
    --events shutter_events.txt \
    --output camera_positions.csv \
    --lever-arm "0.05,-0.02,0.15"
```

### 6. `gui`
Launch the local interactive Web diagnostic server to inspect residuals, satellite skyplots, and solution trajectories:
```bash
gneiss gui --port 8080 --trajectory trajectory.pos
```

### 7. `live`
Stream real-time differential corrections via serial receiver or NTRIP:
```bash
gneiss live \
    --rover /dev/ttyACM0 \
    --mountpoint RTK_VRS
```
