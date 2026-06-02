# Gneiss Workflows & Mental Model

Gneiss is a tightly-coupled GNSS/INS processing engine built in Rust. To make using Gneiss effortless, we have structured the architecture around two core concepts: **Pipelines** and **Engine Modes**.

## The 6 Core Engine Modes

At its heart, Gneiss is a massive math engine capable of fusing different types of sensor data depending on what is available. The Engine Mode dictates what algorithms are used to solve the position.

1. **SPP (Standalone Point Positioning):** Uses only the rover's GNSS receiver and broadcast satellite ephemerides. Expected accuracy: 2-5 meters.
2. **SPP-INS:** Fuses the SPP solution with IMU data in the Extended Kalman Filter (EKF) to bridge short gaps and smooth the trajectory.
3. **RTK (Real-Time Kinematic):** Fuses the rover's GNSS data with a local stationary base station to cancel out atmospheric errors. Uses integer ambiguity resolution (LAMBDA). Expected accuracy: 1-3 centimeters.
4. **RTK-INS:** Tightly couples the RTK solution with IMU measurements, providing incredible resilience in urban canyons and tunnels.
5. **PPP (Precise Point Positioning):** Uses global high-precision ephemerides and clock corrections (instead of a local base station). Takes time to converge. Expected accuracy: 10-30 centimeters.
6. **PPP-INS:** Tightly couples the PPP solution with IMU measurements.

## The 2 Pipelines

The mathematical engine is completely decoupled from how data is fed into it. We provide two "Wrappers" or "Pipelines" to feed data into the engine:

### 1. Real-Time Pipeline (`gneiss-cli live`)

The Live Pipeline is for physical deployment. It streams real-time data directly into the engine and outputs the position instantly.
- Listens to physical serial ports (UART/USB) for the Rover GNSS receiver and IMU.
- Connects to an NTRIP Caster over the internet to stream real-time Base Station corrections.
- Runs purely forward in time.

### 2. Post-Processing Pipeline (`gneiss-cli process`)

The Process Pipeline is for analyzing historical data. It reads from static files on your hard drive and processes them as fast as your CPU allows.
- Reads standard formats: `.obs` / `.ubx` for rover, `.rtcm3` / `.obs` for base, and `.nav` / `.rnx` for ephemerides.
- **Backward Smoothing:** Because all the data is available locally, the Post-Processing pipeline can optionally run the EKF backward in time (Rauch-Tung-Striebel smoothing) at the very end. This fixes errors and bridges gaps that the forward-running filter couldn't see coming.

---

## The Data Fetcher (`gneiss-cli fetch`)

Running RTK or PPP post-processing requires base station data or high-precision global ephemerides. Gathering this data manually is tedious. We provide a built-in `fetch` tool to automate this.

```bash
# Provide your rover observation file, and Gneiss will automatically figure out where and when you were
gneiss-cli fetch --rover-obs ./rover.obs --source all --out-dir ./downloads/
```

- **NOAA CORS:** Automatically finds the nearest physical base station to your rover's trajectory and downloads its Hatanaka-compressed observation file. It even invokes `crx2rnx` to seamlessly convert it to a standard RINEX `.obs` file.
- **CDDIS:** Automatically downloads the precise global broadcast ephemerides (`BRDC`) for your specific day to enable multi-constellation processing and PPP.

## Example E2E Workflow: Urban Canyon Post-Processing

1. **Fetch Data:** Grab the necessary global ephemeris and base station data based on your rover file.
   ```bash
   gneiss-cli fetch --rover-obs rover_ublox.obs --source all --out-dir .
   ```
2. **Process Data:** Feed the rover file, the base file, the ephemeris, and your IMU lever arm into the engine using `rtk-ins` mode with backward smoothing enabled.
   ```bash
   gneiss-cli process \
       -r rover_ublox.obs \
       -b base_trimble.obs \
       -n BRDC00IGS_R_20231000000_01D_MN.rnx \
       -o trajectory.pos \
       --mode rtk-ins \
       --lever-arm "0.1,0.2,-0.1" \
       --enable-backward-smoothing
   ```
3. **Evaluate:** (Optional) Evaluate the trajectory against a ground-truth reference.
   ```bash
   gneiss-cli eval --solution trajectory.pos --truth reference.pos
   ```
