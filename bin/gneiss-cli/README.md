# gneiss-cli

The primary interface for the Gneiss engine. This tool provides a clinical environment for processing and evaluating raw navigation data.

## Overview

`gneiss-cli` is where the math meets the metal. It’s designed for predictable execution and clinical validation.

### Example Workflow

```bash
# Generate a trajectory
gneiss process \
    --rover rover_data.ubx \
    --base base_stream.rtcm3 \
    --output path.pos \
    --enable-imu-fusion
```

---

*Gneiss-cli: A clinical interface for a predictable engine.*
