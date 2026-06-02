# Future Technologies and Algorithmic Analogies for GNSS/RTK Engines

## 1. Marginal Value Add of Satellites
**Question:** What is the marginal value add of each satellite added to the solved subset? Does it diminish after a certain point (maybe 4 or 5)?

**Analysis:**
- **The Mathematical Minimum:** Solving a 3D position and receiver clock bias (x, y, z, cdt) requires exactly 4 satellites. A 5th satellite provides the first degree of over-determination (1 degree of freedom for residuals), allowing the first form of Receiver Autonomous Integrity Monitoring (RAIM) to detect a single outlier.
- **Dilution of Precision (DOP):** The value of each additional satellite heavily depends on its geometric distribution in the sky. If you have 4 satellites clustered together, adding a 5th satellite in a completely different part of the sky has immense marginal value, shrinking the position covariance ellipsoid drastically. Conversely, if you already have 10 well-distributed satellites, adding an 11th provides sharply diminishing returns for geometry.
- **Diminishing Returns:** Generally, the marginal reduction in covariance scales inversely with the square root of the number of satellites ($1/\sqrt{n}$). Therefore, the leap from 4 to 6 satellites is massive in terms of both precision and fault tolerance. The leap from 10 to 12 provides a much smaller precision benefit, though it marginally aids ambiguity resolution in RTK by providing more integer combinations.

## 2. Dynamic Environment Detection (Open Sky vs. Urban Canyon)
**Question:** Is there a way to detect whether we are in an open sky or in an urban canyon? If so, could we fluctuate the number of satellites to use to reduce processing?

**Analysis:**
- **Detection Strategies:**
  1. **C/N0 (SNR) Variance Analysis:** In open sky, signal strengths (C/N0) smoothly follow an elevation-dependent curve. In an urban canyon, C/N0 values are highly erratic and generally lower due to multipath and blockages.
  2. **Elevation Angle Distribution:** In an urban canyon, low-elevation satellites are physically blocked by buildings. If the engine consistently loses lock on satellites below 30-40 degrees elevation but maintains high-elevation satellites, it's a strong indicator of an urban canyon.
  3. **Innovation Variance:** In the EKF, if the measurement innovations (residuals) suddenly spike in variance, we have entered a multipath-heavy environment (urban canyon or tree canopy).
  4. **Vision/LiDAR Fusion:** Using a sky-facing camera (fisheye) to perform semantic segmentation on the sky vs. buildings.
- **Fluctuating the Satellite Subset:**
  Yes, this is an advanced technique known as **Measurement Selection** or **Subset Pruning**.
  - **In Open Sky:** We can stack-rank satellites by a combination of high elevation and high C/N0, selecting the "best" 6-8 satellites to process. This vastly reduces the search space for Integer Ambiguity Resolution (LAMBDA method), saving CPU cycles and battery.
  - **In Urban Canyon:** The strategy flips. You need *every* satellite you can get, but you aggressively de-weight them based on SNR, or run advanced algorithms (like Factor Graphs or particle filters) to handle the non-Gaussian multipath errors.

## 3. State-of-the-Art Analogies to Other CS Branches
**Question:** Are there any analogy bridges to other branches of computer science that can inform our data structures and workflows?

**1. Robotics and SLAM (Factor Graphs)**
- **Analogy:** Traditional RTK uses the Extended Kalman Filter (EKF), which only maintains the *current* state and marginalizes out the past. This is brittle to multipath because an erroneous measurement permanently corrupts the state.
- **Bridge:** Modern SLAM (Simultaneous Localization and Mapping) in robotics uses **Factor Graphs** (e.g., GTSAM or Ceres Solver) to maintain a window of past states. If we port RTK to a Factor Graph, we can perform non-linear smoothing over a sliding time window. If a satellite is later determined to be a multipath reflection, the graph optimizer can retroactively remove it and instantly fix the entire trajectory.

**2. Machine Learning / Graph Neural Networks (GNNs)**
- **Analogy:** Deciding which satellites have line-of-sight vs. multipath is a classification problem based on temporal SNR, elevation, and residual features.
- **Bridge:** A Graph Neural Network (GNN) can treat the constellation of satellites as a dynamic graph (where edges represent geometric correlations or shared atmospheric paths). The GNN can infer the multipath probability of each satellite in real-time, providing highly accurate variance weights to the EKF/Factor Graph, far outperforming traditional empirical models like elevation-based weighting.

**3. Database Query Optimization (Measurement Selection)**
- **Analogy:** A database query optimizer prunes the search space of joins by evaluating statistics.
- **Bridge:** In RTK ambiguity resolution (LAMBDA), the search space grows factorially with the number of satellites. We can use heuristic search algorithms (A* or Monte Carlo Tree Search) to dynamically select the subset of satellites that provides the optimal trade-off between Geometry (PDOP) and Integer Success Rate, rather than blindly attempting to fix all visible satellites.

**4. Networking / Information Theory (Redundancy & Checksums)**
- **Analogy:** RAIM (Receiver Autonomous Integrity Monitoring) is fundamentally an error-detecting/correcting code applied to geometry.
- **Bridge:** Advanced fault detection and exclusion (FDE) can borrow from Byzantine Fault Tolerance (BFT) in distributed systems. When processing measurements from multiple GNSS constellations (GPS, Galileo, BeiDou), we can treat each constellation as an independent "node" proposing a solution. If GPS and Galileo agree but GLONASS disagrees wildly, we can quarantine the GLONASS measurements.
