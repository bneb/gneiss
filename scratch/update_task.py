with open("/Users/kevin/.gemini/antigravity/brain/ca10468e-f8ed-47bf-a69d-e27a6dfafa12/task.md", "r") as f:
    content = f.read()

content = content.replace("- `[/]` Stabilize Tightly Coupled EKF (IMU Fusion)", "- `[x]` Stabilize Tightly Coupled EKF (IMU Fusion)")
content = content.replace("  - `[ ]` Review and fix Measurement Jacobians in `engine/measurement.rs` (especially `H_att` for lever arm).", "  - `[x]` Review and fix Measurement Jacobians in `engine/measurement.rs` (especially `H_att` for lever arm).")
content = content.replace("  - `[ ]` Review and fix Preintegration Jacobians (`factor_graph/imu_factors.rs`).", "  - `[x]` Review and fix Preintegration Jacobians (`factor_graph/imu_factors.rs`).")
content = content.replace("- `[ ]` Overhaul PPP Models", "- `[/]` Overhaul PPP Models")

with open("/Users/kevin/.gemini/antigravity/brain/ca10468e-f8ed-47bf-a69d-e27a6dfafa12/task.md", "w") as f:
    f.write(content)
