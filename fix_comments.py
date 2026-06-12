import os

replacements = {
    "crates/gneiss-rtk/src/engine/ppp_fg.rs": [
        ('// But handle_cycle_slips runs AFTER this block!', '// handle_cycle_slips runs after this block.'),
        ('// Actually, let\'s just let handle_cycle_slips detect cycle slips!', '// Rely on handle_cycle_slips to detect cycle slips.'),
        ('// Let\'s just calculate median phase diff here!', '// Calculate median phase diff here.'),
        ('// Phase did NOT jump, but PR did! So we must compensate ambiguities.', '// Phase did not jump, but PR did, so we compensate ambiguities.'),
        ('// ONLY if ambiguities were compensated!', '// Only execute if ambiguities were compensated.')
    ],
    "crates/gneiss-rtk/src/engine/tests_updater.rs": [
        ('// The worst outlier (index 5) MUST be rejected!', '// The worst outlier (index 5) should be rejected.'),
        ('// NO measurements should be rejected because they are consistent with a clock jump!', '// No measurements should be rejected because they are consistent with a clock jump.')
    ],
    "crates/gneiss-rtk/src/engine/tight_fg.rs": [
        ('// Initialize the Factor Graph search at the SPP position and clock!', '// Initialize the Factor Graph search at the SPP position and clock.')
    ],
    "crates/gneiss-rtk/src/engine/updater.rs": [
        ('// We should reject everything so the filter resets to SPP!', '// Reject everything so the filter resets to SPP.'),
        ('// Recompute S_ii because we need the pre-fit innovation variance for phase!', '// Recompute S_ii to obtain the pre-fit innovation variance for phase.'),
        ('// Phase: NEVER reject based on absolute pre-fit residual, rely on ambiguity resolution!', '// Phase: Do not reject based on absolute pre-fit residual; rely on ambiguity resolution.')
    ],
    "crates/gneiss-rtk/src/factor_graph/gnss_factors.rs": [
        ('let huber_cp = 3.0; // 3-sigma (30cm) to quickly reject cycle slips!', 'let huber_cp = 3.0; // 3-sigma (30cm) to reject cycle slips.'),
        ('Some(3.0) // 3-sigma threshold to quickly reject cycle slips!', 'Some(3.0) // 3-sigma threshold to reject cycle slips.')
    ],
    "crates/gneiss-rtk/src/factor_graph/imu_factors.rs": [
        ('// Finite difference for now to ensure correctness!', '// Finite difference used to ensure correctness.')
    ],
    "crates/gneiss-rtk/src/spp.rs": [
        ('// But the important part is that we run the functions to catch mutations!', '// Ensure functions are run to catch mutations.')
    ]
}

for filepath, reps in replacements.items():
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        continue
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()
    
    for old, new in reps:
        if old in content:
            content = content.replace(old, new)
        else:
            print(f"Warning: could not find target text in {filepath}:\n{old}")
            
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)

print("Comment text replacements complete.")
