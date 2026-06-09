import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

target = """    let mut current_state = if let Some(state) = seed_state {
        // ... previous logic to start from seed state
        // To be safe we should re-compute a clean seed
"""
# That was probably wrong, let's find `compute_spp`
target2 = """pub fn compute_spp"""

# Let's just sed replace the Returns of `Err` in compute_spp
target_err1 = "return Err(SppError::NotEnoughMeasurements);"
replacement_err1 = "tracing::warn!(\"SPP Error: NotEnoughMeasurements\"); return Err(SppError::NotEnoughMeasurements);"

target_err2 = "return Err(SppError::ConvergenceFailed);"
replacement_err2 = "tracing::warn!(\"SPP Error: ConvergenceFailed\"); return Err(SppError::ConvergenceFailed);"

target_err3 = "return Err(SppError::PoorGeometry);"
replacement_err3 = "tracing::warn!(\"SPP Error: PoorGeometry\"); return Err(SppError::PoorGeometry);"

target_err4 = ".ok_or(SppError::MatrixInversionFailed)?;"
replacement_err4 = ".ok_or_else(|| { tracing::warn!(\"SPP Error: MatrixInversionFailed\"); SppError::MatrixInversionFailed })?;"

content = content.replace(target_err1, replacement_err1)
content = content.replace(target_err2, replacement_err2)
content = content.replace(target_err3, replacement_err3)
content = content.replace(target_err4, replacement_err4)

with open(filename, "w") as f:
    f.write(content)
