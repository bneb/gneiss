import os

def replace_in_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    original_content = content
    content = content.replace('299792458.0', 'gneiss_core::constants::SPEED_OF_LIGHT_M_S')
    content = content.replace('7.2921151467e-5', 'gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S')
    content = content.replace('6378137.0', 'gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M')
    
    # We might have introduced gneiss_core::constants:: in a file that is in gneiss_core, which is fine,
    # but inside gneiss_core it might be better to use `crate::constants::`.
    # Let's fix that if the file is in gneiss-core.
    if 'gneiss-core' in filepath:
        content = content.replace('gneiss_core::constants::', 'crate::constants::')
    
    if content != original_content:
        with open(filepath, 'w') as f:
            f.write(content)
        print(f"Updated {filepath}")

for root, dirs, files in os.walk('/Users/kevin/projects/gneiss/crates'):
    for file in files:
        if file.endswith('.rs'):
            filepath = os.path.join(root, file)
            replace_in_file(filepath)
