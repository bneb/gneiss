import re

def fix_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    # In dataset.rs
    if 'dataset.rs' in filepath:
        content = re.sub(r'get_sat_state\(([^,]+),\s*([^,]+),\s*([^,]+),\s*([^,]+)\)', r'get_sat_state(\1, \2, 0.0, \3, \4)', content)
    
    # In gnn_raim.rs
    if 'gnn_raim.rs' in filepath:
        content = re.sub(r'get_sat_state\(([^,]+),\s*([^,]+),\s*([^,]+),\s*([^,]+)\)', r'get_sat_state(\1, \2, 0.0, \3, \4)', content)

    # In measurement.rs
    if 'measurement.rs' in filepath:
        content = re.sub(r'fn compute_sat_state\(\s*eph[^)]+\)\s*->\s*SatState\s*\{', r'fn compute_sat_state(\n    eph: &gneiss_core::ephemeris::Ephemeris,\n    rov_obs: &DdObservation,\n    bas_obs: &DdObservation,\n    time: gneiss_core::time::GpsTime,\n    geom: &EkfGeometryContext,\n    base_time: gneiss_core::time::GpsTime,\n    rcv_clk_bias_m: f64,\n) -> SatState {', content)
        
        # fix get_sat_state calls inside compute_sat_state
        content = re.sub(r'get_sat_state\(eph, rov_obs.pr_l1, time, geom.pos_apc\)', r'get_sat_state(eph, rov_obs.pr_l1, rcv_clk_bias_m, time, geom.pos_apc)', content)
        content = re.sub(r'get_sat_state\(eph, bas_obs.pr_l1, base_time, geom.base_coord_vec\)', r'get_sat_state(eph, bas_obs.pr_l1, 0.0, base_time, geom.base_coord_vec)', content)
        
        # fix test get_sat_state calls in measurement.rs
        content = re.sub(r'get_sat_state\(&eph, pr, time, rx_pos\)', r'get_sat_state(&eph, pr, 0.0, time, rx_pos)', content)
        content = re.sub(r'get_sat_state\(&eph, 0.0, time, rx_pos\)', r'get_sat_state(&eph, 0.0, 0.0, time, rx_pos)', content)
        
        # fix compute_sat_state calls
        content = re.sub(r'compute_sat_state\(\s*sat_eph,\s*rover_sat_orig,\s*base_sat_orig,\s*time,\s*geom,\s*base_time,\s*\)', r'compute_sat_state(sat_eph, rover_sat_orig, base_sat_orig, time, geom, base_time, state.rcv_clk_bias)', content)
        content = re.sub(r'compute_sat_state\(\s*ref_eph,\s*ref_rover_orig,\s*ref_base_orig,\s*time,\s*geom,\s*env\.base_time,\s*\)', r'compute_sat_state(ref_eph, ref_rover_orig, ref_base_orig, time, geom, env.base_time, state.rcv_clk_bias)', content)
        content = re.sub(r'compute_sat_state\(\s*sat_eph,\s*rov_obs,\s*bas_obs,\s*state\.time,\s*geom,\s*env\.base_time,\s*\)', r'compute_sat_state(sat_eph, rov_obs, bas_obs, state.time, geom, env.base_time, state.rcv_clk_bias)', content)
        content = re.sub(r'compute_sat_state\(\s*ref_eph,\s*ref_rov,\s*ref_bas,\s*state\.time,\s*geom,\s*env\.base_time,\s*\)', r'compute_sat_state(ref_eph, ref_rov, ref_bas, state.time, geom, env.base_time, state.rcv_clk_bias)', content)

    # In ambiguity.rs
    if 'ambiguity.rs' in filepath:
        content = re.sub(r'get_sat_state\(\s*([^\n]+),\s*([^\n]+),\s*rover_time,\s*([^\n]+),\s*\)', r'get_sat_state(\1, \2, state.rcv_clk_bias, rover_time, \3)', content)
        content = re.sub(r'get_sat_state\(\s*([^\n]+),\s*([^\n]+),\s*base_time,\s*([^\n]+),\s*\)', r'get_sat_state(\1, \2, 0.0, base_time, \3)', content)

    with open(filepath, 'w') as f:
        f.write(content)

fix_file('crates/gneiss-rtk/src/engine/measurement.rs')
fix_file('crates/gneiss-rtk/src/engine/ambiguity.rs')
fix_file('crates/gneiss-rtk/src/engine/ml/dataset.rs')
fix_file('crates/gneiss-rtk/src/engine/ml/gnn_raim.rs')
