import re

def process_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Update handle_cycle_slips signature
    if 'ppp_fg' in filepath:
        old_sig = 'fn handle_cycle_slips(state: &mut RtkState, obs: &EpochObs, sats: &[ProcessedSat]) -> Vec<gneiss_core::sat::SatelliteId> {'
        new_sig = 'fn handle_cycle_slips(state: &mut RtkState, obs: &EpochObs, sats: &[ProcessedSat], spp_pos: Option<nalgebra::Vector3<f64>>, spp_cdt: Option<f64>, doppler_threshold: f64) -> Vec<gneiss_core::sat::SatelliteId> {'
    else:
        old_sig = 'fn handle_cycle_slips(state: &mut RtkState, obs: &EpochObs, sats: &[crate::engine::processed_sat::ProcessedSat]) -> std::collections::HashSet<gneiss_core::sat::SatelliteId> {'
        new_sig = 'fn handle_cycle_slips(state: &mut RtkState, obs: &EpochObs, sats: &[crate::engine::processed_sat::ProcessedSat], spp_pos: Option<nalgebra::Vector3<f64>>, spp_cdt: Option<f64>, doppler_threshold: f64) -> std::collections::HashSet<gneiss_core::sat::SatelliteId> {'
    
    content = content.replace(old_sig, new_sig)
    
    # Update Doppler threshold logic
    content = content.replace('(diff - median_diff).abs() > 2.0', '(diff - median_diff).abs() > doppler_threshold')
    
    # Update expected_p calculation
    old_expected_p = '''                let dist = (s.sat_pos_rot - state.position.vector).norm();
                let expected_p = dist + state.rcv_clk_bias - s.dt_sat_m + s.tropo_dry + state.zwd * s.map_wet + s.pcv_correction + s.iono_delay;'''
    new_expected_p = '''                let rcv_pos = spp_pos.unwrap_or(state.position.vector);
                let rcv_clk = spp_cdt.unwrap_or(state.rcv_clk_bias);
                let dist = (s.sat_pos_rot - rcv_pos).norm();
                let expected_p = dist + rcv_clk - s.dt_sat_m + s.tropo_dry + state.zwd * s.map_wet + s.pcv_correction + s.iono_delay;'''
    content = content.replace(old_expected_p, new_expected_p)
    
    # Ensure process_ppp uses the correct call
    if 'ppp_fg' in filepath:
        content = content.replace(
            'let initialized_ambiguities = handle_cycle_slips(engine.current_state.as_mut().unwrap(), rover_obs, &sats_to_process);',
            'let initialized_ambiguities = handle_cycle_slips(engine.current_state.as_mut().unwrap(), rover_obs, &sats_to_process, spp_pos, spp_cdt, engine.config.doppler_slip_threshold_cycles);'
        )
        content = content.replace(
            'let mut spp_cdt = None;',
            'let mut spp_cdt = None;\n    let mut spp_pos = None;'
        )
        content = content.replace(
            'spp_cdt = Some(s.cdt);\n            let st',
            'spp_cdt = Some(s.cdt);\n            spp_pos = Some(s.position.vector);\n            let st'
        )
        content = content.replace(
            'spp_cdt = Some(s.cdt);\n        }',
            'spp_cdt = Some(s.cdt);\n            spp_pos = Some(s.position.vector);\n        }'
        )
        content = content.replace(
            'fn process_ppp_fg<\'a>',
            'use nalgebra::Vector3;\n\npub fn process_ppp_fg<\'a>'
        )
    else:
        content = content.replace(
            'let newly_initialized = handle_cycle_slips(engine.current_state.as_mut().unwrap(), rover_obs, &sats_to_process);',
            'let newly_initialized = handle_cycle_slips(engine.current_state.as_mut().unwrap(), rover_obs, &sats_to_process, spp_pos, spp_cdt, engine.config.doppler_slip_threshold_cycles);'
        )
        content = content.replace(
            'let mut spp_cdt = None;',
            'let mut spp_cdt = None;\n    let mut spp_pos = None;'
        )
        content = content.replace(
            'spp_cdt = Some(s.cdt);\n            let st',
            'spp_cdt = Some(s.cdt);\n            spp_pos = Some(s.position.vector);\n            let st'
        )
        content = content.replace(
            'spp_cdt = Some(s.cdt);\n        }',
            'spp_cdt = Some(s.cdt);\n            spp_pos = Some(s.position.vector);\n        }'
        )
    
    with open(filepath, 'w') as f:
        f.write(content)

process_file('crates/gneiss-rtk/src/engine/ppp.rs')
process_file('crates/gneiss-rtk/src/engine/ppp_fg.rs')
