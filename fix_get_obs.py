import re

content = open("crates/gneiss-rtk/src/engine/ppp.rs").read()

old_func = """        let get_obs = |sat_obs: &'a SatObs, obs_type, freq_band| {
            sat_obs.obs.iter()
                .find(|o| o.code.obs_type == obs_type && o.code.signal.freq_band == freq_band)
                .map(|o| (o.value, o.code))
        };"""

new_func = """        let get_obs = |sat_obs: &'a SatObs, obs_type, freq_band| {
            // Priority 1: Exact bias match
            let exact = sat_obs.obs.iter().find(|o| {
                o.code.obs_type == obs_type && o.code.signal.freq_band == freq_band && 
                engine.sinex_bias.as_ref().map_or(false, |s| s.get_exact_bias(sat_obs.sat, o.code, rover_obs.time).is_some())
            });
            if let Some(o) = exact {
                return Some((o.value, o.code));
            }
            
            // Priority 2: Fallback bias match
            let fallback = sat_obs.obs.iter().find(|o| {
                o.code.obs_type == obs_type && o.code.signal.freq_band == freq_band && 
                engine.sinex_bias.as_ref().map_or(false, |s| s.get_bias(sat_obs.sat, o.code, rover_obs.time).is_some())
            });
            if let Some(o) = fallback {
                return Some((o.value, o.code));
            }
            
            // Priority 3: First available
            sat_obs.obs.iter()
                .find(|o| o.code.obs_type == obs_type && o.code.signal.freq_band == freq_band)
                .map(|o| (o.value, o.code))
        };"""

content = content.replace(old_func, new_func)
open("crates/gneiss-rtk/src/engine/ppp.rs", "w").write(content)
