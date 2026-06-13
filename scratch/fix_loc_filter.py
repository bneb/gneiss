import re

with open("crates/gneiss-rtk/src/filter.rs", "r") as f:
    content = f.read()

replacement = """fn filter_by_locktime(
    state: &RtkState, constellations: &[gneiss_core::sat::Constellation], ar_min_lock: u32,
) -> Vec<(usize, usize, u16)> {
    let mut candidates = Vec::new();
    for &constell in constellations {
        if let Some(ref_idx) = find_best_reference_sat(state, constell, ar_min_lock) {
            collect_candidates_for_constellation(state, constell, ref_idx, ar_min_lock, &mut candidates);
        }
    }
    candidates
}

fn find_best_reference_sat(state: &RtkState, constell: gneiss_core::sat::Constellation, ar_min_lock: u32) -> Option<usize> {
    let mut best_ref_idx = None;
    let mut max_lock = 0;
    for i in 0..state.ambiguities.len() {
        let (sat, freq) = state.ambiguity_keys[i];
        if sat.constellation != constell || freq != 1 { continue; }
        let lock = *state.locktimes.get(&(sat, freq)).unwrap_or(&0);
        if lock >= ar_min_lock as u16 && lock > max_lock {
            max_lock = lock; best_ref_idx = Some(i);
        }
    }
    best_ref_idx
}

fn collect_candidates_for_constellation(
    state: &RtkState, constell: gneiss_core::sat::Constellation, ref_idx: usize, ar_min_lock: u32, candidates: &mut Vec<(usize, usize, u16)>
) {
    let ref_sat_id = state.ambiguity_keys[ref_idx].0;
    let l2_ref_idx = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_sat_id && f == 2);
    for i in 0..state.ambiguities.len() {
        if i == ref_idx || Some(i) == l2_ref_idx { continue; }
        let (rov_sat, freq) = state.ambiguity_keys[i];
        if rov_sat.constellation != constell { continue; }
        let lock = *state.locktimes.get(&(rov_sat, freq)).unwrap_or(&0);
        if lock >= ar_min_lock as u16 {
            if freq == 1 { candidates.push((i, ref_idx, lock)); }
            else if let Some(r2_idx) = l2_ref_idx { candidates.push((i, r2_idx, lock)); }
        }
    }
}"""

start_idx = content.find("fn filter_by_locktime(")
if start_idx == -1:
    print("Could not find filter_by_locktime")
    import sys; sys.exit(1)

end_idx = content.find("fn compute_candidate_variance(", start_idx)

new_content = content[:start_idx] + replacement + "\n\n" + content[end_idx:]

with open("crates/gneiss-rtk/src/filter.rs", "w") as f:
    f.write(new_content)

print("Fixed LOC in filter.rs")
