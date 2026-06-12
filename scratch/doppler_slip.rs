    let mut diffs = Vec::new();
    for s in sats {
        if let Some(cp) = s.cp1 {
            let dop = s.doppler;
            let band = s.sat_obs.primary_band();
            if let Some(&(prev_cp, _, prev_time)) = state.phase_history.get(&(s.sat_obs.sat, band)) {
                let dt = obs.time - prev_time;
                if dt > 0.0 && dt < 10.0 {
                    diffs.push(cp - (prev_cp - dop * dt));
                }
            }
        }
    }
    diffs.sort_by(|a: &f64, b: &f64| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median_diff = if diffs.is_empty() { 0.0 } else { diffs[diffs.len() / 2] };
