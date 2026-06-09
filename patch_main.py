import re

with open("bin/gneiss-cli/src/main.rs", "r") as f:
    content = f.read()

interpolation_code = """
                    let b = if let Some(base_measurements) = b_epochs {
                        let r_t = r.time.tow;
                        let mut b1 = None;
                        let mut b2 = None;
                        
                        for base in base_measurements {
                            if base.time.tow <= r_t {
                                b1 = Some(base);
                            } else if b2.is_none() {
                                b2 = Some(base);
                                break;
                            }
                        }
                        
                        if let (Some(base1), Some(base2)) = (b1, b2) {
                            let t1 = base1.time.tow;
                            let t2 = base2.time.tow;
                            if t2 > t1 && r_t - t1 < engine.config.max_base_age_s && t2 - r_t < engine.config.max_base_age_s {
                                let alpha = (r_t - t1) / (t2 - t1);
                                let mut interp_obs = base1.clone();
                                interp_obs.time = r.time;
                                
                                for sat1 in &mut interp_obs.satellites {
                                    if let Some(sat2) = base2.satellites.iter().find(|s| s.sat == sat1.sat) {
                                        for obs1 in &mut sat1.observations {
                                            if let Some(obs2) = sat2.observations.iter().find(|o| o.code == obs1.code) {
                                                obs1.value = obs1.value + alpha * (obs2.value - obs1.value);
                                            }
                                        }
                                    }
                                }
                                Some(interp_obs)
                            } else {
                                // Fallback to nearest
                                let mut closest_base = None;
                                let mut min_diff = engine.config.max_base_age_s;
                                for base in base_measurements {
                                    let diff = (base.time.tow - r_t).abs();
                                    if diff < min_diff {
                                        min_diff = diff;
                                        closest_base = Some(base.clone());
                                    }
                                }
                                closest_base
                            }
                        } else {
                            // Fallback to nearest
                            let mut closest_base = None;
                            let mut min_diff = engine.config.max_base_age_s;
                            for base in base_measurements {
                                let diff = (base.time.tow - r_t).abs();
                                if diff < min_diff {
                                    min_diff = diff;
                                    closest_base = Some(base.clone());
                                }
                            }
                            closest_base
                        }
                    } else { None };
"""

pattern = r"let b = if let Some\(base_measurements\) = b_epochs \{.*?closest_base\n                    \} else \{ None \};"

new_content = re.sub(pattern, interpolation_code.strip(), content, flags=re.DOTALL)

with open("bin/gneiss-cli/src/main.rs", "w") as f:
    f.write(new_content)
