import re

content = open("crates/gneiss-rtk/src/engine/ppp_fg.rs").read()

# We need to change `expected_base` and `h_func` to use `pcv_l1` for L1 and `pcv_l2` for L2!

old_uduc_block = """            if !sat.is_iono_free && sat.cp2.is_some() && sat.p2.is_some() {
                // UDUC dual-frequency
                let p2 = sat.p2.unwrap();
                let p1 = sat.p1;
                let gamma = (sat.f1 * sat.f1) / (sat.f2 * sat.f2);
                let mut i1_est = (p2 - p1) / (gamma - 1.0);
                
                let cp2_cyc = sat.cp2.unwrap();
                let cp1_cyc = sat.cp1.unwrap();
                
                meas.push(FgMeasurement { sat: sat.sat_obs.sat, val: p1, expected: expected_base + i1_est, var: r_p1,
                    h_func: Box::new(move |x| {
                        let isb = Self::extract_isb(x, sat_c);
                        let ztd = if x.len() > 20 && !x[20].is_nan() && x[20] != 0.0 { x[20] } else { state_zwd };
                        let rcv_pos = Vector3::new(x[0], x[1], x[2]) + tide_offset;
                        let dist = (sat_pos - rcv_pos).norm() - pcv_correction;
                        dist + x[15] + isb - dt_sat_m + tropo_dry + ztd * map_wet + x[CORE_STATE_SIZE + n_i1]
                    })
                });

                meas.push(FgMeasurement { sat: sat.sat_obs.sat, val: p2, expected: expected_base + gamma * i1_est, var: r_p2,
                    h_func: Box::new(move |x| {
                        let isb = Self::extract_isb(x, sat_c);
                        let ztd = if x.len() > 20 && !x[20].is_nan() && x[20] != 0.0 { x[20] } else { state_zwd };
                        let rcv_pos = Vector3::new(x[0], x[1], x[2]) + tide_offset;
                        let dist = (sat_pos - rcv_pos).norm() - pcv_correction;
                        dist + x[15] + isb - dt_sat_m + tropo_dry + ztd * map_wet + gamma * x[CORE_STATE_SIZE + n_i1]
                    })
                });

                meas.push(FgMeasurement { sat: sat.sat_obs.sat, val: cp1_cyc * lam1, expected: expected_base - i1_est + x_i[CORE_STATE_SIZE + n1], var: r_cp1,
                    h_func: Box::new(move |x| {
                        let isb = Self::extract_isb(x, sat_c);
                        let ztd = if x.len() > 20 && !x[20].is_nan() && x[20] != 0.0 { x[20] } else { state_zwd };
                        let rcv_pos = Vector3::new(x[0], x[1], x[2]) + tide_offset;
                        let dist = (sat_pos - rcv_pos).norm() - pcv_correction;
                        dist + x[15] + isb - dt_sat_m + tropo_dry + ztd * map_wet - x[CORE_STATE_SIZE + n_i1] + x[CORE_STATE_SIZE + n1]
                    })
                });
                
                let i2_est = gamma * i1_est;
                meas.push(FgMeasurement { sat: sat.sat_obs.sat, val: cp2_cyc * lam2, expected: expected_base - i2_est + x_i[CORE_STATE_SIZE + n2], var: r_cp2,
                    h_func: Box::new(move |x| {
                        let isb = Self::extract_isb(x, sat_c);
                        let ztd = if x.len() > 20 && !x[20].is_nan() && x[20] != 0.0 { x[20] } else { state_zwd };
                        let rcv_pos = Vector3::new(x[0], x[1], x[2]) + tide_offset;
                        let dist = (sat_pos - rcv_pos).norm() - pcv_correction;
                        dist + x[15] + isb - dt_sat_m + tropo_dry + ztd * map_wet - gamma * x[CORE_STATE_SIZE + n_i1] + x[CORE_STATE_SIZE + n2]
                    })
                });"""

new_uduc_block = """            if !sat.is_iono_free && sat.cp2.is_some() && sat.p2.is_some() {
                // UDUC dual-frequency
                let p2 = sat.p2.unwrap();
                let p1 = sat.p1;
                let gamma = (sat.f1 * sat.f1) / (sat.f2 * sat.f2);
                let mut i1_est = (p2 - p1) / (gamma - 1.0);
                
                let cp2_cyc = sat.cp2.unwrap();
                let cp1_cyc = sat.cp1.unwrap();
                
                let pcv_l1 = sat.pcv_l1;
                let pcv_l2 = sat.pcv_l2;

                let expected_base_l1 = dist + x_i[15] + isb - sat.dt_sat_m + sat.tropo_dry + ztd * sat.map_wet - pcv_l1;
                let expected_base_l2 = dist + x_i[15] + isb - sat.dt_sat_m + sat.tropo_dry + ztd * sat.map_wet - pcv_l2;

                meas.push(FgMeasurement { sat: sat.sat_obs.sat, val: p1, expected: expected_base_l1 + i1_est, var: r_p1,
                    h_func: Box::new(move |x| {
                        let isb = Self::extract_isb(x, sat_c);
                        let ztd = if x.len() > 20 && !x[20].is_nan() && x[20] != 0.0 { x[20] } else { state_zwd };
                        let rcv_pos = Vector3::new(x[0], x[1], x[2]) + tide_offset;
                        let dist = (sat_pos - rcv_pos).norm() - pcv_l1;
                        dist + x[15] + isb - dt_sat_m + tropo_dry + ztd * map_wet + x[CORE_STATE_SIZE + n_i1]
                    })
                });

                meas.push(FgMeasurement { sat: sat.sat_obs.sat, val: p2, expected: expected_base_l2 + gamma * i1_est, var: r_p2,
                    h_func: Box::new(move |x| {
                        let isb = Self::extract_isb(x, sat_c);
                        let ztd = if x.len() > 20 && !x[20].is_nan() && x[20] != 0.0 { x[20] } else { state_zwd };
                        let rcv_pos = Vector3::new(x[0], x[1], x[2]) + tide_offset;
                        let dist = (sat_pos - rcv_pos).norm() - pcv_l2;
                        dist + x[15] + isb - dt_sat_m + tropo_dry + ztd * map_wet + gamma * x[CORE_STATE_SIZE + n_i1]
                    })
                });

                meas.push(FgMeasurement { sat: sat.sat_obs.sat, val: cp1_cyc * lam1, expected: expected_base_l1 - i1_est + x_i[CORE_STATE_SIZE + n1], var: r_cp1,
                    h_func: Box::new(move |x| {
                        let isb = Self::extract_isb(x, sat_c);
                        let ztd = if x.len() > 20 && !x[20].is_nan() && x[20] != 0.0 { x[20] } else { state_zwd };
                        let rcv_pos = Vector3::new(x[0], x[1], x[2]) + tide_offset;
                        let dist = (sat_pos - rcv_pos).norm() - pcv_l1;
                        dist + x[15] + isb - dt_sat_m + tropo_dry + ztd * map_wet - x[CORE_STATE_SIZE + n_i1] + x[CORE_STATE_SIZE + n1]
                    })
                });
                
                let i2_est = gamma * i1_est;
                meas.push(FgMeasurement { sat: sat.sat_obs.sat, val: cp2_cyc * lam2, expected: expected_base_l2 - i2_est + x_i[CORE_STATE_SIZE + n2], var: r_cp2,
                    h_func: Box::new(move |x| {
                        let isb = Self::extract_isb(x, sat_c);
                        let ztd = if x.len() > 20 && !x[20].is_nan() && x[20] != 0.0 { x[20] } else { state_zwd };
                        let rcv_pos = Vector3::new(x[0], x[1], x[2]) + tide_offset;
                        let dist = (sat_pos - rcv_pos).norm() - pcv_l2;
                        dist + x[15] + isb - dt_sat_m + tropo_dry + ztd * map_wet - gamma * x[CORE_STATE_SIZE + n_i1] + x[CORE_STATE_SIZE + n2]
                    })
                });"""

if old_uduc_block in content:
    content = content.replace(old_uduc_block, new_uduc_block)
    open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w").write(content)
    print("Fixed!")
else:
    print("Could not find the block to replace!")
