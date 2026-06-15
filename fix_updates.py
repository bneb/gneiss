import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

# Fix EkfUpdates::push
old_push = "pub fn push(&mut self, u: (f64, Vec<f64>, f64, u8, f64), sat: gneiss_core::sat::SatelliteId) {"
new_push = "pub fn push(&mut self, u: SingleUpdate, sat: gneiss_core::sat::SatelliteId) {"
content = content.replace(old_push, new_push)

old_push_body = """        self.z_all.push(u.0);
        self.h_all.push(u.1);
        self.r_all.push(u.2);
        self.type_all.push((sat, u.3, u.4));"""
new_push_body = """        self.z_all.push(u.z);
        self.h_all.push(u.h);
        self.r_all.push(u.r);
        self.type_all.push((sat, u.type_code, u.r_ref));"""
content = content.replace(old_push_body, new_push_body)

# Fix tests
content = content.replace("updates[0].0", "updates[0].z")
content = content.replace("updates[1].0", "updates[1].z")

# Another error mentioned:
# 140 | ) -> SingleUpdate {
# 154 |     (cp_dd - (comp_pr_dd - iono_dd + n_dd), h_cp, r_val, freq_idx, r_ref_val)
# Wait, let's find `(cp_dd - (comp_pr_dd - iono_dd + n_dd), h_cp, r_val, freq_idx, r_ref_val)`
# Ah, I replaced `+ n_dd * lam_sat` but one place has just `+ n_dd`?
content = content.replace("(cp_dd - (comp_pr_dd - iono_dd + n_dd), h_cp, r_val, freq_idx, r_ref_val)", "SingleUpdate { z: cp_dd - (comp_pr_dd - iono_dd + n_dd), h: h_cp, r: r_val, type_code: freq_idx, r_ref: r_ref_val }")

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)

