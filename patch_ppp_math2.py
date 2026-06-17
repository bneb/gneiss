import os

with open("crates/gneiss-parsers/src/sinex_bia.rs", "r") as f:
    content = f.read()

target = """    pub fn get_exact_bias(&self, sat: SatelliteId, obs: ObsCode, t: GpsTime) -> Option<f64> {
        for rec in &self.records {
            if rec.sat == sat && rec.obs1 == obs && t >= rec.start_time && t <= rec.end_time {
                return Some(rec.value);
            }
        }
        None
    }"""
if target not in content:
    target = """    pub fn get_exact_bias(&self, sat: SatelliteId, obs: ObsCode, t: GpsTime) -> Option<f64> {
        None
    }"""
    
replacement = """    pub fn get_exact_bias(&self, sat: SatelliteId, obs: ObsCode, t: GpsTime) -> Option<f64> {
        for rec in &self.records {
            if rec.sat == sat && rec.obs1 == obs && t >= rec.start_time && t <= rec.end_time {
                return Some(rec.value);
            }
        }
        None
    }"""
content = content.replace(target, replacement)
with open("crates/gneiss-parsers/src/sinex_bia.rs", "w") as f:
    f.write(content)

with open("crates/gneiss-rtk/src/engine/ppp_math.rs", "r") as f:
    content = f.read()

content = content.replace("use crate::engine::ProcessingEngine;", "")
content = content.replace("use chrono::TimeZone;", "")
content = content.replace("if l == 0 || l < prev { slip = true; lk = l; } else { lk = lk.min(l); }", "let l32 = l as u32; if l32 == 0 || l32 < prev { slip = true; lk = l32; } else { lk = lk.min(l32); }")

with open("crates/gneiss-rtk/src/engine/ppp_math.rs", "w") as f:
    f.write(content)
