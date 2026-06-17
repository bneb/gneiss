use std::collections::HashMap;
use nalgebra::Vector3;

#[derive(Debug)]
pub struct Sp3Record {
    pub position: Vector3<f64>,
}

#[derive(Debug)]
pub struct Sp3Epoch {
    pub time: f64,
    pub records: HashMap<String, Sp3Record>,
}

pub fn interpolate_orbit_lagrange(points: &[(f64, Vector3<f64>)], target: f64) -> Option<Vector3<f64>> {
    let n = points.len();
    if n == 0 { return None; }
    if n == 1 { return Some(points[0].1); }

    let mut result = Vector3::zeros();
    for i in 0..n {
        let mut term = points[i].1;
        for j in 0..n {
            if i != j {
                let num = target - points[j].0;
                let den = points[i].0 - points[j].0;
                if den == 0.0 { continue; }
                term *= num / den;
            }
        }
        result += term;
    }
    Some(result)
}
