use nalgebra::Vector3;

pub fn ecef_to_llh(ecef: Vector3<f64>) -> Vector3<f64> {
    let a = 6378137.0;
    let b = 6356752.314245;
    let e2 = 1.0 - (b * b) / (a * a);
    let p = f64::sqrt(ecef.x * ecef.x + ecef.y * ecef.y);
    let lon = f64::atan2(ecef.y, ecef.x);
    let mut lat = f64::atan2(ecef.z, p * (1.0 - e2));
    let mut h = 0.0;
    for _ in 0..5 {
        let sin_lat = f64::sin(lat);
        let n = a / f64::sqrt(1.0 - e2 * sin_lat * sin_lat);
        h = p / f64::cos(lat) - n;
        lat = f64::atan2(ecef.z, p * (1.0 - e2 * n / (n + h)));
    }
    Vector3::new(lat, lon, h)
}

pub fn llh_to_ecef(llh: Vector3<f64>) -> Vector3<f64> {
    let a = 6378137.0;
    let b = 6356752.314245;
    let e2 = 1.0 - (b * b) / (a * a);
    let sin_lat = f64::sin(llh.x);
    let cos_lat = f64::cos(llh.x);
    let sin_lon = f64::sin(llh.y);
    let cos_lon = f64::cos(llh.y);
    let n = a / f64::sqrt(1.0 - e2 * sin_lat * sin_lat);
    let x = (n + llh.z) * cos_lat * cos_lon;
    let y = (n + llh.z) * cos_lat * sin_lon;
    let z = (n * (1.0 - e2) + llh.z) * sin_lat;
    Vector3::new(x, y, z)
}

fn main() {
    let sol = Vector3::new(-2694595.7929, -4296531.1950, 3854851.5973);
    let sol_llh = ecef_to_llh(sol);
    println!("Sol LLH: lat={}, lon={}, h={}", sol_llh.x.to_degrees(), sol_llh.y.to_degrees(), sol_llh.z);
    
    let truth_lat = 37.4235759540f64;
    let truth_lon = -122.0941320350f64;
    let truth_h = 33.21;
    let truth_ecef = llh_to_ecef(Vector3::new(truth_lat.to_radians(), truth_lon.to_radians(), truth_h));
    println!("Truth ECEF: x={}, y={}, z={}", truth_ecef.x, truth_ecef.y, truth_ecef.z);
}
