fn main() {
    let rot = nalgebra::Rotation3::from_euler_angles(0.0, 0.0, std::f64::consts::PI / 2.0);
    let v = nalgebra::Vector3::new(1.0, 0.0, 0.0);
    println!("RotZ(90) * [1, 0, 0] = {:?}", rot * v);
}
