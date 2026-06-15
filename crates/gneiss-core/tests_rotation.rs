fn main() {}
#[test]
fn test_rot() {
    let angles = [0.0, 0.0, 3.141592653589793];
    let r_m_v = nalgebra::Rotation3::from_euler_angles(angles[0], angles[1], angles[2]);
    let meas_accel = nalgebra::Vector3::new(-0.341219, -0.060000, -9.800251);
    let result = r_m_v * meas_accel;
    println!("result: {:?}", result);
}
