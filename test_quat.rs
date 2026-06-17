use nalgebra::{UnitQuaternion, Vector3};

fn main() {
    let q = UnitQuaternion::identity();
    let dq = Vector3::new(0.1, 0.2, 0.3);
    let dq_quat = UnitQuaternion::new(dq);
    let q_new = dq_quat * q;
}
