use nalgebra::{Matrix2, Cholesky};
fn main() {
    let a = Matrix2::new(4.0, 1.0, 1.0, 3.0);
    let chol = a.clone().cholesky().unwrap();
    let inv = chol.inverse();
    let inv2 = a.try_inverse().unwrap();
    println!("chol.inverse: {}", inv);
    println!("try_inverse: {}", inv2);
}
