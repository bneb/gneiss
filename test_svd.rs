use nalgebra::{DMatrix, DVector};
fn main() {
    let h = DMatrix::<f64>::zeros(2, 2);
    let b = DVector::<f64>::zeros(2);
    let svd = h.svd(true, true);
    let sol = svd.solve(&b, 1e-9).unwrap();
}
