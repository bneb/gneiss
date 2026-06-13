use nalgebra::DMatrix;

fn main() {
    let mut m = DMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![1.0, 2.0, 3.0]));
    m = m.remove_row(1).remove_column(1);
    assert_eq!(m.nrows(), 2);
    assert_eq!(m[(1, 1)], 3.0);
    
    m = m.insert_row(1, 0.0).insert_column(1, 0.0);
    m[(1, 1)] = 5.0;
    assert_eq!(m.nrows(), 3);
    assert_eq!(m[(1, 1)], 5.0);
    println!("OK");
}
