use nalgebra::DMatrix;

fn main() {
    let dt_abs = 0.2;
    let q_acc = 1.0;
    
    let mut q = DMatrix::<f64>::zeros(15, 15);
    let q_pos = q_acc * dt_abs.powi(3) / 3.0; 
    let q_vel = q_acc * dt_abs;
    let q_pos_vel = q_acc * dt_abs.powi(2) / 2.0;
    for i in 0..3 { 
        q[(i, i)] = q_pos; 
        q[(i+3, i+3)] = q_vel; 
        q[(i, i+3)] = q_pos_vel;
        q[(i+3, i)] = q_pos_vel;
    }
    
    let mut p = DMatrix::<f64>::zeros(15, 15);
    p[(0,0)] = 0.0088;
    p[(3,3)] = 2.082;
    
    let mut phi = DMatrix::<f64>::identity(15, 15);
    for i in 0..3 { phi[(i, 3 + i)] = dt_abs; }
    
    let p_new = &phi * &p * phi.transpose() + &q;
    println!("P_pred pos var: {:.3e}, vel var: {:.3e}", p_new[(0, 0)], p_new[(3, 3)]);
}
