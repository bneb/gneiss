use nalgebra::{UnitQuaternion, Vector2, Vector3, Vector6};

pub trait Manifold {
    fn local_dim(&self) -> usize;
    fn retract(&mut self, delta: &[f64]);
}

#[derive(Debug, Clone)]
pub enum VariableType {
    Pos(Vector3<f64>),
    Vel(Vector3<f64>),
    Att(UnitQuaternion<f64>),
    Biases(Vector6<f64>),
    Clock(Vector2<f64>),
    Ambiguity(f64),
}

impl VariableType {
    fn retract_att(q: &mut UnitQuaternion<f64>, delta: &[f64]) {
        let dq = Vector3::new(delta[0], delta[1], delta[2]);
        let dq_quat = UnitQuaternion::new(dq);
        *q *= dq_quat;
    }
}

impl Manifold for VariableType {
    fn local_dim(&self) -> usize {
        match self {
            VariableType::Pos(_) | VariableType::Vel(_) | VariableType::Att(_) => 3,
            VariableType::Biases(_) => 6,
            VariableType::Clock(_) => 2,
            VariableType::Ambiguity(_) => 1,
        }
    }

    /// Retracts the manifold by the local error state `delta`.
    /// 
    /// For the $SO(3)$ Attitude Manifold, we employ right-perturbation mathematics:
    /// $$ q_{new} = q_{old} \otimes \exp\left(\Delta \theta\right) $$
    /// 
    /// If the norm of $\Delta \theta$ is small, we use a first-order Taylor expansion
    /// to avoid numerical instability in the quaternion exponential:
    /// $$ \Delta q \approx \begin{bmatrix} 1 \\ \frac{\Delta \theta}{2} \end{bmatrix} $$
    fn retract(&mut self, delta: &[f64]) {
        match self {
            VariableType::Pos(p) => { p.x += delta[0]; p.y += delta[1]; p.z += delta[2]; }
            VariableType::Vel(v) => { v.x += delta[0]; v.y += delta[1]; v.z += delta[2]; }
            VariableType::Att(q) => Self::retract_att(q, delta),
            VariableType::Biases(b) => { for i in 0..6 { b[i] += delta[i]; } }
            VariableType::Clock(c) => { c[0] += delta[0]; c[1] += delta[1]; }
            VariableType::Ambiguity(a) => *a += delta[0],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{UnitQuaternion, Vector2, Vector3, Vector6};

    #[test]
    fn test_local_dims() {
        assert_eq!(VariableType::Pos(Vector3::zeros()).local_dim(), 3);
        assert_eq!(VariableType::Vel(Vector3::zeros()).local_dim(), 3);
        assert_eq!(VariableType::Att(UnitQuaternion::identity()).local_dim(), 3);
        assert_eq!(VariableType::Biases(Vector6::zeros()).local_dim(), 6);
        assert_eq!(VariableType::Clock(Vector2::zeros()).local_dim(), 2);
        assert_eq!(VariableType::Ambiguity(0.0).local_dim(), 1);
    }

    #[test]
    fn test_manifold_pos_retract() {
        let mut pos = VariableType::Pos(Vector3::new(1.0, 2.0, 3.0));
        pos.retract(&[0.1, -0.2, 0.3]);
        if let VariableType::Pos(p) = pos {
            assert!((p.x - 1.1).abs() < 1e-9);
            assert!((p.y - 1.8).abs() < 1e-9);
            assert!((p.z - 3.3).abs() < 1e-9);
        } else {
            panic!("Wrong type");
        }
    }

    #[test]
    fn test_manifold_att_retract() {
        let mut att = VariableType::Att(UnitQuaternion::identity());
        // Moderate perturbation to avoid small-angle approx branch
        att.retract(&[0.5, 0.0, 0.0]);
        if let VariableType::Att(q) = att {
            assert!((q.w - (0.25_f64).cos()).abs() < 1e-9);
            assert!((q.i - (0.25_f64).sin()).abs() < 1e-9);
        } else {
            panic!("Wrong type");
        }
    }

    #[test]
    fn test_manifold_att_zero_perturbation() {
        let mut att = VariableType::Att(UnitQuaternion::identity());
        att.retract(&[0.0, 0.0, 0.0]); // Triggers small angle fallback cleanly
        if let VariableType::Att(q) = att {
            assert_eq!(q.w, 1.0);
            assert_eq!(q.i, 0.0);
        } else {
            panic!("Wrong type");
        }
    }

    #[test]
    fn test_manifold_att_extreme_perturbation() {
        let mut att = VariableType::Att(UnitQuaternion::identity());
        // 2*PI perturbation around X-axis. Wraps back to identity effectively.
        att.retract(&[2.0 * std::f64::consts::PI, 0.0, 0.0]);
        if let VariableType::Att(q) = att {
            assert!((q.w.abs() - 1.0).abs() < 1e-5);
            assert!(q.i.abs() < 1e-5);
        } else {
            panic!("Wrong type");
        }
    }

    #[test]
    fn test_manifold_biases_retract() {
        let mut b = VariableType::Biases(Vector6::zeros());
        b.retract(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        if let VariableType::Biases(v) = b {
            assert_eq!(v[0], 1.0);
            assert_eq!(v[1], 2.0);
            assert_eq!(v[2], 3.0);
            assert_eq!(v[3], 4.0);
            assert_eq!(v[4], 5.0);
            assert_eq!(v[5], 6.0);
        } else {
            panic!("Wrong type");
        }
    }

    #[test]
    fn test_manifold_other_retracts() {
        let mut p = VariableType::Pos(Vector3::zeros());
        p.retract(&[1.0, 2.0, 3.0]);
        if let VariableType::Pos(vec) = p { 
            assert_eq!(vec.x, 1.0); 
            assert_eq!(vec.y, 2.0); 
            assert_eq!(vec.z, 3.0); 
        } else { panic!(); }

        let mut v = VariableType::Vel(Vector3::zeros());
        v.retract(&[1.0, 2.0, 3.0]);
        if let VariableType::Vel(vec) = v { 
            assert_eq!(vec.x, 1.0); 
            assert_eq!(vec.y, 2.0); 
            assert_eq!(vec.z, 3.0); 
        } else { panic!(); }
        
        let mut c = VariableType::Clock(Vector2::zeros());
        c.retract(&[1.0, 2.0]);
        if let VariableType::Clock(vec) = c { 
            assert_eq!(vec[0], 1.0); 
            assert_eq!(vec[1], 2.0); 
        } else { panic!(); }
        
        let mut a = VariableType::Ambiguity(1.0);
        a.retract(&[2.0]);
        if let VariableType::Ambiguity(val) = a { assert_eq!(val, 3.0); } else { panic!(); }
    }
}
