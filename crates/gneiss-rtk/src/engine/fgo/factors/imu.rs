use nalgebra::{DMatrix, Matrix3, Vector3};

pub const STATE_SIZE: usize = 15;
pub const POS_IDX: usize = 0;
pub const VEL_IDX: usize = 3;
pub const ROT_IDX: usize = 6;
pub const BA_IDX: usize = 9;
pub const BG_IDX: usize = 12;

pub struct ImuFactor {}

pub struct JacobianInputs {
    pub r_i: Matrix3<f64>,
    pub r_j: Matrix3<f64>,
    pub dt: f64,
    pub omega_ie: Vector3<f64>,
    pub dp_dba: Matrix3<f64>,
    pub dp_dbg: Matrix3<f64>,
    pub dr_dbg: Matrix3<f64>,
    pub jr_inv_er: Matrix3<f64>,
    pub jl_inv_er: Matrix3<f64>,
    pub exp_omega_dt: Matrix3<f64>,
}

pub fn skew_symmetric(v: &Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(0.0, -v.z, v.y, v.z, 0.0, -v.x, -v.y, v.x, 0.0)
}

impl ImuFactor {
    pub fn jacobians(inputs: &JacobianInputs) -> (DMatrix<f64>, DMatrix<f64>) {
        let mut h_i = DMatrix::zeros(STATE_SIZE, STATE_SIZE);
        let mut h_j = DMatrix::zeros(STATE_SIZE, STATE_SIZE);

        let i_mat = Matrix3::identity();
        let omega_skew = skew_symmetric(&inputs.omega_ie);
        let r_i_t = inputs.r_i.transpose();

        // H_i blocks
        h_i.fixed_view_mut::<3, 3>(POS_IDX, POS_IDX)
            .copy_from(&-r_i_t);
        h_i.fixed_view_mut::<3, 3>(POS_IDX, VEL_IDX)
            .copy_from(&(r_i_t * (-i_mat * inputs.dt + omega_skew * inputs.dt.powi(2))));
        h_i.fixed_view_mut::<3, 3>(POS_IDX, BA_IDX)
            .copy_from(&-inputs.dp_dba);
        h_i.fixed_view_mut::<3, 3>(POS_IDX, BG_IDX)
            .copy_from(&-inputs.dp_dbg);

        // dv_dpi is 0
        h_i.fixed_view_mut::<3, 3>(VEL_IDX, VEL_IDX)
            .copy_from(&(r_i_t * (-i_mat + omega_skew * 2.0 * inputs.dt)));

        h_i.fixed_view_mut::<3, 3>(ROT_IDX, ROT_IDX).copy_from(
            &(-inputs.jr_inv_er * inputs.r_j.transpose() * inputs.exp_omega_dt * inputs.r_i),
        );
        h_i.fixed_view_mut::<3, 3>(ROT_IDX, BG_IDX)
            .copy_from(&(-inputs.jl_inv_er * inputs.dr_dbg));

        // H_j blocks
        h_j.fixed_view_mut::<3, 3>(POS_IDX, POS_IDX)
            .copy_from(&r_i_t);
        h_j.fixed_view_mut::<3, 3>(VEL_IDX, VEL_IDX)
            .copy_from(&r_i_t);
        h_j.fixed_view_mut::<3, 3>(ROT_IDX, ROT_IDX)
            .copy_from(&inputs.jr_inv_er);
        h_j.fixed_view_mut::<3, 3>(BA_IDX, BA_IDX).copy_from(&i_mat);
        h_j.fixed_view_mut::<3, 3>(BG_IDX, BG_IDX).copy_from(&i_mat);

        (h_i, h_j)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_mock_inputs() -> JacobianInputs {
        JacobianInputs {
            r_i: Matrix3::new(1.0, 0.1, 0.2, -0.1, 1.0, 0.3, -0.2, -0.3, 1.0),
            r_j: Matrix3::new(1.0, -0.1, 0.2, 0.1, 1.0, 0.3, -0.2, -0.3, 1.0),
            dt: 0.1,
            omega_ie: Vector3::new(0.01, 0.02, 0.03),
            dp_dba: Matrix3::identity() * 0.1,
            dp_dbg: Matrix3::identity() * 0.2,
            dr_dbg: Matrix3::identity() * 0.3,
            jr_inv_er: Matrix3::identity() * 0.9,
            jl_inv_er: Matrix3::identity() * 1.1,
            exp_omega_dt: Matrix3::identity() * 0.99,
        }
    }

    #[test]
    fn test_imu_factor_jacobians_hi_pos() {
        let inputs = create_mock_inputs();
        let (h_i, _) = ImuFactor::jacobians(&inputs);
        let omega_skew = skew_symmetric(&inputs.omega_ie);
        let r_i_t = inputs.r_i.transpose();
        let i_mat = Matrix3::identity();

        let dp_dpi = h_i.fixed_view::<3, 3>(POS_IDX, POS_IDX);
        assert!((dp_dpi - (-r_i_t)).norm() < 1e-10);

        let dp_dvi = h_i.fixed_view::<3, 3>(POS_IDX, VEL_IDX);
        let expected_dp_dvi = r_i_t * (-i_mat * inputs.dt + omega_skew * inputs.dt.powi(2));
        assert!((dp_dvi - expected_dp_dvi).norm() < 1e-10);

        let dp_dbai = h_i.fixed_view::<3, 3>(POS_IDX, BA_IDX);
        assert!((dp_dbai - (-inputs.dp_dba)).norm() < 1e-10);

        let dp_dbgi = h_i.fixed_view::<3, 3>(POS_IDX, BG_IDX);
        assert!((dp_dbgi - (-inputs.dp_dbg)).norm() < 1e-10);
    }

    #[test]
    fn test_imu_factor_jacobians_hi_vel_rot() {
        let inputs = create_mock_inputs();
        let (h_i, _) = ImuFactor::jacobians(&inputs);
        let omega_skew = skew_symmetric(&inputs.omega_ie);
        let r_i_t = inputs.r_i.transpose();
        let i_mat = Matrix3::identity();

        let dv_dpi = h_i.fixed_view::<3, 3>(VEL_IDX, POS_IDX);
        assert!(dv_dpi.norm() < 1e-10);

        let dv_dvi = h_i.fixed_view::<3, 3>(VEL_IDX, VEL_IDX);
        let expected_dv_dvi = r_i_t * (-i_mat + omega_skew * 2.0 * inputs.dt);
        assert!((dv_dvi - expected_dv_dvi).norm() < 1e-10);

        let dr_dthetai = h_i.fixed_view::<3, 3>(ROT_IDX, ROT_IDX);
        let expected_dr_dthetai =
            -inputs.jr_inv_er * inputs.r_j.transpose() * inputs.exp_omega_dt * inputs.r_i;
        assert!((dr_dthetai - expected_dr_dthetai).norm() < 1e-10);

        let dr_dbgi = h_i.fixed_view::<3, 3>(ROT_IDX, BG_IDX);
        assert!((dr_dbgi - (-inputs.jl_inv_er * inputs.dr_dbg)).norm() < 1e-10);
    }

    #[test]
    fn test_imu_factor_jacobians_hj() {
        let inputs = create_mock_inputs();
        let (_, h_j) = ImuFactor::jacobians(&inputs);
        let r_i_t = inputs.r_i.transpose();
        let i_mat = Matrix3::identity();

        let dp_dpj = h_j.fixed_view::<3, 3>(POS_IDX, POS_IDX);
        assert!((dp_dpj - r_i_t).norm() < 1e-10);

        let dv_dvj = h_j.fixed_view::<3, 3>(VEL_IDX, VEL_IDX);
        assert!((dv_dvj - r_i_t).norm() < 1e-10);

        let dr_dthetaj = h_j.fixed_view::<3, 3>(ROT_IDX, ROT_IDX);
        assert!((dr_dthetaj - inputs.jr_inv_er).norm() < 1e-10);

        let dba_dbaj = h_j.fixed_view::<3, 3>(BA_IDX, BA_IDX);
        assert!((dba_dbaj - i_mat).norm() < 1e-10);

        let dbg_dbgj = h_j.fixed_view::<3, 3>(BG_IDX, BG_IDX);
        assert!((dbg_dbgj - i_mat).norm() < 1e-10);
    }
}
