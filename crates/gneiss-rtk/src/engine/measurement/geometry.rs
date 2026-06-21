use crate::engine::measurement::types::MeasurementEnvironment;
use crate::filter::RtkState;
use gneiss_core::ephemeris::Ephemeris;
use nalgebra::Vector3;

pub struct EkfGeometryContext {
    pub pos_apc: Vector3<f64>,
    pub base_coord_vec: Vector3<f64>,
    pub r_b_e: nalgebra::Matrix3<f64>,
    pub lever_arm: Vector3<f64>,
    pub state_size: usize,
}

impl EkfGeometryContext {
    pub fn new(state: &RtkState, env: &MeasurementEnvironment) -> Self {
        let r_b_e = state.attitude.to_rotation_matrix();
        let mut pos_apc = state.position.vector + r_b_e * env.lever_arm;
        let mut base_coord_vec = env.base_coord.vector;

        let set_rov = gneiss_core::tides::solid_earth_tides_ecef(state.time, pos_apc);
        let set_bas = gneiss_core::tides::solid_earth_tides_ecef(state.time, base_coord_vec);

        pos_apc += set_rov;
        base_coord_vec += set_bas;

        let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();

        Self {
            pos_apc,
            base_coord_vec,
            r_b_e: r_b_e.into_inner(),
            lever_arm: env.lever_arm,
            state_size,
        }
    }

    pub fn compute_attitude_jacobian(&self, h_r: &Vector3<f64>) -> Vector3<f64> {
        crate::engine::measurement_math::range_attitude_jacobian(&(self.r_b_e * self.lever_arm), h_r)
    }
}

pub(crate) fn find_ephemeris(
    ephemerides: &[Ephemeris],
    sat: gneiss_core::sat::SatelliteId,
    time_tow: f64,
) -> Option<&Ephemeris> {
    ephemerides
        .iter()
        .filter(|e| e.sat() == sat)
        .min_by(|a, b| {
            let da = (a.toe().tow - time_tow).abs();
            let db = (b.toe().tow - time_tow).abs();
            da.partial_cmp(&db).unwrap()
        })
}
