use nalgebra::Vector3;
use crate::time::GpsTime;

/// Ionospheric delay model parameters (Klobuchar).
#[derive(Debug, Clone, Copy)]
pub struct KlobucharParams {
    pub alpha: [f64; 4],
    pub beta: [f64; 4],
}

impl Default for KlobucharParams {
    fn default() -> Self {
        Self {
            alpha: [0.1118E-07, -0.7451E-08, -0.5960E-07,  0.1192E-06],
            beta:  [0.1167E+06, -0.2294E+06, -0.1311E+06,  0.1049E+07],
        }
    }
}

/// Tropospheric delay model parameters.
#[derive(Debug, Clone, Copy)]
pub struct TropoParams {
    pub temp_k: f64,
    pub press_hpa: f64,
    pub hum_rel: f64,
}

impl Default for TropoParams {
    fn default() -> Self {
        Self {
            temp_k: 288.15, // 15 C
            press_hpa: 1013.25,
            hum_rel: 0.5,
        }
    }
}

/// Models atmospheric delays for a specific satellite-receiver geometry.
pub struct AtmosphereModel;

impl AtmosphereModel {
    /// Computes Ionospheric delay in meters using the Klobuchar model.
    /// `pos`: Receiver ECEF position.
    /// `az`: Satellite azimuth in radians.
    /// `el`: Satellite elevation in radians.
    /// `time`: GPS time of observation.
    pub fn iono_klobuchar(params: &KlobucharParams, pos_llh: Vector3<f64>, _az: f64, el: f64, time: GpsTime) -> f64 {
        // Implementation based on IS-GPS-200
        let f = 1.0 + 16.0 * libm::pow(0.53 - el / core::f64::consts::PI, 3.0);
        let phi_m = pos_llh.x / core::f64::consts::PI + 0.064 * libm::cos(pos_llh.y - 1.617);
        
        let mut t = 43200.0 * phi_m + time.tow;
        t %= 86400.0;
        if t < 0.0 { t += 86400.0; }

        let mut a = params.alpha[0] + params.alpha[1] * phi_m + params.alpha[2] * phi_m * phi_m + params.alpha[3] * phi_m * phi_m * phi_m;
        if a < 0.0 { a = 0.0; }

        let mut p = params.beta[0] + params.beta[1] * phi_m + params.beta[2] * phi_m * phi_m + params.beta[3] * phi_m * phi_m * phi_m;
        if p < 72000.0 { p = 72000.0; }

        let x = 2.0 * core::f64::consts::PI * (t - 50400.0) / p;
        
        let delay = if libm::fabs(x) < 1.57 {
            5e-9 + a * (1.0 - x * x / 2.0 + x * x * x * x / 24.0)
        } else {
            5e-9
        };

        delay * f * 299792458.0 // Return in meters
    }

    /// Computes Tropospheric delay in meters using the Saastamoinen model.
    pub fn tropo_saastamoinen(params: &TropoParams, el: f64, height: f64) -> f64 {
        let z = core::f64::consts::FRAC_PI_2 - el;
        let p = params.press_hpa * libm::pow(1.0 - 0.000022557 * height, 5.2568);
        let t = params.temp_k - 0.0065 * height;
        let e = 6.108 * libm::exp((17.15 * t - 4684.0) / (t - 38.45)) * params.hum_rel;

        
        0.002277 / libm::cos(z) * (p + (1255.0 / t + 0.05) * e - libm::tan(z) * libm::tan(z))
    }

    /// Computes Tropospheric delay in meters using the EXACT RTKLIB Saastamoinen model.
    pub fn tropo_rtklib_saastamoinen(params: &TropoParams, pos_llh: Vector3<f64>, el: f64) -> f64 {
        if pos_llh.z < -100.0 || pos_llh.z > 10000.0 || el <= 0.0 {
            return 0.0;
        }
        let z = core::f64::consts::FRAC_PI_2 - el;
        let height = pos_llh.z;
        let p = params.press_hpa * libm::pow(1.0 - 0.000022557 * height, 5.2568);
        let t = params.temp_k - 0.0065 * height;
        let e = 6.108 * libm::exp((17.15 * t - 4684.0) / (t - 38.45)) * params.hum_rel;

        let trph = 0.0022768 * p / (1.0 - 0.00266 * libm::cos(2.0 * pos_llh.x) - 0.00028 * height / 1000.0);
        let trpw = 0.002277 * (1255.0 / t + 0.05) * e;
        
        (trph + trpw) / libm::cos(z)
    }

    fn nmf_interpc(coef: &[f64; 5], lat: f64) -> f64 {
        let i = (lat / 15.0) as usize;
        if i < 1 {
            return coef[0];
        } else if i > 4 {
            return coef[4];
        }
        let lat_f = lat / 15.0;
        let i_f = i as f64;
        coef[i - 1] * (1.0 - lat_f + i_f) + coef[i] * (lat_f - i_f)
    }

    fn nmf_mapf(el: f64, a: f64, b: f64, c: f64) -> f64 {
        let sinel = libm::sin(el);
        (1.0 + a / (1.0 + b / (1.0 + c))) / (sinel + (a / (sinel + b / (sinel + c))))
    }

    /// Computes Tropospheric delay in meters using the Saastamoinen zenith delay mapped with Niell Mapping Function (NMF).
    /// `pos_llh`: Receiver position (Lat, Lon, Height) in radians and meters
    /// `el`: Elevation angle in radians
    /// `time`: GPS time of observation
    pub fn tropo_nmf(params: &TropoParams, pos_llh: Vector3<f64>, el: f64, time: GpsTime) -> f64 {
        if el <= 0.0 {
            return 0.0;
        }

        let hgt = pos_llh.z;
        let mut lat = pos_llh.x * 180.0 / core::f64::consts::PI;

        let coef = [
            [ 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3 ],
            [ 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3 ],
            [ 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3 ],
            
            [ 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5 ],
            [ 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5 ],
            [ 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5 ],
            
            [ 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4 ],
            [ 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3 ],
            [ 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2 ]
        ];
        let aht = [ 2.53E-5, 5.49E-3, 1.14E-3 ];

        let fractional_year = time.to_fractional_year();
        let doy_fraction = fractional_year - libm::floor(fractional_year);
        let y = doy_fraction - (28.0 / 365.25) + if lat < 0.0 { 0.5 } else { 0.0 };
        
        let cosy = libm::cos(2.0 * core::f64::consts::PI * y);
        lat = libm::fabs(lat);

        let mut ah = [0.0; 3];
        let mut aw = [0.0; 3];

        for i in 0..3 {
            ah[i] = Self::nmf_interpc(&coef[i], lat) - Self::nmf_interpc(&coef[i+3], lat) * cosy;
            aw[i] = Self::nmf_interpc(&coef[i+6], lat);
        }

        let dm = (1.0 / libm::sin(el) - Self::nmf_mapf(el, aht[0], aht[1], aht[2])) * hgt / 1000.0;
        let m_h = Self::nmf_mapf(el, ah[0], ah[1], ah[2]) + dm;
        let m_w = Self::nmf_mapf(el, aw[0], aw[1], aw[2]);

        // Zenith dry and wet delays (simplified Saastamoinen)
        let z_dry = 0.0022768 * params.press_hpa / (1.0 - 0.00266 * libm::cos(2.0 * pos_llh.x) - 0.00028 * hgt / 1000.0);
        
        let e = 6.108 * libm::exp((17.15 * params.temp_k - 4684.0) / (params.temp_k - 38.45)) * params.hum_rel;
        let z_wet = 0.002277 * (1255.0 / params.temp_k + 0.05) * e;

        z_dry * m_h + z_wet * m_w
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    

    #[test]
    fn test_tropo_delay() {
        let params = TropoParams::default();
        let delay = AtmosphereModel::tropo_saastamoinen(&params, 0.5, 100.0);
        assert!(delay > 2.0 && delay < 10.0);
    }
}
