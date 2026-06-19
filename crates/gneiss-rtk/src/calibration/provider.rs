use nalgebra::Vector3;

/// A generic interface for providing hardware-specific sensor calibrations.
/// This allows the engine to decouple from specific hardware while still supporting
/// temperature-calibrated IMU models or rigorous Antenna Phase Center (APC) corrections.
pub trait CalibrationProvider {
    /// Applies IMU calibration (misalignments, scale factors) to raw accelerometer readings.
    fn calibrate_accelerometer(
        &self,
        raw_accel: Vector3<f64>,
        temperature_c: Option<f64>,
    ) -> Vector3<f64>;

    /// Applies IMU calibration (misalignments, scale factors) to raw gyroscope readings.
    fn calibrate_gyroscope(
        &self,
        raw_gyro: Vector3<f64>,
        temperature_c: Option<f64>,
    ) -> Vector3<f64>;

    /// Provides the Antenna Phase Center (APC) offset for a given frequency band.
    /// Usually varies by azimuth and elevation for high-end antennas.
    fn antenna_phase_center_offset(&self, az_rad: f64, el_rad: f64, freq_band: u8) -> Vector3<f64>;
}

/// A default pass-through calibration provider that applies no custom corrections.
/// Useful for when the engine runs in a pure software-defined mode or uses pre-calibrated data.
pub struct DefaultCalibrationProvider;

impl CalibrationProvider for DefaultCalibrationProvider {
    fn calibrate_accelerometer(
        &self,
        raw_accel: Vector3<f64>,
        _temperature_c: Option<f64>,
    ) -> Vector3<f64> {
        raw_accel
    }

    fn calibrate_gyroscope(
        &self,
        raw_gyro: Vector3<f64>,
        _temperature_c: Option<f64>,
    ) -> Vector3<f64> {
        raw_gyro
    }

    fn antenna_phase_center_offset(
        &self,
        _az_rad: f64,
        _el_rad: f64,
        _freq_band: u8,
    ) -> Vector3<f64> {
        Vector3::zeros()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_calibration_provider() {
        let provider = DefaultCalibrationProvider;
        let accel = Vector3::new(1.0, 2.0, 3.0);
        let gyro = Vector3::new(0.1, 0.2, 0.3);

        assert_eq!(provider.calibrate_accelerometer(accel, None), accel);
        assert_eq!(provider.calibrate_gyroscope(gyro, Some(25.0)), gyro);
        assert_eq!(
            provider.antenna_phase_center_offset(0.0, 0.0, 1),
            Vector3::zeros()
        );
    }
}
