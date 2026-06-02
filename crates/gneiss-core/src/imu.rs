use nalgebra::Vector3;

/// A generic 6-axis Inertial Measurement Unit (IMU) reading.
/// Can be natively ingested from a tightly-coupled receiver (like ZED-F9R)
/// or an external sensor on the same bus.
#[derive(Debug, Clone, PartialEq)]
pub struct ImuMeasurement {
    /// The timestamp of the measurement (can be system time or GNSS time depending on sync).
    pub time_tag: u32, 
    /// 3D Acceleration vector (X, Y, Z) in the vehicle/sensor body frame, typically in m/s^2.
    pub accel: Vector3<f64>,
    /// 3D Angular Velocity (Gyroscope) vector (X, Y, Z) in the vehicle/sensor body frame, typically in rad/s.
    pub gyro: Vector3<f64>,
    /// Optional temperature of the sensor in Celsius, useful for thermal drift compensation.
    pub temperature: Option<f64>,
}

impl ImuMeasurement {
    pub fn new(time_tag: u32, accel: Vector3<f64>, gyro: Vector3<f64>) -> Self {
        Self {
            time_tag,
            accel,
            gyro,
            temperature: None,
        }
    }
}
