//! Opaque-box typestates, scale markers, and relational geometry contracts.

use core::marker::PhantomData;
use core::ops::{Add, Sub};
use nalgebra::{Matrix3, UnitQuaternion, Vector3};

pub const WEEK_NANOS: u64 = 604_800_000_000_000;
pub const GPS_BDT_OFFSET_SECONDS: i64 = 14;

pub trait CoordinateFrame: 'static + Send + Sync + Copy + PartialEq + Eq {
    const NAME: &'static str;
}

pub trait ReferenceFrame: 'static + Send + Sync + Copy + PartialEq + Eq {
    const NAME: &'static str;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Itrf2014;
impl ReferenceFrame for Itrf2014 { const NAME: &'static str = "ITRF2014"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Itrf2020;
impl ReferenceFrame for Itrf2020 { const NAME: &'static str = "ITRF2020"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wgs84;
impl ReferenceFrame for Wgs84 { const NAME: &'static str = "WGS84"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Nad83;
impl ReferenceFrame for Nad83 { const NAME: &'static str = "NAD83"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Jgd2011;
impl ReferenceFrame for Jgd2011 { const NAME: &'static str = "JGD2011"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Pz90;
impl ReferenceFrame for Pz90 { const NAME: &'static str = "PZ-90.11"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ecef<R: ReferenceFrame>(PhantomData<R>);
impl<R: ReferenceFrame> CoordinateFrame for Ecef<R> { const NAME: &'static str = "ECEF"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ned;
impl CoordinateFrame for Ned { const NAME: &'static str = "NED"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Enu;
impl CoordinateFrame for Enu { const NAME: &'static str = "ENU"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BodyFrd;
impl CoordinateFrame for BodyFrd { const NAME: &'static str = "Body-FRD"; }

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Point3<F: CoordinateFrame> {
    coords: Vector3<f64>,
    _frame: PhantomData<F>,
}

impl<F: CoordinateFrame> Point3<F> {
    pub const fn from_coords(coords: Vector3<f64>) -> Self {
        Self { coords, _frame: PhantomData }
    }
    pub const fn coords(&self) -> &Vector3<f64> { &self.coords }
    pub const fn into_coords(self) -> Vector3<f64> { self.coords }
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SpatialVector<F: CoordinateFrame> {
    vector: Vector3<f64>,
    _frame: PhantomData<F>,
}

impl<F: CoordinateFrame> SpatialVector<F> {
    pub const fn from_vector(vector: Vector3<f64>) -> Self {
        Self { vector, _frame: PhantomData }
    }
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self::from_vector(Vector3::new(x, y, z))
    }
    pub fn zero() -> Self { Self::from_vector(Vector3::zeros()) }
    pub const fn vector(&self) -> &Vector3<f64> { &self.vector }
    pub const fn into_vector(self) -> Vector3<f64> { self.vector }
    pub fn norm(&self) -> f64 { self.vector.norm() }
}

impl<F: CoordinateFrame> Add<SpatialVector<F>> for Point3<F> {
    type Output = Point3<F>;
    fn add(self, rhs: SpatialVector<F>) -> Self::Output {
        Point3::from_coords(self.coords + rhs.vector)
    }
}

impl<F: CoordinateFrame> Sub for Point3<F> {
    type Output = SpatialVector<F>;
    fn sub(self, rhs: Point3<F>) -> Self::Output {
        SpatialVector::from_vector(self.coords - rhs.coords)
    }
}

impl<F: CoordinateFrame> Add for SpatialVector<F> {
    type Output = SpatialVector<F>;
    fn add(self, rhs: SpatialVector<F>) -> Self::Output {
        SpatialVector::from_vector(self.vector + rhs.vector)
    }
}

pub type EcefPos<R> = Point3<Ecef<R>>;
pub type EcefVector<R> = SpatialVector<Ecef<R>>;

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SpatialVelocity<F: CoordinateFrame> {
    vector: Vector3<f64>,
    _frame: PhantomData<F>,
}

impl<F: CoordinateFrame> SpatialVelocity<F> {
    pub const fn from_vector(vector: Vector3<f64>) -> Self {
        Self { vector, _frame: PhantomData }
    }
    pub const fn vector(&self) -> &Vector3<f64> { &self.vector }
    pub fn norm(&self) -> f64 { self.vector.norm() }
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SpatialCovariance<F: CoordinateFrame> {
    matrix: Matrix3<f64>,
    _frame: PhantomData<F>,
}

impl<F: CoordinateFrame> SpatialCovariance<F> {
    pub const fn from_matrix(matrix: Matrix3<f64>) -> Self {
        Self { matrix, _frame: PhantomData }
    }
    pub const fn matrix(&self) -> &Matrix3<f64> { &self.matrix }
}

pub type NedCovariance = SpatialCovariance<Ned>;
pub type EnuCovariance = SpatialCovariance<Enu>;

impl NedCovariance {
    pub fn std_north(&self) -> f64 { self.matrix[(0, 0)].max(0.0).sqrt() }
    pub fn std_east(&self) -> f64 { self.matrix[(1, 1)].max(0.0).sqrt() }
    pub fn std_down(&self) -> f64 { self.matrix[(2, 2)].max(0.0).sqrt() }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AntennaLeverArm(pub SpatialVector<BodyFrd>);

impl AntennaLeverArm {
    pub fn new(forward_m: f64, right_m: f64, down_m: f64) -> Self {
        Self(SpatialVector::new(forward_m, right_m, down_m))
    }
    pub fn zero() -> Self { Self(SpatialVector::zero()) }
    pub const fn as_body_vector(&self) -> &SpatialVector<BodyFrd> { &self.0 }
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Attitude<From: CoordinateFrame, To: CoordinateFrame> {
    q: UnitQuaternion<f64>,
    _frames: PhantomData<(From, To)>,
}

impl<From: CoordinateFrame, To: CoordinateFrame> Attitude<From, To> {
    pub const fn from_quaternion(q: UnitQuaternion<f64>) -> Self {
        Self { q, _frames: PhantomData }
    }
    pub const fn quaternion(&self) -> &UnitQuaternion<f64> { &self.q }
    pub fn rotate_vector(&self, v: &SpatialVector<From>) -> SpatialVector<To> {
        SpatialVector::from_vector(self.q * v.vector())
    }
    pub fn rotate_velocity(&self, v: &SpatialVelocity<From>) -> SpatialVelocity<To> {
        SpatialVelocity::from_vector(self.q * v.vector())
    }
    pub fn rotate_cov(&self, cov: &SpatialCovariance<From>) -> SpatialCovariance<To> {
        let r = self.q.to_rotation_matrix().into_inner();
        SpatialCovariance::from_matrix(r * cov.matrix() * r.transpose())
    }
    pub fn inverse(&self) -> Attitude<To, From> {
        Attitude::from_quaternion(self.q.inverse())
    }
}

pub trait TimeScale: Copy + Clone + PartialEq + Eq + PartialOrd + Ord + core::fmt::Debug + 'static {
    const NAME: &'static str;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct GpsScale;
impl TimeScale for GpsScale { const NAME: &'static str = "GPST"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct BdtScale;
impl TimeScale for BdtScale { const NAME: &'static str = "BDT"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct GstScale;
impl TimeScale for GstScale { const NAME: &'static str = "GST"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct GlonassScale;
impl TimeScale for GlonassScale { const NAME: &'static str = "GLONASST"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct UtcScale;
impl TimeScale for UtcScale { const NAME: &'static str = "UTC"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct TimeDelta {
    nanos: i64,
}

impl TimeDelta {
    pub const fn from_nanos(nanos: i64) -> Self { Self { nanos } }
    pub const fn from_millis(millis: i64) -> Self { Self { nanos: millis * 1_000_000 } }
    pub fn from_seconds(secs: f64) -> Self { Self { nanos: (secs * 1e9).round() as i64 } }
    pub const fn as_nanos(&self) -> i64 { self.nanos }
    pub fn as_seconds(&self) -> f64 { self.nanos as f64 * 1e-9 }
    pub fn as_millis(&self) -> i64 { self.nanos / 1_000_000 }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EpochKey<S: TimeScale> {
    continuous_ms: u64,
    _scale: PhantomData<S>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Epoch<S: TimeScale> {
    week: u32,
    tow_nanos: u64,
    _scale: PhantomData<S>,
}

impl<S: TimeScale> Epoch<S> {
    pub fn from_week_nanos(week: u32, tow_nanos: u64) -> Self {
        let (extra_weeks, normalized_nanos) = (tow_nanos / WEEK_NANOS, tow_nanos % WEEK_NANOS);
        Self { week: week + extra_weeks as u32, tow_nanos: normalized_nanos, _scale: PhantomData }
    }
    pub fn from_week_tow(week: u32, tow_sec: f64) -> Self {
        let total_nanos = (tow_sec * 1e9).round() as u64;
        Self::from_week_nanos(week, total_nanos)
    }
    pub const fn week(&self) -> u32 { self.week }
    pub const fn tow_nanos(&self) -> u64 { self.tow_nanos }
    pub fn tow_seconds(&self) -> f64 { self.tow_nanos as f64 * 1e-9 }
    pub fn continuous_ms(&self) -> u64 {
        (self.week as u64) * 604_800_000 + (self.tow_nanos / 1_000_000)
    }
    pub fn to_key(&self) -> EpochKey<S> {
        EpochKey { continuous_ms: self.continuous_ms(), _scale: PhantomData }
    }
    pub fn is_within(&self, other: Self, tolerance: TimeDelta) -> bool {
        let diff = if *self >= other { *self - other } else { other - *self };
        diff.as_nanos().abs() <= tolerance.as_nanos().abs()
    }
}

impl<S: TimeScale> PartialOrd for Epoch<S> {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<S: TimeScale> Ord for Epoch<S> {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        match self.week.cmp(&other.week) {
            core::cmp::Ordering::Equal => self.tow_nanos.cmp(&other.tow_nanos),
            ord => ord,
        }
    }
}

impl<S: TimeScale> Sub for Epoch<S> {
    type Output = TimeDelta;
    fn sub(self, rhs: Epoch<S>) -> Self::Output {
        let d_week = self.week as i64 - rhs.week as i64;
        let d_nanos = self.tow_nanos as i64 - rhs.tow_nanos as i64;
        TimeDelta::from_nanos(d_week * (WEEK_NANOS as i64) + d_nanos)
    }
}

impl<S: TimeScale> Add<TimeDelta> for Epoch<S> {
    type Output = Epoch<S>;
    fn add(self, rhs: TimeDelta) -> Self::Output {
        let total_nanos = self.tow_nanos as i64 + rhs.nanos;
        if total_nanos >= 0 {
            let add_weeks = total_nanos / (WEEK_NANOS as i64);
            let rem_nanos = total_nanos % (WEEK_NANOS as i64);
            Epoch::from_week_nanos(self.week + add_weeks as u32, rem_nanos as u64)
        } else {
            let neg_weeks = (-total_nanos + (WEEK_NANOS as i64) - 1) / (WEEK_NANOS as i64);
            let rem_nanos = total_nanos + neg_weeks * (WEEK_NANOS as i64);
            Epoch::from_week_nanos(self.week - neg_weeks as u32, rem_nanos as u64)
        }
    }
}

impl Epoch<BdtScale> {
    pub fn to_gpst(&self) -> Epoch<GpsScale> {
        let dt = TimeDelta::from_seconds(GPS_BDT_OFFSET_SECONDS as f64);
        Epoch::<GpsScale>::from_week_nanos(self.week, self.tow_nanos) + dt
    }
}

impl Epoch<GpsScale> {
    pub fn to_bdt(&self) -> Epoch<BdtScale> {
        let dt = TimeDelta::from_seconds(-(GPS_BDT_OFFSET_SECONDS as f64));
        Epoch::<BdtScale>::from_week_nanos(self.week, self.tow_nanos) + dt
    }
    pub fn to_utc(&self, leap_seconds: i32) -> Epoch<UtcScale> {
        let dt = TimeDelta::from_seconds(-(leap_seconds as f64));
        Epoch::<UtcScale>::from_week_nanos(self.week, self.tow_nanos) + dt
    }
}

impl Epoch<UtcScale> {
    pub fn to_gpst(&self, leap_seconds: i32) -> Epoch<GpsScale> {
        let dt = TimeDelta::from_seconds(leap_seconds as f64);
        Epoch::<GpsScale>::from_week_nanos(self.week, self.tow_nanos) + dt
    }
}
