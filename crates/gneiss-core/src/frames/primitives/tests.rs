//! Golden-vector tests for the frame-tagged spatial primitives.
//!
//! Every expected number below is a hand evaluation of ordinary 3-vector
//! arithmetic on the triple (1, 2, 3) and (4, 5, 6):
//!   dot (1,2,3).(4,5,6) = 4 + 10 + 18            = 32
//!   (1,2,3) x (4,5,6) = (2*6 - 3*5, 3*4 - 1*6, 1*5 - 2*4) = (-3, 6, -3)
//!   |(3, 4, 0)| = 5
//! Distinct components in every input mean that an axis permutation, a sign
//! flip, or a swapped subtraction fails immediately.

use super::*;
use crate::frames::realizations::{Itrf2014, Nad83_2011};
use alloc::format;
use nalgebra::Vector3;

fn v(x: f64, y: f64, z: f64) -> SpatialVector<Ned> {
    SpatialVector::new(x, y, z)
}

fn p(x: f64, y: f64, z: f64) -> Point3<Ned> {
    Point3::new(Vector3::new(x, y, z))
}

#[test]
fn point_accessors_do_not_permute_axes() {
    let q = Point3::<Ecef<Itrf2014>>::from_coords(1.0, -2.0, 3.5);
    assert!((q.x() - 1.0).abs() < f64::EPSILON);
    assert!((q.y() + 2.0).abs() < f64::EPSILON);
    assert!((q.z() - 3.5).abs() < f64::EPSILON);
    assert_eq!(*q.coords(), Vector3::new(1.0, -2.0, 3.5));
    assert_eq!(*q.vector(), Vector3::new(1.0, -2.0, 3.5));
    assert_eq!(q.into_coords(), Vector3::new(1.0, -2.0, 3.5));
    assert_eq!(q.into_vector(), Vector3::new(1.0, -2.0, 3.5));
    // |(1,-2,3.5)| = sqrt(1 + 4 + 12.25) = sqrt(17.25)
    assert!((q.norm() - 17.25f64.sqrt()).abs() < 1e-15);
}

#[test]
fn distance_to_is_the_euclidean_separation() {
    // (1,2,3) -> (4,6,3): (3,4,0), |.| = sqrt(9+16) = 5 exactly.
    let a = p(1.0, 2.0, 3.0);
    assert!((a.distance_to(&p(4.0, 6.0, 3.0)) - 5.0).abs() < 1e-15);
    assert!((a.distance_to(&a)).abs() < 1e-15);
}

#[test]
fn point_minus_point_is_a_displacement_vector() {
    // (4,5,6) - (1,2,3) = (3,3,3)
    let d = p(4.0, 5.0, 6.0) - p(1.0, 2.0, 3.0);
    assert!((d.x() - 3.0).abs() < 1e-15 && (d.y() - 3.0).abs() < 1e-15 && (d.z() - 3.0).abs() < 1e-15);
}

#[test]
fn point_translates_by_a_same_frame_vector() {
    let q = p(10.0, 20.0, 30.0);
    let d = v(1.0, -2.0, 3.0);
    assert_eq!(*(q + d).coords(), Vector3::new(11.0, 18.0, 33.0));
    assert_eq!(*(q - d).coords(), Vector3::new(9.0, 22.0, 27.0));
    let mut a = q;
    a += d;
    assert_eq!(*a.coords(), Vector3::new(11.0, 18.0, 33.0));
    a -= d;
    assert_eq!(*a.coords(), Vector3::new(10.0, 20.0, 30.0));
}

#[test]
fn vector_accessors_do_not_permute_axes() {
    let a = v(1.0, -2.0, 3.5);
    assert!((a.x() - 1.0).abs() < f64::EPSILON);
    assert!((a.y() + 2.0).abs() < f64::EPSILON);
    assert!((a.z() - 3.5).abs() < f64::EPSILON);
    assert_eq!(*a.vector(), Vector3::new(1.0, -2.0, 3.5));
    assert_eq!(*a.coords(), Vector3::new(1.0, -2.0, 3.5));
    assert_eq!(a.into_vector(), Vector3::new(1.0, -2.0, 3.5));
    assert_eq!(a.into_coords(), Vector3::new(1.0, -2.0, 3.5));
    assert_eq!(SpatialVector::<Ned>::zero().into_vector(), Vector3::zeros());
    assert_eq!(SpatialVector::<Ned>::from_vector(Vector3::new(1.0, 2.0, 3.0)), v(1.0, 2.0, 3.0));
    assert_eq!(SpatialVector::<Ned>::from(Vector3::new(1.0, 2.0, 3.0)), v(1.0, 2.0, 3.0));
    assert_eq!(Vector3::<f64>::from(v(1.0, 2.0, 3.0)), Vector3::new(1.0, 2.0, 3.0));
}

#[test]
fn dot_and_norm_match_hand_arithmetic() {
    let a = v(1.0, 2.0, 3.0);
    let b = v(4.0, 5.0, 6.0);
    assert!((a.dot(&b) - 32.0).abs() < 1e-15, "dot = {}", a.dot(&b));
    assert!((a.norm() - 14.0f64.sqrt()).abs() < 1e-15); // sqrt(1+4+9)
    assert!((b.norm() - 77.0f64.sqrt()).abs() < 1e-15); // sqrt(16+25+36)
}

#[test]
fn cross_follows_the_right_hand_rule() {
    let a = v(1.0, 2.0, 3.0);
    let b = v(4.0, 5.0, 6.0);
    let ab = a.cross(&b);
    assert!((ab.x() + 3.0).abs() < 1e-15 && (ab.y() - 6.0).abs() < 1e-15 && (ab.z() + 3.0).abs() < 1e-15,
        "a x b = {:?}", ab.into_vector());
    // Antisymmetry: swapping the operands negates the result exactly.
    let ba = b.cross(&a);
    assert!((ab + ba).norm() < 1e-15);
    // A vector is perpendicular to its own cross product.
    assert!(a.dot(&ab).abs() < 1e-15);
    // Cross with itself vanishes.
    assert!(a.cross(&a).norm() < 1e-15);
}

#[test]
fn vector_arithmetic_and_scaling() {
    let a = v(1.0, 2.0, 3.0);
    let b = v(4.0, 5.0, 6.0);
    assert_eq!(a + b, v(5.0, 7.0, 9.0));
    assert_eq!(b - a, v(3.0, 3.0, 3.0));
    assert_eq!(-a, v(-1.0, -2.0, -3.0));
    assert_eq!(a * 2.0, v(2.0, 4.0, 6.0));
    assert_eq!(2.0 * a, v(2.0, 4.0, 6.0)); // scalar on the left commutes
    assert_eq!(a / 2.0, v(0.5, 1.0, 1.5));

    let mut c = a;
    c += b;
    assert_eq!(c, v(5.0, 7.0, 9.0));
    c -= b;
    assert_eq!(c, a);
    c *= 3.0;
    assert_eq!(c, v(3.0, 6.0, 9.0));
    c /= 3.0;
    assert_eq!(c, a);
}

#[test]
fn velocity_accessors_displacement_and_arithmetic() {
    let u = SpatialVelocity::<Ecef<Itrf2014>>::new(1.0, -2.0, 3.5);
    assert!((u.x() - 1.0).abs() < f64::EPSILON);
    assert!((u.y() + 2.0).abs() < f64::EPSILON);
    assert!((u.z() - 3.5).abs() < f64::EPSILON);
    assert!((u.norm() - 17.25f64.sqrt()).abs() < 1e-15);
    assert_eq!(*u.vector(), Vector3::new(1.0, -2.0, 3.5));
    assert_eq!(u.into_vector(), Vector3::new(1.0, -2.0, 3.5));
    // displacement(dt) = v * dt: at dt = 0.5 s, (1,-2,3.5) -> (0.5,-1,1.75)
    let d = u.displacement(0.5);
    assert!((d.x() - 0.5).abs() < 1e-15 && (d.y() + 1.0).abs() < 1e-15 && (d.z() - 1.75).abs() < 1e-15);

    let w = SpatialVelocity::<Ecef<Itrf2014>>::new(-1.0, 4.0, 0.5);
    assert_eq!(u + w, SpatialVelocity::new(0.0, 2.0, 4.0));
    assert_eq!(u - w, SpatialVelocity::new(2.0, -6.0, 3.0));
    assert_eq!(-u, SpatialVelocity::new(-1.0, 2.0, -3.5));
    assert_eq!(u * 2.0, SpatialVelocity::new(2.0, -4.0, 7.0));
    assert_eq!(u / 2.0, SpatialVelocity::new(0.5, -1.0, 1.75));
    assert_eq!(SpatialVelocity::<Ecef<Itrf2014>>::zero(), SpatialVelocity::new(0.0, 0.0, 0.0));
    assert_eq!(
        SpatialVelocity::<Ecef<Itrf2014>>::from_vector(Vector3::new(1.0, 2.0, 3.0)),
        SpatialVelocity::new(1.0, 2.0, 3.0)
    );
    assert_eq!(
        SpatialVelocity::<Ecef<Itrf2014>>::from(Vector3::new(1.0, 2.0, 3.0)),
        SpatialVelocity::new(1.0, 2.0, 3.0)
    );
    assert_eq!(Vector3::<f64>::from(u), Vector3::new(1.0, -2.0, 3.5));
}

#[test]
fn debug_output_names_the_coordinate_frame() {
    assert_eq!(format!("{:?}", v(1.0, 2.0, 3.0)), "SpatialVector<NED>([1, 2, 3])");
    assert_eq!(format!("{:?}", p(1.0, 2.0, 3.0)), "Point3<NED>([1, 2, 3])");
    assert_eq!(
        format!("{:?}", SpatialVelocity::<Ecef<Itrf2014>>::new(1.0, 2.0, 3.0)),
        "SpatialVelocity<ITRF2014>([1, 2, 3])"
    );
}

#[test]
fn frame_conversion_goes_through_the_helmert_hub() {
    // NAD83(2011) -> ITRF2014 at epoch 2010.0 applies t = (-1.0053, +1.9092, +0.5416) m
    // (the published tx/ty/tz in millimetres) on top of the rotation and scale.
    let pos = EcefPos::<Nad83_2011>::new(Vector3::new(0.0, 0.0, 0.0));
    let got = pos.convert_to::<Itrf2014>(2010.0).into_vector();
    let expect = Vector3::new(-1.0053, 1.9092, 0.5416);
    assert!((got - expect).norm() < 1e-9, "got {got:?} want {expect:?}");
    // And the reverse conversion returns to the original point.
    let back = pos.convert_to::<Itrf2014>(2010.0).convert_to::<Nad83_2011>(2010.0);
    assert!(back.into_vector().norm() < 1e-6);
}

#[test]
fn converting_the_hub_frame_to_itself_is_bit_exact() {
    // ITRF2014 is the hub: its HELMERT link is None, so `params_at` must yield
    // identity_at(t) on both sides of the conversion and the point must survive
    // bit-exactly at any epoch.
    let v = Vector3::new(1_234.5, -6_789.0, 10_111.0);
    let pos = EcefPos::<Itrf2014>::new(v);
    for t in [2010.0_f64, 2025.0, 2040.0] {
        assert_eq!(pos.convert_to::<Itrf2014>(t).into_vector(), v);
    }
}