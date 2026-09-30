use core::f64::consts::FRAC_PI_2;
use gneiss_core::coords::{az_el, ecef_cov_to_enu_std};
use gneiss_core::frames::{
    AntennaLeverArm, Attitude, BodyCovariance, BodyFrd, BodyToEcef, BodyVector, BodyVelocity,
    CoordinateFrame, Ecef, EcefCovariance, EcefPos, EcefVector, EcefVelocity, Enu, EnuCovariance,
    EnuVector, Itrf2014, Itrf2020, Jgd2011, LocalTangentPlane, Nad83, Ned, NedCovariance,
    NedVector, Point3, Pz90, ReferenceFrame, SpatialCovariance, SpatialVector, SpatialVelocity,
    Wgs84,
};
use nalgebra::{Matrix3, UnitQuaternion, Vector3};

#[test]
fn test_zero_cost_memory_layouts() {
    use core::mem::{align_of, size_of};
    assert_eq!(size_of::<Point3<Ecef<Itrf2014>>>(), size_of::<Vector3<f64>>());
    assert_eq!(align_of::<Point3<Ecef<Itrf2014>>>(), align_of::<Vector3<f64>>());

    assert_eq!(size_of::<SpatialVector<BodyFrd>>(), size_of::<Vector3<f64>>());
    assert_eq!(size_of::<SpatialVelocity<BodyFrd>>(), size_of::<Vector3<f64>>());
    assert_eq!(size_of::<SpatialCovariance<Ned>>(), size_of::<Matrix3<f64>>());
    assert_eq!(size_of::<AntennaLeverArm>(), size_of::<Vector3<f64>>());
    assert_eq!(size_of::<BodyToEcef<Itrf2014>>(), size_of::<UnitQuaternion<f64>>());
}

#[test]
fn test_reference_frame_realizations() {
    assert_eq!(<Itrf2014 as ReferenceFrame>::NAME, "ITRF2014");
    assert_eq!(<Itrf2020 as ReferenceFrame>::NAME, "ITRF2020");
    assert_eq!(<Wgs84 as ReferenceFrame>::NAME, "WGS84(Broadcast)");
    assert_eq!(<Nad83 as ReferenceFrame>::NAME, "NAD83(2011)");
    assert_eq!(<Jgd2011 as ReferenceFrame>::NAME, "JGD2011");
    assert_eq!(<Pz90 as ReferenceFrame>::NAME, "PZ-90.11");
    assert!(Pz90::HELMERT_TO_ITRF2014.is_some());

    let ecef_fixed = <Ecef<Itrf2014> as CoordinateFrame>::IS_EARTH_FIXED;
    assert!(ecef_fixed);
    let ned_fixed = <Ned as CoordinateFrame>::IS_EARTH_FIXED;
    assert!(!ned_fixed);
    let enu_fixed = <Enu as CoordinateFrame>::IS_EARTH_FIXED;
    assert!(!enu_fixed);
    let body_fixed = <BodyFrd as CoordinateFrame>::IS_EARTH_FIXED;
    assert!(!body_fixed);
}

#[test]
fn test_pz90_helmert_roundtrip() {
    let p_pz90 = EcefPos::<Pz90>::new(Vector3::new(3_000_000.0, 1_000_000.0, 5_000_000.0));
    let p_itrf: EcefPos<Itrf2014> = p_pz90.convert_to(2020.0);
    let p_back: EcefPos<Pz90> = p_itrf.convert_to(2020.0);
    let diff = (p_back.vector() - p_pz90.vector()).norm();
    assert!(diff < 1e-6, "PZ-90 roundtrip diff was {diff:e} m");
}

#[test]
fn test_spatial_vector_linear_algebra() {
    let v1: BodyVector = SpatialVector::<BodyFrd>::new(1.0, 2.0, 3.0);
    let v2: BodyVector = SpatialVector::<BodyFrd>::new(4.0, 5.0, 6.0);
    assert_eq!((v1 + v2).into_vector(), Vector3::new(5.0, 7.0, 9.0));
    assert_eq!((v2 - v1).into_vector(), Vector3::new(3.0, 3.0, 3.0));
    assert_eq!((-v1).into_vector(), Vector3::new(-1.0, -2.0, -3.0));
    assert_eq!((v1 * 2.0).into_vector(), Vector3::new(2.0, 4.0, 6.0));
    assert_eq!((v1 / 2.0).into_vector(), Vector3::new(0.5, 1.0, 1.5));
    assert_eq!(v1.dot(&v2), 32.0);
    assert_eq!(v1.cross(&v2).into_vector(), Vector3::new(-3.0, 6.0, -3.0));
}

#[test]
fn test_spatial_velocity_and_displacement() {
    let vel = SpatialVelocity::<Ned>::new(10.0, -5.0, 2.0);
    let disp = vel.displacement(0.5);
    assert_eq!(disp.into_vector(), Vector3::new(5.0, -2.5, 1.0));
}

#[test]
fn test_spatial_covariance_std_accessors() {
    let ned_cov = NedCovariance::from_stds(0.02, 0.03, 0.05);
    assert!((ned_cov.std_north() - 0.02).abs() < 1e-9);
    assert!((ned_cov.std_east() - 0.03).abs() < 1e-9);
    assert!((ned_cov.std_down() - 0.05).abs() < 1e-9);

    let enu_cov = EnuCovariance::from_stds(0.03, 0.02, 0.05);
    assert!((enu_cov.std_east() - 0.03).abs() < 1e-9);
    assert!((enu_cov.std_north() - 0.02).abs() < 1e-9);
    assert!((enu_cov.std_up() - 0.05).abs() < 1e-9);
}

#[test]
fn test_antenna_lever_arm_and_attitude_rotation() {
    let arm = AntennaLeverArm::new(0.5, 0.2, -1.0);
    assert_eq!(arm.forward(), 0.5);
    assert_eq!(arm.right(), 0.2);
    assert_eq!(arm.down(), -1.0);

    // 90 degree yaw rotation about down (+Z)
    let q = UnitQuaternion::from_euler_angles(0.0, 0.0, FRAC_PI_2);
    let att = Attitude::<BodyFrd, Ned>::from_unit_quaternion(q);
    let rotated = att.rotate(arm.as_body_vector());

    // In Body: +X=0.5 (fwd), +Y=0.2 (right) -> in NED after 90 deg yaw:
    // +X (fwd) -> North rotated 90 deg -> East = +0.5
    // +Y (right) -> East rotated 90 deg -> South = -0.2 (North = -0.2)
    assert!((rotated.x() - (-0.2)).abs() < 1e-9);
    assert!((rotated.y() - 0.5).abs() < 1e-9);
    assert!((rotated.z() - (-1.0)).abs() < 1e-9);

    let back = att.inverse().rotate(&rotated);
    assert!((back.vector() - arm.as_body_vector().vector()).norm() < 1e-9);
}

#[test]
fn test_attitude_velocity_and_covariance_rotation() {
    let q = UnitQuaternion::from_euler_angles(0.0, 0.0, FRAC_PI_2);
    let att = BodyToEcef::<Itrf2014>::from_unit_quaternion(q);

    let v_b = BodyVelocity::new(10.0, 0.0, 0.0);
    let v_e: EcefVelocity<Itrf2014> = att.rotate_velocity(v_b);
    assert!((v_e.norm() - 10.0).abs() < 1e-9);

    let cov_b = BodyCovariance::from_stds(0.1, 0.2, 0.3);
    let cov_e = att.rotate_cov(&cov_b);
    assert!((cov_e.matrix().trace() - cov_b.matrix().trace()).abs() < 1e-9);
}

#[test]
fn test_point_vector_affine_algebra() {
    let p1 = EcefPos::<Itrf2014>::new(Vector3::new(100.0, 200.0, 300.0));
    let v = EcefVector::<Itrf2014>::new(10.0, -20.0, 30.0);
    let p2 = p1 + v;
    assert_eq!(p2.into_vector(), Vector3::new(110.0, 180.0, 330.0));

    let delta = p2 - p1;
    assert_eq!(delta.into_vector(), v.into_vector());

    let p1_recovered = p2 - v;
    assert_eq!(p1_recovered.into_vector(), p1.into_vector());
}

#[test]
fn test_local_tangent_plane_relational_coupling() {
    let origin = EcefPos::<Itrf2014>::new(Vector3::new(4_027_893.0, 307_041.0, 4_919_475.0));
    let plane = LocalTangentPlane::from_origin(origin);

    let target = EcefPos::<Itrf2014>::new(Vector3::new(4_027_900.0, 307_050.0, 4_919_480.0));
    let enu: EnuVector = plane.to_enu(&target);
    let target_rec = plane.from_enu(&enu);
    assert!((target_rec.vector() - target.vector()).norm() < 1e-9);

    let ned: NedVector = plane.to_ned(&target);
    let target_rec_ned = plane.from_ned(&ned);
    assert!((target_rec_ned.vector() - target.vector()).norm() < 1e-9);

    // Azimuth & Elevation relational agreement
    let sat = EcefPos::<Itrf2014>::new(Vector3::new(15_000_000.0, 10_000_000.0, 20_000_000.0));
    let (az_plane, el_plane) = plane.az_el(&sat);
    let (az_legacy, el_legacy) = az_el(*plane.origin_llh(), *origin.coords(), *sat.coords());
    assert!((az_plane - az_legacy).abs() < 1e-12);
    assert!((el_plane - el_legacy).abs() < 1e-12);
}

#[test]
fn test_local_tangent_plane_covariance_projection() {
    let origin = EcefPos::<Itrf2014>::new(Vector3::new(4_027_893.0, 307_041.0, 4_919_475.0));
    let plane = LocalTangentPlane::from_origin(origin);

    let cov_ecef = EcefCovariance::<Itrf2014>::from_diagonal(Vector3::new(0.04, 0.09, 0.16));
    let enu_cov = plane.project_cov_to_enu(&cov_ecef);

    let (std_e_leg, std_n_leg, std_u_leg) = ecef_cov_to_enu_std(*origin.coords(), *cov_ecef.matrix());
    assert!((enu_cov.std_east() - std_e_leg).abs() < 1e-9);
    assert!((enu_cov.std_north() - std_n_leg).abs() < 1e-9);
    assert!((enu_cov.std_up() - std_u_leg).abs() < 1e-9);
}
