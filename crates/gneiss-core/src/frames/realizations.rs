//! Terrestrial reference frame realization tags and their official Helmert transformation links.

use super::helmert::HelmertParams;

/// Marker trait binding a coordinate realization to its published Helmert link into the ITRF2014 hub.
pub trait ReferenceFrame {
    const NAME: &'static str;
    const HELMERT_TO_ITRF2014: Option<HelmertParams>;
}

pub(crate) const ITRF2020_TO_ITRF2014: HelmertParams = HelmertParams {
    tx_mm: -1.4,
    ty_mm: -0.9,
    tz_mm: 1.4,
    scale_ppb: -0.42,
    rx_mas: 0.0,
    ry_mas: 0.0,
    rz_mas: 0.0,
    ref_epoch_yr: 2015.0,
    tx_rate: 0.0,
    ty_rate: -0.1,
    tz_rate: 0.2,
    scale_rate: 0.0,
    rx_rate: 0.0,
    ry_rate: 0.0,
    rz_rate: 0.0,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Itrf2014;
impl ReferenceFrame for Itrf2014 {
    const NAME: &'static str = "ITRF2014";
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = None;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Itrf2020;
impl ReferenceFrame for Itrf2020 {
    const NAME: &'static str = "ITRF2020";
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(ITRF2020_TO_ITRF2014);
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Igs20;
impl ReferenceFrame for Igs20 {
    const NAME: &'static str = "IGS20";
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(ITRF2020_TO_ITRF2014);
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wgs84Broadcast;
impl ReferenceFrame for Wgs84Broadcast {
    const NAME: &'static str = "WGS84(Broadcast)";
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(ITRF2020_TO_ITRF2014);
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Nad83_2011;
impl ReferenceFrame for Nad83_2011 {
    const NAME: &'static str = "NAD83(2011)";
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(HelmertParams {
        tx_mm: 0.99,
        ty_mm: -1.91,
        tz_mm: -0.51,
        scale_ppb: -1.65,
        rx_mas: 0.0267,
        ry_mas: 0.0005,
        rz_mas: 0.0074,
        ref_epoch_yr: 2010.0,
        tx_rate: -0.067,
        ty_rate: 0.757,
        tz_rate: 0.019,
        rx_rate: 0.0,
        ry_rate: 0.0,
        rz_rate: 0.0,
        scale_rate: 0.102,
    });
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Etrs89;
impl ReferenceFrame for Etrs89 {
    const NAME: &'static str = "ETRS89(ETRF2000)";
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(HelmertParams {
        tx_mm: 52.1,
        ty_mm: 49.3,
        tz_mm: -58.5,
        scale_ppb: -1.04,
        rx_mas: 0.891,
        ry_mas: 5.39,
        rz_mas: -8.71,
        ref_epoch_yr: 1989.0,
        tx_rate: 0.1,
        ty_rate: 0.1,
        tz_rate: -1.8,
        rx_rate: 0.0,
        ry_rate: 0.0,
        rz_rate: 0.0,
        scale_rate: -0.08,
    });
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Gda2020;
impl ReferenceFrame for Gda2020 {
    const NAME: &'static str = "GDA2020";
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(HelmertParams {
        tx_mm: 0.0,
        ty_mm: 0.0,
        tz_mm: 0.0,
        scale_ppb: 0.0,
        rx_mas: 0.0,
        ry_mas: 0.0,
        rz_mas: 0.0,
        ref_epoch_yr: 2020.0,
        tx_rate: 0.0,
        ty_rate: 0.0,
        tz_rate: 0.0,
        rx_rate: 0.0,
        ry_rate: 0.0,
        rz_rate: 0.0,
        scale_rate: 0.0,
    });
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Jgd2011;
impl ReferenceFrame for Jgd2011 {
    const NAME: &'static str = "JGD2011";
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(HelmertParams {
        tx_mm: 0.0,
        ty_mm: 0.0,
        tz_mm: 0.0,
        scale_ppb: 0.0,
        rx_mas: 0.0,
        ry_mas: 0.0,
        rz_mas: 0.0,
        ref_epoch_yr: 2011.0,
        tx_rate: 0.0,
        ty_rate: 0.0,
        tz_rate: 0.0,
        rx_rate: 0.0,
        ry_rate: 0.0,
        rz_rate: 0.0,
        scale_rate: 0.0,
    });
}
