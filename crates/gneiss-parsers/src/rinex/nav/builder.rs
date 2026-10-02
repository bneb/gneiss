//! Ephemeris builder functions mapping raw RINEX broadcast values into typed Ephemeris structs.

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

/// Closed-form range of the broadcast sqrt(A) parameter.
///
/// The propagator takes `a = sqrt_a * sqrt_a` as the semi-major axis **in
/// metres**, so the quoted sqrt(A) is sqrt(a) with a in metres. GNSS broadcast
/// orbits run from 25 510 km (GLONASS) to 42 164 km (IGSO / geostationary),
/// i.e. sqrt(a) from sqrt(25 510 000) = 5051.7 to sqrt(42 164 000) = 6493.4.
/// The bounds below are rounded outward from those, so a record is rejected
/// only when its semi-major axis is unambiguously not an orbit.
const SQRT_A_MIN: f64 = 5000.0;
const SQRT_A_MAX: f64 = 6600.0;

/// Rejects a Keplerian ephemeris whose semi-major axis is not physical.
///
/// A RINEX record whose value columns are unreadable decodes to zeros rather
/// than failing, and a zero semi-major axis propagates to the Earth's centre
/// rather than raising - exactly the "confident garbage" this guards against.
fn sqrt_a_is_physical(sqrt_a: f64) -> bool {
    sqrt_a.is_finite() && (SQRT_A_MIN..=SQRT_A_MAX).contains(&sqrt_a)
}

pub(crate) fn build_ephemeris(
    constellation: Constellation,
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<Ephemeris> {
    let built = build_typed(constellation, sat, toc, af0, af1, af2, vals)?;
    match broadcast_sqrt_a(&built) {
        // GLONASS carries a PZ-90 state vector, not a Keplerian orbit.
        None => Some(built),
        Some(sqrt_a) if sqrt_a_is_physical(sqrt_a) => Some(built),
        Some(_) => None,
    }
}

/// The broadcast sqrt(A) of a Keplerian ephemeris, or `None` for GLONASS.
fn broadcast_sqrt_a(e: &Ephemeris) -> Option<f64> {
    match e {
        Ephemeris::Gps(g) => Some(g.sqrt_a),
        Ephemeris::Galileo(g) => Some(g.sqrt_a),
        Ephemeris::Beidou(b) => Some(b.sqrt_a),
        Ephemeris::Qzss(q) => Some(q.sqrt_a),
        Ephemeris::Glonass(_) => None,
    }
}

fn build_typed(
    constellation: Constellation,
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<Ephemeris> {
    match constellation {
        Constellation::Glonass => build_glonass_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Gps => build_gps_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Galileo => build_galileo_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Beidou => build_beidou_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Qzss => build_qzss_ephemeris(sat, toc, af0, af1, af2, vals),
        Constellation::Navic => build_gps_ephemeris(sat, toc, af0, af1, af2, vals),
        _ => None,
    }
}

fn build_glonass_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<Ephemeris> {
    Some(Ephemeris::Glonass(
        gneiss_core::ephemeris::GlonassEphemeris {
            sat,
            toe: toc,
            freq_num: vals[7] as i8,
            tau_n: af0,
            gamma_n: af1,
            delta_tau_n: af2,
            x: vals[0] * 1000.0,
            y: vals[4] * 1000.0,
            z: vals[8] * 1000.0,
            vx: vals[1] * 1000.0,
            vy: vals[5] * 1000.0,
            vz: vals[9] * 1000.0,
            ax: vals[2] * 1000.0,
            ay: vals[6] * 1000.0,
            az: vals[10] * 1000.0,
        },
    ))
}

fn build_gps_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<Ephemeris> {
    Some(Ephemeris::Gps(
        gneiss_core::ephemeris::GpsEphemeris {
            sat,
            toc,
            toe: GpsTime::new(toc.week, vals[8]),
            af0,
            af1,
            af2,
            iode: vals[0] as u32,
            crs: vals[1],
            delta_n: vals[2],
            m0: vals[3],
            cuc: vals[4],
            e: vals[5],
            cus: vals[6],
            sqrt_a: vals[7],
            cic: vals[9],
            omega0: vals[10],
            cis: vals[11],
            i0: vals[12],
            crc: vals[13],
            omega: vals[14],
            omega_dot: vals[15],
            idot: vals[16],
            tgd: vals[22],
            iodc: vals[23] as u32,
        },
    ))
}

fn build_galileo_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<Ephemeris> {
    Some(Ephemeris::Galileo(
        gneiss_core::ephemeris::GalileoEphemeris {
            sat,
            toc,
            toe: GpsTime::new(toc.week, vals[8]),
            af0,
            af1,
            af2,
            iod_nav: vals[0] as u32,
            crs: vals[1],
            delta_n: vals[2],
            m0: vals[3],
            cuc: vals[4],
            e: vals[5],
            cus: vals[6],
            sqrt_a: vals[7],
            cic: vals[9],
            omega0: vals[10],
            cis: vals[11],
            i0: vals[12],
            crc: vals[13],
            omega: vals[14],
            omega_dot: vals[15],
            idot: vals[16],
            bgd_e1_e5a: vals[22],
            bgd_e1_e5b: vals[23],
        },
    ))
}

fn build_beidou_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    v: &[f64; 32],
) -> Option<Ephemeris> {
    // toc in rinex/nav/mod.rs had 14s added; convert back so toe & toc are consistent BDT
    let toc_bdt = GpsTime::new(toc.week, toc.tow - 14.0);
    let toe_bdt = GpsTime::new(toc_bdt.week, v[8]);
    Some(Ephemeris::Beidou(gneiss_core::ephemeris::BeidouEphemeris {
        sat, toc: toc_bdt, toe: toe_bdt, af0, af1, af2,
        aode: v[0] as u32, crs: v[1], delta_n: v[2], m0: v[3],
        cuc: v[4], e: v[5], cus: v[6], sqrt_a: v[7],
        cic: v[9], omega0: v[10], cis: v[11], i0: v[12],
        crc: v[13], omega: v[14], omega_dot: v[15], idot: v[16],
        tgd1: v[22], tgd2: v[23], aodc: v[25] as u32,
    }))
}

fn build_qzss_ephemeris(
    sat: SatelliteId,
    toc: GpsTime,
    af0: f64,
    af1: f64,
    af2: f64,
    vals: &[f64; 32],
) -> Option<Ephemeris> {
    Some(Ephemeris::Qzss(
        gneiss_core::ephemeris::QzssEphemeris {
            sat,
            toc,
            toe: GpsTime::new(toc.week, vals[8]),
            af0,
            af1,
            af2,
            iode: vals[0] as u32,
            crs: vals[1],
            delta_n: vals[2],
            m0: vals[3],
            cuc: vals[4],
            e: vals[5],
            cus: vals[6],
            sqrt_a: vals[7],
            cic: vals[9],
            omega0: vals[10],
            cis: vals[11],
            i0: vals[12],
            crc: vals[13],
            omega: vals[14],
            omega_dot: vals[15],
            idot: vals[16],
            tgd: vals[22],
            iodc: vals[23] as u32,
        },
    ))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// GLONASS field mapping, checked against the reference decoder
    /// (RTKLIB rinex.cc:1170-1176, `decode_geph`):
    ///
    /// ```c
    /// geph->pos[0]=data[3]*1E3; geph->pos[1]=data[7]*1E3; geph->pos[2]=data[11]*1E3;
    /// geph->vel[0]=data[4]*1E3; geph->vel[1]=data[8]*1E3; geph->vel[2]=data[12]*1E3;
    /// geph->acc[0]=data[5]*1E3; geph->acc[1]=data[9]*1E3; geph->acc[2]=data[13]*1E3;
    /// geph->svh   =(int)data[ 6];
    /// geph->frq   =(int)data[10];
    /// ```
    ///
    /// `data[0..3]` are the three values of the record's first line and
    /// `data[3..15]` the four values of each continuation line, so
    /// `data[3+k] == vals[k]` for k in 0..12. The builder must therefore read
    /// X/Y/Z from vals[0]/vals[4]/vals[8], V from 1/5/9, A from 2/6/10 and the
    /// frequency number from vals[7].
    #[test]
    fn glonass_fields_follow_the_reference_decoder_mapping() {
        //            X           Vx        Y           Vy        Z           Vz        Az        frq
        let mut vals = [0.0; 32];
        vals[0..12].copy_from_slice(&[
            -1.184028808594E+03, -2.155310630798E+00, 0.0, 0.0,
            1.322943017578E+04, -2.060539245605E+00, 0.0, 1.0,
            2.176970068359E+04, 1.135528564453E+00, -1.862645149231E-09, 0.0,
        ]);
        let sat = SatelliteId {
            constellation: Constellation::Glonass,
            prn: 1,
        };
        let toc = GpsTime::new(2137, 0.0);
        let eph = build_ephemeris(Constellation::Glonass, sat, toc, 1.0, 2.0, 3.0, &vals);
        let Ephemeris::Glonass(g) = eph.unwrap() else {
            panic!("expected GLONASS");
        };
        assert!((g.x - -1184028.808594).abs() < 1e-3);
        assert!((g.y - 13229430.17578).abs() < 1e-3);
        assert!((g.z - 21769700.68359).abs() < 1e-3);
        assert!((g.vx - -2155.310630798).abs() < 1e-6);
        assert!((g.vy - -2060.539245605).abs() < 1e-6);
        assert!((g.vz - 1135.528564453).abs() < 1e-6);
        assert!((g.az - -1.862645149231E-06).abs() < 1e-15);
        assert_eq!(g.freq_num, 1);
        assert_eq!((g.tau_n, g.gamma_n, g.delta_tau_n), (1.0, 2.0, 3.0));
    }

    /// GPS broadcast field mapping, checked against the real record at the top
    /// of datasets/wtzr_ppp_1224/BRDC00IGS_R_20203590000_01D_MN.rnx. RINEX 3
    /// puts three values on the first line and four on each of seven
    /// continuations, so vals[0..4] are IODE/Crs/delta-n/M0, vals[4..8] are
    /// Cuc/e/Cus/sqrt(A), vals[8..12] toe/Cic/Omega0/Cis, and the last four
    /// rows carry IDOT (vals[16]), TGD (vals[22]) and IODC (vals[23]).
    #[test]
    fn gps_fields_follow_the_broadcast_record_order() {
        let mut vals = [0.0; 32];
        vals[0] = 51.0;
        vals[5] = 1.020454068203E-02;
        vals[7] = 5153.695047379;
        vals[8] = 345600.0;
        vals[12] = 9.828138365184E-01;
        vals[15] = -7.690320332710E-09;
        vals[16] = 1.207193141583E-10;
        vals[22] = 5.122274160385E-09;
        vals[23] = 51.0;
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let toc = GpsTime::new(2137, 345600.0);
        let Ephemeris::Gps(g) =
            build_ephemeris(Constellation::Gps, sat, toc, 0.0, 0.0, 0.0, &vals).unwrap()
        else {
            panic!("expected GPS");
        };
        assert_eq!(g.iode, 51);
        assert_eq!(g.iodc, 51);
        assert!((g.e - 1.020454068203E-02).abs() < 1e-15);
        assert!((g.sqrt_a - 5153.695047379).abs() < 1e-9);
        assert!((g.i0 - 9.828138365184E-01).abs() < 1e-13);
        assert!((g.omega_dot - -7.690320332710E-09).abs() < 1e-20);
        assert!((g.idot - 1.207193141583E-10).abs() < 1e-22);
        assert!((g.tgd - 5.122274160385E-09).abs() < 1e-20);
        assert_eq!((g.toe.week, g.toe.tow), (2137, 345600.0));
    }
}
