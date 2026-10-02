#![allow(clippy::unwrap_used)]

use super::signals::*;
use super::*;
use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::obs::{ObsType, SignalCode};
use gneiss_core::sat::Constellation;

    /// Pack (num_bits, value) pairs into a byte buffer in Msb0 order.
    fn pack_bits(pairs: &[(usize, u64)]) -> Vec<u8> {
        let total_bits: usize = pairs.iter().map(|(b, _)| *b).sum();
        let mut bytes = vec![0u8; total_bits.div_ceil(8)];
        let mut pos = 0;
        for &(bits, val) in pairs {
            for i in 0..bits {
                if (val >> (bits - 1 - i)) & 1 != 0 {
                    bytes[pos / 8] |= 1 << (7 - (pos % 8));
                }
                pos += 1;
            }
        }
        bytes
    }

    fn payload_for_msm4_1sat_1sig() -> Vec<u8> {
        // MSM4 GPS (msg 1074), 1 satellite (sat 1), 1 signal (sig 1), 1 cell.
        // MSM4 has no DF419 extended sat info (MSM5/7 only): DF397 (8 bits)
        // then DF398 (10 bits) directly, per RTKLIB decode_msm4 / RTCM
        // 10403.3 field order.
        pack_bits(&[
            (12, 1074), (12, 1), (30, 5000), (1, 0), (3, 2), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 2), (10, 100),
            (15, 50), (22, 500), (4, 3), (1, 0), (6, 35),
        ])
    }

    fn payload_for_msm5_1sat_2sig() -> Vec<u8> {
        // Note: parse_signal_data reads field-by-field across all cells, not cell-by-cell.
        // Sat data order: DF397 (int ms) -> DF419 (ext info) -> DF398 (rough range) -> DF399 (rate).
        // Order: fine_pr[0..1], fine_ph[0..1], lock[0..1], half_cycle[0..1], cnr[0..1], rates[0..1]
        pack_bits(&[
            (12, 1075), (12, 1), (30, 5000), (1, 0), (3, 2), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, (1u64 << 31) | (1u64 << 30)), (2, 0b11),
            (8, 2), (4, 5), (10, 100), (14, 25),
            (15, 50), (15, 100),  // fine_pr[0], fine_pr[1]
            (22, 500), (22, 1000), // fine_ph[0], fine_ph[1]
            (4, 3), (4, 5),        // lock[0], lock[1]
            (1, 0), (1, 1),        // half_cycle[0]=false, half_cycle[1]=true
            (6, 35), (6, 40),      // cnr[0], cnr[1]
            (15, 10), (15, 20),    // fine_phase_range_rates[0], [1]
        ])
    }

    fn payload_for_msm6_1sat_1sig() -> Vec<u8> {
        // MSM6 has no DF419 extended sat info either (MSM5/7 only).
        pack_bits(&[
            (12, 1076), (12, 1), (30, 5000), (1, 0), (3, 2), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 2), (10, 100),
            (20, 50), (24, 500), (10, 7), (1, 0), (10, 35),
        ])
    }

    fn payload_for_msm7_1sat_1sig() -> Vec<u8> {
        pack_bits(&[
            (12, 1077), (12, 1), (30, 5000), (1, 0), (3, 2), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 2), (4, 5), (10, 100), (14, 25),
            (20, 50), (24, 500), (10, 7), (1, 0), (10, 35), (15, 10),
        ])
    }

mod part1;
mod part2;
