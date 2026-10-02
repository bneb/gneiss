#![allow(clippy::unwrap_used)]

use super::*;
use gneiss_core::ephemeris::Ephemeris;
use std::io::BufReader;

    /// RINEX 3 navigation records are written 80 columns wide: a field the
    /// producer left empty is blank space, not a missing byte. These fixtures
    /// are hand-written, so pad every line to the record width; without it a
    /// short line looks like a record whose value columns never arrived.
    fn padded(text: &str) -> String {
        let mut out = String::with_capacity(text.len() + text.lines().count() * 20);
        for line in text.lines() {
            out.push_str(line);
            for _ in line.len()..80 {
                out.push(' ');
            }
            out.push('\n');
        }
        out
    }

mod part1;
mod part2;
