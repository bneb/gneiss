//! Diagnostic reporting: per-channel verdicts, one-line summaries and the
//! CSV dump consumed by benchmark reports.

use std::io::{BufWriter, Write};
use std::path::Path;

use super::phase::{BinStat, Diagnostics, SIDEREAL_PERIOD_S};

/// Diagnostic dump for one channel: verdict plus its full bin histogram.
#[derive(Debug, Clone)]
pub struct ChannelDump {
    pub channel: String,
    pub diag: Diagnostics,
    pub corr_rms_m: f64,
    pub active_bins: usize,
    pub bins: Vec<BinStat>,
}

/// Result of [`super::apply_to_trajectory`] for logging and CSV dumps.
///
/// `before`/`after` are both computed over the SECOND half of the session
/// (the only part mitigation touches): pre-correction vs post-correction.
#[derive(Debug, Clone, Default)]
pub struct SiderealReport {
    pub period_s: f64,
    pub n_bins: usize,
    pub split_index: usize,
    pub samples: usize,
    pub channels: Vec<ChannelDump>,
    pub before: super::HalfMetrics,
    pub after: super::HalfMetrics,
}

impl SiderealReport {
    /// Compact one-line-per-channel log summary with before/after metrics.
    #[must_use]
    pub fn summary(&self, label: &str) -> String {
        let head = format!(
            "SIDEREAL [{label}] period={:.0}s bins={} epochs={} split@{}",
            self.period_s, self.n_bins, self.samples, self.split_index
        );
        let chans: String = self.channels.iter().map(fmt_channel).collect();
        let tail = format!(
            " | 2nd-half pre->post(N={}): h_p50 {:.3}->{:.3} h_p95 {:.3}->{:.3} v_p50 {:+.3}->{:+.3} v_rms {:.3}->{:.3}",
            self.after.n,
            self.before.h_p50, self.after.h_p50,
            self.before.h_p95, self.after.h_p95,
            self.before.v_p50, self.after.v_p50,
            self.before.v_rms, self.after.v_rms,
        );
        format!("{head}{chans}{tail}")
    }
}

fn fmt_channel(c: &ChannelDump) -> String {
    let verdict = if c.diag.is_structured() { "STRUCTURED" } else { "UNSTRUCTURED" };
    format!(
        " | {} chi2/dof={:.2} p={:.1e} {} corr={:.0}mm({} bins)",
        c.channel,
        c.diag.chi2 / c.diag.dof.max(1) as f64,
        c.diag.p_value,
        verdict,
        c.corr_rms_m * 1000.0,
        c.active_bins,
    )
}

/// Row counts returned by [`write_diag_csv`] (for cheap sanity checks).
#[derive(Debug, Clone, Copy, Default)]
pub struct CsvWritten {
    pub channel_rows: usize,
    pub bin_rows: usize,
}

/// Dump per-channel verdicts and full per-bin histograms as CSV.
/// Comment lines (#) carry session metadata; the two tables are
/// machine-parseable by splitting on their header lines.
pub fn write_diag_csv(
    path: &Path,
    target: &str,
    dumps: &[ChannelDump],
    n_bins: usize,
    split_index: usize,
) -> std::io::Result<CsvWritten> {
    let f = std::fs::File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(w, "# gneiss sidereal diagnostic")?;
    writeln!(
        w, "# target={target} period_s={:.0} bins={} split_index={}",
        SIDEREAL_PERIOD_S, n_bins, split_index
    )?;
    let mut written = CsvWritten::default();
    writeln!(w, "channel,chi2,dof,p_value,structured,corr_rms_m,active_bins")?;
    written.channel_rows = write_channel_rows(&mut w, dumps)?;
    writeln!(w, "channel,bin,phase_lo,count,mean,std")?;
    written.bin_rows = write_bin_rows(&mut w, dumps, n_bins)?;
    w.flush()?;
    Ok(written)
}

fn write_channel_rows(w: &mut BufWriter<std::fs::File>, dumps: &[ChannelDump]) -> std::io::Result<usize> {
    let mut n = 0;
    for d in dumps {
        writeln!(
            w, "{},{:.6},{},{:.6e},{},{:.6e},{}",
            d.channel, d.diag.chi2, d.diag.dof, d.diag.p_value,
            d.diag.is_structured(), d.corr_rms_m, d.active_bins,
        )?;
        n += 1;
    }
    Ok(n)
}

fn write_bin_rows(
    w: &mut BufWriter<std::fs::File>,
    dumps: &[ChannelDump],
    n_bins: usize,
) -> std::io::Result<usize> {
    let mut n = 0;
    for d in dumps {
        for (i, b) in d.bins.iter().enumerate() {
            writeln!(
                w, "{},{},{:.6},{},{:.6e},{:.6e}",
                d.channel, i, i as f64 / n_bins.max(1) as f64,
                b.count, b.mean, b.std(),
            )?;
            n += 1;
        }
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_written_default_is_zeroed() {
        let w = CsvWritten::default();
        assert_eq!(w.channel_rows, 0);
        assert_eq!(w.bin_rows, 0);
    }

    fn dump(channel: &str, structured: bool) -> ChannelDump {
        // Build the histogram through the public folding API.
        let bins = crate::post_process::sidereal::phase::fold(
            &[
                crate::post_process::sidereal::phase::Sample { tow_s: 10.0, value: 0.01 },
                crate::post_process::sidereal::phase::Sample { tow_s: 20.0, value: 0.02 },
                crate::post_process::sidereal::phase::Sample { tow_s: 30.0, value: 0.00 },
            ],
            1,
        );
        ChannelDump {
            channel: channel.to_string(),
            diag: Diagnostics {
                chi2: if structured { 480.0 } else { 230.0 },
                dof: 239,
                p_value: if structured { 1.0e-6 } else { 0.5 },
            },
            corr_rms_m: 0.0412,
            active_bins: 238,
            bins,
        }
    }

    #[test]
    fn summary_reports_verdict_and_before_after_metrics() {
        let rep = SiderealReport {
            period_s: SIDEREAL_PERIOD_S,
            n_bins: 240,
            split_index: 1440,
            samples: 2880,
            channels: vec![dump("east", true), dump("north", false)],
            before: super::super::HalfMetrics { h_p50: 0.142, h_p95: 0.402, v_p50: 0.010, ..Default::default() },
            after: super::super::HalfMetrics { h_p50: 0.118, h_p95: 0.351, v_p50: -0.002, ..Default::default() },
        };
        let s = rep.summary("Smoothed/P181");
        assert!(s.contains("SIDEREAL [Smoothed/P181]"), "{s}");
        assert!(s.contains("east chi2/dof=2.01"), "{s}");
        assert!(s.contains("STRUCTURED") && s.contains("UNSTRUCTURED"), "{s}");
        assert!(s.contains("corr=41mm"), "{s}");
        assert!(s.contains("h_p50 0.142->0.118"), "{s}");
    }

    #[test]
    fn write_diag_csv_produces_header_channel_rows_and_bins() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("diag.csv");
        let dumps = vec![dump("east", true), dump("up", false)];
        let written =
            write_diag_csv(&path, "P181/Smoothed", &dumps, 240, 360).expect("csv write");
        assert_eq!(written.channel_rows, 2);
        assert_eq!(written.bin_rows, 2);
        let text = std::fs::read_to_string(&path).expect("read back");
        assert!(text.contains("# gneiss sidereal diagnostic"), "{text}");
        assert!(text.contains("target=P181/Smoothed"), "{text}");
        assert!(text.contains("channel,bin,phase_lo,count,mean,std"), "{text}");
        // Per-bin rows carry channel name and phase bin index.
        assert!(text.contains("east,0,0.000000,3,"), "{text}");
    }

    #[test]
    fn fmt_channel_handles_zero_dof_without_nan() {
        let mut c = dump("up", false);
        c.diag.dof = 0;
        c.diag.chi2 = 0.0;
        let s = fmt_channel(&c);
        assert!(!s.contains("NaN"), "{s}");
    }
}
