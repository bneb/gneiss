//! Paired comparison of two evaluation runs over the same epochs.
//!
//! Reads two per-epoch error dumps (written by `GNEISS_ERR_DUMP` in
//! `eval_qinertia_ppk`) and reports:
//!   - the raw CDF of each run, which is the salient result
//!   - a Wilcoxon signed-rank test, since the runs share epochs and are paired
//!   - a Weibull fit of each, to compare tail shape
//!
//! Usage: eval_compare <a.txt> <b.txt> [label_a label_b]

use gneiss_core::stats::{bootstrap_cdf_band_paired, cdf_levels, weibull_with_tail, wilcoxon_signed_rank};

type Row = (u32, u8, f64);

fn load(path: &str) -> Result<Vec<Row>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path))?;
    Ok(text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut it = l.split_whitespace();
            let tow = it.next()?.parse().ok()?;
            let q = it.next()?.parse().ok()?;
            let h = it.next()?.parse().ok()?;
            Some((tow, q, h))
        })
        .collect())
}

fn percentile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let i = ((sorted.len() as f64 * q).floor() as usize).min(sorted.len() - 1);
    sorted[i]
}

fn cdf_line(label: &str, errs: &[f64]) {
    let mut v = errs.to_vec();
    v.sort_by(|a, b| a.total_cmp(b));
    println!(
        "{label:<22} n={:<5} p50={:.3} p90={:.3} p95={:.3} p99={:.3} max={:.3}",
        v.len(),
        percentile(&v, 0.50),
        percentile(&v, 0.90),
        percentile(&v, 0.95),
        percentile(&v, 0.99),
        percentile(&v, 0.99).max(v.last().copied().unwrap_or(0.0)),
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: eval_compare <a.txt> <b.txt> [label_a label_b]");
        std::process::exit(2);
    }
    let (la, lb) = (
        args.get(3).cloned().unwrap_or_else(|| "A".into()),
        args.get(4).cloned().unwrap_or_else(|| "B".into()),
    );
    let a = match load(&args[1]) {
        Ok(v) => v,
        Err(e) => { eprintln!("FAIL: {e}"); std::process::exit(1); }
    };
    let b = match load(&args[2]) {
        Ok(v) => v,
        Err(e) => { eprintln!("FAIL: {e}"); std::process::exit(1); }
    };

    println!("=== raw CDF (the salient result) ===");
    let ea: Vec<f64> = a.iter().map(|r| r.2).collect();
    let eb: Vec<f64> = b.iter().map(|r| r.2).collect();
    cdf_line(&la, &ea);
    cdf_line(&lb, &eb);

    println!("\n=== fix rate ===");
    let fa = a.iter().filter(|r| r.1 == 1).count();
    let fb = b.iter().filter(|r| r.1 == 1).count();
    println!("{la:<22} {}/{} ({:.1}%)", fa, a.len(), 100.0 * fa as f64 / a.len().max(1) as f64);
    println!("{lb:<22} {}/{} ({:.1}%)", fb, b.len(), 100.0 * fb as f64 / b.len().max(1) as f64);

    println!("\n=== Weibull tail fit ===");
    for (label, e) in [(&la, &ea), (&lb, &eb)] {
        match weibull_with_tail(e) {
            Ok(f) => println!(
                "{label:<22} shape k={:.3} scale={:.3} m  p99.9={:.3} m",
                f.shape,
                f.scale,
                f.p999.unwrap_or(f64::NAN)
            ),
            Err(e) => println!("{label:<22} fit failed: {e}"),
        }
    }

    // Paired test: align on epoch key, keep only shared epochs.
    let bm: std::collections::HashMap<u32, f64> = b.iter().map(|r| (r.0, r.2)).collect();
    let mut pa = Vec::new();
    let mut pb = Vec::new();
    for r in &a {
        if let Some(&other) = bm.get(&r.0) {
            pa.push(r.2);
            pb.push(other);
        }
    }

    println!("\n=== bootstrap CDF band (paired, 2000 reps, 95% pointwise) ===");
    println!("negative diff = {la} better; '*' marks a level where the interval excludes zero");
    println!("{:>5} {:>8} {:>8} {:>9} [{:>8}, {:>8}]", "level", "a", "b", "a-b", "lo", "hi");
    if let Ok(band) = bootstrap_cdf_band_paired(&pa, &pb, &cdf_levels(), 2000, 0.95, 0x5EED) {
        for pt in band {
            let verdict = if pt.a_better() { "a better *" }
                else if pt.b_better() { "b better *" } else { "no difference" };
            println!("{:>5.2} {:>8.3} {:>8.3} {:>9.3} [{:>8.3}, {:>8.3}]  {}",
                pt.level, pt.a, pt.b, pt.diff, pt.lo, pt.hi, verdict);
        }
    }


    println!("\n=== Wilcoxon signed-rank (paired, shared epochs, smaller is better) ===");
    println!("shared epochs: {} of {}", pa.len(), a.len());
    match wilcoxon_signed_rank(&pa, &pb) {
        Ok(r) => {
            println!(
                "{la} vs {lb}: n={} (zero-diff dropped {})",
                r.n, r.n_zero_diff
            );
            println!("  W+={:.1}  W-={:.1}", r.w_plus, r.w_minus);
            println!("  two-sided p = {:.3e}", r.p_value);
            println!(
                "  P({la} better) = {:.4}   direction = {}",
                r.probability_of_superiority,
                r.direction()
            );
            let verdict = if r.p_value < 0.05 {
                if r.direction() > 0 { "significant: A better" } else { "significant: B better" }
            } else {
                "not significant"
            };
            println!("  verdict: {verdict}");
        }
        Err(e) => println!("  test failed: {e}"),
    }
}
