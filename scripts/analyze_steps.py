#!/usr/bin/env python3
"""Level-shift ("wrong fix") episode detector for gneiss PPK dump CSVs.

Prototype for branch network-rtk-long-baseline. Detects sustained level
shifts in *proxy* channels — never in the truth-referenced h/v columns.

Proxy channels per epoch (causal, trailing reference window only):
  sep    : Smoothed separation_3d = ||fwd_pos - bwd_pos|| at fuse time
           (post_process/combiner.rs combine_bidirectional_epoch).
  dpos   : cross-pass disagreement hypot(h_fwd-h_smt, v_fwd-v_smt) — h/v used
           ONLY as a between-pass difference, never as absolute error.
  qtrans : binary fix-quality transition (either pass), small score bonus.

Detector: B_MAX-centered robust z-scores -> one-sided CUSUM (drift k, alarm H,
hysteresis release H/4, min duration, gap merge, hard length cap).

Ground truth (evaluation only, never a feature): maximal runs of >=gt_min
consecutive Smoothed epochs with |v| > gt_mm (spec: 100 mm, 5 epochs); stable
stretches: runs of >=20 consecutive epochs with |v| < 50 mm.

Usage: python3 scripts/analyze_steps.py [--dir /tmp/steps] [--selftest]
      [--out-episodes PATH] [--thresholds 9,12,18,25]
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import os
import sys
from dataclasses import dataclass, field

DEFAULT_BASES = ["P181", "OHLN", "CAPO", "P225", "P222", "SLAC"]
DEFAULT_THRESHOLDS = [9.0, 12.0, 18.0, 25.0]
B_MAX = 0.5  # null mean of max() of two unit-scale channels; centers the score

@dataclass
class Row:
    tow: float
    h_mm: float
    v_mm: float
    q: int
    sep: float
    nsat: int

@dataclass
class Epoch:
    """Aligned smoothed+forward epoch. Features: sep, dpos, qtrans ONLY;
    v_mm/h_mm are truth references kept for evaluation, never detection."""
    tow: float
    q_s: int
    q_f: int
    sep: float
    dpos_mm: float
    qtrans: bool
    v_mm: float
    h_mm: float

@dataclass
class Episode:
    start_i: int
    end_i: int  # exclusive
    start_tow: float
    end_tow: float
    peak: float  # peak CUSUM statistic (sigma*epochs)
    sep_shift_m: float
    dpos_shift_mm: float
    qtrans_frac: float
    matched_gt: int = -1  # max-overlap GT index; -1 = detection w/o GT overlap

    @property
    def n(self) -> int:
        return self.end_i - self.start_i

@dataclass
class BaseResult:
    label: str
    n_epochs: int
    gt: list = field(default_factory=list)       # [(i0, i1)] wrong-fix runs
    stable: list = field(default_factory=list)   # [(i0, i1)] stable stretches
    episodes: list = field(default_factory=list)
    tp_gt: set = field(default_factory=set)
    fp_stable: set = field(default_factory=set)
    latencies: list = field(default_factory=list)  # det_start - gt_start
    auc_sep: float = float("nan")  # channel AUC for GT membership, epoch-level
    auc_dpos: float = float("nan")

def load_pass(path: str) -> list[Row]:
    rows = []
    with open(path, newline="") as f:
        for rec in csv.DictReader(f):
            rows.append(Row(float(rec["tow"]), float(rec["h"]) * 1e3,
                            float(rec["v"]) * 1e3, int(rec["q"]),
                            float(rec["sep"]), int(rec["nsat"])))
    return rows

def align_passes(sm: list[Row], fw: list[Row]) -> list[Epoch]:
    """Join passes on integer TOW; build proxy features + eval-only truth."""
    fmap = {int(r.tow): r for r in fw}
    out, prev_s_q, prev_f_q = [], None, None
    for s in sm:
        f = fmap.get(int(s.tow))
        if f is None:  # pass gap resets previous-q memory
            prev_s_q = prev_f_q = None
            continue
        dpos = math.hypot(f.h_mm - s.h_mm, f.v_mm - s.v_mm)
        qt = ((prev_s_q is not None and s.q != prev_s_q)
              or (prev_f_q is not None and f.q != prev_f_q))
        out.append(Epoch(s.tow, s.q, f.q, s.sep, dpos, bool(qt), s.v_mm, s.h_mm))
        prev_s_q, prev_f_q = s.q, f.q
    return out

def runs_of(mask: list[bool], min_len: int) -> list[tuple[int, int]]:
    """Maximal [i0, i1) runs of True with length >= min_len."""
    runs, i = [], 0
    while i < len(mask):
        j = i
        while j < len(mask) and mask[j]:
            j += 1
        if mask[i] and j - i >= min_len:
            runs.append((i, j))
        i = max(j, i + 1)
    return runs

def ground_truth(eps: list[Epoch], gt_mm: float, min_len: int) -> list[tuple[int, int]]:
    return runs_of([abs(e.v_mm) > gt_mm for e in eps], min_len)

def stable_stretches(eps: list[Epoch], st_mm: float, min_len: int) -> list[tuple[int, int]]:
    return runs_of([abs(e.v_mm) < st_mm for e in eps], min_len)

def _med(xs: list[float]) -> float:
    if not xs:
        return 0.0
    s = sorted(xs)
    m = len(s) // 2
    return s[m] if len(s) % 2 else 0.5 * (s[m - 1] + s[m])

def causal_z(xs: list[float], win: int, floor: float) -> list[float]:
    """Trailing robust z: (x - med[past]) / (1.4826*MAD[past] + floor).

    The window excludes the sample itself; until min(30, win//2) samples exist
    z is 0. Sustained shifts contaminate the reference only after ~win epochs.
    """
    zs, min_obs = [], min(30, win // 2)
    for i in range(len(xs)):
        ref = xs[max(0, i - win):i]
        zs.append(0.0 if len(ref) < min_obs
                  else (xs[i] - _med(ref)) / (1.4826 * _med([abs(v - _med(ref))
                                                            for v in ref]) + floor))
    return zs

def score_series(eps: list[Epoch], cfg) -> list[float]:
    """B_MAX-centered combined score; null increments then average -k."""
    z_sep = causal_z([e.sep for e in eps], cfg.window, cfg.sep_floor_m)
    z_dps = causal_z([e.dpos_mm for e in eps], cfg.window, cfg.dpos_floor_mm)
    return [max(a, b) + cfg.c_q * e.qtrans - B_MAX
            for a, b, e in zip(z_sep, z_dps, eps)]

def cusum(scores: list[float], k: float, h_alarm: float, h_rel: float,
          min_dur: int, merge_gap: int, max_dur: int) -> list[Episode]:
    """One-sided CUSUM; statistic resets after every close or max_dur cut so
    sustained elevation must re-earn alarm level."""
    segs, s, i = [], 0.0, 0
    while i < len(scores):
        s = max(0.0, s + scores[i] - k)
        if s < h_alarm:
            i += 1
            continue
        start, peak, j = i, s, i + 1
        while j < len(scores) and j - start < max_dur:
            s = max(0.0, s + scores[j] - k)
            peak = max(peak, s)
            if s < h_rel:
                break
            j += 1
        segs.append([start, j, peak])
        s, i = 0.0, j
    segs = [sg for sg in segs if sg[1] - sg[0] >= min_dur]
    return merge_segments(segs, merge_gap, max_dur)

def merge_segments(segs: list, gap: int, max_dur: int) -> list[Episode]:
    merged: list[list] = []
    for sg in segs:
        if merged and sg[0] - merged[-1][1] <= gap \
                and sg[1] - merged[-1][0] <= max_dur:
            merged[-1][1], merged[-1][2] = sg[1], max(merged[-1][2], sg[2])
        else:
            merged.append(list(sg))
    return [Episode(a, b, 0.0, 0.0, p, 0.0, 0.0, 0.0) for a, b, p in merged]

def annotate(eps: list[Epoch], det: list[Episode], cfg) -> None:
    """Fill tow bounds + native-unit shift magnitudes for each episode."""
    seps, dpss = [e.sep for e in eps], [e.dpos_mm for e in eps]
    for ep in det:
        lo = max(0, ep.start_i - cfg.window)
        ep.start_tow, ep.end_tow = eps[ep.start_i].tow, eps[ep.end_i - 1].tow
        ep.sep_shift_m = _med(seps[ep.start_i:ep.end_i]) - _med(seps[lo:ep.start_i])
        ep.dpos_shift_mm = _med(dpss[ep.start_i:ep.end_i]) - _med(dpss[lo:ep.start_i])
        ep.qtrans_frac = sum(e.qtrans for e in eps[ep.start_i:ep.end_i]) / ep.n

def overlap(a: tuple[int, int], b: tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def epoch_auc(scores: list[float], gt_mask: list[bool]) -> float:
    """Mann-Whitney AUC vs GT membership; 0.5 = uninformative."""
    pos = sorted(s for s, m in zip(scores, gt_mask) if m)
    neg = sorted(s for s, m in zip(scores, gt_mask) if not m)
    if not pos or not neg:
        return float("nan")
    win = sum(bisect.bisect_right(neg, p) for p in pos)
    tie = sum(bisect.bisect_left(neg, p) - bisect.bisect_right(neg, p)
              for p in pos)
    return (win + 0.5 * tie) / (len(pos) * len(neg))

def evaluate(label: str, res: BaseResult) -> BaseResult:
    """TP: some detection covers >=50% of the GT episode (latency: first
    touch); stable-stretch FP: any overlap; detections get best GT match."""
    spans = [(ep.start_i, ep.end_i) for ep in res.episodes]
    for gi, gt in enumerate(res.gt):
        need = max(1, (gt[1] - gt[0]) // 2)
        touches = [s for s in spans if overlap(gt, s) >= need]
        if touches:
            res.tp_gt.add(gi)
            res.latencies.append(min(s for s, _ in touches) - gt[0])
    for si, stch in enumerate(res.stable):
        if any(overlap(stch, s) > 0 for s in spans):
            res.fp_stable.add(si)
    for ep, span in zip(res.episodes, spans):
        ovs = [overlap(gt, span) for gt in res.gt]
        ep.matched_gt = max(range(len(ovs)), key=ovs.__getitem__) if any(ovs) else -1
    res.label = label
    return res

def analyze(paths: tuple[str, str], label: str, cfg) -> BaseResult:
    eps = align_passes(load_pass(paths[0]), load_pass(paths[1]))
    res = BaseResult(label, len(eps))
    res.gt = ground_truth(eps, cfg.gt_mm, cfg.gt_min)
    res.stable = stable_stretches(eps, cfg.stable_mm, cfg.stable_min)
    gt_mask = [any(a <= i < b for a, b in res.gt) for i in range(len(eps))]
    res.auc_sep = epoch_auc(causal_z([e.sep for e in eps], cfg.window,
                                     cfg.sep_floor_m), gt_mask)
    res.auc_dpos = epoch_auc(causal_z([e.dpos_mm for e in eps], cfg.window,
                                      cfg.dpos_floor_mm), gt_mask)
    res.episodes = cusum(score_series(eps, cfg), cfg.k, cfg.threshold,
                         cfg.threshold / 4.0, cfg.min_dur, cfg.merge_gap,
                         cfg.max_dur)
    annotate(eps, res.episodes, cfg)
    return evaluate(label, res)

def fmt_row(cells: list) -> str:
    return "| " + " | ".join(str(c) for c in cells) + " |"

def lat_str(r: BaseResult) -> str:
    if not r.latencies:
        return "-"
    return f"{sum(r.latencies) / len(r.latencies):+.1f}/{_med(r.latencies):+.1f}"

def bucket(n: int) -> str:
    return "5-9" if n < 10 else ("10-29" if n < 30 else "30+")

def bucket_tpr(r: BaseResult) -> dict[str, float]:
    buckets: dict[str, list[bool]] = {}
    for gi, gt in enumerate(r.gt):
        buckets.setdefault(bucket(gt[1] - gt[0]), []).append(gi in r.tp_gt)
    return {b: sum(v) / len(v) for b, v in buckets.items()}

def summarize(results: list[BaseResult]) -> dict:
    n_gt = sum(len(r.gt) for r in results)
    n_tp = sum(len(r.tp_gt) for r in results)
    n_st, n_fp = sum(len(r.stable) for r in results), sum(len(r.fp_stable) for r in results)
    lat = [x for r in results for x in r.latencies]
    hits = sum(sum(1 for ep in r.episodes if ep.matched_gt >= 0) for r in results)
    n_det = sum(len(r.episodes) for r in results)
    return {"tpr": n_tp / max(1, n_gt), "fpr": n_fp / max(1, n_st),
            "n_gt": n_gt, "n_tp": n_tp, "n_st": n_st, "n_fp": n_fp,
            "n_det": n_det, "prec": hits / n_det if n_det else float("nan"),
            "lat_mean": sum(lat) / len(lat) if lat else float("nan"),
            "lat_med": _med(lat) if lat else float("nan")}

def print_base_table(results: list[BaseResult], thr: float, cfg) -> None:
    print(f"\n## Per-base performance @ H={thr:g}, k={cfg.k:g}, W={cfg.window} "
          f"(GT |v|>{cfg.gt_mm:.0f}mm >={cfg.gt_min}ep; "
          f"stable |v|<{cfg.stable_mm:.0f}mm >={cfg.stable_min}ep)\n")
    hdr = ["base", "epochs", "GT_ep", "det_ep", "TPR", "FP_st", "FPR",
           "lat_mean/med", "TPR 5-9", "TPR 10-29", "TPR 30+",
           "AUC_sep", "AUC_dpos", "det_prec"]
    print(fmt_row(hdr))
    print(fmt_row(["---"] * len(hdr)))
    nan = float("nan")
    for r in results:
        bt = bucket_tpr(r)
        print(fmt_row([r.label, r.n_epochs, len(r.gt), len(r.episodes),
                       f"{len(r.tp_gt)}/{len(r.gt)}", len(r.fp_stable),
                       f"{len(r.fp_stable)}/{len(r.stable)}", lat_str(r),
                       f"{bt.get('5-9', nan):.2f}", f"{bt.get('10-29', nan):.2f}",
                       f"{bt.get('30+', nan):.2f}", f"{r.auc_sep:.2f}",
                       f"{r.auc_dpos:.2f}",
                       f"{summarize([r])['prec']:.2f}"]))

def print_sensitivity(by_thr: dict, calib_label: str) -> None:
    print("\n## Threshold sensitivity (pooled over all bases)\n")
    hdr = ["H", "TPR", "FP_st/FPR", "n_det", "prec", "lat_mean", "lat_med"]
    print(fmt_row(hdr))
    print(fmt_row(["---"] * len(hdr)))
    for thr in sorted(by_thr):
        agg, per = by_thr[thr]
        nan = float("nan")
        print(fmt_row([f"{thr:g}", f"{agg['n_tp']}/{agg['n_gt']}={agg['tpr']:.2f}",
                       f"{agg['n_fp']}/{agg['n_st']}={agg['fpr']:.2f}", agg["n_det"],
                       f"{agg['prec']:.2f}", f"{agg['lat_mean']:+.1f}",
                       f"{agg['lat_med']:+.1f}"]))
        for lbl, r in per.items():
            s = summarize([r])
            mark = "*" if lbl == calib_label else " "
            lat = f"{s['lat_mean']:+.1f}" if s["lat_mean"] == s["lat_mean"] else "-"
            print(fmt_row([mark + " " + lbl,
                           f"{s['n_tp']}/{s['n_gt']}={s['tpr']:.2f}",
                           f"{s['n_fp']}/{s['n_st']}={s['fpr']:.2f}", s["n_det"],
                           f"{s['prec']:.2f}", lat, "-"]))

def recommend(by_thr: dict, calib_label: str) -> float:
    """Pick H maximizing TPR-FPR on the calibration base; frozen afterwards."""

    def margin(item):
        _, (_, per) = item
        s = summarize([per[calib_label]])
        return s["tpr"] - s["fpr"]

    return max(by_thr.items(), key=margin)[0]

def print_header(args, thr_rec: float) -> None:
    print("# Step-detector report — proxy: sep (fwd/bwd separation at fuse) "
          "+ fwd-vs-smoothed dpos + q transitions\n")
    print(f"Detector: causal robust-z (trailing W={args.window}) on "
          f"max(z_sep, z_dpos) + {args.c_q:g}*qtrans; CUSUM k={args.k:g}, "
          f"alarm H, release H/4, min_dur={args.min_dur}, "
          f"merge_gap={args.merge_gap}, max_dur={args.max_dur}.")
    print("h/v enter ONLY as cross-pass differences; ground truth below is "
          f"evaluation-only.\nRecommended threshold (max TPR-FPR on "
          f"{args.calib}): H={thr_rec:g}")

def write_episodes(per: dict, path: str) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["base", "start_tow", "end_tow", "n_epochs", "peak_sigma",
                    "sep_shift_m", "dpos_shift_mm", "qtrans_frac", "gt_hit"])
        for lbl, r in per.items():
            for ep in r.episodes:
                w.writerow([lbl, f"{ep.start_tow:.0f}", f"{ep.end_tow:.0f}",
                            ep.n, f"{ep.peak:.1f}", f"{ep.sep_shift_m:+.3f}",
                            f"{ep.dpos_shift_mm:+.1f}", f"{ep.qtrans_frac:.2f}",
                            ep.matched_gt])

def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dir", default="/tmp/steps")
    p.add_argument("--bases", default=",".join(DEFAULT_BASES))
    p.add_argument("--thresholds", default=",".join(map(str, DEFAULT_THRESHOLDS)),
                   help="comma-separated CUSUM alarm levels H")
    p.add_argument("--k", type=float, default=2.0, help="CUSUM drift allowance")
    p.add_argument("--window", type=int, default=180, help="trailing ref window")
    p.add_argument("--c-q", dest="c_q", type=float, default=0.5)
    p.add_argument("--min-dur", dest="min_dur", type=int, default=5)
    p.add_argument("--merge-gap", dest="merge_gap", type=int, default=10)
    p.add_argument("--max-dur", dest="max_dur", type=int, default=180)
    p.add_argument("--gt-mm", dest="gt_mm", type=float, default=100.0)
    p.add_argument("--gt-min", dest="gt_min", type=int, default=5)
    p.add_argument("--stable-mm", dest="stable_mm", type=float, default=50.0)
    p.add_argument("--stable-min", dest="stable_min", type=int, default=20)
    p.add_argument("--sep-floor-m", dest="sep_floor_m", type=float, default=0.005)
    p.add_argument("--dpos-floor-mm", dest="dpos_floor_mm", type=float, default=1.0)
    p.add_argument("--calib", default="P181")
    p.add_argument("--out-episodes", dest="out_episodes", default=None)
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args(argv)
    args.thresholds = [float(x) for x in args.thresholds.split(",")]
    return args

def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.selftest:
        return run_selftest()
    bases = args.bases.split(",")
    paths = {b: tuple(os.path.join(args.dir, f"dump_{b}_{p}.csv")
                      for p in ("Smoothed", "Forward")) for b in bases}
    missing = [p for ps in paths.values() for p in ps if not os.path.exists(p)]
    if missing:
        sys.exit(f"missing dumps: {missing}")
    by_thr: dict[float, tuple] = {}
    for thr in args.thresholds:
        args.threshold = thr
        per = {b: analyze(paths[b], b, args) for b in bases}
        by_thr[thr] = (summarize(list(per.values())), per)
    thr_rec = recommend(by_thr, args.calib)
    print_header(args, thr_rec)
    print_base_table(list(by_thr[thr_rec][1].values()), thr_rec, args)
    print_sensitivity(by_thr, args.calib)
    if args.out_episodes:
        write_episodes(by_thr[thr_rec][1], args.out_episodes)
        print(f"\nEpisodes written: {args.out_episodes}")
    return 0

# ---------------------------------------------------------------- selftests

def _mk_cfg(threshold: float = 9.0):
    import types
    return types.SimpleNamespace(
        window=60, k=1.0, threshold=threshold, c_q=0.5,
        min_dur=5, merge_gap=10, max_dur=180,
        gt_mm=100.0, gt_min=5, stable_mm=50.0, stable_min=20,
        sep_floor_m=0.005, dpos_floor_mm=1.0)

def _synth_eps(n: int, step_at: int = -1, step_amp: float = 0.0,
               seed: int = 7) -> list[Epoch]:
    import random
    rng = random.Random(seed)
    out = []
    for i in range(n):
        bad = step_at >= 0 and i >= step_at
        sep = max(0.07 + (step_amp if bad else 0.0) + rng.gauss(0, 0.01), 0.0)
        v = 30.0 + (150.0 if bad else 0.0) + rng.gauss(0, 20.0)
        out.append(Epoch(86400.0 + 30 * i, 1, 1, sep, abs(rng.gauss(18, 6)),
                         False, v, 25.0))
    return out

def test_causal_z_flags_shift():
    xs = [0.07] * 100 + [0.30] * 60
    zs = causal_z(xs, 180, 0.005)
    assert max(zs[:100]) < 2.0, "null region must stay quiet"
    assert max(zs[105:]) > 4.0, "sustained shift must exceed alarm scale"

def test_causal_z_constant_series_safe():
    zs = causal_z([0.5] * 200, 60, 0.005)
    assert all(abs(z) < 1e-9 for z in zs), "zero MAD must not explode"

def _det(cfg, eps):
    return cusum(score_series(eps, cfg), cfg.k, cfg.threshold,
                 cfg.threshold / 4, cfg.min_dur, cfg.merge_gap, cfg.max_dur)

def test_cusum_detect_latency_and_cap():
    cfg = _mk_cfg()
    det = _det(cfg, _synth_eps(600, step_at=100, step_amp=0.15))
    assert det, "sustained shift must produce an episode"
    assert det[0].start_i <= 110, f"latency too high: {det[0].start_i - 100}"
    assert det[0].n >= 50, f"episode must persist: n={det[0].n}"
    assert all(e.n <= cfg.max_dur for e in det), "max_dur cap violated"

def test_cusum_quiet_series_no_false_alarm():
    assert _det(_mk_cfg(), _synth_eps(600)) == [], \
        "quiet series must yield no episodes"

def test_ground_truth_and_stable_runs():
    eps = [Epoch(30.0 * i, 1, 1, 0.1, 10.0, False,
                 120.0 if 5 <= i < 12 else 30.0, 20.0) for i in range(40)]
    assert ground_truth(eps, 100.0, 5) == [(5, 12)]
    assert stable_stretches(eps, 50.0, 20) == [(12, 40)]

def test_align_skips_unmatched_tow():
    sm = [Row(86400.0 + 30 * i, 10.0, 5.0, 1, 0.1, 8) for i in range(3)]
    eps = align_passes(sm, [sm[0], sm[2]])  # middle tow missing forward
    assert len(eps) == 2 and eps[0].tow == sm[0].tow
    assert eps[0].qtrans is False, "gap resets previous-q memory"

def run_selftest() -> int:
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    for fn in tests:
        fn()
        print(f"ok  {fn.__name__}")
    print(f"{len(tests)} selftests passed")
    return 0

if __name__ == "__main__":
    sys.exit(main())
