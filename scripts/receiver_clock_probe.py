#!/usr/bin/env python3
"""Inter-receiver clock probe: base-minus-rover pseudorange single differences.

For every common epoch/satellite, computes
    tau = ((P_rov - rho_rov) - (P_base - rho_base)) / c
which equals (clk_rov - clk_base) plus common-mode terms that cancel when we
take the per-epoch MEDIAN over satellites. The resulting series is the
relative receiver clock of the two streams, sampled every 30 s. Drift,
steps, or end-of-day excursions in this series are exactly the timing errors
a DD processor silently absorbs into position/ambiguities.
"""

import math
import re
import sys

C = 299792458.0
MU = 3.986005e14
OMEGA_E = 7.2921151467e-5
HALF_WEEK = 302400.0


def parse_nav(path):
    """RINEX2 nav -> {prn: [record dicts]}, records sorted by ToC.

    Uses fixed-width column parsing (teqc emits unspaced D-exponents).
    """
    eph = {}
    float_pat = re.compile(r"[-+]?\d*\.\d+(?:D|E)[-+]?\d+", re.I)

    def d(s):
        s = s.strip()
        return float(s.replace("D", "E").replace("d", "e")) if s else 0.0

    def row_floats(k):
        """All D/E-exponent floats on record line k, in order."""
        return [float(m.group(0).replace("D", "E"))
                for m in float_pat.finditer(lines[i + k])]

    with open(path) as f:
        lines = f.readlines()
    # locate END OF HEADER
    start = 0
    for i, line in enumerate(lines):
        if "END OF HEADER" in line:
            start = i + 1
            break
    i = start
    while i + 8 <= len(lines):
        line = lines[i]
        try:
            prn = int(line[0:3])
            yy = int(line[3:5]); mon = int(line[5:8]); dd = int(line[8:11])
            hh = int(line[11:14]); mi = int(line[14:17])
            ss = float(line[17:22].strip() or "0")
            if not (1980 <= 2000 + yy <= 2099 or yy >= 80):
                raise ValueError("bad year")
            af0 = d(line[22:41]); af1 = d(line[41:60]); af2 = d(line[60:79])

            def row4(k):
                return row_floats(k)

            g = []
            for k in range(1, 6):
                g.extend(row4(k))
            # Week-aligned ToC: nav files carry records from adjacent days
            # (e.g. day-134 evening in a day-135 file). Sorting raw
            # seconds-of-day interleaves those with today's and corrupts
            # broadcast selection. Anchor to day-of-month relative to the
            # file's primary day (mode of all dd values).
            rec = {
                "toc": hh * 3600 + mi * 60 + ss,
                "dd": dd,
                "af0": af0, "af1": af1, "af2": af2,
                "iode": g[0], "crs": g[1], "dn": g[2], "m0": g[3],
                "cuc": g[4], "e": g[5], "cus": g[6], "sqrta": g[7],
                "toe": g[8], "cic": g[9], "om0": g[10],
                "cis": g[11], "i0": g[12], "crc": g[13], "w": g[14],
                "omdot": g[15], "idot": g[16],
                "tgd": g[19] if len(g) > 19 else 0.0,
            }
            eph.setdefault(prn, []).append(rec)
        except Exception:
            pass
        i += 8
    # align ToC to week time using modal day-of-month as the reference day
    from collections import Counter
    for prn in eph:
        days = Counter(r["dd"] for r in eph[prn])
        ref_day = days.most_common(1)[0][0]
        for r in eph[prn]:
            r["wtoc"] = (r["dd"] - ref_day) * 86400.0 + r["toc"]
        eph[prn].sort(key=lambda r: r["wtoc"])
    return eph


def select_eph(eph_prn, tow):
    """Broadcast selection: latest week-aligned ToC <= t, else nearest.

    tow is seconds-of-day on the reference day; records carry wtoc already
    day-offset relative to it.
    """
    best = None
    for r in eph_prn:
        if r["wtoc"] <= tow and (best is None or r["wtoc"] > best["wtoc"]):
            best = r
    if best is None:
        best = min(eph_prn, key=lambda r: abs(r["wtoc"] - tow))
    return best


def sv_position(rec, tow, rx_hint):
    """Return (pos_ecef_m, sv_clock_s_at_tx, t_tx).

    Standard two-step: propagate to transmit time t_tx = t_recv - P/c + dt_sv,
    iterating because P depends on the satellite position at t_tx.
    """
    t_tx = tow
    pos = None
    svdt = 0.0
    for _ in range(3):
        pos, svdt_raw = _kepler(rec, t_tx)
        # rough range from any nearby point (receiver hint) to SV
        rng = math.dist(rx_hint, pos)
        dt = t_tx - rec["toc"]
        if dt > HALF_WEEK:
            dt -= 2 * HALF_WEEK
        elif dt < -HALF_WEEK:
            dt += 2 * HALF_WEEK
        svdt = svdt_raw
        t_tx = tow - rng / C + svdt
        # wrap into week
        if t_tx < rec["toc"] - HALF_WEEK:
            t_tx += 2 * HALF_WEEK
    return pos, svdt, t_tx


def _kepler(rec, t_eval):
    """ECEF position + sv clock at evaluation time t_eval (seconds of week)."""
    dt = t_eval - rec["toe"]
    if dt > HALF_WEEK:
        dt -= 2 * HALF_WEEK
    elif dt < -HALF_WEEK:
        dt += 2 * HALF_WEEK
    a = rec["sqrta"] ** 2
    n0 = math.sqrt(MU / (a ** 3))
    n = n0 + rec["dn"]
    M = rec["m0"] + n * dt
    E = M
    for _ in range(12):
        E -= (E - rec["e"] * math.sin(E) - M) / (1 - rec["e"] * math.cos(E))
    sinE, cosE = math.sin(E), math.cos(E)
    v = math.atan2(math.sqrt(1 - rec["e"] ** 2) * sinE, cosE - rec["e"])
    phi = v + rec["w"]
    du = rec["cuc"] * math.cos(2 * phi) + rec["cus"] * math.sin(2 * phi)
    di = rec["cic"] * math.cos(2 * phi) + rec["cis"] * math.sin(2 * phi)
    r = a * (1 - rec["e"] * cosE)
    phic = phi + du
    inc = rec["i0"] + di + rec["idot"] * dt
    xp = r * math.cos(phic)
    yp = r * math.sin(phic)
    om = (rec["om0"]
          + rec["omdot"] * dt
          - OMEGA_E * (rec["toe"] + dt))
    coso, sino = math.cos(om), math.sin(om)
    cosi, sini = math.cos(inc), math.sin(inc)
    pos = (
        xp * coso - yp * cosi * sino,
        xp * sino + yp * cosi * coso,
        yp * sini,
    )
    svdt = rec["af0"] + rec["af1"] * dt + rec["af2"] * dt * dt \
        - rec["tgd"]
    return pos, svdt


EPOCH_PAT = re.compile(
    r"^\s*(\d{2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+([\d.]+)")
SAT_TOK = re.compile(r"[GRECJS]\s*\d+")


def _parse_types(lines):
    """Return the ordered observables list from RINEX2 header lines.

    Handles wrapped '# / TYPES OF OBSERV' blocks (e.g. 20 types over 3
    lines): n comes from the FIRST line of the block, tokens from ALL of
    them concatenated.
    """
    idxs = [k for k, l in enumerate(lines[:80])
            if "# / TYPES OF OBSERV" in l]
    if not idxs:
        return ["L1", "L2", "C1", "P1", "P2", "S1", "S2"]
    # keep only the contiguous block containing the first marker line
    block = [idxs[0]]
    for k in idxs[1:]:
        if k == block[-1] + 1:
            block.append(k)
        else:
            break
    try:
        n = int(lines[block[0]][:6])
    except ValueError:
        return []
    text = "".join(lines[k] for k in block)
    found = re.findall(r"\b[LCSP]\d[A-Za-z]?\b", text)
    return found[:n] if len(found) >= n else found


def _obs_value(cell):
    """RINEX2 observation cell -> (value_or_None).

    Layout: value right-justified in cols 0-13, LLI col 14, SSI col 15.
    teqc glues digits when fields are full ('95416356.99644' == value .996,
    LLI '4', SSI '4'), so slicing [0:14] is mandatory; float() on the whole
    cell corrupts values or throws.
    """
    head = cell[0:14]
    if not head.strip():
        return None
    try:
        return float(head)
    except ValueError:
        return None


def parse_obs(path):
    """-> {tow: {(sys, num): {'C1': m, ...}}} for ALL constellations listed."""
    epochs = {}
    with open(path) as f:
        lines = f.readlines()
    types = _parse_types(lines)
    ntypes = len(types)
    i = next((k for k, l in enumerate(lines) if "END OF HEADER" in l), -1) + 1
    while i < len(lines):
        m = EPOCH_PAT.match(lines[i])
        if not m:
            i += 1
            continue
        hh, mi = int(m.group(4)), int(m.group(5))
        ss = float(m.group(6))
        tow = hh * 3600 + mi * 60 + round(ss)
        # flag at col 28, nsat right-justified in cols 29-31 (I3). The old
        # single-digit regex group grabbed '2' of '23' and silently broke
        # continuation handling.
        try:
            nsat = int(lines[i][29:32])
        except ValueError:
            i += 1
            continue
        toks = SAT_TOK.findall(lines[i][30:])
        j = i + 1
        while len(toks) < nsat and j < len(lines):
            if EPOCH_PAT.match(lines[j]):
                break
            toks.extend(SAT_TOK.findall(lines[j]))
            j += 1
        obs = {}
        idx = j
        ok_epoch = True
        for tok in toks:
            sysc, num = tok[0], int(tok[1:].strip() or "0")
            vals = []
            while len(vals) < ntypes and idx < len(lines):
                row = lines[idx].rstrip("\n")
                for c in range(0, 80, 16):
                    vals.append(_obs_value(row[c:c + 16]))
                idx += 1
            if len(vals) < ntypes:
                ok_epoch = False
                break
            obs[(sysc, num)] = dict(zip(types, vals))
        if ok_epoch and nsat > 0:
            epochs[tow] = obs
            i = idx
        else:
            i = j  # resync at first data row on malformed epoch
    return epochs


def approx_pos(path):
    with open(path) as f:
        for line in f:
            if "APPROX POSITION XYZ" in line:
                return tuple(float(line[c:c + 14]) for c in (0, 14, 28))
    return None


def main():
    rov_path, bas_path = sys.argv[1], sys.argv[2]
    nav_path = sys.argv[3]
    eph = parse_nav(nav_path)
    print(f"nav: {sum(len(v) for v in eph.values())} records, "
          f"{len(eph)} PRNs", file=sys.stderr)

    rov = parse_obs(rov_path)
    bas = parse_obs(bas_path)
    common = sorted(set(rov) & set(bas))
    print(f"epochs: rov={len(rov)} base={len(bas)} common={len(common)}",
          file=sys.stderr)

    def approx(path):
        with open(path) as f:
            for line in f:
                if "APPROX POSITION XYZ" in line:
                    return tuple(float(line[c:c + 14]) for c in (0, 14, 28))
        return None

    rx = approx(rov_path)
    bx = approx(bas_path)

    print("#tow  relclock_us  nsv  spread_us")
    series = []
    for tow in common:
        taus = []
        for key, entry in rov[tow].items():
            if key not in bas[tow]:
                continue
            c_r = entry.get("C1") or entry.get("P1")
            c_b = bas[tow][key].get("C1") or bas[tow][key].get("P1")
            if not c_r or not c_b:
                continue
            prn = key[1]
            if prn not in eph:
                continue
            rec = select_eph(eph[prn], tow)
            try:
                sp, svdt, t_tx = sv_position(rec, tow, rx)
            except Exception:
                continue
            rho = lambda st, sv: math.dist(st, sv)
            # NOTE: the SV clock cancels exactly in this difference (both
            # receivers observe the same signal); do NOT add svdt here.
            tau = ((c_r - rho(rx, sp)) - (c_b - rho(bx, sp))) / C
            taus.append(tau)
        if len(taus) >= 4:
            taus.sort()
            med = taus[len(taus) // 2]
            spread = (taus[int(len(taus) * .9)] - taus[int(len(taus) * .1)]) * 1e6
            series.append((tow, med * 1e6, len(taus), spread))

    for tow, us, n, spread in series:
        print(f"{tow} {us:.3f} {n} {spread:.3f}")


if __name__ == "__main__":
    main()
