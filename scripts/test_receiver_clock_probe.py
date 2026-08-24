#!/usr/bin/env python3
"""Red-green tests for receiver_clock_probe parsing.

Ground truth values below were transcribed by hand from raw file bytes
(column-ruler dump of datasets/cors_short_baseline/p2241350.20o lines
46-49 and ohln1350.20o first epoch). If one of these fails, the parser
is wrong -- not the expectation.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from receiver_clock_probe import parse_nav, parse_obs  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                    "datasets", "cors_short_baseline")
ROV = os.path.join(DATA, "p2241350.20o")
BAS = os.path.join(DATA, "ohln1350.20o")
NAV = os.path.join(DATA, "brdc1350.20n")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def approx(a, b, tol):
    return a is not None and abs(a - b) <= tol


def main():
    rov = parse_obs(ROV)
    bas = parse_obs(BAS)

    # --- observables lists -------------------------------------------------
    # Rover header declares 20 types across 3 wrapped lines.
    # We expose them via a module-level helper for introspection.
    import receiver_clock_probe as p

    def read_types(path):
        with open(path) as f:
            return p._parse_types(f.readlines())

    check("rover types exact",
          read_types(ROV) == ["L1", "L2", "C1", "P2", "P1", "S1", "S2",
                              "C2", "L5", "C5", "S5", "L6", "C6", "S6",
                              "L7", "C7", "S7", "L8", "C8", "S8"],
          str(read_types(ROV)))
    check("base types exact",
          read_types(BAS) == ["L1", "L2", "C1", "P1", "P2", "S1", "S2"],
          str(read_types(BAS)))

    # --- epoch grid --------------------------------------------------------
    check("rover has 2880 epochs", len(rov) == 2880, str(len(rov)))
    check("base has 2880 epochs", len(bas) == 2880, str(len(bas)))
    grid_ok = all(t in rov for t in range(0, 86400, 30))
    check("rover epochs land exactly on 30s grid", grid_ok)

    # --- first-epoch satellite list ---------------------------------------
    e0 = rov[0]
    check("first epoch has 23 satellites", len(e0) == 23, str(len(e0)))
    check("sat list includes GLONASS R11", ("R", 11) in e0)
    check("sat list includes Galileo E33", ("E", 33) in e0)
    check("first listed SV is G21", ("G", 21) in e0)

    # --- ground-truth observation values (hand-transcribed bytes) ----------
    g21 = e0[("G", 21)]
    check("G21 L1 phase == 122451126.087",
          approx(g21.get("L1"), 122451126.087, 0.001), str(g21.get("L1")))
    check("G21 L2 phase == 95416356.996 (LLI=4 SSI=4 glued after)",
          approx(g21.get("L2"), 95416356.996, 0.001), str(g21.get("L2")))
    check("G21 C1 code == 23301641.148",
          approx(g21.get("C1"), 23301641.148, 0.001), str(g21.get("C1")))
    check("G21 P2 code == 23301647.195",
          approx(g21.get("P2"), 23301647.195, 0.001), str(g21.get("P2")))
    check("G21 P1 absent (blank cell)", g21.get("P1") is None,
          str(g21.get("P1")))
    check("G21 S1 snr == 43.300",
          approx(g21.get("S1"), 43.300, 0.01), str(g21.get("S1")))
    check("G21 S2 snr == 26.000",
          approx(g21.get("S2"), 26.000, 0.01), str(g21.get("S2")))

    b21 = bas[0][("G", 21)]
    # Base epoch-0 row (transcribed): L1=122457018.108 L2=95421056.079
    # C1=23302794.525 P2=23302790.600
    check("base G21 C1 == 23302794.525",
          approx(b21.get("C1"), 23302794.525, 0.001), str(b21.get("C1")))
    check("base G21 P2 == 23302790.600",
          approx(b21.get("P2"), 23302790.600, 0.001), str(b21.get("P2")))

    # --- nav ----------------------------------------------------------------
    eph = parse_nav(NAV)
    check("nav has 32 PRNs", len(eph) == 32, str(len(eph)))
    nrec = sum(len(v) for v in eph.values())
    check("nav record count plausible (>=470)", nrec >= 470, str(nrec))
    rec1 = eph[1][0]
    check("nav PRN1 earliest af0 == -3.792415373027e-4",
          approx(rec1["af0"], -3.792415373027e-4, 1e-12), str(rec1["af0"]))
    check("nav PRN1 earliest toe sane (0 < toe < week)",
          0 < rec1["toe"] < 604800, str(rec1["toe"]))

    # --- physics gate: per-SV clock-term spread must be microsecond-scale --
    import math
    from receiver_clock_probe import select_eph, sv_position, approx_pos

    rx = approx_pos(ROV)
    bx = approx_pos(BAS)
    C = 299792458.0
    worst = 0.0
    checked = 0
    for tow in (21600, 43200, 64800):
        taus = []
        for key in sorted(set(rov[tow]) & set(bas[tow])):
            if key[0] != "G":
                continue
            c_r = rov[tow][key].get("C1")
            c_b = bas[tow][key].get("C1")
            prn = key[1]
            if not c_r or not c_b or prn not in eph:
                continue
            sp, _, _ = sv_position(select_eph(eph[prn], tow), tow, rx)
            rho_r = math.dist(rx, sp)
            rho_b = math.dist(bx, sp)
            taus.append(((c_r - rho_r) - (c_b - rho_b)) / C * 1e6)
        if len(taus) >= 4:
            checked += 1
            spread = max(taus) - min(taus)
            worst = max(worst, spread)
    check(f"per-SV clock-term spread < 500us (worst={worst:.1f}us, "
          f"epochs checked={checked})",
          checked == 3 and worst < 500.0)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        return 1
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
