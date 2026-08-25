"""Fetch multi-GNSS precise orbit/clock (SP3 + RINEX clock) products, no auth.

Verified working open sources (probed live, all anonymous HTTPS):

  FINAL multi-GNSS (latency ~9-13 days):
    AIUB (CODE home server), per calendar year directory:
      https://www.aiub.unibe.ch/download/CODE_MGEX/CODE/{YYYY}/
        COD0MGXFIN_{YYYYDDD}0000_01D_05M_ORB.SP3.gz   (5 min orbits, SP3-d)
        COD0MGXFIN_{YYYYDDD}0000_01D_30S_CLK.CLK.gz   (30 s clocks, RINEX 3)
      Constellations: GPS+GLO+GAL+BDS+QZSS. Best quality for GPS+Galileo PPK.
      Also hosts _OSB.BIA code biases and _ERP.

  RAPID multi-GNSS (latency ~17 h):
    SIO Garner data center, per GPS-week directory:
      https://garner.ucsd.edu/pub/products/{gpsweek}/
        GFZ0MGXRAP_{YYYYDDD}0000_01D_05M_ORB.SP3.gz
        GFZ0MGXRAP_{YYYYDDD}0000_01D_30S_CLK.CLK.gz   (+ OSB.BIA)
      Same week dirs also carry IAX0MGXFIN-style AC finals (JAX0MGXFIN,
      IAC0MGXFIN) and GPS-only combined IGS0OPS* streams.

  FALLBACKS on the same garner week dirs: JAX0MGXFIN (JAXA MADOCA final,
  verified multi-GNSS C/E/G/J/R clocks) and IAC0MGXFIN finals, plus
  GPS-only combined IGS0OPSFIN/RAP (also mirrored by BKG under
  root_ftp/IGS/products/{gpsweek}/; NOTE the IGS "OPS" stream is
  GPS+GLONASS only - no Galileo).

Known-dead or auth-walled (do not rely on): files.igs.org/pub/product
(products removed, see its readme.txt), CDDIS (Earthdata login wall on file
GETs; directory pages return 200 but files return an HTML login page),
BKG MGEX tree root_ftp/MGEX/products (empty) and .../IGS/products/mgex
(frozen at GPS week 2237), WHU igs.gnsswhu.cn, IGN gnss-data-portal.ign.fr
and igs.ensg.ign.fr, ftp.gfz.de, ASI gsc.eur.ac.it (DNS dead),
data.geo.ga.gov.au (unreachable).
Caveat: AIUB redirects www -> download.aiub.unibe.ch, whose DNS record is
intermittently unresolvable; retries hours later usually succeed, and the
garner fallbacks cover the gap meanwhile.

Output layout: {out_dir}/{YYYY}/{DDD}/{STREAM}_{YYYYMMDD}.sp3 / .clk
(gunzipped, header-validated). Re-runs reuse valid cached files.

Usage:
  python scripts/fetch_precise_products.py --date 2025-06-09
  python scripts/fetch_precise_products.py --year 2025 --doy 160 --latency rapid
"""

import argparse
import gzip
import os
import ssl
import sys
import time
import urllib.error
import urllib.request
from datetime import date, datetime, timedelta

GPS_EPOCH = date(1980, 1, 6)
HTTP_TIMEOUT_S = 120
MAX_ATTEMPTS = 3
USER_AGENT = "gneiss-precise-products/1.0 (static PPK engine)"
MIN_VALID_BYTES = 1024

AIUB_FINAL_DIR = "https://www.aiub.unibe.ch/download/CODE_MGEX/CODE/{year}/"
GARNER_WEEK_DIR = "https://garner.ucsd.edu/pub/products/{gpsweek}/"
BKG_WEEK_DIR = "https://igs.bkg.bund.de/root_ftp/IGS/products/{gpsweek}/"

# Ordered fallback chains per latency tier.
# Entry: (name, dir_template, stream, scope). scope "gps" = GPS-only product.
CHAINS = {
    "final": [
        ("aiub_cod_mgxfin", AIUB_FINAL_DIR, "COD0MGXFIN", "multi"),
        ("garner_jax_mgxfin", GARNER_WEEK_DIR, "JAX0MGXFIN", "multi"),
        ("garner_iac_mgxfin", GARNER_WEEK_DIR, "IAC0MGXFIN", "multi"),
        ("bkg_igs_opsfin", BKG_WEEK_DIR, "IGS0OPSFIN", "gps"),
    ],
    # NOTE: there is no operational combined multi-GNSS rapid yet; GFZ's MGX
    # rapid is the highest-quality open multi-GNSS rapid stream.
    "rapid": [
        ("garner_gfz_mgxrap", GARNER_WEEK_DIR, "GFZ0MGXRAP", "multi"),
        ("aiub_cod_mgxfin", AIUB_FINAL_DIR, "COD0MGXFIN", "multi"),
        ("bkg_igs_opsrap", BKG_WEEK_DIR, "IGS0OPSRAP", "gps"),
    ],
}


def gps_week(d):
    """GPS week number (weeks since 1980-01-06; weeks start Sunday)."""
    return (d - GPS_EPOCH).days // 7


def doy_of(d):
    return d.timetuple().tm_yday


def long_stem(stream, d):
    """IGS long-product stem, e.g. COD0MGXFIN_20251600000."""
    return f"{stream}_{d.year}{doy_of(d):03d}0000"


def product_url(dir_tmpl, stream, d, kind):
    """Full URL of one .gz product file ('sp3' orbit or 'clk' clock)."""
    orbit_spec = "05M_ORB.SP3"
    if stream.startswith("IGS0OPS"):
        orbit_spec = "15M_ORB.SP3"
    suffix = orbit_spec if kind == "sp3" else "30S_CLK.CLK"
    fname = f"{long_stem(stream, d)}_01D_{suffix}.gz"
    return dir_tmpl.format(year=d.year, gpsweek=gps_week(d)) + fname


def cached_paths(out_dir, d, stream):
    """(dir, sp3_path, clk_path) for a stream/day under out_dir."""
    base = os.path.join(out_dir, str(d.year), f"{doy_of(d):03d}")
    stem = f"{stream}_{d:%Y%m%d}"
    return base, os.path.join(base, stem + ".sp3"), os.path.join(base, stem + ".clk")


def looks_valid(path, head_marker):
    """File exists, non-trivial size, expected marker in first 4 KB."""
    try:
        if os.path.getsize(path) < MIN_VALID_BYTES:
            return False
        with open(path, "r", errors="replace") as fh:
            return head_marker in fh.read(4096)
    except OSError:
        return False


def tls_verify_failed(exc):
    """True when a URLError wraps an SSL certificate-verification failure."""
    reason = getattr(exc, "reason", None)
    if isinstance(reason, ssl.SSLCertVerificationError):
        return True
    return "CERTIFICATE_VERIFY_FAILED" in str(reason or "")


def http_get(url):
    """Download url bytes with retries/backoff; raises ConnectionError.

    TLS verification is attempted first; on certificate failure (common on
    macOS python.org builds and for mirrors with private CAs) the fetch is
    retried once unverified with a loud warning, since these are public
    read-only scientific mirrors whose payloads are format-validated anyway.
    """
    last_err = None
    for attempt in range(MAX_ATTEMPTS):
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            return _open_verified_or_unverified(req, url)
        except (urllib.error.HTTPError, OSError, TimeoutError) as exc:
            last_err = exc
            time.sleep(2 ** attempt)
    raise ConnectionError(f"{url}: {last_err}")


def _open_verified_or_unverified(req, url):
    try:
        return _urlopen_bytes(req, insecure=False)
    except urllib.error.URLError as exc:
        if not tls_verify_failed(exc):
            raise
        print(f"[warn] TLS verify failed for {url} ({exc.reason}); retrying unverified")
        return _urlopen_bytes(req, insecure=True)


def _urlopen_bytes(req, insecure):
    if insecure:
        ctx = ssl._create_unverified_context()  # noqa: SLF001 - documented fallback
        opener = urllib.request.build_opener(urllib.request.HTTPSHandler(context=ctx))
    else:
        opener = urllib.request.build_opener()
    with opener.open(req, timeout=HTTP_TIMEOUT_S) as resp:
        return resp.read()


def install_gunzip(raw, dest_path, marker):
    """Decompress payload, validate header marker, atomically place file."""
    text = gzip.decompress(raw).decode("utf-8", errors="replace")
    if marker not in text[:4096]:
        raise ValueError(f"payload does not look like expected format ({marker!r})")
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    tmp_path = dest_path + ".part"
    with open(tmp_path, "w") as fh:
        fh.write(text)
    os.replace(tmp_path, dest_path)


def fetch_pair(out_dir, d, source):
    """Fetch one source's orbit+clock pair with caching.

    Returns (paths, notes): paths maps 'sp3'/'clk' -> installed file;
    notes maps 'sp3'/'clk' -> 'cached' or '<n> bytes downloaded'.
    """
    _, sp3_path, clk_path = cached_paths(out_dir, d, source.stream)
    markers = {"sp3": "#", "clk": "RINEX VERSION"}
    paths, notes = {}, {}
    for kind in ("sp3", "clk"):
        dest = sp3_path if kind == "sp3" else clk_path
        if looks_valid(dest, markers[kind]):
            paths[kind] = dest
            notes[kind] = "cached"
            continue
        raw = http_get(product_url(source.dir_tmpl, source.stream, d, kind))
        install_gunzip(raw, dest, markers[kind])
        paths[kind] = dest
        notes[kind] = f"{len(raw)} B gz downloaded"
    return paths, notes


class Source:
    def __init__(self, name, dir_tmpl, stream, scope):
        self.name = name
        self.dir_tmpl = dir_tmpl
        self.stream = stream
        self.scope = scope


def chain_sources(latency):
    return [Source(*entry) for entry in CHAINS[latency]]


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Fetch multi-GNSS precise SP3+CLK products from open mirrors.")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--date", help="calendar date YYYY-MM-DD")
    src.add_argument("--year", type=int, help="year (use with --doy)")
    ap.add_argument("--doy", type=int, help="day of year (with --year)")
    ap.add_argument("--latency", choices=("final", "rapid", "auto"), default="auto",
                    help="product tier; auto picks rapid when recent, else final")
    ap.add_argument("--out-dir",
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         os.pardir, "datasets", "precise"))
    ap.add_argument("--neighbors", type=int, default=0,
                    help="also fetch N days before and after (for sessions "
                         "crossing midnight; interpolation needs edge epochs)")
    return ap.parse_args(argv)


def resolve_date(args):
    if args.date:
        return datetime.strptime(args.date, "%Y-%m-%d").date()
    return date(args.year, 1, 1) + timedelta(days=args.doy - 1)


def pick_tier(latency_arg, d):
    if latency_arg != "auto":
        return latency_arg
    age_days = (date.today() - d).days
    return "rapid" if 0 <= age_days <= 21 else "final"


def fetch_date(args, d):
    """Run the fallback chain for one date. Returns (exit_code, file_paths)."""
    print(f"target {d} (DOY {doy_of(d):03d}, GPS week {gps_week(d)})")
    failures = []
    for source in chain_sources(pick_tier(args.latency, d)):
        print(f"[try] {source.name} ({source.scope})")
        try:
            paths, notes = fetch_pair(args.out_dir, d, source)
            for kind in sorted(paths):
                print(f"[ok ] {kind}: {os.path.abspath(paths[kind])} ({notes[kind]})")
            return 0, [os.path.abspath(p) for p in sorted(paths.values())]
        except Exception as exc:  # noqa: BLE001 - whole point is the fallback chain
            print(f"[err] {source.name}: {exc}")
            failures.append(f"{source.name}: {exc}")
    print("ALL SOURCES FAILED:")
    for line in failures:
        print(f"  - {line}")
    return 1, []


def run(argv=None):
    args = parse_args(argv)
    first = resolve_date(args)
    dates = [first + timedelta(days=off) for off in range(-args.neighbors, args.neighbors + 1)]
    all_files, failed = [], []
    for d in dates:
        code, files = fetch_date(args, d)
        all_files.extend(files)
        if code != 0:
            failed.append(d)
    if failed:
        print(f"DAYS FAILED: {[str(d) for d in failed]}")
        return 1
    print(f"{len(all_files)} product files ready under {os.path.abspath(args.out_dir)}")
    return 0


if __name__ == "__main__":
    sys.exit(run())
