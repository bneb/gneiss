import gzip
import os
import tempfile
import unittest
from datetime import date
from unittest import mock

from fetch_precise_products import (
    CHAINS,
    Source,
    cached_paths,
    chain_sources,
    doy_of,
    gps_week,
    install_gunzip,
    looks_valid,
    long_stem,
    parse_args,
    pick_tier,
    product_url,
    resolve_date,
    run,
)

FAKE_SP3 = "#dP2025  6  9  0  0  0.00000000     289   u+U IGb20 FIT\n" + "x" * 2048
FAKE_CLK = "     3.04                 C                    M   RINEX VERSION / TYPE\n" + "y" * 2048


class TestCalendarMath(unittest.TestCase):
    def test_gps_week_doy160_2025(self):
        self.assertEqual(gps_week(date(2025, 6, 9)), 2370)
        self.assertEqual(doy_of(date(2025, 6, 9)), 160)

    def test_gps_epoch_is_week_zero_sunday(self):
        self.assertEqual((gps_week(date(1980, 1, 6)), date(1980, 1, 6).weekday()), (0, 6))

    def test_long_stem(self):
        self.assertEqual(long_stem("COD0MGXFIN", date(2025, 6, 9)), "COD0MGXFIN_20251600000")


class TestUrlBuilders(unittest.TestCase):
    def test_aiub_final_sp3(self):
        url = product_url(
            "https://www.aiub.unibe.ch/download/CODE_MGEX/CODE/{year}/",
            "COD0MGXFIN", date(2025, 6, 9), "sp3")
        self.assertEqual(url,
            "https://www.aiub.unibe.ch/download/CODE_MGEX/CODE/2025/"
            "COD0MGXFIN_20251600000_01D_05M_ORB.SP3.gz")

    def test_aiub_final_clk(self):
        url = product_url(
            "https://www.aiub.unibe.ch/download/CODE_MGEX/CODE/{year}/",
            "COD0MGXFIN", date(2025, 6, 9), "clk")
        self.assertIn("COD0MGXFIN_20251600000_01D_30S_CLK.CLK.gz", url)

    def test_garner_rapid_uses_gps_week_dir(self):
        url = product_url(
            "https://garner.ucsd.edu/pub/products/{gpsweek}/",
            "GFZ0MGXRAP", date(2025, 6, 9), "sp3")
        self.assertEqual(
            "https://garner.ucsd.edu/pub/products/2370/GFZ0MGXRAP_20251600000_01D_05M_ORB.SP3.gz",
            url)

    def test_bkg_gps_only_orbit_sampling(self):
        url = product_url(
            "https://igs.bkg.bund.de/root_ftp/IGS/products/{gpsweek}/",
            "IGS0OPSFIN", date(2025, 6, 9), "sp3")
        self.assertIn("IGS0OPSFIN_20251600000_01D_15M_ORB.SP3.gz", url)


class TestCacheValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, "x.sp3")

    def tearDown(self):
        self.tmp.cleanup()

    def test_missing_file_invalid(self):
        self.assertFalse(looks_valid(self.path, "#"))

    def test_valid_sp3_header(self):
        with open(self.path, "w") as fh:
            fh.write(FAKE_SP3)
        self.assertTrue(looks_valid(self.path, "#"))

    def test_too_small_invalid(self):
        with open(self.path, "w") as fh:
            fh.write("#")
        self.assertFalse(looks_valid(self.path, "#"))

    def test_install_gunzip_roundtrip_and_rejects_bad_payload(self):
        install_gunzip(gzip.compress(FAKE_SP3.encode()), self.path, "#")
        self.assertTrue(looks_valid(self.path, "#"))
        with open(self.path) as fh:
            self.assertEqual(fh.readline().strip(), FAKE_SP3.splitlines()[0])
        with self.assertRaises(ValueError):
            install_gunzip(gzip.compress(b"not sp3 at all"), self.path, "#")

    def test_cached_paths_layout(self):
        base, sp3, clk = cached_paths("/tmp/out", date(2025, 6, 9), "GFZ0MGXRAP")
        self.assertEqual(base, os.path.join("/tmp/out", "2025", "160"))
        self.assertTrue(sp3.endswith("GFZ0MGXRAP_20250609.sp3"))
        self.assertTrue(clk.endswith("GFZ0MGXRAP_20250609.clk"))


class TestTierSelectionAndArgs(unittest.TestCase):
    def test_pick_tier_explicit(self):
        self.assertEqual(pick_tier("final", date(2025, 6, 9)), "final")

    def test_pick_tier_auto_old_date_prefers_final(self):
        self.assertEqual(pick_tier("auto", date(2000, 1, 1)), "final")

    def test_pick_tier_auto_recent_prefers_rapid(self):
        from datetime import date as _date, timedelta as _td
        recent = _date.today() - _td(days=2)
        self.assertEqual(pick_tier("auto", recent), "rapid")

    def test_resolve_date_from_year_doy(self):
        args = parse_args(["--year", "2025", "--doy", "160"])
        self.assertEqual(resolve_date(args), date(2025, 6, 9))

    def test_resolve_date_from_iso(self):
        args = parse_args(["--date", "2025-06-09"])
        self.assertEqual(resolve_date(args), date(2025, 6, 9))


class TestRunEndToEnd(unittest.TestCase):
    """Full run() with the network layer mocked; verifies caching + fallback."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def _fake_http(self, routes):
        """routes: list of (url_substring, {suffix_key: payload_text})."""
        def fake_get(url):
            for marker, payloads in routes:
                if marker in url:
                    for suffix, payload in payloads.items():
                        if suffix in url:
                            return gzip.compress(payload.encode())
            raise ConnectionError(f"404 for {url}")
        return fake_get

    def test_run_downloads_pair_and_caches(self):
        src = Source("aiub_cod_mgxfin",
                     "https://example.test/CODE_MGEX/CODE/{year}/",
                     "COD0MGXFIN", "multi")
        fake = self._fake_http([("example.test", {"ORB.SP3": FAKE_SP3, "CLK.CLK": FAKE_CLK})])
        with mock.patch("fetch_precise_products.http_get", side_effect=fake), \
             mock.patch("fetch_precise_products.chain_sources", return_value=[src]):
            code = run(["--date", "2025-06-09", "--out-dir", self.out,
                        "--latency", "final"])
        _, sp3, clk = cached_paths(self.out, date(2025, 6, 9), "COD0MGXFIN")
        files = [sp3, clk]
        self.assertEqual(code, 0)
        self.assertTrue(os.path.exists(sp3) and os.path.exists(clk))
        # Second run must be served entirely from cache.
        with mock.patch("fetch_precise_products.http_get",
                        side_effect=AssertionError("network hit on cache path")), \
             mock.patch("fetch_precise_products.chain_sources", return_value=[src]):
            code2 = run(["--date", "2025-06-09", "--out-dir", self.out,
                         "--latency", "final"])
        self.assertEqual(code2, 0)
        self.assertTrue(os.path.exists(sp3) and os.path.exists(clk),
                        "cache re-run must leave both files in place")

    def test_run_falls_back_when_primary_fails(self):
        primary = Source("bad_src", "https://bad.test/{gpsweek}/", "GFZ0MGXRAP", "multi")
        fallback = Source("good_src", "https://good.test/CODE_MGEX/CODE/{year}/",
                          "COD0MGXFIN", "multi")
        fake = self._fake_http([("bad.test", {}),
                                ("good.test", {"ORB.SP3": FAKE_SP3, "CLK.CLK": FAKE_CLK})])
        with mock.patch("fetch_precise_products.http_get", side_effect=fake), \
             mock.patch("fetch_precise_products.chain_sources",
                        return_value=[primary, fallback]), \
             mock.patch("fetch_precise_products.time.sleep"):
            code = run(["--date", "2025-06-09", "--out-dir", self.out,
                        "--latency", "rapid"])
        self.assertEqual(code, 0)
        _, sp3, _ = cached_paths(self.out, date(2025, 6, 9), "COD0MGXFIN")
        self.assertTrue(os.path.exists(sp3), "fallback source file should exist")

    def test_run_returns_error_when_all_fail(self):
        src = Source("dead_src", "https://dead.test/", "GFZ0MGXRAP", "multi")

        def always_fail(url):
            raise ConnectionError("nope")
        with mock.patch("fetch_precise_products.http_get", side_effect=always_fail), \
             mock.patch("fetch_precise_products.chain_sources", return_value=[src]), \
             mock.patch("time.sleep"):
            code = run(["--date", "2025-06-09", "--out-dir", self.out,
                        "--latency", "rapid"])
        self.assertEqual(code, 1)

    def test_neighbors_fetch_adjacent_days(self):
        src = Source("aiub_cod_mgxfin",
                     "https://example.test/CODE_MGEX/CODE/{year}/",
                     "COD0MGXFIN", "multi")
        fake = self._fake_http([("example.test", {"ORB.SP3": FAKE_SP3, "CLK.CLK": FAKE_CLK})])
        with mock.patch("fetch_precise_products.http_get", side_effect=fake), \
             mock.patch("fetch_precise_products.chain_sources", return_value=[src]), \
             mock.patch("fetch_precise_products.time.sleep"):
            code = run(["--date", "2025-06-09", "--out-dir", self.out,
                        "--latency", "final", "--neighbors", "1"])
        self.assertEqual(code, 0)
        for day in (8, 9, 10):
            _, sp3, _ = cached_paths(self.out, date(2025, 6, day), "COD0MGXFIN")
            self.assertTrue(os.path.exists(sp3), f"missing {sp3}")


if __name__ == "__main__":
    unittest.main()
