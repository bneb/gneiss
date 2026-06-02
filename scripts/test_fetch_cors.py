import unittest
import os
import tempfile
from datetime import datetime
from fetch_cors import CorsFetcher

class TestCorsFetcher(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.fetcher = CorsFetcher(output_dir=self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_date_to_doy(self):
        dt = datetime(2020, 5, 14)
        doy = self.fetcher._get_doy(dt)
        self.assertEqual(doy, 135)
        
    def test_build_ngs_url(self):
        dt = datetime(2020, 5, 14)
        url = self.fetcher._build_ngs_url("p222", dt)
        self.assertEqual(url, "https://geodesy.noaa.gov/corsdata/rinex/2020/135/p222/p2221350.20d.gz")

    def test_end_to_end_fetch(self):
        # We test against a real known file but mock the download in a real TDD suite.
        # Since this is an integration test, we will actually try to fetch it.
        dt = datetime(2020, 5, 14)
        obs_file = self.fetcher.fetch("p222", dt)
        
        self.assertTrue(os.path.exists(obs_file), "The uncompressed .obs file should exist")
        self.assertTrue(obs_file.endswith(".20o"), "The file should be a decompressed observation file")
        
        # Check if the file has content and RINEX headers
        with open(obs_file, 'r') as f:
            header = f.readline()
            self.assertIn("RINEX VERSION", header)
            self.assertIn("OBSERVATION DATA", header)

if __name__ == '__main__':
    unittest.main()
