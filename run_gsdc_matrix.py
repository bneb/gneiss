import sys
sys.path.append('scripts')
import run_full_matrix

# override datasets
run_full_matrix.DATASETS = {
    "GSDC (Pixel 4)": {
        "rover": "datasets/gsdc/Pixel4_GnssLog.20o",
        "base": "datasets/gsdc/p2221350.20o",
        "nav": "datasets/gsdc/rover.nav",
        "truth": "datasets/gsdc/reference.csv",
        "conf": None,
        "gneiss_config": "datasets/gsdc/gsdc_config.json",
    }
}
run_full_matrix.OUT_DIR = "benchmarks/gsdc_matrix"
run_full_matrix.main()
