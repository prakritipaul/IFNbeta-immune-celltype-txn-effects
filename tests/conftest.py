from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "realdata: tests that read the Kang 2018 repository data files"
    )
    config.addinivalue_line(
        "markers", "r: tests that run an R script via Rscript"
    )
