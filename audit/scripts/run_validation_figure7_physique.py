"""Wrapper scratchpad : execute livrable/validation_figure7_physique sans
modifier les originaux. Patch du chemin PDF (le PDF est dans 'Article support').
"""
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

WORKSPACE = Path(
    r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents"
    r"\Stage muscle artificièle\modèle\modèle chinois\Espace de travail"
)
LIVRABLE = WORKSPACE / "livrable"
sys.path.insert(0, str(LIVRABLE))
sys.path.insert(0, str(WORKSPACE))

import validation_figure7 as figure7  # noqa: E402

REAL_PDF = WORKSPACE / "Article support" / figure7.PDF_FILENAME
assert REAL_PDF.exists(), REAL_PDF

_orig_extract = figure7.extract_figure7_image
_cache = {}


def patched_extract(pdf_path=REAL_PDF):
    if "img" not in _cache:
        _cache["img"] = _orig_extract(REAL_PDF)
    return _cache["img"]


figure7.extract_figure7_image = patched_extract

import validation_figure7_physique as phys  # noqa: E402

if __name__ == "__main__":
    t0 = time.time()
    rows = phys.run_screening(n_layers=1, n_phi=8, workers=10, names=None)
    phys.print_results(rows, limit=96)
    print(f"\nelapsed: {time.time() - t0:.1f} s")
