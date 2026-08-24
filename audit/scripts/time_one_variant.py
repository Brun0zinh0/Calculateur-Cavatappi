"""Chronometre une seule variante du criblage physique (2 cas eps)."""
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
_orig = figure7.extract_figure7_image
_cache = {}


def patched(pdf_path=REAL_PDF):
    if "img" not in _cache:
        _cache["img"] = _orig(REAL_PDF)
    return _cache["img"]


figure7.extract_figure7_image = patched

import validation_figure7_physique as phys  # noqa: E402

t0 = time.time()
rows = phys.run_screening(n_layers=1, n_phi=8, workers=1, names={"GM-EXP-ETAB-BLIN-SFIX-P0"})
dt = time.time() - t0
phys.print_results(rows, limit=2)
print(f"\nelapsed for 1 variant (2 eps cases): {dt:.1f} s")
