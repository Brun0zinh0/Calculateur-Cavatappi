"""Wrapper: exécute livrable/validation_figure7_physique.py sans modifier l'original.

Corrige uniquement le chemin du PDF (déplacé dans 'Article support').
"""
import sys
import time
from pathlib import Path

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

import validation_figure7_physique as physique  # noqa: E402

if __name__ == "__main__":
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 48
    t0 = time.time()
    rows = physique.run_screening(n_layers=1, n_phi=8, workers=workers)
    physique.print_results(rows, limit=limit)
    print(f"\ntotal time: {time.time() - t0:.1f} s, variants: {len(rows)}")
