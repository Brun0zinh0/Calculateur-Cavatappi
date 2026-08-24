"""Wrapper: exécute validation_figure7.plot_validation sans modifier l'original.

Corrige uniquement le chemin du PDF (déplacé dans 'Article support') et
sauvegarde la figure au lieu de plt.show().
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

WORKSPACE = Path(
    r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents"
    r"\Stage muscle artificièle\modèle\modèle chinois\Espace de travail"
)
SCRATCH = Path(
    r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude"
    r"\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail"
    r"\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
)
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

fig, metrics = figure7.plot_validation(show=False)
out = SCRATCH / "validation_figure7_overlay.png"
fig.savefig(out, dpi=160)
print(f"figure saved: {out}")

print("\n=== METRICS DICT ===")
for eps, rows in metrics.items():
    for key, vals in rows.items():
        print(f"eps={eps} {key}: theory={vals['theory']:.3f} experiment={vals['experiment']:.3f}")
