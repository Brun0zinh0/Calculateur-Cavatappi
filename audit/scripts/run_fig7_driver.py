"""Driver scratchpad : execute validation_figure7.plot_validation sans modifier le projet.

- Ajoute l'espace de travail au sys.path.
- Redirige le chemin du PDF vers 'Article support'.
- Sauvegarde la figure au lieu de plt.show().
"""
import sys
import time
from functools import partial
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

WORKSPACE = Path(r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail")
SCRATCH = Path(r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad")
PDF_REAL = WORKSPACE / "Article support" / "Blocked actuation of twisted and coiled tube-based polymer actuators driven by hydraulic pressure.pdf"

sys.path.insert(0, str(WORKSPACE))

import validation_figure7 as f7  # noqa: E402

# Redirection du PDF (l'argument par defaut est fige a l'import, on repatche la fonction).
_orig_extract = f7.extract_figure7_image
f7.PDF_PATH = PDF_REAL
f7.extract_figure7_image = partial(_orig_extract, PDF_REAL)

if __name__ == "__main__":
    t0 = time.time()
    print(f"Base module used: {f7.Base.__file__}")
    print(f"Base MODEL_VERSION: {getattr(f7.Base, 'MODEL_VERSION', 'unknown')}")
    fig, metrics = f7.plot_validation(show=False)
    out = SCRATCH / "fig7_validation_overlay.png"
    fig.savefig(out, dpi=150)
    print(f"Figure sauvegardee: {out}")
    print(f"Duree totale: {time.time() - t0:.1f} s")
