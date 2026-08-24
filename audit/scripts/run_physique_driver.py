"""Driver scratchpad : execute livrable/validation_figure7_physique.run_screening
sans modifier le projet (PDF redirige vers 'Article support')."""
import sys
import time
from functools import partial
from pathlib import Path

WORKSPACE = Path(r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail")
LIVRABLE = WORKSPACE / "livrable"
PDF_REAL = WORKSPACE / "Article support" / "Blocked actuation of twisted and coiled tube-based polymer actuators driven by hydraulic pressure.pdf"

for p in (str(LIVRABLE), str(WORKSPACE)):
    if p not in sys.path:
        sys.path.insert(0, p)

import validation_figure7 as f7  # noqa: E402

_orig_extract = f7.extract_figure7_image
f7.PDF_PATH = PDF_REAL
f7.extract_figure7_image = partial(_orig_extract, PDF_REAL)

import validation_figure7_physique as phys  # noqa: E402

if __name__ == "__main__":
    t0 = time.time()
    print(f"Base module used: {phys.Base.__file__}")
    print(f"Base MODEL_VERSION: {getattr(phys.Base, 'MODEL_VERSION', 'unknown')}")
    rows = phys.run_screening(n_layers=1, n_phi=8, workers=4)
    phys.print_results(rows, limit=72)
    print(f"total time: {time.time() - t0:.1f} s, variants: {len(rows)}")
