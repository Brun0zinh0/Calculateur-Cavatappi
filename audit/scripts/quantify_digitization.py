"""Quantifie l'attenuation des pics par la numerisation (sans simulation)."""
import sys
from functools import partial
from pathlib import Path

import numpy as np

WORKSPACE = Path(r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail")
LIVRABLE = WORKSPACE / "livrable"
PDF_REAL = WORKSPACE / "Article support" / "Blocked actuation of twisted and coiled tube-based polymer actuators driven by hydraulic pressure.pdf"
for p in (str(LIVRABLE), str(WORKSPACE)):
    if p not in sys.path:
        sys.path.insert(0, p)

import validation_figure7 as f7

_orig = f7.extract_figure7_image
f7.PDF_PATH = PDF_REAL
f7.extract_figure7_image = partial(_orig, PDF_REAL)

targets = f7.digitize_figure7()
theory = f7.digitize_figure7_theory()

PAPER = {
    0.8: {"F_peak": 1824.60, "T_peak": 877.01, "P_peak": (1.41, 1.43)},
    1.0: {"F_peak": 2424.67, "T_peak": 873.41, "P_peak": (1.41, 1.43)},
}

for eps in (0.8, 1.0):
    tp, pp = targets[eps]["pressure"]
    tf, ff = targets[eps]["force"]
    tt, tq = targets[eps]["torque"]
    thf = theory[eps]["force"][1]
    tht = theory[eps]["torque"][1]
    print(f"eps={eps:.1f}")
    print(f"  pression numerisee: max {np.nanmax(pp):.3f} MPa, min {np.nanmin(pp):.3f} (article: pics 1.41-1.43 MPa)")
    print(f"  force exp numerisee: max {np.nanmax(ff):.1f} mN (article 11e cycle: {PAPER[eps]['F_peak']:.1f} mN)")
    print(f"  couple exp numerise: max {np.nanmax(tq):.1f} uNm (article 11e cycle: {PAPER[eps]['T_peak']:.1f} uNm)")
    print(f"  force theorie numerisee: max {np.nanmax(thf):.1f} mN ; couple theorie numerise: max {np.nanmax(tht):.1f} uNm")
    print(f"  nb points: pression {len(pp)}, force {len(ff)}, couple {len(tq)}")
