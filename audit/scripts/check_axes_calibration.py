"""Verifie que les crops de PANEL_SPECS coincident avec les cadres des axes.

Detecte les lignes horizontales/verticales noires longues (cadres de panneaux)
dans l'image de la figure 7 et les compare aux crops du script.
"""
import sys
from pathlib import Path

import numpy as np

WORKSPACE = Path(
    r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents"
    r"\Stage muscle artificièle\modèle\modèle chinois\Espace de travail"
)
sys.path.insert(0, str(WORKSPACE))

import validation_figure7 as figure7  # noqa: E402

REAL_PDF = WORKSPACE / "Article support" / figure7.PDF_FILENAME
img = figure7.extract_figure7_image(REAL_PDF)
arr = np.asarray(img)
dark = (arr[..., 0] < 120) & (arr[..., 1] < 120) & (arr[..., 2] < 120)

H, W = dark.shape
print("image:", W, "x", H)

for label, xa, xb in (("panel a", 130, 690), ("panel b", 910, 1470)):
    rows = dark[:, xa:xb].sum(axis=1)
    frame_rows = [int(i) for i in np.where(rows > 0.85 * (xb - xa))[0]]
    # group consecutive
    groups = []
    for i in frame_rows:
        if groups and i - groups[-1][-1] <= 2:
            groups[-1].append(i)
        else:
            groups.append([i])
    print(f"{label}: horizontal frame lines at rows {[g[0] for g in groups]}")

for label, ya, yb in (("left col", 60, 940),):
    cols = dark[ya:yb, :].sum(axis=0)
    frame_cols = [int(i) for i in np.where(cols > 0.85 * (yb - ya))[0]]
    groups = []
    for i in frame_cols:
        if groups and i - groups[-1][-1] <= 2:
            groups[-1].append(i)
        else:
            groups.append([i])
    print(f"vertical frame lines (rows {ya}-{yb}): {[g[0] for g in groups]}")

print("\ncrops in script:")
for eps, panel in figure7.PANEL_SPECS.items():
    for key in ("force", "torque", "pressure"):
        print(f"  eps={eps} {key}: crop={panel[key]['crop']} ylim={panel[key]['ylim']}")
