"""Detecte les graduations (ticks) des axes Y pour verifier les ylim supposes."""
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


def ticks_rows(x_frame, y_top, y_bot, side="in"):
    if side == "in":
        band = dark[y_top:y_bot, x_frame + 2 : x_frame + 9]
    else:
        band = dark[y_top:y_bot, x_frame - 9 : x_frame - 2]
    counts = band.sum(axis=1)
    rows = np.where(counts >= 5)[0]
    groups = []
    for i in rows:
        if groups and i - groups[-1][-1] <= 2:
            groups[-1].append(int(i))
        else:
            groups.append([int(i)])
    return [g[len(g) // 2] + y_top for g in groups]


# left spine of panel a is around x=120-121, panel b around x=899-900
for label, x_frame in (("panel a", 121), ("panel b", 900)):
    for name, ytop, ybot in (("force", 51, 348), ("torque", 348, 646), ("pressure", 646, 944)):
        for side in ("in", "out"):
            t = ticks_rows(x_frame, ytop, ybot, side)
            if t:
                print(f"{label} {name} ticks ({side}side): rows {t}")
                break
        else:
            print(f"{label} {name}: no ticks found")


def value_at(row, crop, ylim):
    x0, y0, x1, y1 = crop
    height = max(1.0, float(y1 - y0 - 1))
    return ylim[1] - (row - y0) / height * (ylim[1] - ylim[0])


print("\nvaleurs des ticks selon le mapping du script :")
for eps, panel in figure7.PANEL_SPECS.items():
    x_frame = 121 if eps == 0.8 else 900
    for name in ("force", "torque", "pressure"):
        spec = panel[name]
        crop = spec["crop"]
        ytop, ybot = (51, 348) if name == "force" else ((348, 646) if name == "torque" else (646, 944))
        for side in ("in", "out"):
            t = ticks_rows(x_frame, ytop, ybot, side)
            if t:
                break
        vals = [round(value_at(r, crop, spec["ylim"]), 1) for r in t]
        print(f"eps={eps} {name}: tick values -> {vals}")
