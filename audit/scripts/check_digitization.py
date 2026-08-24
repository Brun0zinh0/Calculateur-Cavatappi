"""Verifie l'attenuation des pics par la numerisation (median vs p99)."""
import sys
from pathlib import Path

WORKSPACE = Path(
    r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents"
    r"\Stage muscle artificièle\modèle\modèle chinois\Espace de travail"
)
sys.path.insert(0, str(WORKSPACE))

import numpy as np
import validation_figure7 as figure7

REAL_PDF = WORKSPACE / "Article support" / figure7.PDF_FILENAME
img = figure7.extract_figure7_image(REAL_PDF)
arr = np.asarray(img)

for eps, panel in figure7.PANEL_SPECS.items():
    xlim = panel["xlim"]
    for key in ("force", "torque"):
        spec = panel[key]
        x0, y0, x1, y1 = spec["crop"]
        crop_rgb = arr[y0:y1, x0:x1]
        # experimental colored markers
        mask = figure7.color_mask(crop_rgb, spec["color"])
        ys, xs = np.where(mask)
        w, h = x1 - x0, y1 - y0
        inside = (xs > 8) & (xs < w - 8) & (ys > 4) & (ys < h - 4)
        x, y = figure7.pixels_to_data(xs[inside], ys[inside], spec["crop"], xlim, spec["ylim"])
        ylo, yhi = spec["valid_ylim"]
        v = (y >= ylo) & (y <= yhi)
        bx_med, by_med = figure7.bin_curve(x[v], y[v], 0.45, "median")
        bx_p99, by_p99 = figure7.bin_curve(x[v], y[v], 0.45, "p99")
        print(
            f"eps={eps} {key} EXP: raw max={np.max(y[v]):.1f}  "
            f"median-binned max={np.max(by_med):.1f}  p99-binned max={np.max(by_p99):.1f}"
        )
        # theory black curve
        ysb, xsb = np.where(figure7.color_mask(crop_rgb, "black"))
        insb = (xsb > 8) & (xsb < w - 8) & (ysb > 4) & (ysb < h - 4)
        xt, yt = figure7.pixels_to_data(xsb[insb], ysb[insb], spec["crop"], xlim, spec["ylim"])
        if key == "force":
            tlo, thi = (1380.0, 1900.0) if eps == 0.8 else (1950.0, 2500.0)
        else:
            tlo, thi = (-180.0, 900.0)
        vt = (yt >= tlo) & (yt <= thi)
        btx_med, bty_med = figure7.bin_curve(xt[vt], yt[vt], 0.45, "median")
        btx_p99, bty_p99 = figure7.bin_curve(xt[vt], yt[vt], 0.45, "p99")
        print(
            f"eps={eps} {key} THEORY: raw max={np.max(yt[vt]):.1f}  "
            f"median-binned max={np.max(bty_med):.1f}  p99-binned max={np.max(bty_p99):.1f}  "
            f"raw min={np.min(yt[vt]):.1f}  median-binned min={np.min(bty_med):.1f}"
        )
