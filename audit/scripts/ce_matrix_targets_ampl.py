# -*- coding: utf-8 -*-
"""Amplitudes cibles de la theorie article numerisee (unites vraies)."""
import sys
from pathlib import Path

import numpy as np
from PIL import Image

SCRATCH = Path(
    r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude"
    r"\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-"
    r"artifici-le-mod-le-mod-le-chinois-Espace-de-travail"
    r"\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
)
WORKSPACE = Path(
    r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents"
    r"\Stage muscle artificièle\modèle\modèle chinois\Espace de travail"
)
sys.path.insert(0, str(WORKSPACE))
import validation_figure7 as f7  # noqa: E402

REAL_PDF = WORKSPACE / "Article support" / f7.PDF_FILENAME
_orig = f7.extract_figure7_image
f7.extract_figure7_image = lambda pdf_path=REAL_PDF: _orig(REAL_PDF)

img = np.asarray(Image.open(SCRATCH / "p12_img0.png").convert("RGB")).astype(int)
NONWHITE = img.sum(axis=2) < 3 * 225
PANELS = {
    (0.8, "force"):  dict(crop=(120, 51, 704, 350),  ticks=[2200, 2000, 1800, 1600, 1400, 1200], ps=(1200., 2200.)),
    (0.8, "torque"): dict(crop=(120, 351, 704, 647), ticks=[1650, 1100, 550, 0],                 ps=(-250., 1650.)),
    (1.0, "force"):  dict(crop=(899, 51, 1484, 350), ticks=[2800, 2600, 2400, 2200, 2000],      ps=(2000., 2800.)),
    (1.0, "torque"): dict(crop=(899, 351, 1484, 647), ticks=[1440, 960, 480, 0],                ps=(-250., 1800.)),
}


def label_rows(crop, n_expected):
    x0, y0, x1, y1 = crop
    band = NONWHITE[y0:y1, x0 - 48:x0 - 6]
    rr = np.where(band.any(axis=1))[0]
    groups, start, prev = [], rr[0], rr[0]
    for v in rr[1:]:
        if v - prev <= 3:
            prev = v
        else:
            groups.append((start, prev)); start = prev = v
    groups.append((start, prev))
    groups = [g for g in groups if 8 <= (g[1] - g[0]) <= 34]
    assert len(groups) == n_expected
    return [0.5 * (a + b) for a, b in groups]


def script_to_true(eps, kind, y_script):
    spec = PANELS[(eps, kind)]
    rows_l = np.asarray(label_rows(spec["crop"], len(spec["ticks"])), float)
    A, B = np.polyfit(rows_l, np.asarray(spec["ticks"], float), 1)
    lo, hi = spec["ps"]
    h = spec["crop"][3] - spec["crop"][1]
    height = max(1.0, float(h - 1))
    row = (hi - np.asarray(y_script)) * height / (hi - lo)
    return A * row + B


def cycle_extrema(t, v, t0=185.0, min_gap=8.0):
    t = np.asarray(t, float); v = np.asarray(v, float)
    m = t >= t0
    t, v = t[m], v[m]
    med = np.median(v)

    def _find(sig, ref):
        out = []
        for i in range(1, len(sig) - 1):
            if sig[i] >= sig[i - 1] and sig[i] > sig[i + 1] and sig[i] > ref:
                if not out or t[i] - out[-1][0] > min_gap:
                    out.append((t[i], sig[i]))
                elif sig[i] > out[-1][1]:
                    out[-1] = (t[i], sig[i])
        return out

    pk = np.asarray([p[1] for p in _find(v, med)])
    vl = np.asarray([-p[1] for p in _find(-v, -med)])
    return pk, vl


th = f7.digitize_figure7_theory()
print("Theorie article numerisee (courbe binned mediane, unites vraies) :")
for eps in (0.8, 1.0):
    for kind in ("force", "torque"):
        tx, ty = th[eps][kind]
        yt = script_to_true(eps, kind, ty)
        pk, vl = cycle_extrema(tx, yt)
        unit = "mN" if kind == "force" else "uNm"
        print(f"  eps={eps} {kind:6s}: pics n={pk.size} moy={pk.mean():.1f} "
          f"[{pk.min():.1f}..{pk.max():.1f}]  vallees moy={vl.mean():.1f}  "
          f"amplitude={pk.mean()-vl.mean():.1f} {unit}")
