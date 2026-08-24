"""Analyse rapide de la numérisation seule (pas de run modèle).

- période des cycles de pression numérisés
- pics par cycle des courbes théorie/expérience numérisées
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
_orig = figure7.extract_figure7_image
_cache = {}


def patched(pdf_path=REAL_PDF):
    if "img" not in _cache:
        _cache["img"] = _orig(REAL_PDF)
    return _cache["img"]


figure7.extract_figure7_image = patched

exp = figure7.digitize_figure7()
theo = figure7.digitize_figure7_theory()


def cycle_peaks(t, y, min_gap=8.0):
    """Pics locaux séparés d'au moins min_gap secondes."""
    t = np.asarray(t)
    y = np.asarray(y)
    order = np.argsort(t)
    t, y = t[order], y[order]
    peaks = []
    i = 0
    while i < len(t):
        window = (t >= t[i]) & (t < t[i] + min_gap)
        j = np.argmax(y[window]) + np.flatnonzero(window)[0]
        # étendre si le max est au bord de la fenêtre
        peaks.append((t[j], y[j]))
        i = np.flatnonzero(window)[-1] + 1
    # fusionner pics trop proches
    merged = []
    for tp, yp in peaks:
        if merged and tp - merged[-1][0] < min_gap:
            if yp > merged[-1][1]:
                merged[-1] = (tp, yp)
        else:
            merged.append((tp, yp))
    return merged


for eps in (0.8, 1.0):
    print(f"\n=== eps = {eps} ===")
    pt, pp = exp[eps]["pressure"]
    pk = [(t, p) for t, p in cycle_peaks(pt, pp, min_gap=10.0) if p > 0.8]
    times = np.array([t for t, _ in pk])
    print(f"pression: {len(pk)} pics, t = {np.round(times, 1)}")
    if len(times) > 1:
        print(
            f"  période moyenne = {np.mean(np.diff(times)):.2f} s, "
            f"pics = {np.round([p for _, p in pk], 2)}"
        )
    for name, data in (("exp", exp), ("theorie", theo)):
        for key in ("force", "torque"):
            if key not in data[eps]:
                continue
            tt, yy = data[eps][key]
            pk = cycle_peaks(tt, yy, min_gap=12.0)
            # ne garder que les vrais pics de cycle (au-dessus de la médiane)
            med = np.median(yy)
            span = np.max(yy) - med
            pk = [(t, y) for t, y in pk if y > med + 0.4 * span]
            print(f"{name} {key}: pics par cycle =")
            print("   t: " + " ".join(f"{t:7.1f}" for t, _ in pk))
            print("   y: " + " ".join(f"{y:7.1f}" for _, y in pk))
            vals = np.array([y for _, y in pk])
            if len(vals) > 2:
                print(
                    f"   premier {vals[0]:.1f} -> dernier {vals[-1]:.1f} "
                    f"(decroissance {100 * (vals[0] - vals[-1]) / max(abs(vals[0]), 1e-9):.1f} %)"
                )
