"""Analyse par cycle : pics, lignes de base, derive cycle a cycle.

Compare le modele livrable aux cibles numerisees (exp et theorie, p99 pour les
pics afin d'eviter l'attenuation du binning median).
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


def cycle_stats(t, y, period=18.0, t0=180.0, t1=380.0):
    peaks, bases = [], []
    edges = np.arange(t0, t1 + 1e-9, period)
    for a, b in zip(edges[:-1], edges[1:]):
        m = (t >= a) & (t < b)
        if np.count_nonzero(m) < 3:
            continue
        peaks.append(float(np.nanmax(y[m])))
        bases.append(float(np.nanmin(y[m])))
    return np.array(peaks), np.array(bases)


def raw_theory(eps, key):
    img = figure7.extract_figure7_image()
    arr = np.asarray(img)
    panel = figure7.PANEL_SPECS[eps]
    spec = panel[key]
    x0, y0, x1, y1 = spec["crop"]
    crop_rgb = arr[y0:y1, x0:x1]
    ys, xs = np.where(figure7.color_mask(crop_rgb, "black"))
    w, h = x1 - x0, y1 - y0
    ins = (xs > 8) & (xs < w - 8) & (ys > 4) & (ys < h - 4)
    x, y = figure7.pixels_to_data(xs[ins], ys[ins], spec["crop"], panel["xlim"], spec["ylim"])
    if key == "force":
        lo, hi = (1380.0, 1900.0) if eps == 0.8 else (1950.0, 2500.0)
    else:
        lo, hi = (-180.0, 900.0)
    v = (y >= lo) & (y <= hi)
    bx, by = figure7.bin_curve(x[v], y[v], 0.45, "p99")
    return bx, by


targets = figure7.digitize_figure7()
for eps in (0.8, 1.0):
    cfg, data = figure7.run_model_for_eps(eps, targets=targets)
    win = figure7.window_model_data(data)
    t = win["time"]
    print(f"\n===== eps={eps} =====")
    for key, model_key, unit in (
        ("force", "force_total_mN", "mN"),
        ("torque", "torque_act_microNm", "uNm"),
    ):
        mp, mb = cycle_stats(t, win[model_key])
        tx99, ty99 = raw_theory(eps, key)
        tp, tb = cycle_stats(tx99, ty99)
        ex, ey = targets[eps][key]
        # p99 experimental for peaks
        img = figure7.extract_figure7_image()
        arr = np.asarray(img)
        spec = figure7.PANEL_SPECS[eps][key]
        x0, y0, x1, y1 = spec["crop"]
        crop_rgb = arr[y0:y1, x0:x1]
        m = figure7.color_mask(crop_rgb, spec["color"])
        ys, xs = np.where(m)
        w, h = x1 - x0, y1 - y0
        ins = (xs > 8) & (xs < w - 8) & (ys > 4) & (ys < h - 4)
        x, y = figure7.pixels_to_data(xs[ins], ys[ins], spec["crop"], figure7.PANEL_SPECS[eps]["xlim"], spec["ylim"])
        lo, hi = spec["valid_ylim"]
        v = (y >= lo) & (y <= hi)
        ex99, ey99 = figure7.bin_curve(x[v], y[v], 0.45, "p99")
        ep, eb = cycle_stats(ex99, ey99)
        print(f"-- {key} ({unit}) --")
        print(f"model  peaks: {np.round(mp,1)}")
        print(f"theory peaks (p99): {np.round(tp,1)}")
        print(f"exp    peaks (p99): {np.round(ep,1)}")
        print(f"model  bases: {np.round(mb,1)}")
        print(f"theory bases (p99-binned min): {np.round(tb,1)}")
        print(f"exp    bases: {np.round(eb,1)}")
        n = min(len(mp), len(tp))
        print(
            f"peak deficit model-theory: mean {np.mean(mp[:n]-tp[:n]):+.1f} {unit}, "
            f"drift model (last-first peak) {mp[-1]-mp[0]:+.1f}, "
            f"drift theory {tp[-1]-tp[0]:+.1f}, drift exp {ep[-1]-ep[0]:+.1f}"
        )
