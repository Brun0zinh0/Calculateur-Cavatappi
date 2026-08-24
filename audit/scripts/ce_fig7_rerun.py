# -*- coding: utf-8 -*-
"""Contre-expertise : le biais de numerisation de la pression explique-t-il le deficit ?

1. Baseline : reproduction des runs de validation (pression numerisee biaisee).
2. Pression corrigee : conversion script->vrai via calibration des ticks
   (pics ~1.41 au lieu de 1.32/1.37 ; vallees ~0.01 au lieu de 0.08).
3. Convergence : (5,32), dt=0.5 a eps=1.0 avec pression corrigee.
"""
import sys
import time
from pathlib import Path

import numpy as np

WORKSPACE = Path(
    r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents"
    r"\Stage muscle artificièle\modèle\modèle chinois\Espace de travail"
)
sys.path.insert(0, str(WORKSPACE))

import validation_figure7 as f7  # noqa: E402

REAL_PDF = WORKSPACE / "Article support" / f7.PDF_FILENAME
_orig = f7.extract_figure7_image
_cache = {}


def patched(pdf_path=REAL_PDF):
    if "img" not in _cache:
        _cache["img"] = _orig(REAL_PDF)
    return _cache["img"]


f7.extract_figure7_image = patched

# Calibration vraie (fit sur rows des etiquettes mesurees dans ce_fig7_measure2)
TICK_FITS = {
    0.8: ([21.0, 107.5, 195.0, 282.0], [1.5, 1.0, 0.5, 0.0], (120, 648, 704, 946)),
    1.0: ([15.0, 104.0, 192.5, 282.0], [1.5, 1.0, 0.5, 0.0], (899, 648, 1484, 946)),
}


def script_to_true_pressure(eps, y_script):
    rows_t, vals_t, crop = TICK_FITS[eps]
    A, B = np.polyfit(rows_t, vals_t, 1)
    h = crop[3] - crop[1]
    height = max(1.0, float(h - 1))
    # script: y = 1.5 - row/height*1.5  -> row = (1.5 - y)*height/1.5
    row = (1.5 - np.asarray(y_script)) * height / 1.5
    return A * row + B


def peaks_summary(t, v, label, unit, t0=185.0):
    """Pics par cycle sur la fenetre etablie (>t0 pour ecarter le transitoire)."""
    t = np.asarray(t); v = np.asarray(v)
    m = t >= t0
    t, v = t[m], v[m]
    peaks = []
    for i in range(1, len(v) - 1):
        if v[i] >= v[i - 1] and v[i] > v[i + 1] and v[i] > np.median(v):
            if not peaks or t[i] - peaks[-1][0] > 8.0:
                peaks.append((t[i], v[i]))
            elif v[i] > peaks[-1][1]:
                peaks[-1] = (t[i], v[i])
    pv = np.asarray([p[1] for p in peaks])
    print(f"    {label}: min={v.min():.1f} max={v.max():.1f} ; pics n={pv.size} "
          f"max={pv.max():.1f} moy={pv.mean():.1f} {unit}")
    return pv


def run_case(eps, targets, tag, dt=1.0, n_layers=3, n_phi=16):
    t0 = time.time()
    _, data = f7.run_model_for_eps(eps, targets=targets, dt=dt, n_layers=n_layers, n_phi=n_phi)
    win = f7.window_model_data(data)
    print(f"  [{tag}] eps={eps} dt={dt} layers={n_layers} phi={n_phi} ({time.time()-t0:.1f} s)")
    pf = peaks_summary(win["time"], win["force_total_mN"], "force", "mN")
    pt = peaks_summary(win["time"], win["torque_act_microNm"], "torque", "uNm")
    return pf, pt


print("Base:", f7.Base.__file__, getattr(f7.Base, "MODEL_VERSION", "?"))
targets = f7.digitize_figure7()

for eps in (0.8, 1.0):
    tp, pp = targets[eps]["pressure"]
    print(f"eps={eps}: pression script min={np.min(pp):.3f} max={np.max(pp):.3f}")

print("\n--- 1. BASELINE (pression numerisee biaisee, comme la validation) ---")
base = {}
for eps in (0.8, 1.0):
    base[eps] = run_case(eps, targets, "baseline")

print("\n--- 2. PRESSION CORRIGEE (calibration vraie des ticks) ---")
targets_corr = {e: dict(targets[e]) for e in targets}
for eps in (0.8, 1.0):
    tp, pp = targets[eps]["pressure"]
    pp_true = np.clip(script_to_true_pressure(eps, pp), 0.0, None)
    targets_corr[eps]["pressure"] = (tp, pp_true)
    print(f"eps={eps}: pression corrigee min={pp_true.min():.3f} max={pp_true.max():.3f}")
corr = {}
for eps in (0.8, 1.0):
    corr[eps] = run_case(eps, targets_corr, "pression corrigee")

print("\n--- 3. CONVERGENCE (pression corrigee, eps=1.0) ---")
run_case(1.0, targets_corr, "convergence", dt=0.5, n_layers=5, n_phi=32)
run_case(0.8, targets_corr, "convergence", dt=0.5, n_layers=5, n_phi=32)

print("\n--- Reference article (mesure pixel, unites vraies) ---")
print("  theorie force  (a): pics 1758..1806 (moy 1783) ; (b): 2389..2466 (moy 2427) mN")
print("  theorie torque (a): pics 591..733 (moy 664) ; (b): 702..876 (hors artefact) uNm")
print("  exp force (a): 1762..1843 ; (b): 2346..2440 mN (sommet de marqueur)")
print("  exp torque (a): 837..919 ; (b): 857..928 uNm (sommet de marqueur)")
