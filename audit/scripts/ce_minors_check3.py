# -*- coding: utf-8 -*-
"""Contre-expertise figure 7 : chemin PDF, numerisation (pics, occlusion), jonction t=180 s."""
import sys
from pathlib import Path
import numpy as np

WORKSPACE = Path(r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail")
sys.path.insert(0, str(WORKSPACE))

import validation_figure7 as VF

print("=== chemin PDF ===")
print(f"PDF_PATH attendu par le script : {VF.PDF_PATH}")
print(f"existe ? {VF.PDF_PATH.exists()}")
real_pdf = WORKSPACE / "Article support" / VF.PDF_FILENAME
print(f"chemin reel : {real_pdf} ; existe ? {real_pdf.exists()}")
try:
    VF.extract_figure7_image()
    print("extract_figure7_image() : OK ??")
except FileNotFoundError as exc:
    print(f"extract_figure7_image() -> FileNotFoundError : {exc}")

print()
print("=== quel moteur est importe par validation_figure7 ? ===")
print(f"module Base importe : {VF.Base.__file__}")
print(f"version : {VF.Base.MODEL_VERSION}")

print()
print("=== numerisation avec le PDF reel ===")
image = VF.extract_figure7_image(real_pdf)
print(f"image figure 7 : {image.size}")
targets = VF.digitize_figure7(image)
theory = VF.digitize_figure7_theory(image)
for eps in (0.8, 1.0):
    pt, pp = targets[eps]["pressure"]
    ft, fv = targets[eps]["force"]
    tt, tv = targets[eps]["torque"]
    tf_t, tf_v = theory[eps]["force"]
    ttq_t, ttq_v = theory[eps]["torque"]
    print(f"eps={eps} : P num max = {np.nanmax(pp):.3f} MPa (texte article : 1.41-1.43)")
    print(f"          F exp num max = {np.nanmax(fv):.1f} mN ; T exp num max = {np.nanmax(tv):.1f} uNm")
    print(f"          F theorie num max = {np.nanmax(tf_v):.1f} mN ; T theorie num max = {np.nanmax(ttq_v):.1f} uNm")
    print(f"          (ancrages texte 11e cycle : F 1824.60/2424.67 ; T 877.01/873.41)")

print()
print("=== jonction prefixe genere / pression numerisee a t=180 s ===")
for eps in (0.8, 1.0):
    t, p = VF.digitized_pressure_history(eps, targets)
    i180 = np.searchsorted(t, VF.PAPER_T_MIN)
    print(f"eps={eps} : P(prefixe, t={t[i180-1]:.1f})={p[i180-1]:.3f} MPa -> "
          f"P(figure, t={t[i180]:.1f})={p[i180]:.3f} MPa ; saut={p[i180]-p[i180-1]:+.3f} MPa")
    # espacement des pics numerises
    pt, pp = targets[eps]["pressure"]
    thr = np.nanmin(pp) + 0.65 * np.nanptp(pp)
    peaks = []
    for i in range(1, len(pp) - 1):
        if pp[i] >= pp[i-1] and pp[i] > pp[i+1] and pp[i] >= thr:
            peaks.append(pt[i])
    if len(peaks) > 2:
        print(f"          periode apparente des pics numerises = {np.median(np.diff(peaks)):.2f} s (protocole nominal 18 s)")
