# -*- coding: utf-8 -*-
"""Caracterisation des essais force/pression de Bruno (actionnement bloque)."""
import numpy as np
import pandas as pd

BRUNO = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\expérimentale\Bruno"

FILES = [
    "essai_force_pression_20260730_153353.csv",
    "essai_force_pression_20260730_155135.csv",
]


def plateaus(t, p, tol=0.02, min_dur=3.0):
    """Detecte les paliers de pression (dp/dt ~ 0 sur duree > min_dur)."""
    dp = np.gradient(p, t)
    flat = np.abs(dp) < tol  # MPa/s
    segs = []
    i = 0
    n = len(t)
    while i < n:
        if flat[i]:
            j = i
            while j + 1 < n and flat[j + 1]:
                j += 1
            if t[j] - t[i] >= min_dur:
                segs.append((i, j))
            i = j + 1
        else:
            i += 1
    return segs


for f in FILES:
    df = pd.read_csv(f"{BRUNO}\\{f}")
    t = df["time_s"].to_numpy()
    p_bar = df["pressure_bar"].to_numpy()
    p_mpa = p_bar * 0.1
    F = df["force_mN"].to_numpy()
    Fu = df["force_unfiltered_mN"].to_numpy()
    print("=" * 78)
    print(f)
    print(f"  points={len(t)}  duree={t[-1]-t[0]:.1f} s  dt median={np.median(np.diff(t))*1000:.0f} ms")
    print(f"  pression: min={p_mpa.min():.4f} max={p_mpa.max():.4f} MPa (={p_bar.max():.2f} bar)")
    print(f"  force: init={F[0]:.1f} min={F.min():.1f} max={F.max():.1f} mN | delta max={F.max()-F[0]:.1f} mN")
    print(f"  bruit force (std brut - filtre) = {np.std(Fu - F):.2f} mN")
    # correlation force-pression
    c = np.corrcoef(p_mpa, F)[0, 1]
    print(f"  corr(P, F) = {c:.4f}")
    # paliers
    segs = plateaus(t, p_mpa)
    print(f"  paliers detectes (>3 s, |dP/dt|<0.02 MPa/s): {len(segs)}")
    for (i, j) in segs:
        dur = t[j] - t[i]
        pmean = p_mpa[i:j + 1].mean()
        dF = F[j] - F[i]
        # pente de relaxation en fin de palier (regression sur le segment)
        if j - i > 5:
            slope = np.polyfit(t[i:j + 1], F[i:j + 1], 1)[0]
        else:
            slope = np.nan
        print(f"    t=[{t[i]:7.1f},{t[j]:7.1f}] s ({dur:5.1f} s)  P={pmean:6.3f} MPa  "
              f"F: {F[i]:7.1f}->{F[j]:7.1f} mN  dF={dF:+7.1f}  pente={slope:+.2f} mN/s")
    # cycles ?
    from numpy import ptp
    thr = p_mpa.min() + 0.65 * ptp(p_mpa)
    peaks = np.where((p_mpa[1:-1] >= p_mpa[:-2]) & (p_mpa[1:-1] > p_mpa[2:]) & (p_mpa[1:-1] >= thr))[0] + 1
    print(f"  pics de pression au dessus de 65%: {len(peaks)} -> t = {np.round(t[peaks],1)}")
