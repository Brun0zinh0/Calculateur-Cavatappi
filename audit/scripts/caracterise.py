# -*- coding: utf-8 -*-
"""Caracterisation des essais force/pression de Bruno."""
import numpy as np
import pandas as pd

BRUNO = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\expérimentale\Bruno"

for name in ["essai_force_pression_20260730_153353.csv", "essai_force_pression_20260730_155135.csv"]:
    df = pd.read_csv(BRUNO + "\\" + name)
    t = df["time_s"].to_numpy()
    p_bar = df["pressure_bar"].to_numpy()
    f = df["force_mN"].to_numpy()
    dt = np.diff(t)
    print("=" * 70)
    print(name)
    print(f"  n={len(t)}, duree={t[-1]-t[0]:.1f} s, dt median={np.median(dt)*1000:.0f} ms, dt max={dt.max():.3f} s")
    print(f"  pression: min={p_bar.min():.3f} bar, max={p_bar.max():.3f} bar ({p_bar.max()*0.1:.3f} MPa)")
    print(f"  force: init={f[:20].mean():.1f} mN, min={f.min():.1f}, max={f.max():.1f} mN, gain max={f.max()-f[:20].mean():+.1f} mN")
    # detection de paliers : derivee de pression lissee
    p_s = pd.Series(p_bar).rolling(15, center=True, min_periods=1).mean().to_numpy()
    dpdt = np.gradient(p_s, t)
    plateau = np.abs(dpdt) < 0.02  # bar/s
    # segments de plateau > 4 s
    segs = []
    i = 0
    while i < len(t):
        if plateau[i]:
            j = i
            while j + 1 < len(t) and plateau[j + 1]:
                j += 1
            if t[j] - t[i] > 4.0:
                segs.append((i, j))
            i = j + 1
        else:
            i += 1
    print(f"  paliers detectes (>4 s a dP/dt<0.02 bar/s): {len(segs)}")
    for (i, j) in segs:
        fs = f[i:j + 1]
        ts = t[i:j + 1]
        # pente de relaxation sur le palier (regression lineaire)
        A = np.polyfit(ts, fs, 1)
        print(f"    t=[{t[i]:7.1f},{t[j]:7.1f}] s  P~{p_bar[i:j+1].mean():5.2f} bar  F deb={fs[0]:7.1f} fin={fs[-1]:7.1f} mN  pente={A[0]*1000:+7.2f} mN/ks ({A[0]:+.4f} mN/s)")
    # correlation force-pression
    c = np.corrcoef(p_bar, f)[0, 1]
    print(f"  correlation P-F: {c:.4f}")
    # nb de cycles (montees au dessus de 65% de l'amplitude)
    thr = p_bar.min() + 0.65 * np.ptp(p_bar)
    above = p_bar > thr
    ncross = int(np.sum(~above[:-1] & above[1:]))
    print(f"  montees au-dessus de 65% de Pmax: {ncross}")
