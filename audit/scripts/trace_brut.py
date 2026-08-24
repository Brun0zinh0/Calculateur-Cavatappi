# -*- coding: utf-8 -*-
"""Trace brut des essais de Bruno."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BRUNO = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\expérimentale\Bruno"
OUT = r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"

FILES = [
    "essai_force_pression_20260730_153353.csv",
    "essai_force_pression_20260730_155135.csv",
]

fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex="col")
for k, f in enumerate(FILES):
    df = pd.read_csv(f"{BRUNO}\\{f}")
    t = df["time_s"].to_numpy()
    p = df["pressure_bar"].to_numpy() * 0.1
    F = df["force_mN"].to_numpy()
    ax = axes[0, k]
    ax.plot(t, p, lw=0.8, color="tab:blue")
    ax.set_ylabel("Pression (MPa)")
    ax.set_title(f.replace("essai_force_pression_", "").replace(".csv", ""))
    ax.grid(alpha=0.3)
    ax2 = axes[1, k]
    ax2.plot(t, F, lw=0.8, color="tab:red")
    ax2.set_ylabel("Force (mN)")
    ax2.set_xlabel("Temps (s)")
    ax2.grid(alpha=0.3)
fig.suptitle("Essais bruts Bruno : pression et force bloquée")
fig.tight_layout()
fig.savefig(f"{OUT}\\essais_bruts.png", dpi=130)
print("ok essais_bruts.png")
