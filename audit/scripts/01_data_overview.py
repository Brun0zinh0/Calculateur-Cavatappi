# -*- coding: utf-8 -*-
"""Inspection des 5 fichiers experimentaux de Sacha : trace brut P(t) et F(t)."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/expérimentale/Sacha"
OUT = r"C:/Users/B3DCB~1.PER/AppData/Local/Temp/claude/C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail/c5852e10-6d38-4703-b218-a1649a5b12c9/scratchpad"


def load_xlsx(name, dt_fallback=None):
    df = pd.read_excel(BASE + "/" + name, header=0)
    P = df.iloc[:, 0].astype(float).values
    F = df.iloc[:, 1].astype(float).values
    # periode d echantillonnage stockee en en-tete de colonne 5 si presente
    dt = None
    for c in df.columns:
        if isinstance(c, (int, float)) and not isinstance(c, bool):
            dt = float(c)
    if dt is None:
        dt = dt_fallback
    return P, F, dt


def load_csv(name):
    df = pd.read_csv(BASE + "/" + name, sep=";", decimal=",", usecols=[0, 1])
    return df.iloc[:, 0].astype(float).values, df.iloc[:, 1].astype(float).values


datasets = []
P, F, dt = load_xlsx("Maintient sous pression D 20 mn.xlsx")
datasets.append(("Maintien D 20 min (dt=%.2g s)" % dt, np.arange(len(P)) * dt, P, F))
P, F, dt = load_xlsx("Mise sous pression muscle B 200 g.xlsx", dt_fallback=0.2)
datasets.append(("Mise sous pression B 200 g (dt suppose 0.2 s)", np.arange(len(P)) * 0.2, P, F))
P, F = load_csv("Mise sous pression muscle I 200 g.csv")
datasets.append(("Mise sous pression I 200 g (dt suppose 0.2 s)", np.arange(len(P)) * 0.2, P, F))
P, F = load_csv("Mise sous pression muscle J 200 g.csv")
datasets.append(("Mise sous pression J 200 g (dt suppose 0.2 s)", np.arange(len(P)) * 0.2, P, F))
P, F, dt = load_xlsx("Variation rappide de position (D).xlsx")
datasets.append(("Variation rapide de position D (dt=%.2g s)" % dt, np.arange(len(P)) * dt, P, F))

fig, axes = plt.subplots(len(datasets), 1, figsize=(10, 14), sharex=False)
for ax, (title, t, P, F) in zip(axes, datasets):
    ax2 = ax.twinx()
    l1, = ax.plot(t, F, "b-", lw=1.0, label="Force (N)")
    l2, = ax2.plot(t, P, "r-", lw=0.8, alpha=0.7, label="Pression (bar)")
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("F (N)", color="b")
    ax2.set_ylabel("P (bar)", color="r")
    ax.set_xlabel("t (s)")
    ax.legend(handles=[l1, l2], loc="best", fontsize=8)
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(OUT + "/fig1_data_overview.png", dpi=130)
print("saved fig1_data_overview.png")

poids = 0.200 * 9.81
for title, t, P, F in datasets:
    print("\n===", title)
    print(f"  duree={t[-1]:.1f} s  n={len(t)}")
    print(f"  P: {P.min():.3f}..{P.max():.3f} bar | F: {F.min():.3f}..{F.max():.3f} N")
    if "200 g" in title:
        base_mask = P < 0.5
        Fb = np.median(F[:20])
        print(f"  F baseline (P basse) = {Fb:.3f} N vs poids 200 g = {poids:.3f} N (ecart {100*(Fb-poids)/poids:+.1f} %)")
        print(f"  Delta F pic = {F.max()-Fb:.3f} N pour Delta P = {P.max()-P.min():.3f} bar -> {(F.max()-Fb)/((P.max()-P.min())/10.0):.3f} N/MPa")
