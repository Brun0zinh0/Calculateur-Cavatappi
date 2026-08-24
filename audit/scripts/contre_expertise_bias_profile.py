# -*- coding: utf-8 -*-
"""Contre-expertise du constat 'profil radial d'angle de biais par defaut'.

1) Verification analytique des ecarts d'angle lineaire vs arctan.
2) Run complet du moteur reel (alpha V2/Base.py) avec les deux profils,
   comparaison des sorties force / couple.
NE MODIFIE AUCUN FICHIER DU PROJET.
"""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

# ---------- 1) Ecarts d'angle ----------
theta_f = np.deg2rad(37.91)
Rin, Rout, n = 0.4, 1.0, 4
edges = np.linspace(Rin, Rout, n + 1)
centers = 0.5 * (edges[:-1] + edges[1:])
x = centers / Rout
lin = np.rad2deg(x * theta_f)
arc = np.rad2deg(np.arctan(x * np.tan(theta_f)))
print("centres R/Rout :", np.round(x, 4))
print("theta lineaire :", np.round(lin, 2))
print("theta arctan   :", np.round(arc, 2))
print("ecart (lin-arc):", np.round(lin - arc, 2))

# ---------- 2) Runs comparatifs ----------
def run(profile, n_layers=None, n_cycles=2):
    geom = Base.default_geometry_params(bias_angle_profile=profile)
    kwargs = dict(eps=0.8, n_cycles=n_cycles, Pmax=1.3, geom=geom)
    if n_layers is not None:
        kwargs["n_layers"] = n_layers
    model, arr = Base.run_blocked_actuation(**kwargs)
    return arr

def compare(tag, a_lin, a_arc):
    print(f"\n----- {tag} -----")
    keys = sorted(set(a_lin.keys()) & set(a_arc.keys()))
    for k in keys:
        v1, v2 = np.asarray(a_lin[k], float), np.asarray(a_arc[k], float)
        if v1.shape != v2.shape or v1.size == 0:
            continue
        if any(s in k for s in ("torque", "force", "Force", "Torque")):
            m1, m2 = np.max(np.abs(v1)), np.max(np.abs(v2))
            # amplitude d'actionnement = max - min sur la phase d'actionnement
            r1 = np.max(v1) - np.min(v1)
            r2 = np.max(v2) - np.min(v2)
            rel_max = 100.0 * (m1 - m2) / m2 if m2 else float("nan")
            rel_rng = 100.0 * (r1 - r2) / r2 if r2 else float("nan")
            print(f"{k:42s} max|lin|={m1:12.5g} max|arc|={m2:12.5g} ecart={rel_max:7.2f}%  "
                  f"ampl lin={r1:12.5g} arc={r2:12.5g} ecart={rel_rng:7.2f}%")

for nl, tag in ((4, "n_layers=4 (defaut interface)"), (12, "n_layers=12")):
    a_lin = run("paper_linear", n_layers=nl)
    a_arc = run("uniform_twist", n_layers=nl)
    compare(tag, a_lin, a_arc)
