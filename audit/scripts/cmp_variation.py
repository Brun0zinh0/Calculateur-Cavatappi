# -*- coding: utf-8 -*-
"""Variation rapide de position (D) : trace exp + demo modele d'un echelon de position a P constante."""
import sys, time, json
import numpy as np
import pandas as pd

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRATCH = r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
DATA = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\expérimentale\Sacha"

raw = pd.read_excel(DATA + r"\Variation rappide de position (D).xlsx", header=0)
P_bar = pd.to_numeric(raw.iloc[:, 0], errors="coerce").values
F_N = pd.to_numeric(raw.iloc[:, 1], errors="coerce").values
t_exp = pd.to_numeric(raw.iloc[:, 3], errors="coerce").values
m = np.isfinite(P_bar) & np.isfinite(F_N) & np.isfinite(t_exp)
P_bar, F_N, t_exp = P_bar[m], F_N[m], t_exp[m]
print(f"exp: n={len(t_exp)}, t=[{t_exp[0]},{t_exp[-1]}] s, P=[{P_bar.min():.2f},{P_bar.max():.2f}] bar, F=[{F_N.min():.3f},{F_N.max():.3f}] N")

# ---------------- demo modele : echelons de position a P constante ----------------
# Construire un modele bloque, precontraint eps=0.8, pressurise a P moyenne de l'essai,
# puis imposer des echelons de h_target (etirement +2%, retour, -1%) et observer la force.
cfg = Base.default_simulation_config()
disc = Base.default_discretization(n_layers=cfg.n_layers, n_phi=cfg.n_phi, pre_steps=cfg.pre_steps, dw_bracket=(-0.05, 0.05))
model = Base.TCPAMaxwellBlockedModel(mat=cfg.mat, geom=cfg.geom, disc=disc, integration="exponential")
model.prestretch_to(0.8)
P0 = float(np.mean(P_bar)) / 10.0  # ~0.64 MPa
model.step(P0, 0.0, h_target=model.h_blocked)

h_b = model.h_blocked
dt = 0.25
seq = []  # (h_target, duration)
seq += [(h_b, 2.0)]
seq += [(1.02 * h_b, 4.0)]   # etirement rapide +2 %
seq += [(h_b, 4.0)]          # retour
seq += [(0.99 * h_b, 4.0)]   # relachement -1 %
seq += [(h_b, 6.0)]
t_m, F_m, h_m = [0.0], [model.helix.time], []
t_acc = 0.0
times, forces, targets = [], [], []
t0 = time.time()
for h_tgt, dur in seq:
    n = max(1, int(round(dur / dt)))
    for i in range(n):
        # premier pas de l'echelon : changement de h instantane (dt petit)
        res = model.step(P0, dt, h_target=h_tgt)
        t_acc += dt
        times.append(t_acc); forces.append(res.Ft); targets.append(h_tgt / h_b)
print("demo echelons: %.0f s wall, %d pas" % (time.time() - t0, len(times)))
times = np.array(times); forces = np.array(forces); targets = np.array(targets)

met = {
    "exp_F_range_N": [float(F_N.min()), float(F_N.max())],
    "exp_P_mean_bar": float(np.mean(P_bar)),
    "model_F_after_pressurisation_N": float(forces[0]),
    "model_F_peak_stretch2pct_N": float(forces.max()),
    "model_F_min_N": float(forces.min()),
}
with open(SCRATCH + r"\metrics_variation.json", "w") as fh:
    json.dump(met, fh, indent=2)
print(json.dumps(met, indent=2))

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
ax = axes[0]
ax.plot(t_exp, F_N, "r-", lw=1.2)
ax.set_xlabel("t (s)"); ax.set_ylabel("Force (N)", color="r")
ax2 = ax.twinx(); ax2.plot(t_exp, P_bar, "b:", lw=1, alpha=0.7); ax2.set_ylabel("P (bar)", color="b")
ax.set_title("Essai : variation rapide de position (D), P ~ 6.4 bar")
ax.grid(alpha=0.3)

ax = axes[1]
ax.plot(times, forces, "k-", lw=1.4, label="force modele")
ax2 = ax.twinx(); ax2.plot(times, targets, "g--", lw=1, alpha=0.7); ax2.set_ylabel("h/h_bloque", color="g")
ax.set_xlabel("t (s)"); ax.set_ylabel("Force (N)")
ax.set_title("Modele : echelons de position a P=%.2f MPa (demo qualitative)" % P0)
ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(SCRATCH + r"\fig_variation_rapide_D.png", dpi=130)
print("figure sauvee")
