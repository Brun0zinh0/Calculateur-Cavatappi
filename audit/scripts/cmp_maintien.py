# -*- coding: utf-8 -*-
"""Comparaison : maintien sous pression muscle D (20 min) vs modele bloque Base.py."""
import sys, time, json
import numpy as np
import pandas as pd

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRATCH = r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
DATA = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\expérimentale\Sacha"

# ---------------- donnees experimentales ----------------
raw = pd.read_excel(DATA + r"\Maintient sous pression D 20 mn.xlsx", header=0)
P_bar = pd.to_numeric(raw.iloc[:, 0], errors="coerce").values
F_N = pd.to_numeric(raw.iloc[:, 1], errors="coerce").values
t_exp = pd.to_numeric(raw.iloc[:, 3], errors="coerce").values
m = np.isfinite(P_bar) & np.isfinite(F_N) & np.isfinite(t_exp)
P_bar, F_N, t_exp = P_bar[m], F_N[m], t_exp[m]
P_MPa = P_bar / 10.0
print(f"exp: n={len(t_exp)}, t=[{t_exp[0]},{t_exp[-1]}] s, P=[{P_bar.min():.3f},{P_bar.max():.3f}] bar, F=[{F_N.min():.3f},{F_N.max():.3f}] N")

# ---------------- fit bi-exponentiel de la force exp ----------------
def biexp(t, A0, A1, tau1, A2, tau2):
    return A0 + A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2)

p0 = [F_N[-1], 0.1, 60.0, 0.15, 800.0]
popt_F, _ = curve_fit(biexp, t_exp, F_N, p0=p0, maxfev=40000)
popt_P, _ = curve_fit(biexp, t_exp, P_bar, p0=[P_bar[-1], 1.0, 60.0, 1.5, 800.0], maxfev=40000)
taus_F = sorted([popt_F[2], popt_F[4]])
taus_P = sorted([popt_P[2], popt_P[4]])
print("fit F exp: A0=%.3f  A1=%.3f tau=%.0f s  A2=%.3f tau=%.0f s" % (popt_F[0], popt_F[1], popt_F[2], popt_F[3], popt_F[4]))
print("fit P exp: A0=%.3f  A1=%.3f tau=%.0f s  A2=%.3f tau=%.0f s" % (popt_P[0], popt_P[1], popt_P[2], popt_P[3], popt_P[4]))

# ---------------- run 1 : modele bloque pilote par la pression MESUREE ----------------
DT = 4.0
t_ds = np.arange(0.0, t_exp[-1] + 1e-9, DT)
P_ds = np.interp(t_ds, t_exp, P_MPa)
t0 = time.time()
model1, arr1 = Base.run_blocked_actuation(pressure_time=t_ds, pressure_MPa=P_ds)
print("run mesure-P: %.0f s wall, %d pts" % (time.time() - t0, len(arr1["time"])))

# ---------------- run 2 : modele bloque a pression CONSTANTE (relaxation pure) ----------------
t0 = time.time()
model2, arr2 = Base.run_hold_relaxation(P_hold=float(P_MPa[0]), hold_time=float(t_exp[-1]), ramp_time=8.0, dt=DT)
print("run P-constante: %.0f s wall" % (time.time() - t0))
i0 = int(arr2["hold_start_index"])
th = arr2["time"][i0:] - arr2["time"][i0]
Fh = arr2["force_N"][i0:]
popt_M, _ = curve_fit(biexp, th, Fh, p0=[Fh[-1], 0.05, 60.0, 0.05, 1500.0], maxfev=40000)
taus_M = sorted([popt_M[2], popt_M[4]])
print("fit F modele (P cst): A0=%.3f A1=%.3f tau=%.0f  A2=%.3f tau=%.0f" % tuple(popt_M))

# ---------------- metriques ----------------
mx = Base.default_maxwell_tensile_params()
taus_branches = [mx.eta1 / mx.E1, mx.eta2 / mx.E2, mx.eta3 / mx.E3]

met = {
    "exp": {
        "F0_N": float(F_N[0]), "Fend_N": float(F_N[-1]),
        "decay_frac": float(1 - F_N[-1] / F_N[0]),
        "P0_bar": float(P_bar[0]), "Pend_bar": float(P_bar[-1]),
        "P_decay_frac": float(1 - P_bar[-1] / P_bar[0]),
        "fit_taus_F_s": [float(x) for x in taus_F],
        "fit_taus_P_s": [float(x) for x in taus_P],
        "corr_F_P": float(np.corrcoef(P_bar, F_N)[0, 1]),
    },
    "model_measuredP": {
        "F0_N": float(arr1["force_N"][0]), "Fend_N": float(arr1["force_N"][-1]),
        "decay_frac": float(1 - arr1["force_N"][-1] / arr1["force_N"][0]),
    },
    "model_constP": {
        "F_holdstart_N": float(Fh[0]), "Fend_N": float(Fh[-1]),
        "decay_frac": float(1 - Fh[-1] / Fh[0]),
        "fit_taus_s": [float(x) for x in taus_M],
    },
    "maxwell_branch_taus_s": [float(x) for x in taus_branches],
    "dt_model_s": DT,
}
with open(SCRATCH + r"\metrics_maintien.json", "w") as fh:
    json.dump(met, fh, indent=2)
print(json.dumps(met, indent=2))

# ---------------- figure ----------------
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
ax = axes[0, 0]
ax.plot(t_exp, P_bar, "b-", lw=1)
ax.set_xlabel("t (s)"); ax.set_ylabel("Pression (bar)", color="b")
ax2 = ax.twinx(); ax2.plot(t_exp, F_N, "r-", lw=1)
ax2.set_ylabel("Force (N)", color="r")
ax.set_title("Essai : maintien D (P et F mesures)")

ax = axes[0, 1]
ax.plot(t_exp, F_N / F_N[0], "r-", lw=1.2, label="exp F/F0")
ax.plot(arr1["time"], arr1["force_N"] / arr1["force_N"][0], "k--", lw=1.5, label="modele (P mesuree)")
ax.plot(arr2["time"][i0:] - arr2["time"][i0], Fh / Fh[0], "g-.", lw=1.5, label="modele (P constante)")
ax.set_xlabel("t (s)"); ax.set_ylabel("F/F0"); ax.legend(); ax.grid(alpha=0.3)
ax.set_title("Force normalisee : exp vs modele")

ax = axes[1, 0]
ax.plot(P_bar, F_N, "r.", ms=2, label="exp")
ax.plot(10 * arr1["pressure_MPa"], arr1["force_N"], "k--", lw=1.5, label="modele (P mesuree)")
ax.set_xlabel("P (bar)"); ax.set_ylabel("F (N)"); ax.legend(); ax.grid(alpha=0.3)
ax.set_title("Trajectoire F-P pendant le maintien")

ax = axes[1, 1]
resid_exp = F_N - biexp(t_exp, *popt_F)
ax.semilogx(t_exp[t_exp > 0], (F_N[t_exp > 0] - popt_F[0]) / (F_N[0] - popt_F[0]), "r-", lw=1, label="exp (F-A0)/(F0-A0)")
ax.semilogx(th[th > 0], (Fh[th > 0] - popt_M[0]) / (Fh[0] - popt_M[0]), "g-.", lw=1.5, label="modele P cst")
for tb in taus_branches:
    ax.axvline(tb, color="gray", ls=":", lw=1)
ax.text(taus_branches[0], 1.02, "tau1", fontsize=8); ax.text(taus_branches[1], 1.02, "tau2", fontsize=8); ax.text(taus_branches[2], 1.02, "tau3", fontsize=8)
ax.set_xlabel("t (s, log)"); ax.set_ylabel("decroissance normalisee"); ax.legend(); ax.grid(alpha=0.3)
ax.set_title("Echelles de temps (traits: tau branches Maxwell)")

fig.suptitle("Maintien sous pression D 20 min vs modele Base.py (parametres article)", fontsize=13)
fig.tight_layout()
fig.savefig(SCRATCH + r"\fig_maintien_D.png", dpi=130)
print("figure sauvee")
