# -*- coding: utf-8 -*-
"""Comparaison : mises sous pression B/I/J avec 200 g vs modele bloque (precharge) + info suspendu."""
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
DT_ASSUMED = 0.2  # s, cadence supposee (identique aux fichiers dates: 0.2 s) — non documentee dans les CSV

def load_B():
    raw = pd.read_excel(DATA + r"\Mise sous pression muscle B 200 g.xlsx", header=0)
    p = pd.to_numeric(raw.iloc[:, 0], errors="coerce").values
    f = pd.to_numeric(raw.iloc[:, 1], errors="coerce").values
    m = np.isfinite(p) & np.isfinite(f)
    return p[m], f[m]

def load_csv(name):
    raw = pd.read_csv(DATA + "\\" + name, sep=";", decimal=",", header=0, usecols=[0, 1])
    p = pd.to_numeric(raw.iloc[:, 0], errors="coerce").values
    f = pd.to_numeric(raw.iloc[:, 1], errors="coerce").values
    m = np.isfinite(p) & np.isfinite(f)
    return p[m], f[m]

tests = {
    "B": load_B(),
    "I": load_csv("Mise sous pression muscle I 200 g.csv"),
    "J": load_csv("Mise sous pression muscle J 200 g.csv"),
}

results = {}
model_runs = {}
for name, (P_bar, F_N) in tests.items():
    t = np.arange(len(P_bar)) * DT_ASSUMED
    P_MPa = P_bar / 10.0
    t0 = time.time()
    model, arr = Base.run_blocked_actuation(pressure_time=t, pressure_MPa=P_MPa)
    wall = time.time() - t0
    model_runs[name] = (t, P_bar, F_N, arr)
    # zone de fin utile : pour I et J la fin (masse posee, F->0) n'est pas comparable
    F0_exp = float(np.mean(F_N[:8]))
    iPmax = int(np.argmax(P_bar))
    dF_exp_at_Pmax = float(F_N[iPmax] - F0_exp)
    dF_exp_max = float(np.max(F_N) - F0_exp)
    dF_mod_at_Pmax = float(arr["force_N"][iPmax] - arr["force_N"][0])
    dF_mod_max = float(np.max(arr["force_N"]) - arr["force_N"][0])
    dP = float(P_MPa[iPmax] - P_MPa[0])
    # relaxation pendant le palier : entre indice du max de force et la fin du palier haut
    results[name] = {
        "n_samples": len(P_bar),
        "duree_supposee_s": float(t[-1]),
        "preload_exp_N": F0_exp,
        "P_range_bar": [float(P_bar.min()), float(P_bar.max())],
        "dP_MPa": dP,
        "dF_exp_at_Pmax_N": dF_exp_at_Pmax,
        "dF_exp_max_N": dF_exp_max,
        "dF_model_at_Pmax_N": dF_mod_at_Pmax,
        "dF_model_max_N": dF_mod_max,
        "preload_model_N": float(arr["force_N"][0]),
        "slope_exp_N_per_MPa": dF_exp_max / dP if dP else None,
        "slope_model_N_per_MPa": dF_mod_max / dP if dP else None,
        "wall_s": round(wall, 1),
    }
    print(name, json.dumps(results[name]))

# ---------------- info : reponse libre avec 200 g suspendu (muscle B, meme pression) ----------------
tB, PB_bar, FB_N, arrB = model_runs["B"]
try:
    t0 = time.time()
    model_s, arr_s = Base.run_suspended_actuation(
        load_N=0.2 * 9.81, pressure_time=tB, pressure_MPa=PB_bar / 10.0
    )
    susp = {
        "free_disp_mm_at_end": float(arr_s["free_displacement_mm"][-1]),
        "free_disp_mm_max_abs": float(np.max(np.abs(arr_s["free_displacement_mm"]))),
        "free_actuation_percent_end": float(arr_s["free_actuation_percent"][-1]),
        "wall_s": round(time.time() - t0, 1),
    }
    print("suspendu:", json.dumps(susp))
except Exception as e:  # noqa
    arr_s = None
    susp = {"error": repr(e)}
    print("suspendu ECHEC:", repr(e))

with open(SCRATCH + r"\metrics_200g.json", "w") as fh:
    json.dump({"blocked": results, "suspended_B": susp, "dt_assumed_s": DT_ASSUMED}, fh, indent=2)

# ---------------- figure ----------------
fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
for k, name in enumerate(["B", "I", "J"]):
    t, P_bar, F_N, arr = model_runs[name]
    F0_exp = np.mean(F_N[:8])
    ax = axes[0, k]
    ax.plot(t, F_N - F0_exp, "r-", lw=1.2, label="exp  F - F0")
    ax.plot(arr["time"], arr["force_act_mN"] / 1000.0, "k--", lw=1.4, label="modele F_act")
    ax2 = ax.twinx()
    ax2.plot(t, P_bar, "b:", lw=1, alpha=0.6)
    ax2.set_ylabel("P (bar)", color="b", fontsize=8)
    ax.set_title(f"Muscle {name} — 200 g (dt suppose {DT_ASSUMED} s)")
    ax.set_xlabel("t (s)"); ax.set_ylabel("dF (N)")
    ax.legend(fontsize=8, loc="center right"); ax.grid(alpha=0.3)

    ax = axes[1, k]
    ax.plot(P_bar, F_N - F0_exp, "r.", ms=3, label="exp")
    ax.plot(10 * arr["pressure_MPa"], arr["force_act_mN"] / 1000.0, "k--", lw=1.4, label="modele")
    ax.set_xlabel("P (bar)"); ax.set_ylabel("dF (N)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title("dF vs P")

fig.suptitle("Mises sous pression 200 g (B, I, J) : force mesuree vs modele bloque avec precharge (defauts article)", fontsize=12)
fig.tight_layout()
fig.savefig(SCRATCH + r"\fig_mise_sous_pression_BIJ.png", dpi=130)
print("figure sauvee")

if arr_s is not None:
    fig2, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(arr_s["time"], arr_s["free_displacement_mm"], "m-", lw=1.5)
    ax.set_xlabel("t (s)"); ax.set_ylabel("deplacement libre (mm)")
    ax.set_title("Modele suspendu 200 g, pression du test B (info : signe = allongement>0)")
    ax.grid(alpha=0.3)
    ax2 = ax.twinx(); ax2.plot(tB, PB_bar, "b:", lw=1, alpha=0.6); ax2.set_ylabel("P (bar)", color="b")
    fig2.tight_layout()
    fig2.savefig(SCRATCH + r"\fig_suspendu_B_info.png", dpi=130)
    print("figure suspendu sauvee")
