# -*- coding: utf-8 -*-
"""Relaxation de l'increment de force pendant le palier haut du test B (exp vs modele)."""
import sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

DATA = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\expérimentale\Sacha"
DT = 0.2

raw = pd.read_excel(DATA + r"\Mise sous pression muscle B 200 g.xlsx", header=0)
P_bar = pd.to_numeric(raw.iloc[:, 0], errors="coerce").values
F_N = pd.to_numeric(raw.iloc[:, 1], errors="coerce").values
m = np.isfinite(P_bar) & np.isfinite(F_N)
P_bar, F_N = P_bar[m], F_N[m]
t = np.arange(len(P_bar)) * DT

model, arr = Base.run_blocked_actuation(pressure_time=t, pressure_MPa=P_bar / 10.0)

F0_exp = float(np.mean(F_N[:8]))
dF_exp = F_N - F0_exp
dF_mod = arr["force_act_mN"] / 1000.0

i_peak_exp = int(np.argmax(dF_exp))
i_peak_mod = int(np.argmax(dF_mod))
out = {
    "t_peak_exp_s": float(t[i_peak_exp]),
    "dF_exp_peak_N": float(dF_exp[i_peak_exp]),
    "dF_exp_end_N": float(np.mean(dF_exp[-5:])),
    "perte_exp_frac_du_pic": float(1 - np.mean(dF_exp[-5:]) / dF_exp[i_peak_exp]),
    "t_peak_mod_s": float(arr["time"][i_peak_mod]),
    "dF_mod_peak_N": float(dF_mod[i_peak_mod]),
    "dF_mod_end_N": float(np.mean(dF_mod[-5:])),
    "perte_mod_frac_du_pic": float(1 - np.mean(dF_mod[-5:]) / dF_mod[i_peak_mod]),
    "duree_palier_apres_pic_s": float(t[-1] - t[i_peak_mod]),
}
print(json.dumps(out, indent=2))
