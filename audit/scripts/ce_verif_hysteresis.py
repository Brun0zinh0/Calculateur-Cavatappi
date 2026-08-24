# -*- coding: utf-8 -*-
"""Contre-expertise : hysteresis modele (signe et amplitude) + correlations.

1. Recalcul depuis les sorties sauvegardees essai2 (_sim.csv / _el.csv) :
   delta descente-montee du modele complet et quasi-elastique, residu final.
2. Verification independante : cycle triangulaire 0->0.58->0 MPa en 100 s,
   modele complet grossier (n_layers=2, n_phi=8) : signe de la boucle.
3. Correlations sim/mesure recalculees depuis les fichiers sauvegardes.
"""
import sys
import numpy as np
import pandas as pd

ALPHA = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2"
OUT = r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
sys.path.insert(0, ALPHA)
import Base as modele  # noqa: E402

print("=== 1. Hysteresis essai 2 depuis les sorties sauvegardees ===")
sim = pd.read_csv(OUT + r"\essai2_155135_sim.csv")
el = pd.read_csv(OUT + r"\essai2_155135_el.csv")
t = sim["time_s"].to_numpy(float)
p = sim["pressure_MPa"].to_numpy(float)
fm = sim["force_act_meas_mN"].to_numpy(float)
fv = sim["force_act_sim_mN"].to_numpy(float)
fe = el["force_act_el_mN"].to_numpy(float)
k_v = float(np.dot(fv, fm) / np.dot(fv, fv))
k_e = float(np.dot(fe, fm) / np.dot(fe, fe))
i_pk = int(np.argmax(p))
for p_probe in (0.2, 0.3, 0.4):
    up, dn = slice(0, i_pk + 1), slice(i_pk, None)
    o_up, o_dn = np.argsort(p[up]), np.argsort(p[dn])
    hyst = {}
    for lbl, f, k in (("mesure", fm, 1.0), ("maxwell", fv, k_v), ("quasi-el", fe, k_e)):
        fu = float(np.interp(p_probe, p[up][o_up], (k * f)[up][o_up]))
        fd = float(np.interp(p_probe, p[dn][o_dn], (k * f)[dn][o_dn]))
        hyst[lbl] = fd - fu
    print(f"  P={p_probe:.1f} MPa : mesure {hyst['mesure']:+7.1f} | maxwell(x{k_v:.2f}) {hyst['maxwell']:+7.1f} | quasi-el(x{k_e:.2f}) {hyst['quasi-el']:+7.1f} mN")
print(f"  residu final (P={p[-1]*1000:.0f} kPa) : mesure {fm[-1]:+.1f} | maxwell brut {fv[-1]:+.1f} ({100*fv[-1]/np.max(fv):+.0f}% de son pic) | quasi-el brut {fe[-1]:+.1f} mN")
print(f"  correlations brutes : maxwell {np.corrcoef(fv, fm)[0,1]:.4f} | quasi-el {np.corrcoef(fe, fm)[0,1]:.4f}")

sim1 = pd.read_csv(OUT + r"\essai1_153353_sim.csv")
el1 = pd.read_csv(OUT + r"\essai1_153353_el.csv")
fm1 = sim1["force_act_meas_mN"].to_numpy(float)
fv1 = sim1["force_act_sim_mN"].to_numpy(float)
fe1 = el1["force_act_el_mN"].to_numpy(float)
print(f"  essai 1 correlations : maxwell {np.corrcoef(fv1, fm1)[0,1]:.4f} | quasi-el {np.corrcoef(fe1, fm1)[0,1]:.4f}")

print("\n=== 2. Cycle triangulaire independant (modele complet grossier) ===")
t2 = np.arange(0.0, 101.0, 1.0)
p2 = np.where(t2 <= 50.0, 0.58 * t2 / 50.0, 0.58 * (100.0 - t2) / 50.0)
p2 = np.maximum(p2, 0.0)
_, arr = modele.run_blocked_actuation(pressure_time=t2, pressure_MPa=p2, n_layers=2, n_phi=8, pre_steps=16)
pr, fr, tr = arr["pressure_MPa"], arr["force_act_mN"], arr["time"]
i_pk2 = int(np.argmax(pr))
for p_probe in (0.2, 0.3, 0.4):
    up, dn = slice(0, i_pk2 + 1), slice(i_pk2, None)
    o_up, o_dn = np.argsort(pr[up]), np.argsort(pr[dn])
    fu = float(np.interp(p_probe, pr[up][o_up], fr[up][o_up]))
    fd = float(np.interp(p_probe, pr[dn][o_dn], fr[dn][o_dn]))
    print(f"  P={p_probe:.1f} MPa : montee {fu:+7.2f} descente {fd:+7.2f} delta={fd-fu:+7.2f} mN (brut)")
print(f"  pic {np.max(fr):.2f} mN ; residu final (P=0) {fr[-1]:+.2f} mN = {100*fr[-1]/np.max(fr):+.0f}% du pic")
print("\nTermine.")
