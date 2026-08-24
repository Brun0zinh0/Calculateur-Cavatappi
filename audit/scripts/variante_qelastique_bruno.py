# -*- coding: utf-8 -*-
"""Diagnostic : variante quasi-elastique (eta x 1e9 -> branches Maxwell figees).

Meme raideur instantanee que le modele complet, mais sans memoire visqueuse.
But : tester si la FORME des essais de Bruno correspond mieux au squelette
elastique du modele qu'au modele viscoelastique complet.
"""
import sys
import time as _time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ALPHA = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2"
BRUNO = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\expérimentale\Bruno"
OUT = r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"

sys.path.insert(0, ALPHA)
import Base as modele  # noqa: E402
from pression import measured_pressure_payload  # noqa: E402

ESSAIS = {
    "essai1_153353": "essai_force_pression_20260730_153353.csv",
    "essai2_155135": "essai_force_pression_20260730_155135.csv",
}
PLATEAUX_ESSAI1 = [(112.4, 121.1), (122.2, 129.7), (143.9, 154.6), (158.4, 169.1), (179.6, 186.0), (186.2, 196.4)]

mw = modele.default_maxwell_tensile_params()
mw_frozen = modele.default_maxwell_tensile_params(eta1=mw.eta1 * 1e9, eta2=mw.eta2 * 1e9, eta3=mw.eta3 * 1e9)
mat_frozen = modele.default_material_params(maxwell=mw_frozen)


def slope(t, y):
    return float(np.polyfit(t, y, 1)[0])


for tag, name in ESSAIS.items():
    df = pd.read_csv(BRUNO + "\\" + name)
    cols = {"t": df["time_s"].to_numpy(float), "p": df["pressure_bar"].to_numpy(float)}
    payload = measured_pressure_payload(cols, "t", "p", unit="bar", subtract_initial=False)
    t_meas, p_meas = payload["time"], payload["pressure_MPa"]
    step = max(1, int(round(0.5 / np.median(np.diff(t_meas)))))
    idx = np.arange(0, len(t_meas), step)
    if idx[-1] != len(t_meas) - 1:
        idx = np.append(idx, len(t_meas) - 1)

    t0 = _time.perf_counter()
    _, arr = modele.run_blocked_actuation(mat=mat_frozen, pressure_time=t_meas[idx], pressure_MPa=p_meas[idx])
    print(f"[{tag}] variante quasi-elastique simulee en {_time.perf_counter()-t0:.1f} s")

    t_sim = arr["time"]
    f_el = arr["force_act_mN"]
    p_sim = arr["pressure_MPa"]

    ref = pd.read_csv(OUT + f"\\{tag}_sim.csv")
    f_meas = ref["force_act_meas_mN"].to_numpy(float)
    f_ve = ref["force_act_sim_mN"].to_numpy(float)

    k_el = float(np.dot(f_el, f_meas) / np.dot(f_el, f_el))
    rmse_el = float(np.sqrt(np.mean((k_el * f_el - f_meas) ** 2)))
    corr_el = float(np.corrcoef(f_el, f_meas)[0, 1])
    peak = float(np.max(f_meas))
    print(f"[{tag}] quasi-elastique : pic={np.max(f_el):.1f} mN, k={k_el:.2f}, RMSE={rmse_el:.1f} mN ({100*rmse_el/peak:.1f}% du pic), corr={corr_el:.4f}")

    if tag == "essai1_153353":
        for (ta, tb) in PLATEAUX_ESSAI1:
            m = (t_sim >= ta) & (t_sim <= tb)
            print(f"    palier [{ta:6.1f},{tb:6.1f}] pente mesure={slope(t_sim[m], f_meas[m]):+.4f}  quasi-el(x{k_el:.1f})={slope(t_sim[m], k_el*f_el[m]):+.4f} mN/s")
    if tag == "essai2_155135":
        i_pk = int(np.argmax(p_sim))
        for p_probe in (0.2, 0.3, 0.4):
            fu = np.interp(p_probe, p_sim[:i_pk], k_el * f_el[:i_pk])
            fd = np.interp(p_probe, p_sim[i_pk:][::-1], k_el * f_el[i_pk:][::-1])
            print(f"    hysteresis a P={p_probe:.1f} MPa : quasi-el(x{k_el:.1f}) descente-montee={fd-fu:+.1f} mN")
        print(f"    fin de descente (P~0) : quasi-el F_act={f_el[-1]:.1f} mN")

    fig, ax = plt.subplots(figsize=(10, 5))
    k_ve = float(np.dot(f_ve, f_meas) / np.dot(f_ve, f_ve))
    ax.plot(t_sim, f_meas, color="tab:green", lw=0.9, label="mesure (F - precharge)")
    ax.plot(ref["time_s"], k_ve * f_ve, color="tab:red", lw=1.2, label=f"modele complet x {k_ve:.2f}")
    ax.plot(t_sim, k_el * f_el, color="tab:blue", lw=1.2, ls="--", label=f"quasi-elastique x {k_el:.2f}")
    ax.set_xlabel("Temps (s)"); ax.set_ylabel("Force d'actionnement (mN)")
    ax.set_title(f"{name} — mesure vs modele complet vs quasi-elastique (echelles ajustees)")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT + f"\\{tag}_el.png", dpi=150)
    plt.close(fig)

    pd.DataFrame({"time_s": t_sim, "pressure_MPa": p_sim, "force_act_el_mN": f_el}).to_csv(OUT + f"\\{tag}_el.csv", index=False)

print("Termine.")
