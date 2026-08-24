# -*- coding: utf-8 -*-
"""Confrontation modele bloque (Base.py, specimen article) vs essais de Bruno.

Chemin identique a l'interface : pression mesuree (bar -> MPa) injectee dans
run_blocked_actuation(pressure_time=..., pressure_MPa=...), comparaison sur
force_act_mN (variation de force par rapport a l'etat bloque a P=0), la
precharge mesuree (~570-590 mN) etant soustraite cote essai.
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

# fenetres de paliers detectees pour l'essai 1 (voir caracterise.py)
PLATEAUX_ESSAI1 = [(112.4, 121.1), (122.2, 129.7), (143.9, 154.6), (158.4, 169.1), (179.6, 186.0), (186.2, 196.4)]


def baseline_force(p_mpa, f_mn):
    """Precharge mesuree : mediane de la force quand P est au plancher (comme interface.py)."""
    span = float(np.ptp(p_mpa))
    lim = float(np.min(p_mpa)) + max(0.02 * span, 1e-9)
    mask = p_mpa <= lim
    if int(np.count_nonzero(mask)) < 3:
        mask = np.zeros(len(f_mn), dtype=bool)
        mask[: max(1, len(f_mn) // 20)] = True
    return float(np.median(f_mn[mask]))


def slope(t, y):
    return float(np.polyfit(t, y, 1)[0])


results = {}
for tag, name in ESSAIS.items():
    df = pd.read_csv(BRUNO + "\\" + name)
    t_full = df["time_s"].to_numpy(float)
    p_bar_full = df["pressure_bar"].to_numpy(float)
    f_full = df["force_mN"].to_numpy(float)

    # --- chemin pression.py : bar -> MPa, temps recale a 0, clamp >= 0 ---
    cols = {"t": t_full, "p": p_bar_full}
    payload = measured_pressure_payload(cols, "t", "p", unit="bar", subtract_initial=False)
    t_meas = payload["time"]
    p_meas = payload["pressure_MPa"]
    f_meas = np.interp(t_meas, t_full - t_full[0], f_full)

    # --- decimation a ~0.5 s pour le cout de calcul (forme conservee) ---
    step = max(1, int(round(0.5 / np.median(np.diff(t_meas)))))
    idx = np.arange(0, len(t_meas), step)
    if idx[-1] != len(t_meas) - 1:
        idx = np.append(idx, len(t_meas) - 1)
    t_sim_in = t_meas[idx]
    p_sim_in = p_meas[idx]

    print(f"[{tag}] {len(t_meas)} pts mesures -> {len(idx)} pts simules (dt~{np.median(np.diff(t_sim_in)):.2f} s)")
    t0 = _time.perf_counter()
    _, arr = modele.run_blocked_actuation(pressure_time=t_sim_in, pressure_MPa=p_sim_in)
    print(f"[{tag}] simulation en {_time.perf_counter()-t0:.1f} s")

    t_sim = arr["time"]
    f_sim_act = arr["force_act_mN"]
    p_sim = arr["pressure_MPa"]
    f_sim_total0 = float(arr["force_total_mN"][0])

    # --- mesure : force d'actionnement = force - precharge ---
    f0 = baseline_force(p_meas, f_meas)
    f_meas_act = f_meas - f0
    # mesure interpolee sur la grille simulation (elle-meme sous-ensemble des temps mesures)
    f_meas_act_g = np.interp(t_sim, t_meas, f_meas_act)

    # --- metriques ---
    rmse_raw = float(np.sqrt(np.mean((f_sim_act - f_meas_act_g) ** 2)))
    denom = float(np.dot(f_sim_act, f_sim_act))
    k = float(np.dot(f_sim_act, f_meas_act_g) / denom) if denom > 0 else np.nan
    rmse_scaled = float(np.sqrt(np.mean((k * f_sim_act - f_meas_act_g) ** 2)))
    peak_meas = float(np.max(f_meas_act_g))
    peak_sim = float(np.max(f_sim_act))
    sens_meas = slope(p_sim[p_sim > 0.02], f_meas_act_g[p_sim > 0.02])
    sens_sim = slope(p_sim[p_sim > 0.02], f_sim_act[p_sim > 0.02])
    corr = float(np.corrcoef(f_sim_act, f_meas_act_g)[0, 1])
    nrmse_scaled = rmse_scaled / max(peak_meas, 1e-12)

    res = dict(k=k, rmse_raw=rmse_raw, rmse_scaled=rmse_scaled, nrmse_scaled=nrmse_scaled,
               peak_meas=peak_meas, peak_sim=peak_sim, ratio=peak_meas / peak_sim,
               sens_meas=sens_meas, sens_sim=sens_sim, corr=corr,
               f0_meas=f0, f0_sim=f_sim_total0)
    print(f"[{tag}] precharge mesuree={f0:.1f} mN | precharge bloquee modele={f_sim_total0:.1f} mN")
    print(f"[{tag}] pic mesure={peak_meas:.1f} mN, pic simule={peak_sim:.1f} mN, rapport={peak_meas/peak_sim:.2f}")
    print(f"[{tag}] sensibilite dF/dP: mesure={sens_meas:.0f} mN/MPa, modele={sens_sim:.0f} mN/MPa")
    print(f"[{tag}] correlation sim/mesure={corr:.4f}")
    print(f"[{tag}] RMSE brut={rmse_raw:.1f} mN | facteur d'echelle k={k:.2f} | RMSE apres echelle={rmse_scaled:.1f} mN ({100*nrmse_scaled:.1f}% du pic)")

    # --- pentes de relaxation aux paliers (essai 1) ---
    if tag == "essai1_153353":
        print(f"[{tag}] pentes aux paliers (mN/s), mesure vs k*simulation :")
        pl = []
        for (ta, tb) in PLATEAUX_ESSAI1:
            m = (t_sim >= ta) & (t_sim <= tb)
            if np.count_nonzero(m) < 4:
                continue
            sm = slope(t_sim[m], f_meas_act_g[m])
            ss = slope(t_sim[m], k * f_sim_act[m])
            sp = slope(t_sim[m], p_sim[m])
            pl.append((ta, tb, sm, ss, sp))
            print(f"    [{ta:6.1f},{tb:6.1f}] s  dP/dt={sp*1000:+6.2f} kPa/s  mesure={sm:+.4f}  modele(x{k:.1f})={ss:+.4f} mN/s")
        res["plateaux"] = pl

    # --- hysteresis (essai 2 : montee puis descente) ---
    if tag == "essai2_155135":
        i_pk = int(np.argmax(p_sim))
        for p_probe in (0.2, 0.3, 0.4):
            fu_m = np.interp(p_probe, p_sim[:i_pk], f_meas_act_g[:i_pk])
            fd_m = np.interp(p_probe, p_sim[i_pk:][::-1], f_meas_act_g[i_pk:][::-1])
            fu_s = np.interp(p_probe, p_sim[:i_pk], k * f_sim_act[:i_pk])
            fd_s = np.interp(p_probe, p_sim[i_pk:][::-1], k * f_sim_act[i_pk:][::-1])
            print(f"[{tag}] hysteresis a P={p_probe:.1f} MPa : mesure descente-montee={fd_m-fu_m:+.1f} mN | modele(x{k:.1f})={fd_s-fu_s:+.1f} mN")

    results[tag] = (res, t_sim, p_sim, f_meas_act_g, f_sim_act)

    # --- figure temps ---
    fig, ax = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    ax[0].plot(t_meas, 10.0 * p_meas, color="tab:blue", lw=0.8, label="P mesuree")
    ax[0].plot(t_sim, 10.0 * p_sim, color="k", lw=1.0, ls="--", label="P injectee (decimee)")
    ax[0].set_ylabel("Pression (bar)")
    ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[0].set_title(f"{ESSAIS[tag]} — pression mesuree injectee dans run_blocked_actuation (specimen article, eps=0.8)")
    ax[1].plot(t_meas, f_meas_act, color="tab:green", lw=0.8, label=f"F mesuree - precharge ({f0:.0f} mN)")
    ax[1].plot(t_sim, f_sim_act, color="tab:red", lw=1.4, label="F modele (force_act_mN), brute")
    ax[1].plot(t_sim, k * f_sim_act, color="tab:orange", lw=1.4, ls="--", label=f"F modele x {k:.2f} (echelle ajustee)")
    ax[1].set_ylabel("Force d'actionnement (mN)")
    ax[1].legend(); ax[1].grid(alpha=0.3)
    ax[2].plot(t_sim, f_meas_act_g - k * f_sim_act, color="tab:purple", lw=1.0)
    ax[2].axhline(0.0, color="k", lw=0.6)
    ax[2].set_ylabel("Residu mesure - k*modele (mN)")
    ax[2].set_xlabel("Temps (s)")
    ax[2].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT + f"\\{tag}_t.png", dpi=150)
    plt.close(fig)

    # --- figure force-pression ---
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.plot(10.0 * p_sim, f_meas_act_g, color="tab:green", lw=1.0, label="mesure")
    ax.plot(10.0 * p_sim, f_sim_act, color="tab:red", lw=1.2, label="modele brut")
    ax.plot(10.0 * p_sim, k * f_sim_act, color="tab:orange", lw=1.2, ls="--", label=f"modele x {k:.2f}")
    ax.set_xlabel("Pression (bar)"); ax.set_ylabel("Force d'actionnement (mN)")
    ax.set_title(f"{ESSAIS[tag]} — force vs pression")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT + f"\\{tag}_fp.png", dpi=150)
    plt.close(fig)

    pd.DataFrame({"time_s": t_sim, "pressure_MPa": p_sim,
                  "force_act_meas_mN": f_meas_act_g,
                  "force_act_sim_mN": f_sim_act}).to_csv(OUT + f"\\{tag}_sim.csv", index=False)

print("Termine.")
