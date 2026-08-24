# -*- coding: utf-8 -*-
"""Variante 'equilibree' : 7200 s a P=0 (bloque) avant l'historique mesure.

But : retirer le transitoire de relaxation de la precontrainte (protocole article
= pompe demarree immediatement apres la precontrainte) pour comparer a l'essai
de Bruno ou le muscle etait equilibre (force stable a t=0).
"""
import sys
import json
import time as _time

import numpy as np

ALPHA = r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/alpha V2"
DATA = r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/expérimentale/Bruno"
OUT = r"C:/Users/B3DCB~1.PER/AppData/Local/Temp/claude/C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail/c5852e10-6d38-4703-b218-a1649a5b12c9/scratchpad"

sys.path.insert(0, ALPHA)
import Base  # noqa: E402
import pression  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

T_SETTLE = 7200.0
DT_SETTLE = 50.0

TRIALS = {
    "essai_153353": {
        "file": DATA + "/essai_force_pression_20260730_153353.csv",
        "holds": [(18.0, 26.0), (29.5, 36.0), (110.0, 176.0)],
    },
    "essai_155135": {
        "file": DATA + "/essai_force_pression_20260730_155135.csv",
        "holds": [(31.5, 56.0), (38.0, 56.0)],
    },
}


def linfit(t, y):
    A = np.vstack([t, np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(coef[0]), float(coef[1])


results = {}
for name, spec in TRIALS.items():
    print("=" * 78, flush=True)
    print(name, "(variante equilibree)", flush=True)
    raw = open(spec["file"], "rb").read()
    cols = pression.parse_uploaded_numeric_csv(raw)
    payload = pression.measured_pressure_payload(
        cols, "time_s", "pressure_bar", unit="bar", subtract_initial=False
    )
    t_meas = payload["time"]
    p_meas = payload["pressure_MPa"]
    exp = pression.experimental_force_pressure_payload(
        cols, "time_s", "pressure_bar", "force_mN",
        pressure_unit="bar", force_unit="mN",
    )
    f_meas = exp["force_mN"]

    # historique etendu : palier P=0 de 7200 s puis historique mesure
    t_pre = np.arange(0.0, T_SETTLE, DT_SETTLE)
    t_ext = np.concatenate([t_pre, T_SETTLE + t_meas])
    p_ext = np.concatenate([np.zeros_like(t_pre), p_meas])

    t0 = _time.perf_counter()
    model, arr = Base.run_blocked_actuation(pressure_time=t_ext, pressure_MPa=p_ext)
    elapsed = _time.perf_counter() - t0
    print(f"  simulation {len(t_ext)} pas en {elapsed:.0f} s", flush=True)

    t_sim = arr["time"]
    f_tot = arr["force_total_mN"]
    # reference = etat equilibre juste avant le demarrage de la pression mesuree
    i_ref = int(np.searchsorted(t_sim, T_SETTLE, side="left"))
    f_ref = float(f_tot[i_ref])
    f_sim_act_i = np.interp(T_SETTLE + t_meas, t_sim, f_tot) - f_ref
    df_meas = f_meas - f_meas[0]

    rmse_raw = float(np.sqrt(np.mean((f_sim_act_i - df_meas) ** 2)))
    corr = float(np.corrcoef(f_sim_act_i, df_meas)[0, 1])
    denom = float(np.dot(f_sim_act_i, f_sim_act_i))
    a_scale = float(np.dot(f_sim_act_i, df_meas) / denom)
    rmse_scaled = float(np.sqrt(np.mean((a_scale * f_sim_act_i - df_meas) ** 2)))
    peak_meas = float(np.max(df_meas))
    peak_sim = float(np.max(f_sim_act_i))

    holds = []
    for (ta, tb) in spec["holds"]:
        m = (t_meas >= ta) & (t_meas <= tb)
        sp, _ = linfit(t_meas[m], p_meas[m] * 10.0)
        sf_m, _ = linfit(t_meas[m], f_meas[m])
        sf_s, _ = linfit(t_meas[m], f_sim_act_i[m])
        holds.append(dict(t=(ta, tb), P_bar=float(np.mean(p_meas[m]) * 10.0),
                          dPdt_bar_s=sp, dFdt_meas=sf_m, dFdt_sim=sf_s))

    res = dict(
        preload_sim_settled_mN=f_ref,
        relax_during_settle_mN=float(f_tot[i_ref] - f_tot[0]),
        dF_peak_meas_mN=peak_meas, dF_peak_sim_mN=peak_sim,
        peak_ratio_sim_over_meas=peak_sim / peak_meas,
        rmse_raw_mN=rmse_raw, corr=corr,
        scale_factor_fit=a_scale, rmse_scaled_mN=rmse_scaled,
        holds=holds, sim_seconds=elapsed,
    )
    results[name] = res
    print(json.dumps(res, indent=2), flush=True)

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(t_meas, df_meas, "k", lw=1.0, label="mesure ΔF")
    ax.plot(t_meas, f_sim_act_i, "tab:red", lw=1.2, label="modèle équilibré ΔF")
    ax.plot(t_meas, a_scale * f_sim_act_i, "tab:orange", ls="--", lw=1.2,
            label=f"modèle équilibré × {a_scale:.2f}")
    ax2 = ax.twinx()
    ax2.plot(t_meas, p_meas * 10, "tab:blue", lw=0.8, alpha=0.5)
    ax2.set_ylabel("P [bar]", color="tab:blue")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("ΔF [mN]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    ax.set_title(f"{name} — modèle équilibré 7200 s avant essai | corr {corr:.3f} | "
                 f"pic sim/mes {peak_sim/peak_meas:.2f} | RMSE éch. {rmse_scaled:.1f} mN")
    fig.tight_layout()
    fig.savefig(f"{OUT}/confrontation_{name}_equilibre.png", dpi=130)
    plt.close(fig)

    np.savez(f"{OUT}/sim_{name}_settled.npz",
             t=t_meas, p_mpa=p_meas, df_meas=df_meas, f_sim_act=f_sim_act_i)

with open(f"{OUT}/metrics_bruno_settled.json", "w", encoding="utf-8") as fh:
    json.dump(results, fh, indent=2)
print("TERMINE", flush=True)
