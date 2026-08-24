# -*- coding: utf-8 -*-
"""Confrontation modele bloque (Base.run_blocked_actuation) vs essais de Bruno.

Chemin de donnees identique a l'interface Streamlit :
  pression.parse_uploaded_numeric_csv -> pression.measured_pressure_payload
  -> Base.run_blocked_actuation(config par defaut = specimen article, eps=0.8)

Sorties : PNG + metriques imprimees (RMSE, rapport des pics, facteur d'echelle,
pentes de relaxation aux paliers).
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

TRIALS = {
    "essai_153353": {
        "file": DATA + "/essai_force_pression_20260730_153353.csv",
        # fenetres de palier (pression quasi constante) reperees lors de la caracterisation
        "holds": [(18.0, 26.0), (29.5, 36.0), (57.0, 64.0), (110.0, 176.0)],
    },
    "essai_155135": {
        "file": DATA + "/essai_force_pression_20260730_155135.csv",
        "holds": [(31.5, 56.0), (67.0, 75.0)],
    },
}


def linfit(t, y):
    """pente (unite/s) et intercept par moindres carres."""
    A = np.vstack([t, np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(coef[0]), float(coef[1])


results = {}
for name, spec in TRIALS.items():
    print("=" * 78, flush=True)
    print(name, flush=True)
    raw = open(spec["file"], "rb").read()
    cols = pression.parse_uploaded_numeric_csv(raw)

    # --- meme chemin que l'interface : payload pression mesuree en MPa ---
    payload = pression.measured_pressure_payload(
        cols, "time_s", "pressure_bar", unit="bar", subtract_initial=False
    )
    t = payload["time"]
    p_mpa = payload["pressure_MPa"]

    # --- force mesuree alignee sur la meme base de temps ---
    exp = pression.experimental_force_pressure_payload(
        cols, "time_s", "pressure_bar", "force_mN",
        pressure_unit="bar", force_unit="mN",
        unfiltered_force_column="force_unfiltered_mN",
    )
    f_meas = exp["force_mN"]
    assert len(f_meas) == len(t)

    # --- simulation, parametres par defaut (specimen de l'article, eps=0.8) ---
    t0 = _time.perf_counter()
    model, arr = Base.run_blocked_actuation(pressure_time=t, pressure_MPa=p_mpa)
    elapsed = _time.perf_counter() - t0
    print(f"  simulation {len(t)} pas en {elapsed:.0f} s", flush=True)

    t_sim = arr["time"]
    f_sim_tot = arr["force_total_mN"]
    f_sim_act = arr["force_act_mN"]

    # interpolation sim -> temps mesures (les grilles coincident a 1 pas pres)
    f_sim_act_i = np.interp(t, t_sim, f_sim_act)
    df_meas = f_meas - f_meas[0]

    # --- metriques globales ---
    rmse_raw = float(np.sqrt(np.mean((f_sim_act_i - df_meas) ** 2)))
    peak_meas = float(np.max(df_meas))
    peak_sim = float(np.max(f_sim_act_i))
    denom = float(np.dot(f_sim_act_i, f_sim_act_i))
    a_scale = float(np.dot(f_sim_act_i, df_meas) / denom) if denom > 0 else np.nan
    rmse_scaled = float(np.sqrt(np.mean((a_scale * f_sim_act_i - df_meas) ** 2)))
    corr = float(np.corrcoef(f_sim_act_i, df_meas)[0, 1])

    # sensibilite quasi statique dF/dP (regression lineaire F vs P)
    sens_meas, _ = linfit(p_mpa, df_meas)   # mN/MPa
    sens_sim, _ = linfit(p_mpa, f_sim_act_i)

    # --- paliers : pentes de relaxation ---
    holds = []
    for (ta, tb) in spec["holds"]:
        m = (t >= ta) & (t <= tb)
        if m.sum() < 5:
            continue
        sp, _ = linfit(t[m], p_mpa[m] * 10.0)       # bar/s
        sf_meas, _ = linfit(t[m], f_meas[m])        # mN/s
        sf_sim, _ = linfit(t[m], f_sim_act_i[m])    # mN/s
        holds.append(dict(
            t=(ta, tb), P_bar=float(np.mean(p_mpa[m]) * 10.0),
            dPdt_bar_s=sp, dFdt_meas=sf_meas, dFdt_sim=sf_sim,
            dF_meas_win=float(f_meas[m][-1] - f_meas[m][0]),
            dF_sim_win=float(f_sim_act_i[m][-1] - f_sim_act_i[m][0]),
        ))

    res = dict(
        n=len(t), duration_s=float(t[-1]),
        P_max_MPa=float(p_mpa.max()),
        preload_meas_mN=float(f_meas[0]),
        preload_sim_mN=float(f_sim_tot[0]),
        dF_peak_meas_mN=peak_meas, dF_peak_sim_mN=peak_sim,
        peak_ratio_sim_over_meas=peak_sim / peak_meas,
        rmse_raw_mN=rmse_raw, corr=corr,
        scale_factor_fit=a_scale, rmse_scaled_mN=rmse_scaled,
        sens_meas_mN_per_MPa=sens_meas, sens_sim_mN_per_MPa=sens_sim,
        holds=holds, sim_seconds=elapsed,
    )
    results[name] = res
    print(json.dumps(res, indent=2), flush=True)

    # --- figure ---
    fig, axes = plt.subplots(3, 1, figsize=(11, 11), sharex=False)
    ax = axes[0]
    ax.plot(t, p_mpa * 10.0, color="tab:blue", lw=1.0)
    ax.set_ylabel("Pression [bar]")
    ax.set_xlabel("t [s]")
    ax.set_title(f"{name} — profil de pression mesuré (max {p_mpa.max()*10:.2f} bar)")
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(t, df_meas, color="k", lw=1.0, label="mesure ΔF (F - précharge)")
    ax.plot(t, f_sim_act_i, color="tab:red", lw=1.2, label="modèle ΔF (défauts article, eps=0.8)")
    ax.plot(t, a_scale * f_sim_act_i, color="tab:orange", lw=1.2, ls="--",
            label=f"modèle × {a_scale:.2f} (échelle ajustée)")
    for (ta, tb) in spec["holds"]:
        ax.axvspan(ta, tb, color="tab:blue", alpha=0.07)
    ax.set_ylabel("ΔF bloquée [mN]")
    ax.set_xlabel("t [s]")
    ax.legend(loc="best", fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_title(
        f"RMSE brut {rmse_raw:.1f} mN | RMSE après échelle {rmse_scaled:.1f} mN | "
        f"corr {corr:.3f} | pic sim/mes {peak_sim/peak_meas:.2f}"
    )

    ax = axes[2]
    ax.plot(p_mpa * 10.0, df_meas, ".", color="k", ms=2, label="mesure")
    ax.plot(p_mpa * 10.0, f_sim_act_i, ".", color="tab:red", ms=2, label="modèle")
    ax.plot(p_mpa * 10.0, a_scale * f_sim_act_i, ".", color="tab:orange", ms=2,
            label=f"modèle × {a_scale:.2f}")
    ax.set_xlabel("Pression [bar]")
    ax.set_ylabel("ΔF [mN]")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_title("Caractéristique ΔF–P (hystérésis)")

    fig.tight_layout()
    png = f"{OUT}/confrontation_{name}.png"
    fig.savefig(png, dpi=130)
    plt.close(fig)
    print("  figure ->", png, flush=True)

    np.savez(
        f"{OUT}/sim_{name}.npz",
        t=t, p_mpa=p_mpa, f_meas=f_meas, df_meas=df_meas,
        f_sim_act=f_sim_act_i, f_sim_tot=np.interp(t, t_sim, f_sim_tot),
    )

with open(f"{OUT}/metrics_bruno.json", "w", encoding="utf-8") as fh:
    json.dump(results, fh, indent=2)
print("TERMINE", flush=True)
