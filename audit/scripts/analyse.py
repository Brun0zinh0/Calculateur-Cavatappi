# -*- coding: utf-8 -*-
"""Analyse et figures : simulation vs mesures Bruno."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"

PLATEAUX = {
    # fenetres [t0, t1] de palier pour la pente de relaxation
    "20260730_153353": [(110.0, 200.0)],
    "20260730_155135": [(35.0, 50.0)],
}


def power_fit(P, dF, pmin=0.05):
    """Ajuste dF = a * P^n sur la branche de charge (P croissant, P > pmin)."""
    m = (P > pmin) & (dF > 0.5)
    if m.sum() < 10:
        return np.nan, np.nan
    lp, lf = np.log(P[m]), np.log(dF[m])
    n, la = np.polyfit(lp, lf, 1)
    return np.exp(la), n


for tag in ["20260730_153353", "20260730_155135"]:
    d = np.load(f"{OUT}\\sim_{tag}.npz")
    t = d["t_meas"]
    P = d["p_meas"]
    Fm = d["F_meas"]
    dFm = Fm - Fm[0]
    ts = d["t_sim"]
    dFs_full = d["F_sim_act"]
    F0s = d["F_sim_total"][0]
    # aligne : la simulation retourne un pas de plus (etat initial a t=0 double ?)
    dFs = np.interp(t, ts, dFs_full)

    # facteur d'echelle optimal (moindre carres, sans decalage : les deux partent de 0)
    k = float(np.dot(dFm, dFs) / np.dot(dFs, dFs))
    rmse_raw = float(np.sqrt(np.mean((dFm - dFs) ** 2)))
    rmse_scaled = float(np.sqrt(np.mean((dFm - k * dFs) ** 2)))
    corr = float(np.corrcoef(dFm, dFs)[0, 1])
    peak_ratio = float(dFm.max() / dFs.max())

    # pentes de relaxation aux plateaux
    rel = []
    for (t0, t1) in PLATEAUX[tag]:
        m = (t >= t0) & (t <= t1)
        sm = np.polyfit(t[m], dFm[m], 1)[0]
        ss = np.polyfit(t[m], dFs[m], 1)[0]
        dP = np.polyfit(t[m], P[m], 1)[0]
        rel.append((t0, t1, sm, ss, dP))

    # exposant de non-linearite sur la charge (jusqu'au max de pression)
    i_pk = int(np.argmax(P))
    a_m, n_m = power_fit(P[:i_pk], dFm[:i_pk])
    a_s, n_s = power_fit(P[:i_pk], dFs[:i_pk])

    print("=" * 72)
    print(f"{tag}  (Pmax={P.max():.3f} MPa)")
    print(f"  precharge:      mesure {Fm[0]:7.1f} mN | simulation {F0s:7.1f} mN")
    print(f"  pic dF:         mesure {dFm.max():7.1f} mN | simulation {dFs.max():7.1f} mN | rapport meas/sim = {peak_ratio:.2f}")
    print(f"  RMSE brut  = {rmse_raw:7.1f} mN")
    print(f"  facteur d'echelle optimal k = {k:.2f}  ->  RMSE apres echelle = {rmse_scaled:.1f} mN "
          f"({100*rmse_scaled/max(dFm.max(),1e-9):.1f} % du pic mesure)")
    print(f"  correlation temporelle dF_meas / dF_sim = {corr:.4f}")
    print(f"  exposant n de dF ~ a P^n (charge) : mesure n={n_m:.2f} | simulation n={n_s:.2f}")
    for (t0, t1, sm, ss, dP) in rel:
        print(f"  palier [{t0:.0f},{t1:.0f}] s : pente mesure {sm:+.3f} mN/s | simulation {ss:+.3f} mN/s "
              f"| derive pression {dP*1000:+.3f} kPa/s")

    # figure
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.4, 1.4])
    ax0 = fig.add_subplot(gs[0, :])
    ax0.plot(t, P, color="tab:blue", lw=0.9)
    ax0.set_ylabel("Pression (MPa)")
    ax0.set_title(f"Essai {tag} : entrée pression mesurée (historique injecté dans le modèle)")
    ax0.grid(alpha=0.3)

    ax1 = fig.add_subplot(gs[1, :])
    ax1.plot(t, dFm, color="tab:red", lw=1.0, label="Mesure ΔF (mN)")
    ax1.plot(t, dFs, color="tab:green", lw=1.0, label="Modèle ΔF (défauts article)")
    ax1.plot(t, k * dFs, color="tab:green", lw=1.0, ls="--", label=f"Modèle × {k:.2f} (échelle ajustée)")
    ax1.set_ylabel("ΔF = F − F(0)  (mN)")
    ax1.legend(loc="upper left", fontsize=9)
    ax1.grid(alpha=0.3)
    ax1.set_title(f"Force d'actionnement bloquée — RMSE échelle ajustée {rmse_scaled:.1f} mN, corr {corr:.3f}")

    i_pk = int(np.argmax(P))
    ax2 = fig.add_subplot(gs[2, 0])
    ax2.plot(P[:i_pk], dFm[:i_pk], color="tab:red", lw=1.0, label="Mesure (charge)")
    ax2.plot(P[i_pk:], dFm[i_pk:], color="tab:red", lw=0.8, ls=":", label="Mesure (décharge)")
    ax2.plot(P[:i_pk], k * dFs[:i_pk], color="tab:green", lw=1.0, label=f"Modèle × {k:.2f} (charge)")
    ax2.plot(P[i_pk:], k * dFs[i_pk:], color="tab:green", lw=0.8, ls=":", label=f"Modèle × {k:.2f} (décharge)")
    ax2.set_xlabel("Pression (MPa)")
    ax2.set_ylabel("ΔF (mN)")
    ax2.set_title(f"ΔF(P) — exposants : mesure {n_m:.2f}, modèle {n_s:.2f}")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    ax3 = fig.add_subplot(gs[2, 1])
    (t0, t1, sm, ss, dP) = rel[0]
    m = (t >= t0) & (t <= t1)
    ax3.plot(t[m], dFm[m] - dFm[m][0], color="tab:red", lw=1.0, label=f"Mesure ({sm:+.3f} mN/s)")
    ax3.plot(t[m], (k * dFs[m]) - (k * dFs[m][0]), color="tab:green", lw=1.0,
             label=f"Modèle × {k:.2f} ({k*ss:+.3f} mN/s)")
    ax3.set_xlabel("Temps (s)")
    ax3.set_ylabel("ΔF − ΔF(début palier)  (mN)")
    ax3.set_title(f"Palier [{t0:.0f},{t1:.0f}] s : évolution pendant le maintien")
    ax3.legend(fontsize=8)
    ax3.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(f"{OUT}\\confrontation_{tag}.png", dpi=130)
    plt.close(fig)
    print(f"  figure -> confrontation_{tag}.png")
