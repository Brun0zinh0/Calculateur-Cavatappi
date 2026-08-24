# -*- coding: utf-8 -*-
"""Contre-expertise : verifications independantes sur les CSV de Bruno.

1. Precharge + derive du capteur a P=0 (debut et fin d'essai)
2. Detection independante des paliers (dP/dt ~ 0) et pentes de force mesurees
3. Courbure de F(P) sur la montee (fit a*P + b*P^2, b/a)
4. Hysteresis essai 2 (montee vs descente)
5. Artefacts capteur : quantification pression (ADC), filtre force
"""
import numpy as np
import pandas as pd

BRUNO = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\expérimentale\Bruno"
ESSAIS = {
    "essai1_153353": "essai_force_pression_20260730_153353.csv",
    "essai2_155135": "essai_force_pression_20260730_155135.csv",
}


def slope(t, y):
    if len(t) < 3:
        return np.nan
    return float(np.polyfit(t, y, 1)[0])


for tag, name in ESSAIS.items():
    df = pd.read_csv(BRUNO + "\\" + name)
    t = df["time_s"].to_numpy(float)
    p_bar = df["pressure_bar"].to_numpy(float)
    p = 0.1 * p_bar  # MPa
    f = df["force_mN"].to_numpy(float)
    f_raw = df["force_unfiltered_mN"].to_numpy(float)

    print(f"\n================= {tag} ({name}) =================")
    print(f"duree {t[-1]-t[0]:.1f} s, {len(t)} pts, dt median {np.median(np.diff(t))*1000:.0f} ms")
    print(f"P max {p.max():.3f} MPa ; F min/max {f.min():.1f}/{f.max():.1f} mN")

    # --- 1. precharge et derive a P=0 en debut d'essai ---
    m0 = p <= 1e-9
    # premier bloc contigu a P=0
    i_first_p = np.argmax(p > 1e-9)  # premier indice sous pression
    t_pre, f_pre = t[:i_first_p], f[:i_first_p]
    print(f"precharge (mediane avant 1re pression, {i_first_p} pts, {t_pre[-1]-t_pre[0]:.1f} s) = {np.median(f_pre):.1f} mN")
    print(f"derive capteur a P=0 en debut d'essai : pente = {slope(t_pre, f_pre)*1000:+.3f} mN/ks = {slope(t_pre, f_pre):+.4f} mN/s (etendue {np.ptp(f_pre):.2f} mN)")
    # bruit du capteur : ecart-type residuel autour de la droite
    if len(t_pre) > 5:
        res = f_pre - np.polyval(np.polyfit(t_pre, f_pre, 1), t_pre)
        print(f"bruit force (sd residuel a P=0) = {np.std(res):.3f} mN ; filtre vs brut sd = {np.std(f_pre - f_raw[:i_first_p]):.3f} mN")

    # fin d'essai : P revenu ~0 ?
    m_end = (t > t[-1] - 5.0)
    print(f"fin d'essai : P mediane {np.median(p[m_end])*1000:.1f} kPa, F mediane {np.median(f[m_end]):.1f} mN (delta precharge {np.median(f[m_end]) - np.median(f_pre):+.1f} mN)")

    # --- 5. quantification du capteur de pression ---
    v = df["voltage"].to_numpy(float)
    dv = np.unique(np.round(np.diff(np.unique(v)), 6))
    p_steps = np.unique(np.round(np.diff(np.unique(p_bar[p_bar > 0])), 6))
    print(f"quantification : pas de tension min {dv[dv>0].min():.4f} V ; pas de pression min {p_steps[p_steps>0].min()*100:.2f} kPa")

    # --- 2. paliers independants : fenetres glissantes ou |dP/dt| petit ---
    f0 = np.median(f_pre)
    fa = f - f0
    # lissage pression pour derivee
    win = 21
    ker = np.ones(win) / win
    p_s = np.convolve(p, ker, mode="same")
    dpdt = np.gradient(p_s, t)
    is_plat = (np.abs(dpdt) < 0.002) & (p > 0.05)  # <2 kPa/s et sous pression
    # segments contigus > 5 s
    print("paliers detectes (|dP/dt|<2 kPa/s, P>0.05 MPa, duree>5 s) :")
    i = 0
    n = len(t)
    while i < n:
        if is_plat[i]:
            j = i
            while j + 1 < n and is_plat[j + 1]:
                j += 1
            if t[j] - t[i] >= 5.0:
                sl_f = slope(t[i:j], fa[i:j])
                sl_p = slope(t[i:j], p[i:j])
                # sensibilite locale dF/dP (secante autour du niveau du palier)
                print(f"  [{t[i]:7.1f},{t[j]:7.1f}] s ({t[j]-t[i]:5.1f} s)  P~{np.mean(p[i:j]):.3f} MPa  dP/dt={sl_p*1000:+7.3f} kPa/s  dF/dt mesure={sl_f:+.4f} mN/s")
            i = j + 1
        else:
            i += 1

    # --- 3. courbure F(P) sur la montee ---
    i_pk = int(np.argmax(p))
    # montee : du premier point sous pression au pic
    mm = np.arange(i_first_p, i_pk + 1)
    P_up, F_up = p[mm], fa[mm]
    # fit F = a P + b P^2 (contraint F(0)=0)
    A = np.vstack([P_up, P_up**2]).T
    coef, *_ = np.linalg.lstsq(A, F_up, rcond=None)
    a2, b2 = coef
    pred = A @ coef
    rmse2 = np.sqrt(np.mean((pred - F_up) ** 2))
    # fit lineaire seul
    a1 = float(np.dot(P_up, F_up) / np.dot(P_up, P_up))
    rmse1 = np.sqrt(np.mean((a1 * P_up - F_up) ** 2))
    # fit puissance F = c P^n
    mpos = (P_up > 0.02) & (F_up > 0.5)
    n_exp, logc = np.polyfit(np.log(P_up[mpos]), np.log(F_up[mpos]), 1)
    print(f"montee ({t[mm[0]]:.1f}->{t[i_pk]:.1f} s, {len(mm)} pts, Pmax {P_up.max():.3f} MPa) :")
    print(f"  fit F=aP     : a={a1:.1f} mN/MPa, RMSE={rmse1:.1f} mN")
    print(f"  fit F=aP+bP2 : a={a2:.1f}, b={b2:.1f}, b/a={b2/a2:+.1f} /MPa, RMSE={rmse2:.1f} mN")
    print(f"  fit F=cP^n   : n={n_exp:.2f}")
    # verification brute a mi-pression
    P_half = 0.5 * P_up.max()
    F_half = float(np.interp(P_half, P_up, F_up))  # approx (montee non monotone possible)
    F_lin = float(np.interp(P_up.max(), P_up, F_up)) * 0.5
    print(f"  F(P_max/2)={F_half:.1f} mN contre lineaire {F_lin:.1f} mN (rapport {F_half/max(F_lin,1e-9):.2f})")

    # --- 4. hysteresis (essai 2) : montee vs descente ---
    if tag == "essai2_155135":
        dn = np.arange(i_pk, n)
        P_dn, F_dn = p[dn], fa[dn]
        for p_probe in (0.1, 0.2, 0.3, 0.4, 0.5):
            # montee : interpolation sur portion monotone (tri par P)
            o_up = np.argsort(P_up)
            o_dn = np.argsort(P_dn)
            fu = float(np.interp(p_probe, P_up[o_up], F_up[o_up]))
            fd = float(np.interp(p_probe, P_dn[o_dn], F_dn[o_dn]))
            print(f"  hysteresis P={p_probe:.1f} MPa : montee {fu:7.1f} mN, descente {fd:7.1f} mN, delta={fd-fu:+6.1f} mN")
        print(f"  fin de descente : P={p[-1]*1000:.1f} kPa, F_act={fa[-1]:+.1f} mN ({100*fa[-1]/F_up.max():+.1f}% du pic)")

print("\nTermine.")
