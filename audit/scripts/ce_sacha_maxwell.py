# -*- coding: utf-8 -*-
"""Contre-expertise independante des 4 constats (3 Sacha + 1 stiffness_sanity_report)."""
import sys, time, json
import numpy as np
import pandas as pd

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base
from scipy.optimize import curve_fit

DATA = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\expérimentale\Sacha"
out = {}

# =====================================================================
# 1) MAINTIEN 20 min (muscle D) : la pression chute-t-elle ? F suit-elle P ?
# =====================================================================
raw = pd.read_excel(DATA + r"\Maintient sous pression D 20 mn.xlsx", header=0)
print("colonnes maintien:", list(raw.columns))
P = pd.to_numeric(raw.iloc[:, 0], errors="coerce").values
F = pd.to_numeric(raw.iloc[:, 1], errors="coerce").values
T = pd.to_numeric(raw.iloc[:, 3], errors="coerce").values
m = np.isfinite(P) & np.isfinite(F) & np.isfinite(T)
P, F, T = P[m], F[m], T[m]

def biexp(t, A0, A1, tau1, A2, tau2):
    return A0 + A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2)

pF, _ = curve_fit(biexp, T, F, p0=[F[-1], 0.1, 60.0, 0.15, 800.0], maxfev=60000)
pP, _ = curve_fit(biexp, T, P, p0=[P[-1], 1.0, 60.0, 1.5, 800.0], maxfev=60000)

# regression F = a*P + b sur tout l'essai (F suit-elle P lineairement ?)
A = np.vstack([P, np.ones_like(P)]).T
coef, res, _, _ = np.linalg.lstsq(A, F, rcond=None)
F_pred = A @ coef
r2 = 1 - np.sum((F - F_pred) ** 2) / np.sum((F - np.mean(F)) ** 2)

out["maintien"] = {
    "n": int(len(T)), "t_end_s": float(T[-1]),
    "P0_bar": float(P[0]), "Pend_bar": float(P[-1]),
    "P_drop_frac": float(1 - P[-1] / P[0]),
    "F0_N": float(F[0]), "Fend_N": float(F[-1]),
    "F_drop_frac": float(1 - F[-1] / F[0]),
    "corr_F_P": float(np.corrcoef(P, F)[0, 1]),
    "taus_F_s": sorted([float(pF[2]), float(pF[4])]),
    "taus_P_s": sorted([float(pP[2]), float(pP[4])]),
    "R2_F_vs_P_lineaire": float(r2),
    "pente_F_P_N_per_bar": float(coef[0]),
}
print(json.dumps(out["maintien"], indent=1))

# =====================================================================
# 2) MISES SOUS PRESSION 200 g B/I/J : facteur dF/dP
#    - rechargement independant
#    - estimateur pente par REGRESSION (pas seulement endpoints)
# =====================================================================
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

tests = {"B": load_B(),
         "I": load_csv("Mise sous pression muscle I 200 g.csv"),
         "J": load_csv("Mise sous pression muscle J 200 g.csv")}

DT = 0.2
out["200g"] = {}
for name, (Pb, Fn) in tests.items():
    t = np.arange(len(Pb)) * DT
    Pmpa = Pb / 10.0
    t0 = time.time()
    model, arr = Base.run_blocked_actuation(pressure_time=t, pressure_MPa=Pmpa)
    wall = time.time() - t0
    Fmod = arr["force_N"]

    F0e = float(np.mean(Fn[:8]))
    dFe = Fn - F0e
    dFm = Fmod - Fmod[0]

    # zone montee : jusqu'au max de pression (exclut fin depressurisee I/J)
    iPmax = int(np.argmax(Pb))
    # regression dF vs P sur la montee (points ou P > P0 + 0.2 bar)
    sel = np.arange(iPmax + 1)
    sel = sel[Pb[sel] > Pb[0] + 0.2]
    if len(sel) >= 5:
        Ae = np.vstack([Pmpa[sel], np.ones(len(sel))]).T
        ce, _, _, _ = np.linalg.lstsq(Ae, dFe[sel], rcond=None)
        cm, _, _, _ = np.linalg.lstsq(Ae, dFm[sel], rcond=None)
        slope_reg_exp, slope_reg_mod = float(ce[0]), float(cm[0])
    else:
        slope_reg_exp = slope_reg_mod = None

    dP = float(Pmpa[iPmax] - Pmpa[0])
    se = float(np.max(dFe[: iPmax + 5]) - 0.0) / dP
    sm = float(np.max(dFm[: iPmax + 5]) - 0.0) / dP
    out["200g"][name] = {
        "n": int(len(Pb)), "wall_s": round(wall, 1),
        "preload_exp_N": F0e, "preload_model_N": float(Fmod[0]),
        "dP_MPa": dP,
        "slope_endpoint_exp": se, "slope_endpoint_mod": sm,
        "ratio_endpoint": se / sm,
        "slope_reg_exp": slope_reg_exp, "slope_reg_mod": slope_reg_mod,
        "ratio_reg": (slope_reg_exp / slope_reg_mod) if slope_reg_mod else None,
    }
    print(name, json.dumps(out["200g"][name]))

    # ---- plateau (test B uniquement) ----
    if name == "B":
        ipe = int(np.argmax(dFe)); ipm = int(np.argmax(dFm))
        # verification : P est-elle constante pendant le palier ?
        plateau = slice(ipm, len(Pb))
        out["plateau_B"] = {
            "t_total_s": float(t[-1]),
            "t_peak_exp_s": float(t[ipe]), "t_peak_mod_s": float(t[ipm]),
            "P_plateau_min_bar": float(np.min(Pb[plateau])),
            "P_plateau_max_bar": float(np.max(Pb[plateau])),
            "dF_exp_peak_N": float(dFe[ipe]),
            "dF_exp_end_N": float(np.mean(dFe[-5:])),
            "perte_exp_frac": float(1 - np.mean(dFe[-5:]) / dFe[ipe]),
            "dF_mod_peak_N": float(dFm[ipm]),
            "dF_mod_end_N": float(np.mean(dFm[-5:])),
            "perte_mod_frac": float(1 - np.mean(dFm[-5:]) / dFm[ipm]),
        }
        # sensibilite a l'hypothese dt : refaire avec dt=0.1 et dt=0.5
        for dt_alt in (0.1, 0.5):
            t_alt = np.arange(len(Pb)) * dt_alt
            _, arr_a = Base.run_blocked_actuation(pressure_time=t_alt, pressure_MPa=Pmpa)
            dFa = arr_a["force_N"] - arr_a["force_N"][0]
            ipa = int(np.argmax(dFa))
            out["plateau_B"][f"perte_mod_frac_dt{dt_alt}"] = float(
                1 - np.mean(dFa[-5:]) / dFa[ipa]
            )
        print("plateau_B", json.dumps(out["plateau_B"], indent=1))

# =====================================================================
# 3) stiffness_sanity_report : 31.24 vs 37.76
# =====================================================================
rep = Base.stiffness_sanity_report()
C3124 = Base.ti_stiffness_from_paper(31.24, 8.82, 7.24, 0.205, 0.422)
C3776 = Base.ti_stiffness_from_paper(37.76, 8.82, 7.24, 0.205, 0.422)
Crot_3776 = Base.rotate_stiffness_bias(C3776, 0.7)
out["sanity"] = {
    "report_asis": {"iso_rot_err": rep.isotropic_rotation_error,
                    "sym_err": rep.rotated_stiffness_symmetry_error,
                    "min_eig": rep.local_stiffness_min_eigenvalue},
    "min_eig_31_24": float(np.linalg.eigvalsh(C3124).min()),
    "min_eig_37_76": float(np.linalg.eigvalsh(C3776).min()),
    "sym_err_37_76": float(np.max(np.abs(Crot_3776 - Crot_3776.T))),
    "somme_maxwell_sans_E2": 6.36 + 20.67 + 4.75,
    "somme_maxwell_totale": 6.36 + 20.67 + 5.98 + 4.75,
    "E_axial_defaut_moteur": Base.default_material_params().E_axial,
}
print("sanity", json.dumps(out["sanity"], indent=1))

with open(r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad\ce_sacha_maxwell.json", "w") as fh:
    json.dump(out, fh, indent=2)
print("OK")
