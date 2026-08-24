# -*- coding: utf-8 -*-
"""Contre-expertise cote modele.

A. Refit des sorties sauvegardees (essai1/2 _sim.csv et _el.csv) :
   - courbure b/a du modele complet et de la variante quasi-elastique
   - chute de force du modele complet a P quasi constant (t=57 -> 128 s)
B. Simulations independantes grossieres (n_layers=2, n_phi=8) :
   - rampe 0->0.39 MPa en 55 s puis maintien 75 s : relaxation du modele complet
   - rampe lente 0->0.6 MPa : courbure F(P) en mode fixed / updated / quasi-elastique
"""
import sys
import time as _time
import numpy as np
import pandas as pd

ALPHA = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2"
OUT = r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
sys.path.insert(0, ALPHA)
import Base as modele  # noqa: E402


def fit_ba(P, F):
    A = np.vstack([P, P**2]).T
    coef, *_ = np.linalg.lstsq(A, F, rcond=None)
    return coef[0], coef[1]


print("=== A. Refit des sorties sauvegardees ===")
for tag in ("essai1_153353", "essai2_155135"):
    sim = pd.read_csv(OUT + f"\\{tag}_sim.csv")
    el = pd.read_csv(OUT + f"\\{tag}_el.csv")
    t = sim["time_s"].to_numpy(float)
    p = sim["pressure_MPa"].to_numpy(float)
    fm = sim["force_act_meas_mN"].to_numpy(float)
    fv = sim["force_act_sim_mN"].to_numpy(float)
    fe = el["force_act_el_mN"].to_numpy(float)
    i_pk = int(np.argmax(p))
    up = slice(0, i_pk + 1)
    a_m, b_m = fit_ba(p[up], fm[up])
    a_v, b_v = fit_ba(p[up], fv[up])
    a_e, b_e = fit_ba(p[up], fe[up])
    print(f"[{tag}] montee jusqu'au pic (t={t[i_pk]:.1f} s, P={p[i_pk]:.3f} MPa)")
    print(f"   mesure       : a={a_m:7.1f} b={b_m:7.1f} b/a={b_m/a_m:+6.1f} /MPa")
    print(f"   modele Maxwell: a={a_v:7.1f} b={b_v:7.1f} b/a={b_v/a_v:+6.2f} /MPa")
    print(f"   quasi-elast. : a={a_e:7.1f} b={b_e:7.1f} b/a={b_e/a_e:+6.2f} /MPa")
    if tag == "essai1_153353":
        for tt in (57.0, 128.0):
            i = int(np.argmin(np.abs(t - tt)))
            print(f"   t={t[i]:6.1f} s  P={p[i]:.3f} MPa  F_maxwell={fv[i]:6.1f} mN  F_mesure={fm[i]:6.1f} mN  F_el={fe[i]:6.1f} mN")

print("\n=== B. Simulations independantes grossieres ===")
common = dict(n_layers=2, n_phi=8, pre_steps=16)

# B1 : rampe + maintien, modele complet
t1 = np.concatenate([np.arange(0.0, 55.0, 1.0), np.arange(55.0, 130.0, 1.0)])
p1 = np.where(t1 < 55.0, 0.39 * t1 / 55.0, 0.39)
t0 = _time.perf_counter()
_, arr = modele.run_blocked_actuation(pressure_time=t1, pressure_MPa=p1, **common)
print(f"[B1 rampe+maintien, Maxwell complet] simule en {_time.perf_counter()-t0:.1f} s")
tt, ff, pp = arr["time"], arr["force_act_mN"], arr["pressure_MPa"]
i55 = int(np.argmin(np.abs(tt - 55.0)))
i128 = int(np.argmin(np.abs(tt - 128.0)))
print(f"   F(fin rampe t=55 s) = {ff[i55]:.2f} mN ; F(t=128 s, P constant) = {ff[i128]:.2f} mN ; variation {100*(ff[i128]-ff[i55])/ff[i55]:+.1f}%")
sl = np.polyfit(tt[i55 + 5 : i128], ff[i55 + 5 : i128], 1)[0]
print(f"   pente moyenne au maintien : {sl:+.4f} mN/s (brute, avant tout facteur d'echelle)")

# B2 : rampe lente 0->0.6 MPa en 60 s, trois variantes
t2 = np.arange(0.0, 61.0, 1.0)
p2 = 0.6 * t2 / 60.0
mw = modele.default_maxwell_tensile_params()
mw_frozen = modele.default_maxwell_tensile_params(eta1=mw.eta1 * 1e9, eta2=mw.eta2 * 1e9, eta3=mw.eta3 * 1e9)

variants = {
    "Maxwell complet, section fixe": dict(),
    "quasi-elastique, section fixe": dict(mat=modele.default_material_params(maxwell=mw_frozen)),
    "quasi-elastique, section updated": dict(
        mat=modele.default_material_params(maxwell=mw_frozen),
        geom=modele.default_geometry_params(section_update_mode="updated"),
    ),
}
for label, kw in variants.items():
    t0 = _time.perf_counter()
    _, arr = modele.run_blocked_actuation(pressure_time=t2, pressure_MPa=p2, **common, **kw)
    p_r, f_r = arr["pressure_MPa"], arr["force_act_mN"]
    m = p_r > 1e-6
    a, b = fit_ba(p_r[m], f_r[m])
    n_exp = np.polyfit(np.log(p_r[(p_r > 0.02) & (f_r > 0.01)]), np.log(np.maximum(f_r[(p_r > 0.02) & (f_r > 0.01)], 1e-9)), 1)[0]
    print(f"[B2 {label}] a={a:.1f} b={b:+.1f} b/a={b/a:+.3f} /MPa, exposant n={n_exp:.3f}, F(0.3)={np.interp(0.3, p_r, f_r):.1f} F(0.6)={np.interp(0.6, p_r, f_r):.1f} mN ({_time.perf_counter()-t0:.1f} s)")

print("\nTermine.")
