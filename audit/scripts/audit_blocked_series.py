# -*- coding: utf-8 -*-
"""Audit numerique independant : solveur bloque + extremites series (Base.py alpha V2)."""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

np.set_printoptions(precision=6, suppress=True)

print("MODEL_VERSION =", Base.MODEL_VERSION)

# ---------------------------------------------------------------------------
# 1) Run bloque de reference (ends rigides), verification de la fermeture Eq. (3)
# ---------------------------------------------------------------------------
model, arr = Base.run_blocked_actuation(eps=0.8, n_cycles=1, Pmax=1.3, dt=0.5,
                                        n_layers=4, n_phi=16, pre_steps=24)
res = np.abs(arr["residual"])
print("\n[1] Bloque rigide : max|residu moment| =", res.max(), "N.mm")
# fermeture manuelle Eq. (3) sur chaque pas
alpha = arr["alpha_rad"]; rho = arr["rho_mm"]
Ft = arr["force_N"]; Tt = arr["torque_Nmm"]
Ftube = arr["Ftube_axis_N"]; Fny = arr["Fnylon_axis_N"]
Mtube = arr["Mtube_Nmm"]; Mny = arr["Mnylon_Nmm"]
Ttube = arr["Ttube_axis_Nmm"]; Tny = arr["Tnylon_axis_Nmm"]
e1 = np.abs(Ftube + Fny - Ft*np.sin(alpha))
e2 = np.abs(Mtube + Mny - (Tt*np.cos(alpha) - Ft*rho*np.sin(alpha)))
e3 = np.abs(Ttube + Tny - (Tt*np.sin(alpha) + Ft*rho*np.cos(alpha)))
print("    fermeture Eq.(3) : max e1 =", e1.max(), " max e2 =", e2.max(), " max e3 =", e3.max())
print("    h constant (Eq.9) : ecart max rho.tan.alpha =", np.abs(rho*np.tan(alpha) - rho[0]*np.tan(alpha[0])).max(), "mm")
print("    F crete =", arr["force_mN"].max(), "mN ; T crete =", arr["torque_microNm"].max(), "uN.m")

# ---------------------------------------------------------------------------
# 2) Extremites compliantes : compatibilite serie + reduction de force
# ---------------------------------------------------------------------------
geom_c = Base.default_geometry_params(uncoiled_length=10.0)  # 2 x 5 mm
model2, arr2 = Base.run_blocked_actuation(eps=0.8, n_cycles=1, Pmax=1.3, dt=0.5,
                                          n_layers=4, n_phi=16, pre_steps=24, geom=geom_c)
print("\n[2] Serie (L_unc=10 mm, tangent_beam)")
print("    compliance =", model2.uncoiled_compliance_mm_per_N, "mm/N ; k_ext =", model2.uncoiled_stiffness_N_per_mm, "N/mm")
print("    max|residu moment| =", np.abs(arr2["residual"]).max())
print("    max|residu compat serie| =", np.abs(arr2["series_compatibility_residual_mm"]).max(), "mm")
# compat recomputee a la main : 2 pi N (h - h_ref) + c (F - F_ref)
h = arr2["h_mm_per_rad"]
manual = 2*np.pi*model2.turns*(h - model2.series_reference_h) + \
         model2.uncoiled_compliance_mm_per_N*(arr2["force_N"] - model2.series_reference_force_N)
print("    compat recomputee max =", np.abs(manual).max(), "mm")
print("    F crete rigide =", arr["force_mN"].max(), " vs serie =", arr2["force_mN"].max(), "mN")
print("    delta force actionnement (rigide vs serie) :",
      (arr["force_mN"].max()-arr["force_mN"][0]), "vs", (arr2["force_mN"].max()-arr2["force_mN"][0]), "mN")

# ---------------------------------------------------------------------------
# 3) Compliance tangent_beam recomputee independamment
# ---------------------------------------------------------------------------
g = model2.geom; m = model2.mat
EA = m.E_axial*np.pi*(g.Rout**2 - g.Rin**2) + m.E_nylon*np.pi*g.r_nylon**2
EI = m.E_axial*0.25*np.pi*(g.Rout**4 - g.Rin**4) + m.E_nylon*0.25*np.pi*g.r_nylon**4
a0 = np.deg2rad(g.alpha0_deg)
L2 = 0.5*model2.uncoiled_length
c_manual = 2.0*(L2*np.sin(a0)**2/EA + L2**3*np.cos(a0)**2/(3.0*EI))
print("\n[3] EA =", EA, "N ; EI =", EI, "N.mm2")
print("    compliance manuelle =", c_manual, " code =", model2.uncoiled_compliance_mm_per_N)
# variante avec alpha_tk (apres prestretch 0.8) pour mesurer l'effet du choix alpha0
a_tk = np.arctan((1.8)*np.tan(a0))
c_tk = 2.0*(L2*np.sin(a_tk)**2/EA + L2**3*np.cos(a_tk)**2/(3.0*EI))
print("    compliance si alpha_tk =", c_tk, " ecart relatif =", (c_tk-c_manual)/c_manual*100, "%")
# variante axial_rod
geom_r = Base.default_geometry_params(uncoiled_length=10.0, uncoiled_compliance_mode="axial_rod")
m3 = Base.TCPAMaxwellBlockedModel(geom=geom_r)
print("    axial_rod : compliance =", m3.uncoiled_compliance_mm_per_N, " manuelle =", 10.0/EA)

# effet du module effectif du tube torsade vs E_axial dans EA/EI des extremites
EL_layers = [v[0] for v in model2.vbar]
print("    E_L effectif par couche (tube torsade) =", np.array(EL_layers), " vs E_axial =", m.E_axial)

# ---------------------------------------------------------------------------
# 4) Etat Maxwell a la fin de la precontrainte (doit etre nul en mode elastic_tk)
# ---------------------------------------------------------------------------
m4 = Base.TCPAMaxwellBlockedModel(disc=Base.default_discretization(n_layers=4, n_phi=16, pre_steps=24, dw_bracket=(-0.05, 0.05)))
m4.prestretch_to(0.8)
print("\n[4] fin precontrainte : max|sigma_i| =", np.abs(m4.sigma_i).max(),
      " max|sigma_reference| =", np.abs(m4.sigma_reference).max(),
      " max|sigma0| =", np.abs(m4.sigma0).max())
print("    Fnylon apres prestretch =", m4.Fnylon, "N ; temps =", m4.helix.time, "s")

# ---------------------------------------------------------------------------
# 5) Multi-racines / continuite du residu en dw (P premier pas d'actionnement)
# ---------------------------------------------------------------------------
dP = 0.1; dt = 0.5
grid = np.linspace(-0.05, 0.05, 201)
vals = []
for x in grid:
    try:
        vals.append(float(m4._trial_state(x, dP, dt, m4.h_blocked)["residual"]))
    except Exception:
        vals.append(np.nan)
vals = np.array(vals)
sign_changes = np.sum(np.isfinite(vals[:-1]) & np.isfinite(vals[1:]) & (vals[:-1]*vals[1:] < 0))
print("\n[5] residu(dw) sur [-0.05,0.05] : nb changements de signe =", int(sign_changes),
      " nb NaN =", int(np.sum(~np.isfinite(vals))))
print("    residu extremes :", np.nanmin(vals), np.nanmax(vals))

# ---------------------------------------------------------------------------
# 6) dt=0 pas initial et re-zero du temps
# ---------------------------------------------------------------------------
print("\n[6] arr time[0] =", arr["time"][0], " arr2 time[0] =", arr2["time"][0],
      " force_act[0] =", arr["force_act_mN"][0])
print("    pression[0] =", arr["pressure_MPa"][0])
