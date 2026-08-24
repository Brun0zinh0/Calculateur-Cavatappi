# -*- coding: utf-8 -*-
"""Cas limites : P(0)>0, saut de pression, series compliance extreme."""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

# [A] historique mesure demarrant a P=1.2 MPa (saut instantane dt=0)
t = np.array([0.0, 1.0, 2.0, 3.0])
p = np.array([1.2, 1.3, 1.2, 1.1])
try:
    model, arr = Base.run_blocked_actuation(eps=0.8, n_layers=4, n_phi=16, pre_steps=24,
                                            pressure_time=t, pressure_MPa=p)
    print("[A] P(0)=1.2 : OK, F[0] =", arr["force_mN"][0], "mN, residu max =", np.abs(arr["residual"]).max())
    print("    force_act[0] =", arr["force_act_mN"][0], "(reference inclut le saut P(0))")
except Exception as e:
    print("[A] ECHEC :", type(e).__name__, e)

# [B] series compliance avec P(0)>0 : la reference serie est verrouillee APRES le saut rigide
geom_c = Base.default_geometry_params(uncoiled_length=10.0)
try:
    m2, arr2 = Base.run_blocked_actuation(eps=0.8, n_layers=4, n_phi=16, pre_steps=24,
                                          geom=geom_c, pressure_time=t, pressure_MPa=p)
    print("[B] serie + P(0)=1.2 : OK ; F_ref serie =", m2.series_reference_force_N,
          "N ; uncoiled_ext[0] =", arr2["uncoiled_extension_mm"][0], "mm")
except Exception as e:
    print("[B] ECHEC :", type(e).__name__, e)

# [C] compliance tres grande (extremites 60 mm) : le solveur serie converge-t-il ?
geom_x = Base.default_geometry_params(uncoiled_length=60.0)
try:
    m3, arr3 = Base.run_blocked_actuation(eps=0.8, n_cycles=1, Pmax=1.3, dt=0.5,
                                          n_layers=4, n_phi=16, pre_steps=24, geom=geom_x)
    print("[C] L_unc=60 mm : OK ; c =", m3.uncoiled_compliance_mm_per_N, "mm/N ;",
          "F crete =", arr3["force_mN"].max(), "mN ; compat max =",
          np.abs(arr3["series_compatibility_residual_mm"]).max())
    print("    variation h :", arr3["h_mm_per_rad"].min(), "-", arr3["h_mm_per_rad"].max(),
          " (h_blocked =", m3.h_blocked, ")")
except Exception as e:
    print("[C] ECHEC :", type(e).__name__, e)

# [D] saut de pression brutal 0 -> 1.5 MPa en un pas dt=0.5 (dw reste-t-il dans le bracket ?)
m4 = Base.TCPAMaxwellBlockedModel(disc=Base.default_discretization(n_layers=4, n_phi=16, pre_steps=24, dw_bracket=(-0.05, 0.05)))
m4.prestretch_to(0.8)
try:
    out = m4.step(1.5, 0.5)
    print("[D] saut 0->1.5 MPa : OK ; dw =", out.dw, " residu =", out.residual)
except Exception as e:
    print("[D] ECHEC :", type(e).__name__, e)

# [E] coherence du couplage nylon en mode serie : _nylon_axial_coupling renvoie
#     prestrain_coupling des que h_target != h_blocked (etiquette trompeuse, valeurs=1)
m5 = Base.TCPAMaxwellBlockedModel(geom=geom_c,
                                  disc=Base.default_discretization(n_layers=4, n_phi=16, pre_steps=24, dw_bracket=(-0.05, 0.05)))
print("[E] coupling(h_blocked) =", m5._nylon_axial_coupling(m5.h_blocked),
      " coupling(h_blocked*0.99) =", m5._nylon_axial_coupling(m5.h_blocked * 0.99))

# [F] periode du profil par defaut vs article
tt, pp = Base.cyclic_pressure_history(n_cycles=1, Pmax=1.3, flow_rate_mL_min=10.0, volume_mL=1.5, dt=0.5, nonlinear=False)
print("[F] duree cycle profil defaut =", tt[-1], "s (article Fig.7 : ~22 s apparent, 18 s theorique volume/debit)")
