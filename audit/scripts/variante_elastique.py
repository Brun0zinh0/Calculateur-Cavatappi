# -*- coding: utf-8 -*-
"""Test de forme : variante quasi elastique du modele (branches de Maxwell
annulees, E0 porte au module instantane 37.76 MPa) sur l'essai 2 de Bruno.

Ce N'EST PAS le modele de l'article : c'est un test diagnostique pour verifier
que l'ecart de forme vient de la viscoelasticite et non de la cinematique.
"""
import sys
import time
import numpy as np

ALPHA = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2"
BRUNO = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\expérimentale\Bruno"
OUT = r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"

sys.path.insert(0, ALPHA)
import Base  # noqa: E402
import pression  # noqa: E402

with open(f"{BRUNO}\\essai_force_pression_20260730_155135.csv", "rb") as fh:
    cols = pression.parse_uploaded_numeric_csv(fh.read())
payload = pression.measured_pressure_payload(cols, "time_s", "pressure_bar", "bar", subtract_initial=True)

maxwell_el = Base.default_maxwell_tensile_params(E0=37.76, E1=1e-9, E2=1e-9, E3=1e-9)
mat_el = Base.default_material_params(maxwell=maxwell_el)
t0 = time.perf_counter()
model, arr = Base.run_blocked_actuation(
    mat=mat_el, pressure_time=payload["time"], pressure_MPa=payload["pressure_MPa"]
)
print(f"elastique: {time.perf_counter()-t0:.0f} s")
np.savez(
    f"{OUT}\\sim_elastique_155135.npz",
    t_sim=arr["time"], p_sim=arr["pressure_MPa"],
    F_sim_total=arr["force_total_mN"], F_sim_act=arr["force_act_mN"],
)
print(f"F0={arr['force_total_mN'][0]:.1f} mN  dFmax={arr['force_act_mN'].max():.1f} mN")
print("FINI_EL")
