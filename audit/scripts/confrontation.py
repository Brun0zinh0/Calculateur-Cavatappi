# -*- coding: utf-8 -*-
"""Confrontation modele (run_blocked_actuation, specimen article) vs essais Bruno.

Chemin identique a l'interface : pression.parse_uploaded_numeric_csv +
pression.measured_pressure_payload (bar -> MPa, soustraction du zero,
clamp >= 0), puis Base.run_blocked_actuation(pressure_time=..., pressure_MPa=...)
avec la configuration par defaut (= specimen de l'article BLOCKED).
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

FILES = [
    "essai_force_pression_20260730_153353.csv",
    "essai_force_pression_20260730_155135.csv",
]

for fname in FILES:
    tag = fname.replace("essai_force_pression_", "").replace(".csv", "")
    with open(f"{BRUNO}\\{fname}", "rb") as fh:
        raw = fh.read()
    cols = pression.parse_uploaded_numeric_csv(raw)
    # historique de pression : meme chemin que l'interface (bar -> MPa, zero soustrait)
    payload = pression.measured_pressure_payload(
        cols, "time_s", "pressure_bar", "bar", subtract_initial=True
    )
    t_meas = payload["time"]
    p_meas = payload["pressure_MPa"]
    # force mesuree (mN), realignee sur la meme base de temps dedupliquee
    exp = pression.experimental_force_pressure_payload(
        cols, "time_s", "pressure_bar", "force_mN", "bar", "mN"
    )
    F_meas = exp["force_mN"]
    assert len(F_meas) == len(t_meas), (len(F_meas), len(t_meas))

    print(f"=== {tag} : {len(t_meas)} pas, duree {t_meas[-1]:.1f} s, Pmax {p_meas.max():.3f} MPa", flush=True)
    t0 = time.perf_counter()
    model, arr = Base.run_blocked_actuation(pressure_time=t_meas, pressure_MPa=p_meas)
    el = time.perf_counter() - t0
    print(f"    simulation terminee en {el:.0f} s", flush=True)
    np.savez(
        f"{OUT}\\sim_{tag}.npz",
        t_meas=t_meas, p_meas=p_meas, F_meas=F_meas,
        t_sim=arr["time"], p_sim=arr["pressure_MPa"],
        F_sim_total=arr["force_total_mN"], F_sim_act=arr["force_act_mN"],
        F_tube=arr.get("force_tube_mN", np.array([])),
        F_nylon=arr.get("force_nylon_mN", np.array([])),
        residual=arr.get("residual", np.array([])),
    )
    print(f"    F0_sim={arr['force_total_mN'][0]:.1f} mN  dFmax_sim={arr['force_act_mN'].max():.1f} mN  "
          f"F0_meas={F_meas[0]:.1f} mN  dFmax_meas={(F_meas - F_meas[0]).max():.1f} mN", flush=True)
print("FINI")
