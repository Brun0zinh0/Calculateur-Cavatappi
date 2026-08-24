# -*- coding: utf-8 -*-
"""Phase 3.1c — ordre de l'erreur de fermeture de cycle (quasi élastique, updated).

Un cycle P: 0 -> 1.4 -> 0, trois pas de temps (h, h/2, h/4), avec et sans le
correcteur de point milieu. Mesure : écarts nets de Rin, Rout, aire de paroi,
rho, alpha, force après le cycle fermé.
  - erreur ~ O(h)  -> biais d'intégration explicite (corrigible par schéma)
  - erreur ~ O(1)  -> incohérence de formulation (bug d'équation)
"""
import sys
from pathlib import Path

import numpy as np

APP = Path(__file__).resolve().parent / "alpha V3"
sys.path.insert(0, str(APP))

import Base  # noqa: E402
import parametres  # noqa: E402


def one_cycle(dt, use_corrector):
    s = dict(parametres.DEFAULT_SETTINGS)
    s.update(
        {
            "eps": 0.8,
            "n_cycles": 1,
            "use_fixed_duration": False,
            "p_max_mpa": 1.4,
            "dt": dt,
            "n_layers": 3,
            "n_phi": 16,
            "pre_steps": 24,
            "nonlinear_pressure": False,
            "section_update_mode": "updated",
            "maxwell_eta1_mpa_s": 154.57e9,
            "maxwell_eta2_mpa_s": 977.79e9,
            "maxwell_eta3_mpa_s": 11044.83e9,
        }
    )
    cfg = parametres.build_config(s)
    disc = Base.default_discretization(
        n_layers=cfg.n_layers, n_phi=cfg.n_phi, pre_steps=cfg.pre_steps, dw_bracket=(-0.05, 0.05)
    )
    model = Base.TCPAMaxwellBlockedModel(
        mat=cfg.mat, geom=cfg.geom, disc=disc, integration=cfg.integration,
        prestrain_reference_mode=cfg.prestrain_reference_mode,
    )
    if not use_corrector:
        model._midpoint_corrected_coeffs = (
            lambda coeffs, dw, dv, dP, dt_, C_alg, hist: (coeffs, None)
        )
    model.prestretch_to(cfg.eps)
    t, p = Base.cyclic_pressure_history(
        n_cycles=1, Pmax=1.4, flow_rate_mL_min=10.0, volume_mL=1.5, dt=dt, nonlinear=False
    )
    model.step(float(p[0]), 0.0, h_target=model.h_blocked)
    model.lock_blocked_series_reference()
    ref = dict(
        Rin=float(model.R_edges[0]),
        Rout=float(model.R_edges[-1]),
        area=float(np.pi * (model.R_edges[-1] ** 2 - model.R_edges[0] ** 2)),
        rho=float(model.helix.rho),
        alpha=float(model.helix.alpha),
        F=float(model.history[-1].Ft),
        umax=0.0,
    )
    for k in range(1, len(t)):
        prev = model.R_edges.copy()
        model.step(float(p[k]), float(t[k] - t[k - 1]), h_target=model.h_blocked)
        ref["umax"] = max(ref["umax"], float(np.max(np.abs(model.R_edges - prev))))
    return dict(
        dRin=float(model.R_edges[0]) - ref["Rin"],
        dRout=float(model.R_edges[-1]) - ref["Rout"],
        darea=float(np.pi * (model.R_edges[-1] ** 2 - model.R_edges[0] ** 2)) - ref["area"],
        drho=float(model.helix.rho) - ref["rho"],
        dalpha_deg=float(np.rad2deg(model.helix.alpha - ref["alpha"])),
        dF_mN=1000.0 * (float(model.history[-1].Ft) - ref["F"]),
        umax=ref["umax"],
        n_steps=len(t) - 1,
    )


print(f"{'corr':>5} {'dt':>5} {'pas':>4} | {'dRin':>11} {'dRout':>11} {'dAire':>11} {'drho':>11} {'dalpha':>9} {'dF(mN)':>9} | {'u_max/pas':>10}")
for use_corr in (False, True):
    prev = None
    for dt in (1.0, 0.5, 0.25):
        r = one_cycle(dt, use_corr)
        ratio = "" if prev is None else f"  (x{r['darea'] / prev:.2f})" if prev != 0 else ""
        print(
            f"{str(use_corr):>5} {dt:>5.2f} {r['n_steps']:>4} | {r['dRin']:>11.3e} {r['dRout']:>11.3e} "
            f"{r['darea']:>11.3e} {r['drho']:>11.3e} {r['dalpha_deg']:>9.2e} {r['dF_mN']:>9.2f} | {r['umax']:>10.3e}{ratio}",
            flush=True,
        )
        prev = r["darea"]
print("\nlecture : si dAire est divisee par ~2 quand dt est divise par 2 -> O(h) ;")
print("si dAire reste ~constante -> incoherence de formulation (bug).")
