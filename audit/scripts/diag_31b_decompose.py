# -*- coding: utf-8 -*-
"""Phase 3.1b — décomposition du ratchet : quel terme de u ne se referme pas ?

Un seul cycle P: 0 -> 1.4 -> 0 (linéaire, 18 pas), quasi élastique, updated.
À chaque pas, on décompose le déplacement radial appliqué aux bords :
  u = [C1 R^mu + C2 R^-mu]  (réponse homogène, pilotée par la pression/CL)
    + [A dv R^2]            (couplage vrille)
    + [B dw R]              (couplage axial)
et on somme par demi-cycle (charge / décharge).
"""
import sys
from pathlib import Path

import numpy as np

APP = Path(__file__).resolve().parent / "alpha V3"
sys.path.insert(0, str(APP))

import Base  # noqa: E402
import parametres  # noqa: E402

s = dict(parametres.DEFAULT_SETTINGS)
s.update(
    {
        "eps": 0.8,
        "n_cycles": 1,
        "use_fixed_duration": False,
        "p_max_mpa": 1.4,
        "dt": 0.5,
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
disc = Base.default_discretization(n_layers=cfg.n_layers, n_phi=cfg.n_phi, pre_steps=cfg.pre_steps, dw_bracket=(-0.05, 0.05))
model = Base.TCPAMaxwellBlockedModel(
    mat=cfg.mat, geom=cfg.geom, disc=disc, integration=cfg.integration,
    prestrain_reference_mode=cfg.prestrain_reference_mode,
)
model.prestretch_to(cfg.eps)
t, p = Base.cyclic_pressure_history(n_cycles=1, Pmax=1.4, flow_rate_mL_min=10.0, volume_mL=1.5, dt=0.5, nonlinear=False)
model.step(float(p[0]), 0.0, h_target=model.h_blocked)
model.lock_blocked_series_reference()

print(f"{'k':>3} {'P':>6} {'dP':>7} {'dw':>10} {'dv':>10} | {'u_in':>10} {'hom_in':>10} {'Adv_in':>10} {'Bdw_in':>10} | {'u_out':>10}")

sums = {"load": np.zeros(4), "unload": np.zeros(4)}
Rin0 = float(model.R_edges[0])

for k in range(1, len(t)):
    dt_k = float(t[k] - t[k - 1])
    P_new = float(p[k])
    dP = P_new - model.helix.pressure
    # reproduire la resolution du pas pour decomposer u AVANT de committer
    dw = model._find_dw(dP, dt_k, model.h_blocked)
    rho_new, alpha_new, dw2, dv, dkappa = model._kinematic_increments(dw, model.h_blocked)
    C_alg, hist = model._algorithmic_data(dt_k)
    coeffs = model._solve_radial_constants(dw2, dv, dP, dt_k, C_alg, hist)

    def decompose(edge_index, layer):
        R = float(model.R_edges[edge_index])
        A, B, mu = model._layer_AB_mu(layer, C_alg)
        C1, C2 = coeffs[layer]
        hom = C1 * R**mu + C2 * R ** (-mu)
        adv = A * dv * R * R
        bdw = B * dw2 * R
        return hom + adv + bdw, hom, adv, bdw

    u_in, hom_in, adv_in, bdw_in = decompose(0, 0)
    u_out, _, _, _ = decompose(len(model.R_edges) - 1, model.disc.n_layers - 1)
    phase = "load" if dP > 0 else "unload"
    sums[phase] += np.array([u_in, hom_in, adv_in, bdw_in])
    print(f"{k:>3} {P_new:>6.3f} {dP:>7.3f} {dw2:>10.3e} {dv:>10.3e} | {u_in:>10.3e} {hom_in:>10.3e} {adv_in:>10.3e} {bdw_in:>10.3e} | {u_out:>10.3e}")

    model.step(P_new, dt_k, h_target=model.h_blocked)

print("\n=== bilan par demi-cycle (bord interne) ===")
print(f"{'':>10} {'u_total':>12} {'homogene':>12} {'A.dv.R2':>12} {'B.dw.R':>12}")
for phase in ("load", "unload"):
    v = sums[phase]
    print(f"{phase:>10} {v[0]:>12.4e} {v[1]:>12.4e} {v[2]:>12.4e} {v[3]:>12.4e}")
net = sums["load"] + sums["unload"]
print(f"{'NET/cycle':>10} {net[0]:>12.4e} {net[1]:>12.4e} {net[2]:>12.4e} {net[3]:>12.4e}")
print(f"\nRin initial {Rin0:.5f} -> final {float(model.R_edges[0]):.5f} mm "
      f"(net mesure {float(model.R_edges[0]) - Rin0:+.5e} mm)")
print(f"helice : rho {model.helix.rho:.6f} mm, alpha {np.rad2deg(model.helix.alpha):.4f} deg "
      f"(retour attendu vers l'etat post-precontrainte)")
print(f"somme dw sur le cycle : charge+decharge = {sums['load'][3]/(  -0.9737*Rin0):.3e} (via Bdw) — voir net ci-dessus")
