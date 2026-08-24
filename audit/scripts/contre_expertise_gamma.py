# -*- coding: utf-8 -*-
"""Contre-expertise finding 5 : rampe gamma=3.5 (defaut Base.run_hold_relaxation)
vs rampe lineaire (defaut interface run_relaxation_model)."""
import sys
import numpy as np

sys.path.insert(0, r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/alpha V2")
import Base  # noqa: E402

# Chemin 1 : Base.run_hold_relaxation (defaut -> nonlinear_ramp=True, gamma=3.5)
t1, p1 = Base.ramp_hold_pressure_history(P_hold=1.0, ramp_time=9.0, hold_time=200.0, dt=2.0)
# Chemin 2 : interface (nonlinear_ramp=False par defaut)
t2, p2 = Base.ramp_hold_pressure_history(P_hold=1.0, ramp_time=9.0, hold_time=200.0, dt=2.0,
                                         nonlinear_ramp=False)
print("Profils identiques ?", np.allclose(p1, p2))
print("P a t=4.5 s : gamma=3.5 ->", p1[np.searchsorted(t1, 4.5)], "; lineaire ->", p2[np.searchsorted(t2, 4.5)])

cfg = Base.default_simulation_config(eps=0.5, Pmax=1.0, dt=2.0, n_layers=3, n_phi=8, pre_steps=4)
_, a1 = Base.run_blocked_actuation(cfg, pressure_time=t1, pressure_MPa=p1)
_, a2 = Base.run_blocked_actuation(cfg, pressure_time=t2, pressure_MPa=p2)
i0 = int(np.searchsorted(a1["time"], 9.0, side="left"))
F1, F2 = a1["force_total_mN"], a2["force_total_mN"]
print(f"Force au debut du maintien : gamma=3.5 -> {F1[i0]:.2f} mN ; lineaire -> {F2[i0]:.2f} mN ; ecart {F1[i0]-F2[i0]:+.2f} mN")
print(f"Force fin de maintien       : gamma=3.5 -> {F1[-1]:.2f} mN ; lineaire -> {F2[-1]:.2f} mN ; ecart {F1[-1]-F2[-1]:+.2f} mN")
rel1 = F1[-1] - F1[i0]
rel2 = F2[-1] - F2[i0]
print(f"Relaxation mesuree pendant le maintien : gamma=3.5 -> {rel1:+.2f} mN ; lineaire -> {rel2:+.2f} mN")
