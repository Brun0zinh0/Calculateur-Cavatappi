# -*- coding: utf-8 -*-
"""Contre-expertise : effet du mode section_update_mode fixed vs updated."""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base


def run_case(bias_profile, section_mode, label, eps=0.8, Pmax=1.3, n_cycles=1):
    geom = Base.default_geometry_params(
        bias_angle_profile=bias_profile, section_update_mode=section_mode
    )
    cfg = Base.default_simulation_config(
        geom=geom, n_cycles=n_cycles, Pmax=Pmax, dt=0.25, n_layers=4, n_phi=24,
        pre_steps=24, nonlinear_pressure=False,
    )
    model, arr = Base.run_blocked_actuation(cfg, eps=eps)
    fpk = float(np.max(arr["force_mN"]))
    tpk = float(np.max(arr["torque_microNm"]))
    f_act = float(np.max(arr["force_act_mN"]))
    t_act = float(np.max(arr["torque_act_microNm"]))
    # etat final de la section
    print(f"  {label:44s} Fmax={fpk:9.2f} mN Tmax={tpk:8.2f} uNm dF_act={f_act:8.2f} dT_act={t_act:8.2f}")
    print(f"      R_edges finaux = {np.round(model.R_edges, 4)}  theta(deg) = {np.round(np.rad2deg(model.theta_layers), 2)}")
    return fpk, tpk, f_act, t_act, model


print("=== fixed vs updated (defaut paper_linear, eps=0.8, Pmax=1.3, 1 cycle) ===")
base = run_case("paper_linear", "fixed", "defaut : paper_linear + fixed")
upd = run_case("paper_linear", "updated", "paper_linear + updated")
alt = run_case("uniform_twist", "fixed", "uniform_twist + fixed")
upd2 = run_case("uniform_twist", "updated", "uniform_twist + updated")

print("\nEcarts relatifs vs defaut (fixed) :")
for name, case in (("linear/updated", upd), ("arctan/fixed", alt), ("arctan/updated", upd2)):
    dF = 100.0 * (case[0] - base[0]) / abs(base[0])
    dT = 100.0 * (case[1] - base[1]) / abs(base[1])
    dFa = 100.0 * (case[2] - base[2]) / abs(base[2])
    dTa = 100.0 * (case[3] - base[3]) / abs(base[3])
    print(f"  {name:16s}: Fmax {dF:+6.2f}%  Tmax {dT:+6.2f}%  dF_act {dFa:+6.2f}%  dT_act {dTa:+6.2f}%")

# Amplitude du deplacement radial en un seul pas de pression (ordre de grandeur)
print("\n=== Ordre de grandeur du deplacement radial du a P=1.3 MPa (un increment) ===")
geom = Base.default_geometry_params()
mat = Base.default_material_params()
disc = Base.default_discretization(n_layers=4, n_phi=24, pre_steps=24)
model = Base.TCPAMaxwellBlockedModel(mat=mat, geom=geom, disc=disc, integration="exponential")
dP = 1.3
dw, dv = 0.0, 0.0
C_alg, hist = model._algorithmic_data(0.25)
coeffs = model._solve_radial_constants(dw, dv, dP, 0.25, C_alg, hist)
for edge, R in ((0, model.R_edges[0]), (-1, model.R_edges[-1])):
    j = 0 if edge == 0 else model.disc.n_layers - 1
    u, du = model._u_du_layer(j, R, coeffs, dv, dw, C_alg)
    print(f"  R={R:.2f} mm : u={u:+.4f} mm  (u/R = {u/R*100:+.1f} %)  du/dR={du:+.4f}")
